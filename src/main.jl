### main.jl
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
import ArchGDAL
using Dates
using Rasters
using DataFrames
using CSV
using Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows

include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis

datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"
mkpath(output_dir)

runs = [
    (
        name = "2021_July",
        surface_path = joinpath(datadir_WWFS_input, "surface_2021_cr.tif"),
        thickness_path = joinpath(datadir_WWFS_input, "ice_thickness_2021_july.tif")
    ),
    (
        name = "2024_June",
        surface_path = joinpath(datadir_WWFS_input, "surface_2024_june.tif"),
        thickness_path = joinpath(datadir_WWFS_input, "ice_thickness_2024_june.tif")
    ),
    (
        name = "2024_October",
        surface_path = joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"),
        thickness_path = joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif")
    )
]

# TO CHANGE
min_depths = 2.#[0.0, 2.0]           # meters
smooth_coeffs = 0.1#[0.0, 0.1]        # as fraction of thickness
#filling_fractions = 0.#[0.0, 0.75]   # fraction of supraglacial filling
filling_volume = [0.,100000.0] # m3, volume to fill the largest supraglacial lake

# Summary
summaries = DataFrame()

for min_depth in min_depths
    for smooth_coeff in smooth_coeffs
        for fill_vol in filling_volume
            for run in runs[2:2]

                println("\n🔷 Processing run: $(run.name), fill_vol=$(fill_vol)")

                # Load Rasters
                surface = clean_raster(Raster(run.surface_path))
                thickness = clean_raster(Raster(run.thickness_path))
                bedrock = surface - thickness

                # Grid spacing
                x, y = dims(surface)
                dx = step(x)
                println("  Grid spacing dx = ", dx, " m")

                # Smoothing surface
                if smooth_coeff > 0
                    println("  Smoothing surface with coefficient: ", smooth_coeff)
                    mask = bedrock .> 0
                    smooth_half_window = smooth_coeff / 2
                    y = x # WWFS expects square grid
                    surface = WWFS.smooth_surface(x, y, surface, bedrock, smooth_half_window, mask)
                else
                    println("  No smoothing applied.")
                end

                surface_fill = copy(surface)  # start with a full copy
                
                if fill_vol > 0
                    # Supraglacial lake filling
                    (_, _, dir, _, _, sinks, _, _, _) = WWF.waterflows(surface)
                    surf_fill = fill_dem(surface, sinks, dir)
                    lake_surf = surf_fill .- surface

                    # Analyze initial lake
                    analysis_supra = analyze_lakes(lake_surf, thickness)
                    supralake_volume_m3 = analysis_supra.LargestLake.volume
                    supralake_area_m2 = analysis_supra.LargestLake.area
                    supralake_mask = analysis_supra.LargestLake.mask

                    # Adjust lake surface until target volume is reached
                    c = 0
                    while supralake_volume_m3 > fill_vol
                        println("  Largest supraglacial lake volume: ", supralake_volume_m3, " m³") 
                        println("  Largest supraglacial lake area: ", supralake_area_m2, " m²")
                        lake_surf .-= 0.1   # lower by 10 cm
                        lake_surf = max.(lake_surf, 0.0)  # avoid negative values
                        # or should we control the lake area instead?
                        analysis_supra = analyze_lakes(lake_surf, thickness)
                        supralake_volume_m3 = analysis_supra.LargestLake.volume
                        supralake_mask = analysis_supra.LargestLake.mask
                        c = c + 0.1
                    end
                    println("lake lowering =", c ," m" )
                    supralake_area_m2 = sum(supralake_mask)
                    println("  Largest new supraglacial lake volume: ", supralake_volume_m3, " m³") 
                    println("  Largest new supraglacial lake area: ", supralake_area_m2, " m²")

                    # Initialize surface_fill first
                    surface_fill = copy(surface)

                    # convert mask in bollean

                    # Update only lake pixels with hydro-converted water height
                    surface_fill[supralake_mask] .= surface[supralake_mask] .+ lake_surf[supralake_mask] ./ 0.9
                end

                # Run WWFS
                ((areas, slen, dir, nout, nin, sinks, pits, c, bnds),
                 (sc_locs, kappas, diro, phi),
                 (lakes, lakes_free_surf),
                 sinkout) = WWFS.waterflows_subglacial(
                    surface_fill, bedrock, dx;
                    gamma = [0, WWFS.GAMMA][1],
                    bnd_as_sink = true,
                    drain_pits = true
                )

                # Output naming
                fill_id = "_fill$(string(round(fill_vol)))"
                run_id = run.name * "_md$(Int(min_depth))_sm$(replace(string(smooth_coeff), "." => ""))" * fill_id
                out_prefix = joinpath(output_dir, run_id)


                # Save outputs
                if fill_vol > 0 # otherwise lake_surf not defined
                    supra_lake = surface_fill .- surface
                    write(joinpath(output_dir, run_id * "_filledsupralake.tif"), supra_lake; force=true)
                end
                #write(out_prefix * "_lakes_free.tif", lakes_free_surf; force=true)
                #write(out_prefix * "_phi.tif", phi; force=true)
                #write(out_prefix * "_area.tif", areas[1]; force=true)
                #println("  ➤ Saved 'lakes_free_surf' and 'phi' rasters.")

                # Analyze lakes
                println("  Analyzing lake_free_surface...")
                analysis = analyze_lakes(lakes_free_surf, thickness; min_depth=min_depth)

                df = DataFrame(
                    run = run.name,
                    smooth_surface_ice_fraction = smooth_coeff,
                    min_depth_m = analysis.min_depth,
                    supragl_fill_volume_m3 = fill_vol,
                    n_lakes = nrow(analysis.stats),
                    total_volume_m3 = sum(analysis.stats.volume),
                    mean_area_m2 = mean(analysis.stats.area_m2),
                    mean_depth_m = mean(collect(lakes_free_surf)[analysis.labels .> 0]),
                    max_depth_m = maximum(collect(lakes_free_surf)[analysis.labels .> 0]),
                    largest_single_volume_m3 = analysis.LargestLake.volume
                )

                # Plotting

                plot_lake_depth(lakes_free_surf,thickness,analysis,
                    phi,out_prefix * "_lakes.png";min_depth = analysis.min_depth,show_all_lakes = true, area = areas[1])

                #plot_hydraulic_head(phi, out_prefix * "_phi.png")

                # plot upslope area
                #plot_hydraulic_head_and_flux(phi,areas[1],thickness,out_prefix * "_phi_flux.png";min_threshold = 1e5,max_threshold = 1e6)

                                
                # Append to summary
                append!(summaries, df)

                println("  ✅ Done with run: ", run.name)
            end
        end
    end
end

# Save summary
CSV.write(joinpath(output_dir, "WWFS_lake_summary.csv"), summaries)
println("\n✅ All runs complete. Summary saved to: ", output_dir)


#heatmap_mindepth_vs_smoothing(summaries, "2024_October", :total_volume_m3, joinpath(output_dir, "2024_October_heatmap_vol.png"))
#heatmap_mindepth_vs_smoothing(summaries, "2024_October", :n_lakes, joinpath(output_dir, "2024_October_heatmap_nlakes.png"))
