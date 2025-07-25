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
filling_fractions = 0.# [0.0, 0.75]   # fraction of supraglacial filling

# Summary
summaries = DataFrame()

for min_depth in min_depths
    for smooth_coeff in smooth_coeffs
        for fill_frac in filling_fractions
            for run in runs

                println("\n🔷 Processing run: $(run.name), fill_frac=$(fill_frac)")

                # Load Rasters
                surface = Raster(run.surface_path)
                thickness = Raster(run.thickness_path)
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

                # Supraglacial lake filling
                (_, _, dir, _, _, sinks, _, _, _) = WWF.waterflows(surface)
                surf_fill = fill_dem(surface, sinks, dir)
                lake_surf = surf_fill .- surface
                if fill_frac == 0.0
                    surface_fill = surface
                else
                    surface_fill = surface .+ fill_frac * lake_surf
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
                fill_id = "_fill$(replace(string(Int(round(fill_frac * 100))), "." => ""))"
                run_id = run.name * "_md$(Int(min_depth))_sm$(replace(string(smooth_coeff), "." => ""))" * fill_id
                out_prefix = joinpath(output_dir, run_id)


                # Save outputs
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
                    supragl_fill_fraction = fill_frac,
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
