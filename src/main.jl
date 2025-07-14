### main.jl
# Activate your Julia project environment
using Pkg
Pkg.activate("/scratch-1/cogier/hydraulic_barriers/")
import ArchGDAL
using Dates
using Rasters
using DataFrames
using CSV
using Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows

# Include lake analysis module
include("LakeAnalysis.jl")
include("plot.jl")
using .LakeAnalysis

# Define input and output directories
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"
mkpath(output_dir) # create and do nothing if it already exists



# Define runs to process
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

    # here possibly add surface dems filled, or smooth, etc...
    # TODO add runs with variation in smoothing, minimal lake depth treshold, 
    # and write this in csv file and find a good visualization

]

# Initialize summary DataFrame
summaries = DataFrame()

# Loop over each run
for run in runs[1:1]

    println("\n🔷 Processing run: ", run.name)

    # Load Rasters
    surface = Raster(run.surface_path)
    thickness = Raster(run.thickness_path)
    bedrock = surface - thickness

    # Grid spacing
    dx = step(dims(surface)[1])
    println("  Grid spacing dx = ", dx, " m")

    # Run WWFS and capture all outputs
    ((areas, slen, dir, nout, nin, sinks, pits, c, bnds),
     (sc_locs, kappas, diro, phi),
     (lakes, lakes_free_surf),
     sinkout) = WWFS.waterflows_subglacial(
        surface, bedrock, dx;
        gamma=[0,WWFS.GAMMA][1],
        bnd_as_sink = true,
        drain_pits = true
    )

    # Define output prefix
    out_prefix = joinpath(output_dir, run.name)

    # Save selected WWFS outputs
    write(joinpath(output_dir, run.name * "_lakes_free.tif"), lakes_free_surf; force=true)
    write(joinpath(output_dir, run.name * "_phi.tif"), phi; force=true)

    println("  ➤ Write 'lakes_free_surf' and 'phi' rasters.")

    # Analyze lakes
    println("  Analyzing lake_free_surface...")
    analysis = analyze_lakes(lakes_free_surf, surface, thickness, bedrock)

    # Compute metrics
    n_lakes = maximum(analysis.labels)
    mean_area = mean(analysis.stats.area_m2)
    total_volume = sum(analysis.stats.volume)
    depth = analysis.stats.volume ./ analysis.stats.area_m2
    median_depth = median(depth)
    max_depth = maximum(depth)

    # Add to summary table
    df = DataFrame(
        run = run.name,
        n_lakes = n_lakes,
        mean_area_m2 = mean_area,
        total_volume_m3 = total_volume,
        median_depth_m = median_depth,
        max_depth_m = max_depth
    )
    append!(summaries, df)

    # Save labeled lakes (optional)
    #write(out_prefix * "_lake_labels.tif", Raster(analysis.labels); force=true)

    #plot
    plot_lake_depth(lakes_free_surf, thickness, joinpath(output_dir, run.name * "_lakes_free.png"))
    plot_hydraulic_head(phi, joinpath(output_dir, run.name * "_phi.png"))



    println("  ✅ Done with run: ", run.name)
end

# Save summary
CSV.write(joinpath(output_dir, "WWFS_lake_summary.csv"), summaries)
println("\n✅ All runs complete. Summary saved to: ", output_dir)





