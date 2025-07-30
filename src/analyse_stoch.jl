### analyse stochastic lake depth and plot results
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
using Rasters
using DataFrames
using CSV
using Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows
using Serialization

include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis

# --- Paths ---
name = "2024_October"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# --- Load aggregated Monte Carlo result ---
# all unc.
aggr1 = deserialize(joinpath(output_dir, "aggr1_2024_October.jls"))
# only bed unc.
aggr2 = deserialize(joinpath(output_dir, "aggr2_2024_October.jls"))
# only surf unc.
aggr3 = deserialize(joinpath(output_dir, "aggr3_2024_October.jls"))
# only floatfrac unc.
aggr4 = deserialize(joinpath(output_dir, "aggr4_2024_October.jls"))

# Reconstruct rasters from arrays
thickness = clean_raster(Raster(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif")))
lake_depth_mean = Raster(aggr1.lakes_depth_fs, dims(thickness))
area_stoch = Raster(aggr1.areas,  dims(thickness))



# --- Analyze ---
println("  Analyzing stochastic lake depth (free surface)...")

analysis = analyze_lakes(lake_depth_mean, thickness; min_depth=2.0)

df = DataFrame(
    run = name,
    smooth_surface_ice_fraction = NaN,  # not used here
    min_depth_m = analysis.min_depth,
    supragl_fill_fraction = NaN,        # not used here
    n_lakes = nrow(analysis.stats),
    total_volume_m3 = sum(analysis.stats.volume),
    mean_area_m2 = mean(analysis.stats.area_m2),
    mean_depth_m = mean(collect(lake_depth_mean)[analysis.labels .> 0]),
    max_depth_m = maximum(collect(lake_depth_mean)[analysis.labels .> 0]),
    largest_single_volume_m3 = analysis.LargestLake.volume
)

CSV.write(joinpath(output_dir, "WWFS_stoch_lake_summary.csv"), df)

# --- Plot ---
# mean lake depth and mean area
plot_lake_depth(
    lake_depth_mean,
    thickness,
    analysis,
    nothing, # phi
    joinpath(output_dir, "stochastic_lake_depth.png");
    min_depth = analysis.min_depth,
    show_all_lakes = true,
    area = area_stoch,
    area_threshold = 1e4
)

# Distribution of lake volumes
boxplot_lake_vol_stoch(aggr1;
    aggr2 = aggr2,
    aggr3 = aggr3,
    aggr4 = aggr4,
    savepath = joinpath(output_dir, "lake_volume_comparison.png")
)


# plot boxplot with the 4 sub boxplot: all contributions, surf, bed, f


