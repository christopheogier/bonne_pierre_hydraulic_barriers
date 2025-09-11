### analyse stochastic lake depth and plot results
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
using Rasters, DataFrames, CSV, Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF  = WhereTheWaterFlows
using Serialization

include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis

# --- Params / Paths ---
name = "2024_June"   # or "2024_October"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir         = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# Map run name -> thickness file
thickness_file = Dict(
    "2024_October" => "ice_thickness_2024_oct.tif",
    "2024_June"    => "ice_thickness_2024_june.tif",
)[name]

# --- Load aggregated Monte Carlo results (parametric in `name`) ---
aggr1 = deserialize(joinpath(output_dir, "aggr1_$(name).jls")) # all uncertainties
aggr2 = deserialize(joinpath(output_dir, "aggr2_$(name).jls")) # bed only
aggr3 = deserialize(joinpath(output_dir, "aggr3_$(name).jls")) # surface only
aggr4 = deserialize(joinpath(output_dir, "aggr4_$(name).jls")) # flotation only
aggr5 = deserialize(joinpath(output_dir, "aggr5_$(name).jls")) # no uncertainties (deterministic)

# --- Reconstruct rasters from arrays (use thickness grid) ---
thickness        = clean_raster(Raster(joinpath(datadir_WWFS_input, thickness_file)))
lake_depth_mean  = Raster(aggr1.lakes_depth_fs, dims(thickness))
area_stoch       = Raster(aggr1.areas,          dims(thickness))

# --- Analyze ---
println("  Analyzing stochastic lake depth (free surface) for $(name)...")
analysis = analyze_lakes(lake_depth_mean, thickness; min_depth=2.0)

df = DataFrame(
    run = name,
    smooth_surface_ice_fraction = NaN,  # should be taken form the file name?
    min_depth_m = analysis.min_depth,
    supragl_fill_fraction = NaN,        
    n_lakes = nrow(analysis.stats),
    total_volume_m3 = sum(analysis.stats.volume),
    mean_area_m2 = mean(analysis.stats.area_m2),
    mean_depth_m = mean(collect(lake_depth_mean)[analysis.labels .> 0]),
    max_depth_m = maximum(collect(lake_depth_mean)[analysis.labels .> 0]),
    largest_single_volume_m3 = analysis.LargestLake.volume
)

CSV.write(joinpath(output_dir, "WWFS_stoch_lake_summary_$(name).csv"), df)

# --- Plot mean lake depth + area ---
plot_lake_depth(
    lake_depth_mean,
    thickness,
    analysis,
    nothing, # phi
    joinpath(output_dir, "stochastic_lake_depth_$(name).png");
    min_depth = analysis.min_depth,
    show_all_lakes = true,
    area = area_stoch,
    area_threshold = 1e4
    #depressions_path = "/scratch-3/cogier/data/BonnePierre_input/depressions_BP_20240830.shp"
)

# --- Boxplot: contributions (all, surface, bed, flotation) ---
boxplot_lake_vol_stoch(
    aggr1; aggr2=aggr2, aggr3=aggr3, aggr4=aggr4,aggr5=aggr5,
    savepath = joinpath(output_dir, "lake_volume_comparison_$(name).png")
)

# --- Boxplot: flotation correlation-length sensitivity (optional, filenames fixed) ---
aggr_f10   = deserialize(joinpath(output_dir, "aggr_flot_L10.jls"))
aggr_f100  = deserialize(joinpath(output_dir, "aggr_flot_L100.jls"))
aggr_f1000 = deserialize(joinpath(output_dir, "aggr_flot_L1000.jls"))

boxplot_lake_vol_stoch(
    aggr_f10; aggr2=aggr_f100, aggr3=aggr_f1000,
    labels = ["L=10m", "L=100m", "L=1000m"],
    savepath = joinpath(output_dir, "boxplot_lake_vol_flotation_$(name).png")
)
