# ================
# Analyse stochastic lake depth and plot results (Birch)
# ================
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using Rasters
using DataFrames, CSV, Serialization
using Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF  = WhereTheWaterFlows

# local includes
include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis

# --- helpers: safe reducers for possibly-empty selections ---
safe_mean(v) = isempty(v) ? NaN : mean(v)
safe_max(v)  = isempty(v) ? NaN : maximum(v)

# --- Paths (same convention as birch_workflow) ---
datadir_out = "/scratch-3/cogier/data/Birch_output"
analysis_outdir = datadir_out

# --- Reload already saved rasters ---
surf_23        = Raster(joinpath(datadir_out, "Birch_surface_2023_cr.tif"))
surf_25_resamp = Raster(joinpath(datadir_out, "Birch_surface_2025_resamp.tif"))
bed_resamp     = Raster(joinpath(datadir_out, "Birch_bedrock_resamp.tif"))
thickness_23   = Raster(joinpath(datadir_out, "Birch_thickness_2023.tif"))
thickness_25   = Raster(joinpath(datadir_out, "Birch_thickness_2025.tif"))

# --- Analysis params ---
analysis_min_depth  = 0.0      # meters
area_plot_threshold = 1e3      # m² for overlaying upslope area

# Map run -> thickness
thk_by_run = Dict(
    "2023" => thickness_23,
    "2025" => thickness_25,
)

summaries = DataFrame()

for run_name in ["2023", "2025"]
    println("\n📊 Analyzing stochastic lakes for run: $run_name")

    thk = thk_by_run[run_name]

    # ---- Load aggregated Monte Carlo results ----
    aggr1_path = joinpath(analysis_outdir, "aggr_$(run_name)_aggr1") # all uncertainties
    aggr2_path = joinpath(analysis_outdir, "aggr_$(run_name)_aggr2") # bed-only

    aggr1 = deserialize(aggr1_path)
    aggr2 = deserialize(aggr2_path)

    # ---- Reconstruct rasters from arrays ----
    lake_depth_mean = Raster(aggr1.lakes_depth_fs, dims(thk))
    area_stoch      = Raster(aggr1.areas,          dims(thk))


    # ---- Analyze ----
    println("  ▸ Analyzing stochastic lake depth (free surface)...")
    analysis = analyze_lakes(lake_depth_mean, thk; min_depth=analysis_min_depth)

    labels_pos = analysis.labels .> 0
    depth_vals = collect(lake_depth_mean)[labels_pos]  # may be empty

    n_lakes    = nrow(analysis.stats)
    total_vol  = n_lakes == 0 ? 0.0 : sum(analysis.stats.volume)
    mean_area  = n_lakes == 0 ? NaN  : mean(analysis.stats.area_m2)
    mean_depth = safe_mean(depth_vals)
    max_depth  = safe_max(depth_vals)

    df = DataFrame(
        run = run_name,
        smooth_surface_ice_fraction = NaN,
        min_depth_m = analysis.min_depth,
        supragl_fill_fraction = NaN,
        n_lakes = n_lakes,
        total_volume_m3 = total_vol,
        mean_area_m2 = mean_area,
        mean_depth_m = mean_depth,
        max_depth_m = max_depth,
        largest_single_volume_m3 = n_lakes == 0 ? 0.0 : analysis.LargestLake.volume
    )
    append!(summaries, df)

    # ---- Plot lake depth + upslope area ----
    out_png = joinpath(analysis_outdir, "stochastic_lake_depth_$(run_name).png")
    plot_lake_depth(
        lake_depth_mean,
        thk,
        analysis,
        nothing,
        out_png;
        min_depth = analysis.min_depth,
        show_all_lakes = true,
        area = area_stoch,
        area_threshold = area_plot_threshold
    )

    # ---- Boxplot: all uncertainties vs bed-only ----
    out_box = joinpath(analysis_outdir, "lake_volume_comparison_$(run_name).png")
    boxplot_lake_vol_stoch(aggr1;
        aggr2 = aggr2,
        labels = ["all uncertainties", "bedrock only"],
        savepath = out_box
    )

    println("  ✅ Saved plots for run: $run_name")
end

# ---- Save CSV summary ----
csv_path = joinpath(analysis_outdir, "WWFS_stoch_lake_summary_Birch.csv")
CSV.write(csv_path, summaries)
println("📝 Wrote summary: $csv_path")
