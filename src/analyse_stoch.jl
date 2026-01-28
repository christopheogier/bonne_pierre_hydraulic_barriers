### analyse stochastic lake depth and plot results
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
using Rasters, DataFrames, CSV, Statistics
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF  = WhereTheWaterFlows
using Serialization
using Base: summarysize

include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis


# ======================================================================================
# --- Params / Paths ---
# ======================================================================================
N = 1000                          # number of realizations
name = "2024_June"              # or "2024_October"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir         = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# --- which case do we use for MAPS/PLOTTING? (important) ---
# NOTE: CSV summary uses aggr1–aggr8 anyway.
case_for_maps = 1             # 1=all unc., 2=bed only,.. 4= Lf 100m only ..., 8=f(L=1000m)

# --- plot constraint: plot only big WPs to avoid killing plotting ---
V_thr_plot = 1000.0             # m^3

# --- lake analysis threshold (used to build connected components)
min_depth = 0.1               # m  (WARNING: plotting ALL components may kill)


# ======================================================================================
# --- Input rasters (grid definition)
# ======================================================================================
thickness_file = Dict(
    "2024_October" => "ice_thickness_2024_oct.tif",
    "2024_June"    => "ice_thickness_2024_june.tif",
)[name]

thickness = clean_raster(Raster(joinpath(datadir_WWFS_input, thickness_file)))


# ======================================================================================
# --- Load aggregated Monte Carlo results (parametric in `name`)
# ======================================================================================
aggr1 = deserialize(joinpath(output_dir, "aggr1_n$(N)_$(name).jls")) # all uncertainties
aggr2 = deserialize(joinpath(output_dir, "aggr2_n$(N)_$(name).jls")) # bed only
aggr3 = deserialize(joinpath(output_dir, "aggr3_n$(N)_$(name).jls")) # surface only
aggr4 = deserialize(joinpath(output_dir, "aggr4_n$(N)_$(name).jls")) # flotation only l = 100m
aggr5 = deserialize(joinpath(output_dir, "aggr5_n$(N)_$(name).jls")) # no uncertainties (deterministic)
aggr6 = deserialize(joinpath(output_dir, "aggr6_n$(N)_$(name).jls")) # flotation only, l = 10m
aggr7 = deserialize(joinpath(output_dir, "aggr7_n$(N)_$(name).jls")) # flotation only, l = 50m
aggr8 = deserialize(joinpath(output_dir, "aggr8_n$(N)_$(name).jls")) # flotation only, l = 1000m

aggr_labels = [
    "all unc.",
    "bed. unc.",
    "surf. unc.",
    "f. unc. (L=100 m)",
    "no unc.",
    "f. unc. (L=10 m)",
    "f. unc. (L=50 m)",
    "f. unc. (L=1000 m)",
]

aggr_list = [aggr1, aggr2, aggr3, aggr4, aggr5, aggr6, aggr7, aggr8]

# pick the case for maps/plotting
aggr_map = aggr_list[case_for_maps]
println("📌 Map/plot case: aggr$(case_for_maps) = $(aggr_labels[case_for_maps])")


# ======================================================================================
# --- Reconstruct rasters from arrays (use thickness grid)
# ======================================================================================
lake_depth_mean = Raster(aggr_map.lakes_depth_fs, dims(thickness))
area_stoch      = Raster(aggr_map.areas,         dims(thickness))


# ======================================================================================
# --- Analyze mean lake-depth map (for plotting mask only)
# ======================================================================================
println("  Analyzing lake depth map (free surface) for $(name)...")
analysis = analyze_lakes(lake_depth_mean, thickness; min_depth=min_depth)


# ======================================================================================
# --- CSV 1: one-line summary for aggr1 (all uncertainties) 
# ======================================================================================
df = DataFrame(
    run = name,
    N_runs = N,
    smooth_surface_ice_fraction = NaN,
    min_depth_m = analysis.min_depth,
    supragl_fill_fraction = NaN,
    n_lakes = nrow(analysis.stats),

    total_volume_m3 = mean(aggr1.lake_fs_vol),
    total_volume_std = std(aggr1.lake_fs_vol),

    mean_area_m2 = mean(analysis.stats.area_m2),
    mean_depth_m = mean(collect(lake_depth_mean)[analysis.labels .> 0]),
    max_depth_m = maximum(collect(lake_depth_mean)[analysis.labels .> 0]),

    mean_largest_lake_m3 = mean(aggr1.largest_lake_fs_vol),
    mean_largest_lake_std = std(aggr1.largest_lake_fs_vol)
)

csv_path = joinpath(output_dir, "WWFS_stoch_lake_summary_all_unc_$(name).csv")
CSV.write(csv_path, df)
println("✅ Saved volume all unc summary to: ", csv_path)


# ======================================================================================
# --- CSV 2: summary table across aggr1–aggr8 
# ======================================================================================
means         = [mean(a.lake_fs_vol)          for a in aggr_list]
stds          = [std(a.lake_fs_vol)           for a in aggr_list]
means_largest = [mean(a.largest_lake_fs_vol)  for a in aggr_list]
stds_largest  = [std(a.largest_lake_fs_vol)   for a in aggr_list]
means_gt1000  = [mean(a.lake_fs_vol_gt1000)   for a in aggr_list]
stds_gt1000   = [std(a.lake_fs_vol_gt1000)    for a in aggr_list]
means_n_gt1000 = [mean(a.n_lakes_gt1000)      for a in aggr_list]
stds_n_gt1000  = [std(a.n_lakes_gt1000)       for a in aggr_list]

df_summary = DataFrame(
    label = aggr_labels,
    N_runs = fill(N, length(aggr_labels)),
    min_depth_m = fill(min_depth, length(aggr_labels)),
    smooth_surface_ice_fraction = fill(0.1, length(aggr_labels)),  # manual
    supragl_fill_fraction = fill(0, length(aggr_labels)),          # manual

    mean_n_wp_gt1000 = means_n_gt1000,
    std_n_wp_gt1000  = stds_n_gt1000,

    mean_total_volume_m3 = means,
    std_total_volume_m3  = stds,
    mean_total_volume_gt1000_m3 = means_gt1000,
    std_total_volume_gt1000_m3  = stds_gt1000,

    mean_largest_lake_m3 = means_largest,
    std_largest_lake_m3  = stds_largest
)

csv_path = joinpath(output_dir, "WWFS_stoch_volume_summary_$(name).csv")
CSV.write(csv_path, df_summary)
println("✅ Saved volume summary to: ", csv_path)


# ======================================================================================
# --- Plot: only big WPs (individual volume > 1000 m³), otherwise plotting gets killed
# ======================================================================================
println("  Building plotting raster: keep only WPs with individual volume > $(V_thr_plot) m³ (min_depth=$(min_depth) m)")

keep_labels = analysis.stats.label[analysis.stats.volume .> V_thr_plot]

keep_mask = falses(size(analysis.labels))
for lab in keep_labels
    keep_mask .|= (analysis.labels .== lab)
end

lake_depth_mean_big = copy(lake_depth_mean)
lake_depth_mean_big[.!keep_mask] .= 0.0

analysis_big = analyze_lakes(lake_depth_mean_big, thickness; min_depth=min_depth)

plot_lake_depth(
    lake_depth_mean_big,
    thickness,
    analysis_big,
    nothing,
    joinpath(output_dir, "stochastic_lake_depth_WP1000_$(name)_aggr$(case_for_maps).png");
    min_depth = analysis_big.min_depth,
    show_all_lakes = true,
    area = area_stoch,
    area_threshold = 1e4,
    stochastic = true
)


# ======================================================================================
# --- Boxplots (kept as-is)
# ======================================================================================
labels = ["none", "bedrock", "surface", "f", "all"]

lake_vols = [
    aggr5.lake_fs_vol,
    aggr2.lake_fs_vol,
    aggr3.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr1.lake_fs_vol,
]

lake_vols_gt1000 = [
    aggr5.lake_fs_vol_gt1000,
    aggr2.lake_fs_vol_gt1000,
    aggr3.lake_fs_vol_gt1000,
    aggr4.lake_fs_vol_gt1000,
    aggr1.lake_fs_vol_gt1000,
]

largest_lake_vols = [
    aggr5.largest_lake_fs_vol,
    aggr2.largest_lake_fs_vol,
    aggr3.largest_lake_fs_vol,
    aggr4.largest_lake_fs_vol,
    aggr1.largest_lake_fs_vol,
]

plot_lake_volume_boxplot(
    lake_vols, lake_vols_gt1000,
    labels, joinpath(output_dir, "boxplot_total_vs_largest_$(name).png");
    plot_largest = false, # here true to plot second
    logscale = false
)

Ls_labels = ["L=50 m", "L=100 m", "L=1000 m"]
Ls_vols = [
    aggr7.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr8.lake_fs_vol,
]

fig2 = Figure(size = (200, 420))
ax2  = Axis(fig2[1, 1];
    title  = "",
    ylabel = "",
    xticks  = (1:3, Ls_labels),
)
# --- enforce identical y-scale than total boxplots ---
ylims!(ax2, 0.0, 8e5)

# ticks at 0,2,4,6 (×10⁵)
ytick_vals   = [0.0, 2e5, 4e5, 6e5, 8e5]

cols = [:lightsteelblue, :dodgerblue, :royalblue]
for (i, v) in enumerate(Ls_vols)
    xi = fill(i, length(v))
    boxplot!(ax2, xi, v; color=cols[i], width=0.98, show_outliers=true)
    scatter!(ax2, [i], [mean(v)]; color=:black, marker=:cross, markersize=8)
end
xlims!(ax2, 0.51, 4.49)

save(joinpath(output_dir, "boxplot_flotation_totals_Ls_gt1000_$(name).png"), fig2)
println("✅ Saved: ", joinpath(output_dir, "boxplot_flotation_totals_Ls_gt1000_$(name).png"))

# --- Spaghetti plot of hydraulic head φ along transect A→B ---
#plot_transect_spaghetti(...)
