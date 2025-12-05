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
N = 1000  # number of realizations 
name = "2024_June"   # or "2024_October"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir         = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# Map run name -> thickness file
thickness_file = Dict(
    "2024_October" => "ice_thickness_2024_oct.tif",
    "2024_June"    => "ice_thickness_2024_june.tif",
)[name]

# --- Load aggregated Monte Carlo results (parametric in `name`) ---
aggr1 = deserialize(joinpath(output_dir, "aggr1_n$(N)_$(name).jls")) # all uncertainties
aggr2 = deserialize(joinpath(output_dir, "aggr2_n$(N)_$(name).jls")) # bed only
aggr3 = deserialize(joinpath(output_dir, "aggr3_n$(N)_$(name).jls")) # surface only
aggr4 = deserialize(joinpath(output_dir, "aggr4_n$(N)_$(name).jls")) # flotation only l = 100m
aggr5 = deserialize(joinpath(output_dir, "aggr5_n$(N)_$(name).jls")) # no uncertainties (deterministic)
aggr6 = deserialize(joinpath(output_dir, "aggr6_n$(N)_$(name).jls")) # flotation only, l = 10m
aggr7 = deserialize(joinpath(output_dir, "aggr7_n$(N)_$(name).jls")) # flotation only, l = 50m
aggr8 = deserialize(joinpath(output_dir, "aggr8_n$(N)_$(name).jls")) # flotation only, l = 1000m


# --- Reconstruct rasters from arrays (use thickness grid) ---
thickness        = clean_raster(Raster(joinpath(datadir_WWFS_input, thickness_file)))
lake_depth_mean  = Raster(aggr1.lakes_depth_fs, dims(thickness))
area_stoch       = Raster(aggr1.areas,          dims(thickness))

# write area_stoch raster
#write(joinpath(output_dir, "stochastic_upslope_area_$(name).tif"), area_stoch,force=true)
#println("✅ Saved: ", joinpath(output_dir, "stochastic_upslope_area_$(name).tif"))

# --- Analyze ---
println("  Analyzing stochastic lake depth (free surface) for $(name)...")
analysis = analyze_lakes(lake_depth_mean, thickness; min_depth=2.0)

df = DataFrame(
    run = name,
    smooth_surface_ice_fraction = NaN,  # should be taken fom the file name?
    min_depth_m = analysis.min_depth,
    supragl_fill_fraction = NaN,        
    n_lakes = nrow(analysis.stats),
    total_volume_m3 = mean(aggr1.lake_fs_vol), #and not sum(analysis.stats.volume)! 
    total_volume_std = std(aggr1.lake_fs_vol),   
    #sum(analysis.stats.volume) = volume computed from the mean depth map (lakes_depth_fs) after thresholding at 2 m and relabeling components.
    # which is different that mean(aggr1.lake_fs_vol) in principle, because the mean of lake volumes over all stochastic runs is not equal to the volume computed from the mean lake depth map.
    #total_volume_m3 = sum(analysis.stats.volume),
    mean_area_m2 = mean(analysis.stats.area_m2),
    mean_depth_m = mean(collect(lake_depth_mean)[analysis.labels .> 0]),
    max_depth_m = maximum(collect(lake_depth_mean)[analysis.labels .> 0]),
    mean_largest_lake_m3 = mean(aggr1.largest_lake_fs_vol), #mean of the 20 realizations
    mean_largest_lake_std = std(aggr1.largest_lake_fs_vol)
    #largest_single_volume_m3 = analysis.LargestLake.volume #is the largest lake from the mean depth field, not the mean of the largest lake across realizations.
)

CSV.write(joinpath(output_dir, "WWFS_stoch_lake_summary_$(name).csv"), df)

# --- Save summary table of mean and std for all stochastic runs (aggr1–aggr8) ---

# Define labels corresponding to your 8 runs
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

# Collect all aggr datasets in order
aggr_list = [aggr1, aggr2, aggr3, aggr4, aggr5, aggr6, aggr7, aggr8]

# Compute summary statistics (mean + std) for total volumes
means = [mean(a.lake_fs_vol) for a in aggr_list]
stds  = [std(a.lake_fs_vol)  for a in aggr_list]
means_largest = [mean(a.largest_lake_fs_vol) for a in aggr_list]
stds_largest  = [std(a.largest_lake_fs_vol)  for a in aggr_list]

# Combine into a DataFrame
df_summary = DataFrame(
    label = aggr_labels,
    mean_total_volume_m3 = means,
    std_total_volume_m3 = stds,
    mean_largest_lake_m3 = means_largest,
    std_largest_lake_m3 = stds_largest
)

# Write CSV
csv_path = joinpath(output_dir, "WWFS_stoch_volume_summary_$(name).csv")
CSV.write(csv_path, df_summary)
println("✅ Saved volume summary to: ", csv_path)



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
    area_threshold = 1e4,
    stochastic = true
    #depressions_path = "/scratch-3/cogier/data/BonnePierre_input/depressions_BP_20240830.shp"
)


# --- Boxplot: total lake volume and largest-lake volume across aggr1 to aggr5 ---

labels = ["none", "bedrock", "surface", "f", "all"]

lake_vols = [
    aggr5.lake_fs_vol,
    aggr2.lake_fs_vol,
    aggr3.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr1.lake_fs_vol,
]

largest_lake_vols = [
    aggr5.largest_lake_fs_vol,
    aggr2.largest_lake_fs_vol,
    aggr3.largest_lake_fs_vol,
    aggr4.largest_lake_fs_vol,
    aggr1.largest_lake_fs_vol,
]


# --- Plot only total lake volume
plot_lake_volume_boxplot(
    lake_vols, largest_lake_vols,
    labels, joinpath(output_dir, "boxplot_total_vs_largest_$(name).png");
    plot_largest = false,  # true to plot largest lake volume in addition
    logscale = false
)

# --- Boxplot: flotation uncertainty correlation lengths (linear y) ---
Ls_labels = ["L=10 m", "L=100 m", "L=1000 m"]

# Use the total-volume vectors (NOT the whole aggr structs).
# aggr6 → L=10 m, aggr7 → L=50 m, aggr4 → L=100 m, aggr8 → L=1000 m
Ls_vols = [
    aggr6.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr8.lake_fs_vol,
]

fig2 = Figure(size = (200, 420))
ax2  = Axis(fig2[1, 1];
    title  = "",
    ylabel = "",
    xticks  = (1:3, Ls_labels),
)

# Tight spacing: very wide boxes and tight x-limits
cols = [:lightsteelblue, :dodgerblue, :royalblue]  # keep L=100 as dodgerblue
for (i, v) in enumerate(Ls_vols)
    xi = fill(i, length(v))
    boxplot!(ax2, xi, v; color=cols[i], width=0.98, show_outliers=true)  # width ~1 so boxes touch
    scatter!(ax2, [i], [mean(v)]; color=:black, marker=:cross, markersize=8)
end
xlims!(ax2, 0.51, 4.49)

save(joinpath(output_dir, "boxplot_flotation_totals_Ls_$(name).png"), fig2)
println("✅ Saved: ", joinpath(output_dir, "boxplot_flotation_totals_Ls_$(name).png"))

# --- Spaghetti plot of hydraulic head φ along transect A→B ---

#plot_transect_spaghetti(
    #aggr2, aggr5;  # spagethi versus mean (no unc.)
    #case_label_spag = "bed. unc. (aggr2)",
    #case_label_mean = "no unc. (aggr5)",
    #plot_phi  = true,
    #plot_lake = true,
    #savepath  = joinpath(output_dir, "transect_spaghetti_aggr2_vs_aggr5_$(name).png")
#)


