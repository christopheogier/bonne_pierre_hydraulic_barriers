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
N = 10  # number of realizations 
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

labels = ["no unc.", "bed. unc.", "surf. unc.", "f. unc.", "all unc"]

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

fig = Figure(size = (420, 420))
ax  = Axis(fig[1, 1];
    ylabel = "Water pocket volume (m³)",
    xticks = (1:length(labels), labels)
)

n       = length(labels)
centers = 1:n
offset  = 0.12
w       = 0.3

pos_tot = centers .- offset
pos_lrg = centers .+ offset

x_tot = vcat([fill(pos_tot[i], length(v)) for (i,v) in enumerate(lake_vols)]...)
y_tot = vcat(lake_vols...)
x_lrg = vcat([fill(pos_lrg[i], length(v)) for (i,v) in enumerate(largest_lake_vols)]...)
y_lrg = vcat(largest_lake_vols...)

boxplot!(ax, x_tot, y_tot; color=:dodgerblue, width=w)
boxplot!(ax, x_lrg, y_lrg; color=:orange,     width=w)

scatter!(ax, pos_tot, [mean(v) for v in lake_vols];
         color=:black, marker=:cross, markersize=9)
scatter!(ax, pos_lrg, [mean(v) for v in largest_lake_vols];
         color=:black, marker=:cross, markersize=9)

lines!(ax, [NaN], [NaN]; color=:dodgerblue, label="Total pockets volume")
lines!(ax, [NaN], [NaN]; color=:orange,     label="Largest pocket volume")
#Legend(fig[2, 1], ax; framevisible=false)

save(joinpath(output_dir, "boxplot_total_vs_largest_$(name).png"), fig)
println("✅ Saved: ", joinpath(output_dir, "boxplot_total_vs_largest_$(name).png"))



# --- Boxplot: flotation uncertainty correlation lengths (linear y) ---
Ls_labels = ["L=10 m", "L=50 m", "L=100 m", "L=1000 m"]

# Use the total-volume vectors (NOT the whole aggr structs).
# aggr6 → L=10 m, aggr7 → L=50 m, aggr4 → L=100 m, aggr8 → L=1000 m
Ls_vols = [
    aggr6.lake_fs_vol,
    aggr7.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr8.lake_fs_vol,
]

fig2 = Figure(size = (200, 420))
ax2  = Axis(fig2[1, 1];
    title  = "",
    ylabel = "",
    xticks  = (1:4, Ls_labels),
)

# Tight spacing: very wide boxes and tight x-limits
cols = [:lightsteelblue, :steelblue, :dodgerblue, :royalblue]  # keep L=100 as dodgerblue
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

# catchment probability at a given point 
catch_prob   = aggr1.catchments[:, :, 1]              # Float16
#catch_rast16 = Raster(catch_prob, dims(thickness))    # same grid as thickness
#catch_rast = clean_raster(catch_rast16)               # → Float32, NaN instead of missing
#write(joinpath(output_dir, "catchment_prob_point_$(name).tif"),catch_rast; force = true)
#println("✅ Saved",joinpath(output_dir, "catchment_prob_point_$(name).tif"))

# Quick visualization
fig = Figure()
ax  = Axis(fig[1,1], title="Catchment probability to GPR point")

x,y,_ = get_axes_and_matrix(thickness)

heatmap!(ax, x, y, catch_prob; colormap=:viridis, colorrange=(0,1))
x0,y0 = 962468.722,6431539.971
#i0, j0 = coord_to_index(surface_smooth, x0, y0)
scatter!(ax, [x0], [y0]; color=:red, markersize=12)

save(joinpath(output_dir, "catchment_prob_point_$(name).png"), fig)
println("✅ Saved: ", joinpath(output_dir, "catchment_prob_point_$(name).png"))
