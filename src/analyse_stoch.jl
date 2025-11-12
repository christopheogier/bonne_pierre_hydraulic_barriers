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
aggr4 = deserialize(joinpath(output_dir, "aggr4_$(name).jls")) # flotation only l = 100m
aggr5 = deserialize(joinpath(output_dir, "aggr5_$(name).jls")) # no uncertainties (deterministic)
aggr6 = deserialize(joinpath(output_dir, "aggr6_$(name).jls")) # flotation only, l = 10m
aggr7 = deserialize(joinpath(output_dir, "aggr7_$(name).jls")) # flotation only, l = 50m
aggr8 = deserialize(joinpath(output_dir, "aggr8_$(name).jls")) # flotation only, l = 1000m


# --- Reconstruct rasters from arrays (use thickness grid) ---
thickness        = clean_raster(Raster(joinpath(datadir_WWFS_input, thickness_file)))
lake_depth_mean  = Raster(aggr1.lakes_depth_fs, dims(thickness))
area_stoch       = Raster(aggr1.areas,          dims(thickness))

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


# --- Boxplot: total lake volume and largest-lake volume across aggr1 to aggr5 ---

labels = ["all unc.", "bed. unc.", "surf. unc.", "f. unc.", "No unc"]

lake_vols = [
    aggr1.lake_fs_vol,
    aggr2.lake_fs_vol,
    aggr3.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr5.lake_fs_vol,
]

largest_lake_vols = [
    aggr1.largest_lake_fs_vol,
    aggr2.largest_lake_fs_vol,
    aggr3.largest_lake_fs_vol,
    aggr4.largest_lake_fs_vol,
    aggr5.largest_lake_fs_vol,
]

fig = Figure(size = (700, 420))
ax  = Axis(fig[1, 1];
    ylabel = "Water-pocket volume (m³)",
    xticks = (1:length(labels), labels)
)

n       = length(labels)
centers = 1:n
offset  = 0.12
w       = 0.22

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
Legend(fig[2, 1], ax; framevisible=false)

save(joinpath(output_dir, "boxplot_total_vs_largest_$(name).png"), fig)
println("✅ Saved: ", joinpath(output_dir, "boxplot_total_vs_largest_$(name).png"))


# --- Boxplot: flotation uncertainty correlation lengths ---

Ls_labels = ["L=10 m", "L=50 m", "L=100 m", "L=1000 m"]

# Use the total-volume vectors (NOT the whole aggr structs).
# aggr6 → L=10 m, aggr7 → L=50 m, aggr4 → L=100 m, aggr8 → L=1000 m
Ls_vols = [
    aggr6.lake_fs_vol,
    aggr7.lake_fs_vol,
    aggr4.lake_fs_vol,
    aggr8.lake_fs_vol,
]

fig2 = Figure(size = (520, 420))
ax2  = Axis(fig2[1, 1];
    title = "Flotation uncertainty — correlation length",
    ylabel = "Total water-pocket volume (m³)",
    xticks = (1:4, Ls_labels)
)

# Tight spacing: contiguous positions, relatively wide boxes, minimal x padding
for (i, v) in enumerate(Ls_vols)
    xi = fill(i, length(v))
    boxplot!(ax2, xi, v; color=:steelblue, width=0.65)
    scatter!(ax2, [i], [mean(v)]; color=:black, marker=:cross, markersize=9)
end
xlims!(ax2, 0.5, 4.5)

save(joinpath(output_dir, "boxplot_flotation_totals_Ls_$(name).png"), fig2)
println("✅ Saved: ", joinpath(output_dir, "boxplot_flotation_totals_Ls_$(name).png"))
