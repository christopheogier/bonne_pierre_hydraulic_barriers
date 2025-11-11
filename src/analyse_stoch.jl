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

# --- Boxplot: total vs largest lake volumes across uncertainties ---

# ------- Left panel: total vs largest per uncertainty -------
labels = ["No unc", "all unc.", "bed. unc.", "surf. unc.", "f. unc."]

lake_vols = [
    aggr5.lake_fs_vol,
    aggr1.lake_fs_vol,
    aggr2.lake_fs_vol,
    aggr3.lake_fs_vol,
    aggr4.lake_fs_vol,
    
]

largest_lake_vols = [
    aggr5.largest_lake_fs_vol,
    aggr1.largest_lake_fs_vol,
    aggr2.largest_lake_fs_vol,
    aggr3.largest_lake_fs_vol,
    aggr4.largest_lake_fs_vol
]

# ------- Right panel: flotation-only, TOTAL volume with different L -------
flot_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis_flot"
aggr_f50  = deserialize(joinpath(flot_dir, "aggr_flot_L50_2024_June.jls"))
aggr_f100 = deserialize(joinpath(flot_dir, "aggr_flot_L100_2024_June.jls"))
aggr_f500 = deserialize(joinpath(flot_dir, "aggr_flot_L500_2024_June.jls"))

Ls_labels = ["L=50 m", "L=100 m", "L=500 m"]
Ls_vols   = [aggr_f50.lake_fs_vol,   # <-- TOTAL volume 
             aggr_f100.lake_fs_vol,
             aggr_f500.lake_fs_vol]

# ------- Figure -------
fig = Figure(size = (950, 420))

# Left subplot
axL = Axis(fig[1, 1];
    ylabel = "Water-pocket volume (m³)",
    xticks = (1:length(labels), labels)
)

n        = length(labels)
centers  = 1:n
offset   = 0.12
w        = 0.22

pos_tot  = centers .- offset
pos_lrg  = centers .+ offset

x_tot = vcat([fill(pos_tot[i], length(v)) for (i,v) in enumerate(lake_vols)]...)
y_tot = vcat(lake_vols...)
x_lrg = vcat([fill(pos_lrg[i], length(v)) for (i,v) in enumerate(largest_lake_vols)]...)
y_lrg = vcat(largest_lake_vols...)

boxplot!(axL, x_tot, y_tot; color=:dodgerblue, width=w)
boxplot!(axL, x_lrg, y_lrg; color=:orange,     width=w)

scatter!(axL, pos_tot, [mean(v) for v in lake_vols];
         color=:black, marker=:cross, markersize=9)
scatter!(axL, pos_lrg, [mean(v) for v in largest_lake_vols];
         color=:black, marker=:cross, markersize=9)

lines!(axL, [NaN], [NaN]; color=:dodgerblue, label="Total pockets volume")
lines!(axL, [NaN], [NaN]; color=:orange,     label="Largest pocket volume")
Legend(fig[2, 1], axL; framevisible=false)

# Right subplot (TOTAL volume, flotation correlation length)
axR = Axis(fig[1, 2];
    title = "Flotation uncertainty — correlation length",
    ylabel = "Total water-pocket volume (m³)",
    xticks = (1:3, Ls_labels)
)

cols_R = [:steelblue, :dodgerblue, :royalblue]  # keep L=100 as dodgerblue

for (i, v) in enumerate(Ls_vols)
    xi = fill(i, length(v))
    boxplot!(axR, xi, v; color=cols_R[i], width=0.35)
    scatter!(axR, [i], [mean(v)]; color=:black, marker=:cross, markersize=9)
end

lines!(axR, [NaN], [NaN]; color=:steelblue,  label="L = 50 m")
lines!(axR, [NaN], [NaN]; color=:dodgerblue, label="L = 100 m")
lines!(axR, [NaN], [NaN]; color=:royalblue,  label="L = 500 m")
Legend(fig[2, 2], axR; framevisible=false)

save(joinpath(output_dir, "lake_volume_total_vs_largest_plus_flotL_total_$(name).png"), fig)
println("✅ Saved figure: ",
        joinpath(output_dir, "lake_volume_total_vs_largest_plus_flotL_total_$(name).png"))


# --- Boxplot: contributions (all, surface, bed, flotation) ---
#boxplot_lake_vol_stoch(
    #aggr1; aggr2=aggr2, aggr3=aggr3, aggr4=aggr4,aggr5=aggr5,
    #savepath = joinpath(output_dir, "lake_volume_comparison_$(name).png")
#)

# --- Quick plot: boxplot of largest-lake volumes across realizations ---

#vols = aggr1.largest_lake_fs_vol
#fig = Figure(size=(400, 500))
#ax  = Axis(fig[1,1];
    #title = "Largest water pocket volume — $(name)",
    #ylabel = "Volume (m³)",
    #xticks = ([1], ["Largest lake"])   # <- vector, not scalar
#)

#boxplot!(ax, fill(1, length(vols)), vols; color=:lightblue)
#scatter!(ax, [1], [mean(vols)]; color=:red, marker=:cross, markersize=14, label="Mean")

#save(joinpath(output_dir, "largest_lake_volume_boxplot_$(name).png"), fig; px_per_unit=3)
#println("✅ Saved boxplot to ", joinpath(output_dir, "largest_lake_volume_boxplot_$(name).png"))

# --- Boxplot: flotation correlation-length sensitivity (optional, filenames fixed) ---
#aggr_f10   = deserialize(joinpath(output_dir, "aggr_flot_L10.jls"))
#aggr_f100  = deserialize(joinpath(output_dir, "aggr_flot_L100.jls"))
#aggr_f1000 = deserialize(joinpath(output_dir, "aggr_flot_L1000.jls"))

#boxplot_lake_vol_stoch(
    #aggr_f10; aggr2=aggr_f100, aggr3=aggr_f1000,
    #labels = ["L=10m", "L=100m", "L=1000m"],
    #savepath = joinpath(output_dir, "boxplot_lake_vol_flotation_$(name).png")
#)
