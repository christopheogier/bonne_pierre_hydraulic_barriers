### plot
using Pkg
Pkg.activate("/scratch-1/cogier/hydraulic_barriers/")
using CSV, DataFrames
include("plots_makie.jl")  # or "plot_summary_overview.jl" if separated
using Rasters

output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

df = CSV.read("/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis/WWFS_lake_summary.csv", DataFrame)

#plot_selected_scenarios(df, output_dir)

