### plot_summaries.jl
using Pkg
Pkg.activate("/scratch-1/cogier/hydraulic_barriers/")
using CSV, DataFrames
include("plot.jl")  # or "plot_summary_overview.jl" if separated

output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

df = CSV.read("/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis/WWFS_lake_summary.csv", DataFrame)

plot_selected_scenarios(df, output_dir)
