using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
import ArchGDAL
using Rasters
using DataFrames
using CSV
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows
using Statistics



include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis
using .UncUtils

datadir_input = "/scratch-3/cogier/data/BonnePierre_input"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# Load surface and thickness for 2024 October
name = "2024_October"
surface = Raster(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"))
thickness = Raster(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif"))
bedrock = Raster(joinpath(datadir_WWFS_input,"bedrock_resamp_1m.tif"))

# Load GPR - ±5 m bedrocks
bed_gpr_plus = Raster(joinpath(datadir_WWFS_input,"bedrock_gpr_plus_1m.tif"))
bed_gpr_minus = Raster(joinpath(datadir_WWFS_input,"bedrock_gpr_minus_1m.tif"))

#### Compute GPR uncertainty maps 
u_plus_gpr = bedrock .- bed_gpr_plus
u_minus_gpr = bedrock .- bed_gpr_minus

#### Compute Glate uncertainty maps

# Extract glacier outline points directly from the thickness raster
outline_polylines = extract_outline_from_thickness(thickness)

# Choose the longest outline (main glacier polygon)
outline_points = isempty(outline_polylines) ? [] : reduce(vcat, outline_polylines)

# Load and extend GPR dataset with outline points (h = 0)
gpr_df = load_and_extend_gpr(
    joinpath(datadir_input, "BonnePierre_ice_thickness_GPR_resampled_50cmtxt_reprojectLambert93.txt"),
    outline_points
)

# Optional: save for reproducibility
#CSV.write(joinpath(output_dir, "gpr_plusoutline_df.csv"), gpr_df)

# Compute distance raster to nearest GPR point
distance_raster = compute_distance_to_gpr(gpr_df, bedrock)

# Compute empirical uncertainty bounds based on mean ice thickness
h_mean_2024 = mean(thickness[thickness .> 0]) 
println("Mean ice thickness for 2024: ", h_mean_2024)
u_minus_interp, u_plus_interp = unc_propagate(distance_raster, h_mean_2024)

# Mask uncertainty outside glacier domain
glacier_mask = thickness .> 0
u_minus_interp[.!glacier_mask] .= NaN
u_plus_interp[.!glacier_mask] .= NaN

#### Total uncertainty bedrock
u_plus_bed = u_plus_gpr .+ u_plus_interp
u_minus_bed = u_minus_gpr .+ u_minus_interp

# plot all Uncertainties
plot_uncertainty_bed(u_plus_gpr, u_minus_gpr,
    "Bedrock: GPR Uncertainty", "GPR +5 m", "GPR -5 m",
    joinpath(output_dir, "bedrock_gpr_uncertainty.png"))

plot_uncertainty_bed(u_plus_interp, u_minus_interp,
    "Bedrock: GLATE interpolation uncertainty", "Interpolation +", "Interpolation -",
    joinpath(output_dir, "glate_interp_uncertainty.png"))

plot_uncertainty_bed(u_plus_bed, u_minus_bed,
    "Cumulative bedrock uncertainty (GPR + GLATE-Interpolation)", "+ sigma", "- sigma",
    joinpath(output_dir, "bedrock_all_uncertainty.png"))

### WWFS stochastic

# one need to define a symmetric uncertainty for the bedrock (plus minus sigma, the standard deviation)
bed_std = 0.5 .* (u_plus_bed .- u_minus_bed)
write(joinpath(datadir_WWFS_input, "bedrock_std_1m.tif"), bed_std, force=true)

#VAriogram and correlation_length




# define error in flotation factor + corr length.

# run WWFS stochastic 
# TODO: add lake analyse INSIDE WWFS so the output is a spread for each paramters (3: zs, zb and f)

#Plot this nicely: boxplot ? boxplot for 4 case: the total, and the 3 contribution. Which is the largest ?