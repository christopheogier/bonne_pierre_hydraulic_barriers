using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
using ArchGDAL
using Rasters
using DataFrames
using CSV
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows



include("LakeAnalysis.jl")
include("plots_makie.jl")
include("functions.jl")
using .LakeAnalysis


datadir_input = "/scratch-3/cogier/data/BonnePierre_input"
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"

# Load surface and thickness for 2024 October
name = "2024_October"
surface = clean_raster(Raster(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif")))
thickness = clean_raster(Raster(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif")))
bedrock = clean_raster(Raster(joinpath(datadir_WWFS_input,"bedrock_resamp_1m.tif")))

#load uncertainties
bed_err_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif")))


################################ WWFS stochastic ########################################


# --- Define uncertainty models ---
#kernel = "gauss"
cov_fn = WWFS.GRF.gaussian_kernel #or WWFS.GRF.exponential_kernel
range_bed = 50.0 #m, see XDEM variograms outputs
range_surf = 100.0 # ARBITRARY FOR NOW
corr_length_f = 50 # ARBITRARY FOR NOW

# correlation lengths
corr_length_bed = range_bed / sqrt(3)  # m 
#chatgpt: For models where the variogram approaches the sill asymptotically, 
#the practical range is defined as the distance at which the variogram reaches 95% of the sill. 
#This practical range relates to the correlation length as follows:​
# Gaussian Model: Practical range ≈ sqrt(3) x ℓ​ = 1.73 x l
# Exponential Model: Practical range ≈ 3 x ℓ​
corr_length_surf = 10.0 / sqrt(3)      # placeholder for DEM error corr. length


# Input fields (already loaded), but also convert in float for WWFS
surfdem = surface
beddem = bedrock
rmask     = thickness .> 0
floatfrac = 1 .* ones(size(surfdem))
source    = ones(size(surfdem)) # what is "source" ?

# Uncertainties
surfdem_uc   = Uncertainty(absuc=0.5, reluc=0.0, correlation_length=corr_length_surf, covariance_fn=cov_fn )  # e.g. DEM smoothing
beddem_uc    = Uncertainty(absuc=bed_err_std, reluc=0.0, correlation_length=corr_length_bed, covariance_fn=cov_fn)
floatfrac_uc = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=100)  # example value
source_uc    = Uncertainty()  

# Extract raster grid
x, y = dims(surface)
println("step(x): ", step(x))
println("step(y): ", step(y)) # -1 !

# Sink definitions
sink_areas = (
    outlet = [CartesianIndices((1:10, 1:length(y)))[:],
              CartesianIndices((1:10, 1:(length(y)÷2)))[:]]
)

# Stochastic model
model, get_sample, aggregate = WWFS.make_fns(step(x), 
                                             surfdem, surfdem_uc,
                                             beddem, beddem_uc,
                                             floatfrac, floatfrac_uc,
                                             source, source_uc,
                                             sink_areas,
                                             rmask)

# Single realization
input, output = model(get_sample()...);

# Monte Carlo sampling
aggr = map_mc(model, get_sample, aggregate, 1) #20


# RESULTS 
# Convert lake depth array to Raster
lake_depth_mean = Raster(aggr.lakes_depth_fs, dims(surface))

# Fake a LakeAnalysisResult to satisfy plotting interface
dummy_result = LakeAnalysisResult([], [], nothing, DataFrame(:volume => [0.0]))  # empty placeholder
# PUt real lake analyse using `LakeAnalysis` here 

# Plot lake depth with contours and volume text
plot_lake_depth(
    lake_depth_mean,
    thickness,
    dummy_result,
    nothing,  # phi
    joinpath(output_dir, "stochastic_lake_depth.png");
    min_depth=2.0,
    show_all_lakes=false
)
