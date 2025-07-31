using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
using ArchGDAL
using Rasters
using DataFrames
using CSV
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF = WhereTheWaterFlows
using Serialization



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
surface_smooth_avg = clean_raster(Raster(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_avg_01smooth.tif"))) # average of resampled and smoothed DEM
thickness = clean_raster(Raster(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif")))
bedrock = clean_raster(Raster(joinpath(datadir_WWFS_input,"bedrock_resamp_1m.tif")))

#load uncertainties
bed_err_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif")))
#surface unc
surface_2024_oct_smooth_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "surface_2024_oct_err_std_smooth01.tif")))


################################ WWFS stochastic ########################################


# --- Define uncertainty models ---
#kernel = "gauss"
cov_fn = WWFS.GRF.gaussian_kernel #or WWFS.GRF.exponential_kernel
range_bed = 2900 #m, see XDEM variograms outputs
range_surf = 2900 #m it seems correlated all over the dem area !
corr_length_f = 500 # ARBITRARY FOR NOW, otherwise mae a sensitivity analysis

# the longer the correlation length the smaller the spread in the stochastic runs

# correlation lengths
corr_length_bed = range_bed / sqrt(3)  # m 
#chatgpt: For models where the variogram approaches the sill asymptotically, 
#the practical range is defined as the distance at which the variogram reaches 95% of the sill. 
#This practical range relates to the correlation length as follows:​
# Gaussian Model: Practical range ≈ sqrt(3) x ℓ​ = 1.73 x l
# Exponential Model: Practical range ≈ 3 x ℓ​
corr_length_surf = range_surf / sqrt(3)      # placeholder for DEM error corr. length


# Input fields (already loaded), but also convert in float for WWFS
surfdem = surface_smooth_avg   
beddem = bedrock
rmask     = thickness .> 0
floatfrac = 1 .* ones(size(surfdem))
source    = ones(size(surfdem)) # what is "source" ?

# Uncertainties
surfdem_uc   = Uncertainty(absuc=surface_2024_oct_smooth_std, reluc=0.0, correlation_length=corr_length_surf, covariance_fn=cov_fn )  
beddem_uc    = Uncertainty(absuc=bed_err_std, reluc=0.0, correlation_length=corr_length_bed, covariance_fn=cov_fn)
floatfrac_uc = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=corr_length_f,covariance_fn=cov_fn) 
#  f = 0.6 to 1.11 in Chu et aL 2016 (greenland)
# f = 0.8 to 1.1 in Bowling et al 2015 (greenland)
source_uc    = Uncertainty()  

# Loop over 4 uncertainty cases
for (i, (surf_uc, bed_uc, float_uc)) in enumerate([
    # aggr1: all uncertainties
    (surfdem_uc, beddem_uc, floatfrac_uc),
    # aggr2: only bedrock uncertain
    (Uncertainty(absuc=0.0, reluc=0.0), beddem_uc, Uncertainty(absuc=0.0, reluc=0.0)),
    # aggr3: only surface uncertain
    (surfdem_uc, Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0)),
    # aggr4: only flotation uncertain
    (Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0), floatfrac_uc)
])

    # Extract raster grid
    x, y = dims(surface)
    #println("step(x): ", step(x))
    #println("step(y): ", step(y)) # -1 !

    # Sink definitions
    sink_areas = (
    outlet = [CartesianIndices((1:10, 1:length(y)))[:],
              CartesianIndices((1:10, 1:(length(y)÷2)))[:]])
              
    println("🔄 Running WWFS stochastic for aggr$i...")
    model, get_sample, aggregate = WWFS.make_fns(step(x),
                                                 surfdem, surf_uc,
                                                 beddem, bed_uc,
                                                 floatfrac, float_uc,
                                                 source, source_uc,
                                                 sink_areas,
                                                 rmask)

    aggr = map_mc(model, get_sample, aggregate, 20)
    #(:areas, :areas_extra, :melt_freeze, :lakes_depth, :lakes_mask, :lakes_depth_fs, :lakes_mask_fs, :sc_locs, :kappas, :catchments, :catchment_fluxes, :n_samples)

    serialize(joinpath(output_dir, "aggr$(i)_2024_October.jls"), aggr)
    println("✅ Saved aggr$i to disk.")

    # save TIF results
    #write(joinpath(output_dir, "lake_depth_stoch_mean.tif"), Raster(aggr.lakes_depth_fs, dims(surface)), force=true)
    #write(joinpath(output_dir, "area_mean.tif"), Raster(aggr.areas, dims(surface))  , force=true)
end



