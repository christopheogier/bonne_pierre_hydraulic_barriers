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


# --- Choose which run to execute ---
run_name = "2024_June"       # or "2024_October"
# Name used in all outputs below

# Build file paths per run (matching what you wrote out in prepare.jl)
paths = if run_name == "2024_October"
    Dict(
        :surface_raw   => joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"),
        :surface_smooth=> joinpath(datadir_WWFS_input, "surface_2024_oct_smooth_01.tif"),
        :surface_err   => joinpath(datadir_WWFS_input, "surface_2024_oct_err_avg_01smooth.tif"),
        :thickness     => joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif"),
        :bedrock       => joinpath(datadir_WWFS_input, "bedrock_resamp_1m.tif")
    )
elseif run_name == "2024_June"
    Dict(
        :surface_raw   => joinpath(datadir_WWFS_input, "surface_2024_june.tif"),
        :surface_smooth=> joinpath(datadir_WWFS_input, "surface_2024_june_smooth_01.tif"),
        :surface_smooth_filled => "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis/2024_June_smoothed_surface_with_lakefilled.tif",
        :surface_err   => joinpath(datadir_WWFS_input, "surface_2024_june_err_avg_01smooth.tif"),
        :thickness     => joinpath(datadir_WWFS_input, "ice_thickness_2024_june.tif"),
        :bedrock       => joinpath(datadir_WWFS_input, "bedrock_resamp_1m.tif")
    )
end

surface         = clean_raster(Raster(paths[:surface_raw]))
surface_smooth  = clean_raster(Raster(paths[:surface_smooth]))
surface_err     = clean_raster(Raster(paths[:surface_err]))
thickness       = clean_raster(Raster(paths[:thickness]))
beddem          = clean_raster(Raster(paths[:bedrock]))


#load uncertainties
bed_err_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif")))
# import lus and minus sigma if we can force WWF within two assymetric bound?



################################ WWFS stochastic ########################################


# --- Define uncertainty models ---
#kernel = "gauss"
cov_fn = WWFS.GRF.gaussian_kernel #or WWFS.GRF.exponential_kernel
range_bed = 200 #m, see XDEM variograms outputs
range_surf = 10 #m #variogram indicate glacier-size length, i expect it to be equal to the smoothing length scale
corr_length_f = 100 #[10,100,1000] # ARBITRARY FOR NOW, otherwise mae a sensitivity analysis



#A longer spatial correlation length means that the Gaussian Random Field (GRF) has more smoothly varying, spatially coherent patterns. This causes neighboring pixels to vary together — leading to:

#So, in Monte Carlo simulations: The pixel-wise variability decreases, and The realizations look smoother, with fewer high-frequency perturbations

# correlation lengths
corr_length_bed = 200 #range_bed / sqrt(3)  # m 
#chatgpt: For models where the variogram approaches the sill asymptotically, 
#the practical range is defined as the distance at which the variogram reaches 95% of the sill. 
#This practical range relates to the correlation length as follows:​
# Gaussian Model: Practical range ≈ sqrt(3) x ℓ​ = 1.73 x l
# Exponential Model: Practical range ≈ 3 x ℓ​
# Spherical Model> range ≈ 0.66 x l
corr_length_surf = 10 #range_surf / sqrt(3)      # placeholder for DEM error corr. length



# Input fields (already loaded), but also convert in float for WWFS
surfdem = surface_smooth
rmask     = thickness .> 0
floatfrac = 1 .* ones(size(surfdem))
source    = ones(size(surfdem)) # what is "source" ?

# Uncertainties
surfdem_uc   = Uncertainty(absuc=surface_err, reluc=0.0, correlation_length=corr_length_surf, covariance_fn=cov_fn )  
beddem_uc    = Uncertainty(absuc=bed_err_std, reluc=0.0, correlation_length=corr_length_bed, covariance_fn=cov_fn)
floatfrac_uc = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=corr_length_f,covariance_fn=cov_fn) 
#  f = 0.6 to 1.11 in Chu et aL 2016 (greenland)
# f = 0.8 to 1.1 in Bowling et al 2015 (greenland)
source_uc    = Uncertainty()  

# Loop over 5 uncertainty cases
for (i, (surf_uc, bed_uc, float_uc)) in enumerate([
    # aggr1: all uncertainties
    (surfdem_uc, beddem_uc, floatfrac_uc),
    # aggr2: only bedrock uncertain
    (Uncertainty(absuc=0.0, reluc=0.0), beddem_uc, Uncertainty(absuc=0.0, reluc=0.0)),
    # aggr3: only surface uncertain
    (surfdem_uc, Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0)),
    # aggr4: only flotation uncertain
    (Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0), floatfrac_uc),
    # aggr5: no uncertainties (deterministic)
    (Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0), Uncertainty(absuc=0.0, reluc=0.0))   
])

    # Extract raster grid
    x, y = dims(surface)
    #println("step(x): ", step(x))
    #println("step(y): ", step(y)) # -1 !

    # Sink definitions
    sink_areas = (
    outlet = [CartesianIndices((1:10, 1:length(y)))[:],
              CartesianIndices((1:10, 1:(length(y)÷2)))[:]])
              
    println("🔄 Running WWFS stochastic for aggr$i on $run_name...")
    model, get_sample, aggregate = WWFS.make_fns(step(x),
                                                 surfdem, surf_uc,
                                                 beddem, bed_uc,
                                                 floatfrac, float_uc,
                                                 source, source_uc,
                                                 sink_areas,
                                                 rmask)

    aggr = map_mc(model, get_sample, aggregate, 20)
    #(:areas, :areas_extra, :melt_freeze, :lakes_depth, :lakes_mask, :lakes_depth_fs, :lakes_mask_fs, :sc_locs, :kappas, :catchments, :catchment_fluxes, :n_samples)

    serialize(joinpath(output_dir, "aggr$(i)_$(run_name).jls"), aggr)
    println("✅ Saved aggr$i to disk for $run_name.")

    # save TIF results
    #write(joinpath(output_dir, "lake_depth_stoch_mean.tif"), Raster(aggr.lakes_depth_fs, dims(surface)), force=true)
    #write(joinpath(output_dir, "area_mean.tif"), Raster(aggr.areas, dims(surface))  , force=true)
end



