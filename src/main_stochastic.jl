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

surface_raw         = clean_raster(Raster(paths[:surface_raw]))
surface_smooth  = clean_raster(Raster(paths[:surface_smooth])) # change to :surface_smooth_filled if needed
surface_err     = clean_raster(Raster(paths[:surface_err]))
thickness       = clean_raster(Raster(paths[:thickness]))
beddem          = clean_raster(Raster(paths[:bedrock]))


#load uncertainties
bed_err_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif")))
# import lus and minus sigma if we can force WWF within two assymetric bound?

# --- Transect A→B for profile spaghetti (same as in plot_profiles) ---
surface_raw = clean_raster(Raster(paths[:surface_raw]))

A = (962632.47, 6431521.78)  # upstream
B = (962420.11, 6431549.94)  # downstream
dx_profile = 2.0             # sampling step (m)

x1, y1 = B
x2, y2 = A
L = hypot(x2 - x1, y2 - y1)
n_profile = max(1, floor(Int, L/dx_profile)) + 1
ts_profile = range(0.0, 1.0; length = n_profile)

pts_profile  = [(x1 + t*(x2 - x1), y1 + t*(y2 - y1)) for t in ts_profile]
dist_profile = collect(range(0.0, L; length = n_profile))

# Static bedrock & surface profiles (same for all cases and runs)
z_bed_profile  = _profile_vals(beddem,      pts_profile)
z_surf_profile = _profile_vals(surface_raw, pts_profile)


################################ WWFS stochastic ########################################


# --- Define uncertainty models ---
N = 1000 # number of realization
#kernel = "gauss"
cov_fn = WWFS.GRF.gaussian_kernel #or WWFS.GRF.exponential_kernel
range_bed = 247 #m, see XDEM variograms outputs
#range_surf = 10 #m #variogram indicate glacier-size length, i expect it to be equal to the smoothing length scale
corr_length_f = 100 #[10,100,1000] # ARBITRARY FOR NOW, otherwise mae a sensitivity analysis



#A longer spatial correlation length means that the Gaussian Random Field (GRF) has more smoothly varying, spatially coherent patterns. This causes neighboring pixels to vary together — leading to:

#So, in Monte Carlo simulations: The pixel-wise variability decreases, and The realizations look smoother, with fewer high-frequency perturbations

# correlation lengths
corr_length_bed = range_bed / sqrt(2)  # m = 175m
#chatgpt: For models where the variogram approaches the sill asymptotically, 
#the practical range is defined as the distance at which the variogram reaches 95% of the sill. 
#This practical range relates to the correlation length as follows:​
# Gaussian Model: Practical range ≈ sqrt(3) x ℓ​ = 1.73 x l
# Exponential Model: Practical range ≈ 3 x ℓ​
# Spherical Model> range ≈ 0.66 x l
corr_length_surf = 5 #range_surf / sqrt(3)      # placeholder for DEM error corr. length



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

# --- Helper for quickly defining flotation-uncertainty with a given corr. length ---
float_uc(L) = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=L, covariance_fn=cov_fn)
zero_uc()   = Uncertainty(absuc=0.0, reluc=0.0)

# Keep your existing definitions:
#   surfdem_uc   = Uncertainty(absuc=surface_err, reluc=0.0, correlation_length=corr_length_surf, covariance_fn=cov_fn)
#   beddem_uc    = Uncertainty(absuc=bed_err_std,  reluc=0.0, correlation_length=corr_length_bed,  covariance_fn=cov_fn)
#   floatfrac_uc = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=corr_length_f, covariance_fn=cov_fn)  # L = 100 m (your aggr4)
#   source_uc    = Uncertainty()

# --- Define all cases (indexed order controls aggr#) ---
cases = [
    # aggr1: all uncertainties
    (surfdem_uc, beddem_uc, floatfrac_uc),      # floatfrac_uc here is your L=100 m
    # aggr2: only bedrock uncertain
    (zero_uc(),  beddem_uc,  zero_uc()),
    # aggr3: only surface uncertain
    (surfdem_uc, zero_uc(),  zero_uc()),
    # aggr4: only flotation uncertain (L = 100 m, as defined above)
    (zero_uc(),  zero_uc(),  floatfrac_uc),
    # aggr5: no uncertainties (deterministic)
    (zero_uc(),  zero_uc(),  zero_uc()),
    # aggr6: only flotation uncertain (L = 10 m)
    (zero_uc(),  zero_uc(),  float_uc(10.0)),
    # aggr7: only flotation uncertain (L = 50 m)
    (zero_uc(),  zero_uc(),  float_uc(50.0)),
    # aggr8: only flotation uncertain (L = 1000 m)
    (zero_uc(),  zero_uc(),  float_uc(1000.0))
]

# --- Run all cases and save as aggr1..aggr8 ---
for (i, (surf_uc, bed_uc, float_uc_i)) in enumerate(cases)
#i=8
#(surf_uc, bed_uc, float_uc_i)=cases[8]  # for testing a single case

    # Extract raster grid
    xdim, ydim = dims(surfdem)
    dx = step(xdim)

    # --- SINKS : Define special catchment point  ---
    x0,y0 = 962468.722,6431539.971
    i0, j0 = coord_to_index(surface_smooth, x0, y0)
    println("Pixel index: ", (i0, j0))

    sink_areas = (point = [CartesianIndex(i0, j0)])  # single-pixel sink
    #sink_areas = (outlet = [CartesianIndices((1:10, 1:length(ydim)))[:],CartesianIndices((1:10, 1:(length(ydim)÷2)))[:]])
    


    println("🔄 Running WWFS stochastic for aggr$(i) on $run_name...")
    model, get_sample, aggregate = WWFS.make_fns(
        dx,
        surfdem,  surf_uc,
        beddem,   bed_uc,
        floatfrac, float_uc_i,
        source,   source_uc,
        sink_areas,
        rmask
    )
    
    # Largest-lake volume per realization
    largest_vols = Float64[]
    #  Per-realization transect profiles (for all cases)  # this is to check the mean versus stochastic ensembles 
    phi_profiles  = Vector{Vector{Float32}}()
    lake_profiles = Vector{Vector{Float32}}()

    for _ in 1:N
        s = get_sample()
        _, output = model(s...)
        lakes_free_surf = output[3][2]
        phi             = output[2][4]   # (sc_locs, kappas, diro, phi)

        # 1) largest lake volume
        analysis = analyze_lakes(lakes_free_surf, thickness; min_depth=2.0)
        push!(largest_vols, analysis.LargestLake.volume)

        # 2) profiles along A→B
        phi_r  = Raster(phi,             dims(thickness))
        lake_r = Raster(lakes_free_surf, dims(thickness))

        push!(phi_profiles,  Float32.(_profile_vals(phi_r,  pts_profile)))
        push!(lake_profiles, Float32.(_profile_vals(lake_r, pts_profile)))
    end


    # Aggregate maps/statistics over N runs
    aggr = map_mc(model, get_sample, aggregate, N)

    # Attach largest-lake ensemble + transect spaghetti + static profiles
    aggr = merge(aggr, (
        largest_lake_fs_vol = Float32.(largest_vols),
        phi_profiles        = phi_profiles,              # Vector{Vector{Float32}}
        lake_profiles       = lake_profiles,             # Vector{Vector{Float32}}
        dist_profile        = Float32.(dist_profile),    # 1D distance along A→B
        z_bed_profile       = Float32.(z_bed_profile),
        z_surf_profile      = Float32.(z_surf_profile),
    ))

    # Save with index-matched file name
    outfile = joinpath(output_dir, "aggr$(i)_n$(N)_$(run_name).jls")
    serialize(outfile, aggr)
    println("✅ Saved $(outfile)")

   

end



