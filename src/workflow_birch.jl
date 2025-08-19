# WWFS Birch prelim results

using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")
# Pkg.instantiate()  # Uncomment if you want to install packages
using ArchGDAL
using Dates
using Printf
using Rasters
using Statistics
include("plots_makie.jl")
include("functions.jl")
using .UncUtils
using WhereTheWaterFlowsSubglacially
const WWFS = WhereTheWaterFlowsSubglacially
using Serialization


datadir_in = "/scratch-3/cogier/data/Birch_input"
datadir_out = "/scratch-3/cogier/data/Birch_output"

### BED from grab et al 2021
bed = load_bedrock(joinpath(datadir_in, "B32-06_GlacierBed.tif"))


### SURFACE
#20230823
surf_23 = load_surface(joinpath(datadir_in, "B32-06_birch_20230823_dsm_is-stand_2m_lv95_ln02.tif"),bed) #2m res
#20250523
surf_25 =  load_surface(joinpath(datadir_in, "dronecam-2025-05-23T10_28_00-dem_200mmPP.tif"),bed) # 20 cm res
#write(joinpath(datadir_out, "Birch_surface_2025_cr.tif"), surf_25)

# resample surf 25 to surf 23 resolution (2m)
surf_25_resamp = clean_raster(resample(surf_25; to=surf_23, method=:bilinear))

### BEDROCK RESAMP and uncertainties

# resamp bedrock to surf 23 resolution
bed_resamp = clean_raster(resample(bed; to=surf_23, method=:bilinear))  # downscaling bedrock to 2m resolution
# need to resamp bed over surf_25 too to match surf_25 outline (see differnce in plots)
### Bedrock uncertainties

bedrock_err_plus = clean_raster(
    abs.(resample(Raster(joinpath(datadir_in, "B32-06_IceThicknessUncertaintyMinus.tif")); to=bed_resamp,method=:bilinear))) #minus ice lead to higher bedrock elevation
bed_plus_resamp = bed_resamp .+ bedrock_err_plus # to avoid the union missing
bedrock_err_minus = clean_raster(
    -1 * resample(Raster(joinpath(datadir_in, "B32-06_IceThicknessUncertaintyPlus.tif")); to=bed_resamp,method=:bilinear))
bed_minus_resamp = bed_resamp .+ bedrock_err_minus

# plot
plot_bedrock(bed_resamp; savepath=joinpath(datadir_out, "Birch_bedrock_elevation_resamp.png"))
plot_uncertainty_bed(bedrock_err_plus, bedrock_err_minus,
    "Bedrock: Uncertainty", "Uncertainty +", "Uncertainty -",
    joinpath(datadir_out, "bedrock_uncertainty.png"))

### THICKNESS

thickness = load_ice_thickness(joinpath(datadir_in, "B32-06_IceThickness.tif"))
# thickness 2023 and 2025
method_inter = :bilinear # Bilinear, near (weird) or cubic
thickness_25 = compute_ice_thickness(surf_25_resamp, bed_resamp, method_inter)
thickness_23 = compute_ice_thickness(surf_23, bed_resamp, method_inter)
#plot 
plot_ice_thickness(thickness;savepath=joinpath(datadir_out, "Birch_thickness_grab21"))
plot_ice_thickness(thickness_25;savepath=joinpath(datadir_out, "Birch_thickness_2025"))
plot_ice_thickness(thickness_23;savepath=joinpath(datadir_out, "Birch_thickness_2023"))
plot_ice_thickness(thickness_25 .- thickness_23; savepath=joinpath(datadir_out, "Birch_thickness_diff_2025-2023"))

### write rasters for later analysis ###
write(joinpath(datadir_out, "Birch_surface_2023_cr.tif"), surf_23; force=true)
write(joinpath(datadir_out, "Birch_surface_2025_cr.tif"), surf_25; force=true)
write(joinpath(datadir_out, "Birch_surface_2025_resamp.tif"), surf_25_resamp; force=true)
write(joinpath(datadir_out, "Birch_bedrock_resamp.tif"), bed_resamp; force=true)
write(joinpath(datadir_out, "Birch_thickness_2023.tif"), thickness_23; force=true)
write(joinpath(datadir_out, "Birch_thickness_2025.tif"), thickness_25; force=true)


#### WWFS #####

runs = [
    (
        name = "2023",
        surface = surf_23,
        thickness = thickness_23,
        bed = bed_resamp
    ),
    (
        name = "2025",
        surface = surf_25_resamp,
        thickness = thickness_25,
        bed = bed_resamp
    )
]

# Uncertainties
cov_fn = WWFS.GRF.gaussian_kernel # or WWFS.GRF.exponential_kernel
surfdem_uc   = Uncertainty(absuc=0.0, reluc=0.0)  
beddem_uc    = Uncertainty(absuc=bedrock_err_plus, reluc=0.0, correlation_length=1000, covariance_fn=cov_fn)
floatfrac_uc = Uncertainty(absuc=0.0, reluc=0.1, correlation_length=100, covariance_fn=cov_fn) 
source_uc    = Uncertainty()  

### WWFS deterministic run on 2023 (for φ) ###
println("\n🔷 Deterministic WWFS on 2023 surface")

# Grid spacing
x, _ = dims(surf_23)
dx = step(x)

# Run subglacial solver (no stochasticity)
((areas, slen, dir, nout, nin, sinks, pits, c, bnds),
 (sc_locs, kappas, diro, phi),
 (lakes, lakes_free_surf),
 sinkout) = WWFS.waterflows_subglacial(
    surf_23, bed_resamp, dx;
    gamma = [0, WWFS.GAMMA][1],   
    bnd_as_sink = true,
    drain_pits  = true
)

# Save φ raster + plot
#phi_tif  = joinpath(datadir_out, "Birch_phi_2023.tif")
#write(phi_tif, phi; force=true)

plot_hydraulic_head_and_flux(
    phi,
    areas[1],           # upslope area from waterflows_subglacial
    thickness_23,
    joinpath(datadir_out, "Birch_phi_flux_2023.png");
    min_threshold = 1e4,
    max_threshold = 1e6
)


# Stochastic runs


for run in runs
    println("\n🔷 Processing run: $(run.name)")
    println("eltypes: surface=$(eltype(run.surface)), thickness=$(eltype(run.thickness)), bed=$(eltype(run.bed))")
    println("has missing: surface=$(any(ismissing, run.surface)), thickness=$(any(ismissing, run.thickness)), bed=$(any(ismissing, run.bed))")

    # Input fields
    surfdem   = run.surface
    bed_resamp = run.bed
    thickness = run.thickness

    rmask     = thickness .> 0
    floatfrac = 1 .* ones(size(surfdem))
    source    = ones(size(surfdem)) 

    # Extract raster grid
    x, y = dims(surfdem)

    # Sink definitions (example)
    sink_areas = (
        outlet = [CartesianIndices((1:10, 1:length(y)))[:],
                  CartesianIndices((1:10, 1:(length(y)÷2)))[:]]
    )

    # Loop over 2 uncertainty cases
    for (i, (surf_uc, bed_uc, float_uc)) in enumerate([
        # aggr1: all uncertainties
        (surfdem_uc, beddem_uc, floatfrac_uc),
        # aggr2: only bedrock
        (Uncertainty(absuc=0.0, reluc=0.0), beddem_uc, Uncertainty(absuc=0.0, reluc=0.0))
    ])

        println("🔄 Running WWFS stochastic for $(run.name), aggr$i...")

        model, get_sample, aggregate = WWFS.make_fns(step(x),
                                                     surfdem, surf_uc,
                                                     bed_resamp, bed_uc,
                                                     floatfrac, float_uc,
                                                     source, source_uc,
                                                     sink_areas,
                                                     rmask)

        aggr = map_mc(model, get_sample, aggregate, 20)

        # Save results
        run_id = "$(run.name)_aggr$(i)"
        serialize(joinpath(datadir_out, "aggr_" * run_id), aggr)
        println("✅ Saved $run_id to disk.")

        # Example: save raster results
        # write(joinpath(output_dir, run_id * "_lake_depth_mean.tif"),
        #       Raster(aggr.lakes_depth_fs, dims(surfdem)); force=true)

    end
end




