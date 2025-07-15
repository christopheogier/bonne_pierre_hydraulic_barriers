using Pkg
Pkg.activate("/scratch-1/cogier/hydraulic_barriers/")
# Pkg.instantiate()  # Uncomment if you want to install packages

using Dates
using Printf
using Rasters
include("plot.jl")
include("functions.jl")

"""
prepare.jl

Prepare script to load data, plot ice thickness and bedrock elevation, surface DEMs,
resample data to common resolution, and optionally plot GPR points.
"""

datadir_in = "/scratch-3/cogier/data/BonnePierre_input"
datadir_out = "/scratch-3/cogier/data/BonnePierre_output"
plots_dir = joinpath(datadir_out, "plots")
datadir_WWFS_input = joinpath(datadir_in, "WWFS_input")

ENV["GKSwstype"] = "100"  # Use "file output" mode for GR when plotting


# Load bedrock and ice thickness (original resolution)
bedrock = load_bedrock(joinpath(datadir_in, "BED_10m_l93_20smooth.tif"))
ice_thickness = load_ice_thickness(joinpath(datadir_in, "Hall_10m_l93_GlaTE20.asc"))  # ice thickness from air eth and glate on 6.11.2024

# load bedrock plus and minus from GPR


# Load and clean surface DEM cropped to bedrock extent
surface_2021_cr = load_surface(joinpath(datadir_in, "Berarde_Res1.0_CompElevmean_merged.tif"), bedrock) #1m
surface_2024_june = load_surface(joinpath(datadir_in, "Lidar_juin2024_Bonne_Pierre_1m.tif"), bedrock) #1m
surface_2024_oct = load_surface(joinpath(datadir_in, "Lidar_oct2024_Bonne_Pierre_glacier.tif"), bedrock) # 50cm

# Resample October DEM to 2021 DEM resolution: from 50cm to 1m
surface_2024_oct_resamp = resample(surface_2024_oct; to=surface_2021_cr, method=:bilinear)

# Resample bedrock and ice thickness to surface 2021 resolution
bed_resamp = resample(bedrock; to=surface_2021_cr, method=:bilinear)
ice_thickness_resamp = resample(ice_thickness; to=surface_2021_cr, method=:bilinear)

mask = bed_resamp .> 0

# Load GPR points
gpr_file = joinpath(datadir_in, "BonnePierre_ice_thickness_GPR_resampled_50cmtxt_reprojectLambert93.txt")
gpr_points = load_gpr_points(gpr_file)

# Compute ice thickness for 2021, June 2024, and October 2024
ice_thickness_2021 = compute_ice_thickness(surface_2021_cr, bedrock)
ice_thickness_2024_june = compute_ice_thickness(surface_2024_june, bedrock)
ice_thickness_2024_oct = compute_ice_thickness(surface_2024_oct_resamp, bedrock)


# plotting
plot_ice_thickness(
    ice_thickness_resamp;
    gpr_points = gpr_points,
    savepath = joinpath(plots_dir, "ice_thickness_november2024.png"),
    year = "November 2024"
)
plot_ice_thickness(ice_thickness_2021; gpr_points = nothing, savepath = joinpath(plots_dir, "ice_thickness_july2021.png"), year="July 2021")
plot_ice_thickness(ice_thickness_2024_june; gpr_points = nothing, savepath = joinpath(plots_dir, "ice_thickness_june2024.png"), year="June 2024")
plot_ice_thickness(ice_thickness_2024_oct; gpr_points = nothing, savepath = joinpath(plots_dir, "ice_thickness_oct2024.png"), year="October 2024")

plot_bedrock(
    bed_resamp;
    savepath = joinpath(plots_dir, "bedrock_elevation_resamp.png")
)

# save data for WWFS input:
write(joinpath(datadir_WWFS_input, "ice_thickness_2021_july.tif"), ice_thickness_2021,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_june.tif"), ice_thickness_2024_june,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif"), ice_thickness_2024_oct,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_november_glate.tif"), ice_thickness_resamp,force=true)
write(joinpath(datadir_WWFS_input, "bedrock_resamp_1m.tif"), bed_resamp,force=true)
write(joinpath(datadir_WWFS_input, "surface_2021_cr.tif"), surface_2021_cr,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"), surface_2024_oct_resamp,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_june.tif"), surface_2024_june,force=true)
