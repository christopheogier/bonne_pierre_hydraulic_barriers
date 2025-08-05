# prepare.jl

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


"""
prepare.jl

Prepare script to load data, GPR points, plot ice thickness and bedrock elevation, surface DEMs, uncertainties,
resample data to common resolution.
"""

datadir_in = "/scratch-3/cogier/data/BonnePierre_input"
datadir_out = "/scratch-3/cogier/data/BonnePierre_output"
plots_dir = joinpath(datadir_out, "plots")
datadir_WWFS_input = joinpath(datadir_in, "WWFS_input")

#ENV["GKSwstype"] = "100"  # Use "file output" mode for GR when plotting


# Load bedrock and ice thickness (original resolution)
bedrock = load_bedrock(joinpath(datadir_in, "BED_10m_l93_20smooth.tif"))
ice_thickness_2024_nov = load_ice_thickness(joinpath(datadir_in, "Hall_10m_l93_GlaTE20.asc"))  # ice thickness from air eth and glate on 6.11.2024

# load bedrock plus and minus from GPR
bedrock_gpr_plus_10m = load_bedrock(joinpath(datadir_in, "BED_10m_l93_20smooth_+5m.tif")) 
bedrock_gpr_minus_10m = load_bedrock(joinpath(datadir_in, "BED_10m_l93_20smooth_-5m.tif"))


# Load and clean surface DEM cropped to bedrock extent
surface_2021_cr = load_surface(joinpath(datadir_in, "Berarde_Res1.0_CompElevmean_merged.tif"), bedrock) #1m
surface_2024_june = load_surface(joinpath(datadir_in, "Lidar_juin2024_Bonne_Pierre_1m.tif"), bedrock) #1m
surface_2024_oct = load_surface(joinpath(datadir_in, "Lidar_oct2024_Bonne_Pierre_glacier.tif"), bedrock) # 50cm

# Resample October DEM to 2021 DEM resolution: from 50cm to 1m
surface_2024_oct_resamp = resample(surface_2024_oct; to=surface_2021_cr, method=:bilinear)

# Resample bedrock and ice thickness to surface 2021 resolution
bed_resamp = resample(bedrock; to=surface_2021_cr, method=:bilinear)
ice_thickness_2024_nov_resamp = resample(ice_thickness_2024_nov; to=surface_2021_cr, method=:bilinear)
bedrock_gpr_plus = resample(bedrock_gpr_plus_10m; to=bed_resamp, method=:bilinear)
bedrock_gpr_minus = resample(bedrock_gpr_minus_10m; to=bed_resamp, method=:bilinear)

# mask
mask = bed_resamp .> 0

### Compute surface uncertainties for October 2024 (due to smoothing)
#smoothing
smooth_coeff = 0.1 # as fraction of thickness
smooth_half_window = smooth_coeff / 2
x, y = dims(surface_2024_oct_resamp)
# below y = x is a trick as WWFS.smooth_surface test: @assert dy==dx and here dy = -1 (dx=1)
surface_2024_oct_smooth = WWFS.smooth_surface(x, x, surface_2024_oct_resamp, bed_resamp, smooth_half_window, mask)
# Determinstic error field:
# one sigma standard deviation:
surface_2024_oct_std = abs.(surface_2024_oct_resamp .- surface_2024_oct_smooth)
# if considered the surface in between the two surfaces, then the error is half of the difference
surface_2024_oct_resamp_avg = (surface_2024_oct_resamp .+ surface_2024_oct_smooth) ./ 2  # surface_2024_oct_resamp is the ground truth
surface_2024_oct_std_bis = abs.(surface_2024_oct_resamp_avg .- surface_2024_oct_smooth) ./ 1 # divided by one because hypothesis: err = ±1σ ≈ 68%

### Compute GPR uncertainty maps 
u_plus_gpr  = bedrock_gpr_plus .- bed_resamp    # ≥ 0
u_minus_gpr = bedrock_gpr_minus .- bed_resamp   # ≤ 0

### Compute Glate uncertainty maps

# Extract glacier outline points directly from the thickness raster
thickness = ice_thickness_2024_nov_resamp
outline_polylines = extract_outline_from_thickness(thickness)

# Choose the longest outline (main glacier polygon)
outline_points = isempty(outline_polylines) ? [] : reduce(vcat, outline_polylines)

# Load and extend GPR dataset with outline points (h = 0)
gpr_df = load_and_extend_gpr(
    joinpath(datadir_in, "BonnePierre_ice_thickness_GPR_resampled_50cmtxt_reprojectLambert93.txt"),
    outline_points
)

# Compute distance raster to nearest GPR point
distance_raster = compute_distance_to_gpr(gpr_df, bed_resamp)

# Compute empirical uncertainty bounds based on mean ice thickness
h_mean_2024 = mean(thickness[thickness .> 0]) 
println("Mean ice thickness for November 2024: ", h_mean_2024)
u_minus_interp, u_plus_interp = unc_propagate(distance_raster, h_mean_2024)
# write
write(joinpath(datadir_WWFS_input, "bed_err_plus_interpolation_1m.tif"), u_plus_interp, force=true)

# Mask uncertainty outside glacier domain
u_minus_interp[.!mask] .= NaN
u_plus_interp[.!mask] .= NaN

#### Total uncertainty bedrock

# Load GPR uncertainty maps = assumed to be symetric from now own, based on the minus (larger and thus mor conservative)
u_plus_bed = sqrt.(u_plus_gpr.^2 .+ u_plus_interp.^2)
u_minus_bed = -sqrt.(u_minus_gpr.^2 .+ u_minus_interp.^2)

#one need to define a symmetric uncertainty for the bedrock (plus minus sigma, the standard deviation)
bed_err_std = u_minus_bed


# Load GPR points
gpr_file = joinpath(datadir_in, "BonnePierre_ice_thickness_GPR_resampled_50cmtxt_reprojectLambert93.txt")
gpr_points = load_gpr_points(gpr_file)
# make a filter to remove gpr points that are outside the bedrock extent
geometries = [(x, y) for (x, y) in eachrow(gpr_points)]
vals = extract(ice_thickness_2024_nov, geometries)
# (assumes single-layer raster — extract returns NamedTuple with a single key besides `geometry`)
raster_field = first(keys(vals[1]))  # e.g. :geometry or :ice_thickness
raster_name = filter(k -> k != :geometry, keys(vals[1]))[1]  # get actual raster name key
mask = [!ismissing(v[raster_name]) && v[raster_name] > 0 for v in vals]
# Step 5: Filter GPR points
gpr_points = gpr_points[mask, :]


# Compute ice thickness for 2021, June 2024, and October 2024
ice_thickness_2021 = compute_ice_thickness(surface_2021_cr, bedrock)
ice_thickness_2024_june = compute_ice_thickness(surface_2024_june, bedrock)
ice_thickness_2024_oct = compute_ice_thickness(surface_2024_oct_resamp, bedrock)


# plotting

# Plot all four ice thickness maps
ice_thickness_rasters = [
    (ice_thickness_2024_nov_resamp, "ice_thickness_2024_nov_resamp"),
    (ice_thickness_2021, "ice_thickness_2021"),
    (ice_thickness_2024_june, "ice_thickness_2024_june"),
    (ice_thickness_2024_oct, "ice_thickness_2024_oct")
]

for (raster, name) in ice_thickness_rasters
    savepath = joinpath(plots_dir, name)
    plot_ice_thickness(raster;gpr_points=gpr_points,savepath=savepath)
end

# bedrock
plot_bedrock(bed_resamp; gpr_points=gpr_points, savepath=joinpath(plots_dir, "bedrock_elevation_resamp.png"), glacier_outline_raster = ice_thickness_2024_nov_resamp)

# plot bed Uncertainties
plot_uncertainty_bed(u_plus_gpr, u_minus_gpr,
    "Bedrock: GPR Uncertainty", "GPR +5 m", "GPR -5 m",
    joinpath(plots_dir, "bedrock_gpr_uncertainty.png"))
  
plot_uncertainty_bed(u_plus_interp, u_minus_interp,
    "Bedrock: GLATE interpolation uncertainty", "Interpolation +", "Interpolation -",
    joinpath(plots_dir, "glate_interp_uncertainty.png"))

plot_uncertainty_bed(u_plus_bed,u_minus_bed,
    "Cumulative bedrock uncertainty (GPR + GLATE-Interpolation)", "+ sigma", "- sigma",
    joinpath(plots_dir, "bedrock_all_uncertainty.png"))


# save data for WWFS input:
#ice thickness
write(joinpath(datadir_WWFS_input, "ice_thickness_2021_july.tif"), ice_thickness_2021,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_june.tif"), ice_thickness_2024_june,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif"), ice_thickness_2024_oct,force=true)
write(joinpath(datadir_WWFS_input, "ice_thickness_2024_november_glate.tif"), ice_thickness_2024_nov_resamp,force=true) 
#bedrock
write(joinpath(datadir_WWFS_input, "bedrock_resamp_1m.tif"), bed_resamp,force=true)
write(joinpath(datadir_WWFS_input, "bedrock_gpr_plus5m_1m.tif"),bedrock_gpr_plus, force=true)
write(joinpath(datadir_WWFS_input, "bedrock_gpr_minus5m_1m.tif"), bedrock_gpr_minus, force=true)
write(joinpath(datadir_WWFS_input, "bed_err_plus_gpr5m_1m.tif"), u_plus_gpr, force=true)
write(joinpath(datadir_WWFS_input, "bed_err_minus_gpr5m_1m.tif"), u_minus_gpr, force=true)
write(joinpath(datadir_WWFS_input, "bedrock_err_plus_1m.tif"), u_plus_bed, force=true)
write(joinpath(datadir_WWFS_input, "bedrock_err_minus_1m.tif"), u_minus_bed, force=true)
write(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif"), bed_err_std, force=true)
#surface
write(joinpath(datadir_WWFS_input, "surface_2021_cr.tif"), surface_2021_cr,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"), surface_2024_oct_resamp,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_june.tif"), surface_2024_june,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_oct_err_std_smooth01.tif"), surface_2024_oct_std, force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_oct_smooth_01.tif"), surface_2024_oct_smooth, force=true)
# average of resampled and smoothed DEM
write(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_avg_01smooth.tif"), surface_2024_oct_resamp_avg, force=true)
