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
method_inter =:bilinear #Bilinear, near (weird) or cubic 
surface_2024_oct_resamp = resample(surface_2024_oct; to=surface_2021_cr, method=method_inter)
#resample June too (to match bed and surface for later smoothing)
surface_2024_june = clean_raster(resample(surface_2024_june; to=surface_2021_cr, method=method_inter))
println("⚠️  Warning: resample introduced Missing values — run clean_raster() before using.")

# Resample bedrock and ice thickness to surface 2021 resolution
bed_resamp = resample(bedrock; to=surface_2021_cr, method=method_inter)
# smooth the resample bedrock as the interopolation makes noise at the meters scale
#x, y = dims(bed_resamp)
#dx = x[2] - x[1]                  # grid spacing (m)
#pixrad = max(1, round(Int, 10.0/dx))  # 10 m moving window for smoothing. 
#Wconst = fill(pixrad, size(bed_resamp))
#bed_resamp = WWFS.boxcar(bed_resamp, Wconst, trues(size(bed_resamp)))  # no mask, WARNING THAT ENLARGE the bedrock outline, need to mask

ice_thickness_2024_nov_resamp = resample(ice_thickness_2024_nov; to=surface_2021_cr, method=method_inter)
bedrock_gpr_plus = resample(bedrock_gpr_plus_10m; to=bed_resamp, method=method_inter)
bedrock_gpr_minus = resample(bedrock_gpr_minus_10m; to=bed_resamp, method=method_inter)

# Compute ice thickness for 2021, June 2024, and October 2024
ice_thickness_2021 = compute_ice_thickness(surface_2021_cr, bedrock,method_inter)
ice_thickness_2024_june = compute_ice_thickness(surface_2024_june, bedrock,method_inter)
ice_thickness_2024_oct = compute_ice_thickness(surface_2024_oct_resamp, bedrock,method_inter)

# mask
mask = bed_resamp .> 0

### Compute surface uncertainties due to smoothing 
#(Determinstic error field)# one sigma standard deviation


### Compute/store June surface ensemble for stochastic surface sampling
# raw (0.0) already exists as surface_2024_june.tif

smooth_coeffs = collect(0.1:0.1:1.0)  # 0.1, 0.2, ..., 1.0

jun_surfaces = Dict{Float64, Raster}()

for sc in smooth_coeffs
    println("Computing June smoothed surface for smoothing coefficient = ", sc)
    res = surface_uncertainty_from_smoothing(surface_2024_june, bed_resamp, sc, mask)
    jun_surfaces[sc] = res.smooth
end


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

# Compute empirical uncertainty bounds based on distance to GPR points
h_mean_2024 = mean(thickness[thickness .> 0]) 
println("Mean ice thickness for November 2024: ", h_mean_2024)
u_std_interp = unc_propagate(distance_raster) # symetric
# write
write(joinpath(datadir_WWFS_input, "bed_err_std_interpolation_1m.tif"), u_std_interp, force=true)

# Mask uncertainty outside glacier domain
u_std_interp[.!mask] .= NaN

#### Total uncertainty bedrock

# Load GPR uncertainty maps = assumed to be symetric from now own, based on the minus (larger and thus mor conservative)
u_plus_bed = sqrt.(u_plus_gpr.^2 .+ u_std_interp.^2)
u_minus_bed = -sqrt.(u_minus_gpr.^2 .+ u_std_interp.^2)

#one need to define a symmetric uncertainty for the bedrock (plus minus sigma, the standard deviation)
bed_err_std = abs.(u_plus_bed)


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
    plot_ice_thickness(raster;gpr_points=nothing,savepath=savepath)
end

# plot surface hillshade
plot_surface_hillshade(surface_2024_oct_resamp; outline_raster=ice_thickness_2024_oct, savepath=joinpath(plots_dir, "surface_hillshade.png"))

# bedrock
plot_bedrock(bed_resamp; gpr_points=gpr_points, savepath=joinpath(plots_dir, "bedrock_elevation_resamp.png"), glacier_outline_raster = ice_thickness_2024_nov_resamp)

# plot bed Uncertainties
plot_uncertainty_bed(u_plus_gpr; r2 = u_minus_gpr,
    title    = "Bedrock: GPR uncertainty",
    subtitle1 = "GPR +5 m",
    subtitle2 = "GPR -5 m",
    savepath = joinpath(plots_dir, "bedrock_gpr_uncertainty.png"),
)
  
plot_uncertainty_bed(u_std_interp;  # symetric error field
    title="Bedrock: GLATE interpolation uncertainty", subtitle1="Interpolation std +-σ",
    savepath=joinpath(plots_dir, "glate_interp_uncertainty.png"))

plot_uncertainty_bed(bed_err_std;
    title="" , subtitle1="",
    savepath=joinpath(plots_dir, "bedrock_all_uncertainty.png"))


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
write(joinpath(datadir_WWFS_input, "bedrock_err_interp_std_1m.tif"),u_std_interp , force=true)
#surface
write(joinpath(datadir_WWFS_input, "surface_2021_cr.tif"), surface_2021_cr,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_oct_resamp_1m.tif"), surface_2024_oct_resamp,force=true)
write(joinpath(datadir_WWFS_input, "surface_2024_june.tif"), surface_2024_june,force=true)
# June smoothed surface ensemble for stochastic surface sampling
for sc in smooth_coeffs
    tag = replace(@sprintf("%.1f", sc), "." => "")   # 0.1 -> "01", 1.0 -> "10"
    write(
        joinpath(datadir_WWFS_input, "surface_2024_june_smooth_$(tag).tif"),
        jun_surfaces[sc],
        force=true
    )
end
