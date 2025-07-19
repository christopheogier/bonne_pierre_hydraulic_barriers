

using Rasters
using ArchGDAL
using DelimitedFiles

"""
    load_ice_thickness(filepath::String) -> Raster

Load the ice thickness raster from the given file path.
Replace values < 0 (nodata) by 0.
"""
function load_ice_thickness(filepath::String)
    rt = Raster(filepath)
    rt[rt .< 0] .= 0
    return rt
end

"""
    load_bedrock(filepath::String) -> Raster

Load the bedrock elevation raster from the given file path.
Replace values < 0 (nodata) by NaN.
"""
function load_bedrock(filepath::String)
    rt = Raster(filepath)
    rt[rt .< 0] .= NaN
    return rt
end

function load_surface(filepath::String, bed::Raster)
    # Load raster
    surface = Raster(filepath)
    # Crop to bed extent
    surface_cr = crop(surface; to=bed)
    # Replace negative values with NaN
    surface_cr[surface_cr .< 0] .= NaN
    return surface_cr
end

"""
    load_gpr_points(filepath::String) -> Matrix{Float64}

Load GPR point coordinates from a text file with columns: x, y, thickness.
Returns Nx2 matrix.
"""
function load_gpr_points(filepath::String)
    data = readdlm(filepath, '\t', skipstart=1)
    return data[:, 1:2]
end

"""
    compute_ice_thickness(surface::Raster, bed::Raster) -> Raster

Computes ice thickness as the difference between a glacier surface raster and a bedrock raster.
The bed raster is resampled to match the surface raster grid using bilinear interpolation.
Negative thickness values are set to 0, and invalid surface elevations (e.g., < 0) are masked as 0.

# Arguments
- `surface::Raster`: Glacier surface elevation raster (e.g. LiDAR or photogrammetry).
- `bed::Raster`: Bedrock elevation raster (may have different resolution or grid).

# Returns
- `Raster`: Ice thickness raster aligned with `surface`.
"""
function compute_ice_thickness(surface::Raster, bed::Raster)
    # Crop surface raster to bed extent
    #surface_cropped = crop(surface; to=bed)

    # Optionally mask invalid surface elevations (uncomment if needed)
    # surface_cropped[surface_cropped .< 0] .= NaN

    # Resample bedrock to match surface grid/resolution
    bed_resamp = resample(bed; to=surface, method=:bilinear)

    # Compute thickness (surface minus bed)
    ice_thickness = surface .- bed_resamp

    # Set negative thickness to zero
    ice_thickness[ice_thickness .< 0] .= 0

    # Replace NaN with zero (if any)
    ice_thickness[isnan.(ice_thickness)] .= 0

    return ice_thickness
end


function compute_glacier_outline(thickness::Raster)
    # Convert Raster to raw matrix (assume Band 1)
    thickness_array = Matrix(thickness)
    cs = contours(thickness_array, levels=[0.5])  # 0.5m ice thickness = glacier margin

    # Extract contour paths and map them back to raster coordinates
    xdim, ydim = dims(thickness)
    outline_coords = []

    for c in cs[1].lines
        path = [(xdim[x], ydim[y]) for (x, y) in zip(c.coordinates[:, 1], c.coordinates[:, 2])]
        push!(outline_coords, path)
    end

    return outline_coords  # Vector of vectors of (x, y)
end
