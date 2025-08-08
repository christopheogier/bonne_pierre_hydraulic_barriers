# function.jl

using Rasters
using ArchGDAL
using CSV, DataFrames, NearestNeighbors

"""
    clean_raster(r::Raster) -> Raster{Float32}

Convert a Raster{Union{Missing, Float32}} to Raster{Float32}, replacing `missing` with `NaN` 
and removing internal `missingval` metadata.
"""
function clean_raster(r::Raster)
    # Replace missing values with NaN32 and ensure Float32 array
    clean_data = Float32[ismissing(v) ? NaN32 : Float32(v) for v in r]

    # Rebuild raster with same spatial dims, no missingval metadata
    return Raster(reshape(clean_data, size(r)), dims(r))
end


function load_ice_thickness(filepath::String)
    rt = Raster(filepath)
    rt = clean_raster(rt)
    rt[rt .< 0] .= 0
    return rt
end

function load_bedrock(filepath::String)
    rt = Raster(filepath)
    rt = clean_raster(rt)
    rt[rt .< 0] .= NaN
    return rt
end

function load_surface(filepath::String, bed::Raster)
    surface = Raster(filepath)
    surface_cr = crop(surface; to=bed)
    return clean_raster(surface_cr)
end

"""
    load_gpr_points(filepath::String) -> Matrix{Float64}

Load GPR point coordinates from a text file with columns: x, y, thickness.
Returns Nx2 matrix.
"""

function load_gpr_points(filepath::String)
    df = CSV.read(filepath, DataFrame; header=false)
    return Matrix(df[:, 1:2])
end


"""
    compute_ice_thickness(surface::Raster, bed::Raster) -> Raster

Computes ice thickness as the difference between a glacier surface raster and a bedrock raster.
The bed raster is resampled to match the surface raster grid using bilinear interpolation.

- If either surface or bed value is `missing`, the result is set to `NaN`.
- If the computed thickness is negative, it is clamped to 0.0.

# Arguments
- `surface::Raster`: Glacier surface elevation raster (e.g. LiDAR or photogrammetry).
- `bed::Raster`: Bedrock elevation raster (may have different resolution or grid).

# Returns
- `Raster`: Ice thickness raster aligned with `surface`.
"""
function compute_ice_thickness(surface::Raster, bed_resamp::Raster)
    bed_resamp = resample(bed; to=surface, method=:bilinear)

    # Compute thickness: surface - bed
    thickness_array = map((s, b) -> begin
        if ismissing(s) || ismissing(b)
            NaN
        else
            val = s - b
            val < 0 ? 0.0 : val
        end
    end, surface, bed_resamp)

    return Raster(Float32.(thickness_array), dims(surface))

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


module UncUtils

using Rasters, CSV, DataFrames, NearestNeighbors
import Contour: contour, lines, coordinates


export load_and_extend_gpr, compute_distance_to_gpr, unc_propagate, extract_outline_from_thickness

function extract_outline_from_thickness(thickness::Raster)
    x, y = collect.(dims(thickness))
    Z = Matrix(thickness) .> 0  # binary mask

    contourset = contour(x, y, Z, 0.5)  # get isolines at 0.5 to outline glacier

    outline_points = []

    for line in lines(contourset)
        xs, ys = coordinates(line)
        append!(outline_points, zip(xs, ys))
    end

    return outline_points  # Vector of (x, y) tuples
end

function load_and_extend_gpr(gpr_path::String, outline_points::Vector{Tuple{Float64, Float64}})
    gpr_df = CSV.read(gpr_path, DataFrame; header=false)
    rename!(gpr_df, [:x, :y, :h, :h1, :h2])
    for (x, y) in outline_points
        push!(gpr_df, (; x, y, h = 0.0, h1 = 0.0, h2 = 0.0))
    end
    return gpr_df
end

function compute_distance_to_gpr(gpr_df::DataFrame, template_raster::Raster)
    x_dim, y_dim = dims(template_raster)
    xs, ys = x_dim.val, y_dim.val
    gpr_points = hcat(gpr_df.x, gpr_df.y)'
    tree = KDTree(gpr_points)
    distance_map = Array{Float64}(undef, length(xs), length(ys))
    for (i, x) in enumerate(xs)
        for (j, y) in enumerate(ys)
            point = [x, y]
            _, dists = knn(tree, point, 1)
            distance_map[i, j] = dists[1]
        end
    end
    return Raster(distance_map, dims(template_raster))
end

function unc_propagate(dist_raster::Raster, h_mean::Real) # so that float or integer works
    u_minus  = (-0.27 .- 0.00077 .* dist_raster) .* h_mean
    u_plus = (0.08  .+ 0.00083 .* dist_raster) .* h_mean  # note this is inverted from Grabe et al because here u is for the bedrock, not the ice thickness
    # set u = 0 if d = 0 (a case not accounted above. Seems like that is what Grab et al did but it was not epxlicit in the paper)
   
    mask_zero = dist_raster .== 0
    u_minus[mask_zero] .= 0.0  # Set u = 0 where distance is zero
    u_plus[mask_zero]  .= 0.0 
    
    return u_minus, u_plus
end

end # module
