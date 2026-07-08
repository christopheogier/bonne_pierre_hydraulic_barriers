# function.jl

using Rasters
using ArchGDAL
using CSV, DataFrames, NearestNeighbors
using WhereTheWaterFlows
const WWF = WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially
const WWFR = WhereTheWaterFlows.Randomly
import Contour: contour, lines, coordinates

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
function compute_ice_thickness(surface::Raster, bed::Raster, method_inter::Union{Symbol,String})
    bed_resamp = resample(bed; to=surface, method=method_inter)

    # Compute thickness: surface - bed
    thickness_array = map((s, b) -> begin
        if ismissing(s) || ismissing(b)
            NaN
        else
            val = s - b
            val < 0 ? 0.0 : val  # negative value (bed above surface) are clipped to 0
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

function unc_propagate(dist_raster::Raster) # so that float or integer works
    # Error 1-sigma = 0.0411365519449864 * distance (mètres) + 13.34869312056046
    u_std  = 0.04.*dist_raster .+ 13.35 # in meters
   
    #set err = 0 where GPR points exist (distance = 0)
    #mask_zero = dist_raster .== 0
    #u_std[mask_zero] .= 0.0  # Set u = 0 where distance is zero
    
    
    return u_std
end

function surface_uncertainty_from_smoothing(surface, bed, smooth_coeff, mask)

    smooth_half_window = smooth_coeff / 2

    x, y = dims(surface)

    surface_smooth = WWFS.smooth_surface(x, x, surface, bed, smooth_half_window, mask)

    # full difference between raw and smoothed
    diff = surface .- surface_smooth

    # half difference = plausible amplitude around midpoint
    avg = diff ./ 2

    # midpoint between raw and smoothed
    mid = (surface .+ surface_smooth) ./ 2

    # full difference kept if needed for diagnostics only
    std = abs.(diff)

    return (smooth = surface_smooth, std = std, avg = avg, mid = mid)
end

"""
    coord_to_index(raster, x, y)

Return pixel indices (i, j) of the raster cell closest to coordinate (x, y).

Works for rasters whose x/y axes are ascending or descending.

Throws an error if (x,y) lies outside the raster extent.
"""
function coord_to_index(raster, x, y)
    xd, yd = dims(raster)
    xs = collect(xd)
    ys = collect(yd)

    # helper: get index for ascending or descending axis
    function find_idx(axis, val)
        if axis[1] <= axis[end]
            # ascending
            return searchsortedfirst(axis, val)
        else
            # descending
            k = searchsortedfirst(reverse(axis), val)
            return length(axis) - k + 1
        end
    end

    i = find_idx(xs, x)
    j = find_idx(ys, y)

    # bounds checks
    if !(1 ≤ i ≤ length(xs))
        error("x=$x outside raster x-range $(extrema(xs))")
    end
    if !(1 ≤ j ≤ length(ys))
        error("y=$y outside raster y-range $(extrema(ys))")
    end

    return i, j
end

"""
    _profile_vals(r::Raster, pts; method=:bilinear, fill=NaN)

Sample raster `r` along a polyline given by `pts = [(x,y), ...]`.

- `method = :nearest` or `:bilinear`
- `fill` returned for points outside the raster domain
"""
function _profile_vals(r::Raster, pts::AbstractVector{<:Tuple};
                       method::Symbol = :bilinear,
                       fill = NaN32)

    xdim, ydim = dims(r)
    xs = collect(xdim)
    ys = collect(ydim)
    nx = length(xs)
    ny = length(ys)

    # assume regular spacing
    dx = xs[2] - xs[1]
    dy = ys[2] - ys[1]

    # handle possibly decreasing axes (common in rasters)
    x0 = xs[1]; y0 = ys[1]
    x_increasing = xs[end] > xs[1]
    y_increasing = ys[end] > ys[1]

    out = Vector{Float32}(undef, length(pts))

    @inbounds for k in eachindex(pts)
        x, y = pts[k]

        # fractional index in 1-based coordinates
        fx = (x - x0)/dx + 1
        fy = (y - y0)/dy + 1

        # if axis is decreasing, flip fractional coordinate
        if !x_increasing
            fx = (x0 - x)/dx + 1
        end
        if !y_increasing
            fy = (y0 - y)/dy + 1
        end

        if method === :nearest
            i = round(Int, fx)
            j = round(Int, fy)
            if 1 ≤ i ≤ nx && 1 ≤ j ≤ ny
                v = r[i, j]
                out[k] = ismissing(v) ? fill : Float32(v)
            else
                out[k] = fill
            end

        elseif method === :bilinear
            i0 = floor(Int, fx)
            j0 = floor(Int, fy)
            i1 = i0 + 1
            j1 = j0 + 1

            if 1 ≤ i0 ≤ nx && 1 ≤ i1 ≤ nx && 1 ≤ j0 ≤ ny && 1 ≤ j1 ≤ ny
                tx = Float32(fx - i0)
                ty = Float32(fy - j0)

                v00 = r[i0, j0]
                v10 = r[i1, j0]
                v01 = r[i0, j1]
                v11 = r[i1, j1]

                # if any corner is missing, fall back to fill
                if any(ismissing, (v00, v10, v01, v11))
                    out[k] = fill
                else
                    v00 = Float32(v00); v10 = Float32(v10)
                    v01 = Float32(v01); v11 = Float32(v11)
                    out[k] = (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
                end
            else
                out[k] = fill
            end

        else
            error("Unknown method = $method. Use :nearest or :bilinear.")
        end
    end

    return out
end

using Rasters

"""
Export lakes with volume > vol_thresh to GeoTIFFs for QGIS.

Writes:
- <out_prefix>_lakes_vol_gt_<thresh>_mask.tif   (UInt8, 0/1)

`template` should be a Raster aligned with the masks (e.g. thickness or lakes_free_surf).
"""
function export_big_lake_masks!(
    template::Raster,
    analysis::LakeAnalysisResult,
    out_prefix::String;
    vol_thresh::Float64 = 1000.0
)
    big_inds = findall(v -> v > vol_thresh, analysis.stats.volume)
    n_big = length(big_inds)
    # Preallocate arrays on the template grid
    combined = falses(size(template))
    id = Int32(1)
    for label in big_inds
        mask = analysis.lake_masks[label]
        combined .|= mask
        id += 1
    end
    out_mask   = rebuild(template; data = UInt8.(combined))
    tif_mask   = out_prefix * "_lakes_free_vol_gt_$(Int(vol_thresh)).tif"
    write(tif_mask, out_mask; force=true)
    println("  ✅ Exported big-lake masks ($(n_big) lakes) to:")
    println("     - $tif_mask")
    return (tif_mask=tif_mask, n_big=n_big)
end

function make_fns_surface_ensemble(
    dx,
    surface_ensemble,
    beddem, beddem_uc,
    floatfrac, floatfrac_uc,
    waterinput, waterinput_uc,
    ctch_sinks,
    rmask;
    gamma=[0, WWFS.GAMMA][1],
    min_lake_depth=0.1,
    rhow=WWFS.RHOW,
    rhoi=WWFS.RHOI
)
    base_surface = first(surface_ensemble)
    zero_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.0)

    model(surf, bed, floatfrac, waterinput) = (
        (; surf, bed, dx, floatfrac, waterinput),
        WWFS.waterflows_subglacial(
            surf, bed, dx, floatfrac, waterinput, rmask;
            gamma,
            drain_pits=true,
            bnd_as_sink=true,
            nan_as_sink=true,
            rhow,
            rhoi,
            ctch_sinks
        )
    )

    get_sample = let
        beddem_grf_sampler = WWFR.make_sampler(dx, beddem, beddem_uc)
        floatfrac_grf_sampler = WWFR.make_sampler(dx, floatfrac, floatfrac_uc)
        waterinput_grf_sampler = WWFR.make_sampler(dx, waterinput, waterinput_uc)

        function ()
            surf = rand(surface_ensemble)
            bed = WWFR.make_field_realization(beddem, beddem_grf_sampler, beddem_uc)
            float = WWFR.make_field_realization(floatfrac, floatfrac_grf_sampler, floatfrac_uc)
            water = WWFR.make_field_realization(waterinput, waterinput_grf_sampler, waterinput_uc)
            return surf, bed, float, water
        end
    end

    function reduce!()
        return (
            areas = zeros(Float32, size(base_surface)),
            lakes_depth_fs = zeros(Float32, size(base_surface)),
            n_samples = Ref(0),
            lake_vol_free_surface = Float32[]
        )
    end

    function reduce!(aggr, res)
        input, output = res
        lake_depth_free_surface = output.lakes.depth_free_surface
        aggr.areas .+= output.routing.area.total
        aggr.lakes_depth_fs .+= lake_depth_free_surface
        push!(aggr.lake_vol_free_surface, sum(Float32.(lake_depth_free_surface[lake_depth_free_surface .> min_lake_depth])))
        aggr.n_samples[] += 1
        return aggr
    end

    function reduce!(aggr)
        if aggr.n_samples[] > 0
            aggr.areas ./= aggr.n_samples[]
            aggr.lakes_depth_fs ./= aggr.n_samples[]
        end
        return aggr
    end

    return model, get_sample, reduce!
end


