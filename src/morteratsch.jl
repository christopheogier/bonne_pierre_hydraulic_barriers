### main.jl

using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using ArchGDAL
using Dates
using Printf
using Rasters
using Statistics
using Serialization
using CairoMakie

include("plots_makie.jl")
include("functions.jl")

using WhereTheWaterFlowsSubglacially
const WWFS = WhereTheWaterFlowsSubglacially

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
datadir_in  = "/scratch-3/cogier/data/Morteratsch"
datadir_out = "/scratch-3/cogier/data/Morteratsch/output"
mkpath(datadir_out)

bedrock_path = joinpath(datadir_in, "bedrock_langhammer.tif")
surface_path = joinpath(datadir_in, "DEM_katarina_2018.tif")

# ------------------------------------------------------------------
# Load rasters
# ------------------------------------------------------------------
bed  = load_bedrock(bedrock_path)
surf = load_surface(surface_path, bed)

# ------------------------------------------------------------------
# Resample bedrock to surface grid
# ------------------------------------------------------------------
bed_resamp = clean_raster(resample(bed; to=surf, method=:bilinear))

# ------------------------------------------------------------------
# Compute ice thickness
# ------------------------------------------------------------------
method_inter = :bilinear
thickness = compute_ice_thickness(surf, bed_resamp, method_inter)

# optional clipping of negative values
thickness = map(x -> isfinite(x) ? max(x, 0.0f0) : NaN32, thickness)

# ------------------------------------------------------------------
# Save rasters
# ------------------------------------------------------------------
write(joinpath(datadir_out, "Morteratsch_surface_2018.tif"), surf; force=true)
write(joinpath(datadir_out, "Morteratsch_bedrock_resamp.tif"), bed_resamp; force=true)
write(joinpath(datadir_out, "Morteratsch_thickness_2018.tif"), thickness; force=true)

# ------------------------------------------------------------------
# Plot thickness (this function is generic enough)
# ------------------------------------------------------------------
plot_ice_thickness(thickness; savepath=joinpath(datadir_out, "Morteratsch_thickness_2018"))

function print_raster_stats(name, rt)
    _, _, Z = get_axes_and_matrix(rt)
    vals = filter(isfinite, vec(Z))

    println("\n--- $name ---")
    println("size = ", size(Z))
    println("finite count = ", length(vals))
    println("nan count = ", count(x -> !isfinite(x), vec(Z)))

    if !isempty(vals)
        println("min = ", minimum(vals))
        println("max = ", maximum(vals))
        println("mean = ", mean(vals))
        println("n zeros = ", count(==(0), vals))
        println("n negative = ", count(<(0), vals))
    end
end

# ------------------------------------------------------------------
# Generic bedrock plot for Morteratsch
# avoids hard-coded Bonne-Pierre overlays in plot_bedrock()
# ------------------------------------------------------------------
function plot_bedrock_simple(rt::Raster, savepath::String)
    x, y, Z = get_axes_and_matrix(rt)
    vmin, vmax = finite_minmax(Z)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "East (m)",
        ylabel = "North (m)"
    )

    ticks_vals, ticks_labels, cbmin, cbmax = colorbar_ticks_with_step(vmin, vmax, 300)

    hm = heatmap!(ax, x, y, Z;
        colormap   = :heat,
        colorrange = (cbmin, cbmax)
    )

    contour!(ax, x, y, Z;
        levels    = range(cbmin, stop=cbmax, step=20),
        linewidth = 0.5,
        color     = :black
    )

    Colorbar(fig[1, 2], hm;
        ticks = (ticks_vals, ticks_labels),
        label = "Elevation (m a.s.l.)"
    )

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved generic bedrock plot to: $savepath")
    return fig
end

plot_bedrock_simple(bed_resamp, joinpath(datadir_out, "Morteratsch_bedrock_resamp.png"))

# ------------------------------------------------------------------
# Deterministic WWFS
# ------------------------------------------------------------------
println("\n🔷 Running deterministic WWFS on Morteratsch")

x, _ = dims(surf)
dx = step(x)

((areas, slen, dir, nout, nin, sinks, pits, c, bnds),
 (sc_locs, kappas, diro, phi),
 (lakes, lakes_free_surf),
 sinkout) = WWFS.waterflows_subglacial(
    surf, bed_resamp, dx;
    gamma = [0, WWFS.GAMMA][1],
    bnd_as_sink = true,
    drain_pits  = true
)

## ANALYS

#analysis = analyze_lakes(lakes, thickness)


# ------------------------------------------------------------------
# Generic phi + upslope area plot
# avoids the unwanted red contour from plot_hydraulic_head_and_flux()
# ------------------------------------------------------------------
function plot_phi_and_area_simple(
    phi::Raster,
    upslope_area::Raster,
    thickness::Raster,
    savepath::String;
    min_threshold::Float64 = 1e4,
    max_threshold::Float64 = 1e6
)
    x, y, Z_phi  = get_axes_and_matrix(phi)
    _, _, Z_area = get_axes_and_matrix(upslope_area)
    _, _, Z_thk  = get_axes_and_matrix(thickness)

    glacier_mask = .!isnan.(Z_thk) .& (Z_thk .> 0)

    # mask outside glacier
    Z_phi[.!glacier_mask]  .= NaN
    Z_area[.!glacier_mask] .= NaN
    Z_area[Z_area .<= 0]   .= NaN

    # clip for display
    Z_area_plot = clamp.(Z_area, min_threshold, max_threshold)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "East (m)",
        ylabel = "North (m)"
    )

    hm = heatmap!(ax, x, y, Z_area_plot;
        colormap   = :blues,
        colorrange = (min_threshold, max_threshold),
        lowclip    = :white,
        highclip   = :black
    )

    # hydraulic head contours only
    vals_phi = filter(isfinite, vec(Z_phi))
    if !isempty(vals_phi)
        vmin_phi = floor(minimum(vals_phi))
        vmax_phi = ceil(maximum(vals_phi))
        levels_phi = collect(vmin_phi:20:vmax_phi)

        contour!(ax, x, y, Z_phi;
            levels    = levels_phi,
            linewidth = 0.6,
            color     = :black
        )
    end

    # glacier outline
    contour!(ax, x, y, Int.(glacier_mask);
        levels    = [0.5],
        color     = :black,
        linewidth = 1.2
    )

    Colorbar(fig[1, 2], hm;
        label = "Upslope area (m²)",
        ticks = [min_threshold, max_threshold]
    )

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved phi + upslope area plot to: $savepath")
    return fig
end

function plot_lake_depth_simple(
    lakes_free_surf::Raster,
    thickness::Raster,
    phi::Union{Raster,Nothing},
    savepath::String
)
    x, y, Z_thk = get_axes_and_matrix(thickness)
    _, _, Z_fs_raw = get_axes_and_matrix(lakes_free_surf)

    glacier_mask = .!isnan.(Z_thk) .& (Z_thk .> 0)

    fig = Figure(size=(800, 600))

    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "East (m)",
        ylabel = "North (m)"
    )

    # keep only glacier cells with finite lake free-surface values
    Z_fs = copy(Z_fs_raw)
    Z_fs[.!glacier_mask] .= NaN32
    Z_fs[.!isfinite.(Z_fs)] .= NaN32
    Z_fs[Z_fs .== 0] .= NaN32

    finite_fs = filter(isfinite, vec(Z_fs))
    hm_fs = nothing

    if !isempty(finite_fs)
        vmin = minimum(finite_fs)
        vmax = maximum(finite_fs)

        hm_fs = heatmap!(ax, x, y, Z_fs;
            colormap   = Reverse(:viridis),
            colorrange = (vmin, vmax)
        )
    end

    # hydraulic head contours
    if phi !== nothing
        _, _, Z_phi = get_axes_and_matrix(phi)
        Z_phi[.!glacier_mask] .= NaN32

        vals_phi = filter(isfinite, vec(Z_phi))
        if !isempty(vals_phi)
            vmin_phi = floor(minimum(vals_phi))
            vmax_phi = ceil(maximum(vals_phi))
            levels_phi = collect(vmin_phi:10:vmax_phi)

            contour!(ax, x, y, Z_phi;
                levels    = levels_phi,
                linewidth = 0.7,
                color     = :grey
            )
        end
    end

    # glacier outline
    contour!(ax, x, y, Float32.(glacier_mask);
        levels    = [0.5],
        color     = :black,
        linewidth = 1.2
    )

    if hm_fs !== nothing
        vmin = minimum(finite_fs)
        vmax = maximum(finite_fs)

        ticks_vals, ticks_labels, cbmin, cbmax =
            colorbar_ticks_with_step(vmin, vmax, 10)

        hm_fs.colorrange[] = (cbmin, cbmax)

        cb = Colorbar(fig[1, 2], hm_fs;
            ticks = (ticks_vals, ticks_labels),
            label = "Lake free-surface elevation (m a.s.l.)"
        )
        cb.height[] = 350
    end

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved lake free-surface plot to: $savepath")
    return fig
end
# ------------------------------------------------------------------
#lake depth plot
# ------------------------------------------------------------------

plot_lake_depth_simple(
    lakes_free_surf,
    thickness,
    phi,
    joinpath(datadir_out, "Morteratsch_lake_depth_2018.png");
   
)
# save phi
write(joinpath(datadir_out, "Morteratsch_phi_2018.tif"), phi; force=true)





println("✅ Morteratsch deterministic workflow completed.")