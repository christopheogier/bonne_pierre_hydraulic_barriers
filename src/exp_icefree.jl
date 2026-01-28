### main_bp_stoch_bedonly_withice_and_icefree_oct2024_routes.jl
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using ArchGDAL
using Rasters
using Serialization
using ProgressMeter
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF  = WhereTheWaterFlows

using CairoMakie

include("functions.jl")    # clean_raster, coord_to_index, etc.
include("plots_makie.jl")  # get_axes_and_matrix

# ---------------- paths ----------------
datadir_WWFS_input = "/scratch-3/cogier/data/BonnePierre_input/WWFS_input"
output_dir         = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis"
mkpath(output_dir)

run_name = "2024_October"

paths = Dict(
    :surface_smooth => joinpath(datadir_WWFS_input, "surface_2024_oct_smooth_01.tif"),
    :thickness      => joinpath(datadir_WWFS_input, "ice_thickness_2024_oct.tif"),
    :bedrock        => joinpath(datadir_WWFS_input, "bedrock_resamp_1m.tif"),
)

# ---------------- load rasters ----------------
surfdem   = clean_raster(Raster(paths[:surface_smooth]))
thickness = clean_raster(Raster(paths[:thickness]))
beddem    = clean_raster(Raster(paths[:bedrock]))

bed_err_std = clean_raster(Raster(joinpath(datadir_WWFS_input, "bedrock_err_std_1m.tif")))

# ---------------- settings ----------------
N = 1000  # Monte Carlo samples


gamma       = [0, WWFS.GAMMA][1]
bnd_as_sink = true
drain_pits  = true

# ---------------- uncertainty model: ONLY bedrock uncertain ----------------
cov_fn = WWFS.GRF.gaussian_kernel   # or WWFS.GRF.exponential_kernel
range_bed = 247.0
corr_length_bed = range_bed / sqrt(2)  # 175 m

zero_uc() = Uncertainty(absuc=0.0, reluc=0.0)

beddem_uc = Uncertainty(
    absuc = bed_err_std,
    reluc = 0.0,
    correlation_length = corr_length_bed,
    covariance_fn = cov_fn
)

source_uc = Uncertainty()

# ---------------- sink definition ----------------
x0, y0 = 962468.722, 6431539.971
i0, j0 = coord_to_index(surfdem, x0, y0)
println("Sink pixel index: ", (i0, j0))
sink_areas = (point = [CartesianIndex(i0, j0)],)

# ---------------- grid spacing ----------------
xdim, ydim = dims(surfdem)
dx = step(xdim)
println("dx = ", dx, " m")

# =============================================================================
# 1) STOCHASTIC SUBGLACIAL ROUTING (WITH ICE MASK) — bedrock uncertainty only
# =============================================================================
rmask_ice     = thickness .> 0
floatfrac_ice = 1.0 .* ones(size(surfdem))
source_ice    = ones(size(surfdem))

println("🔄 Running WWFS stochastic SUBGLACIAL (bedrock uncertainty only, WITH ICE) for $run_name | N=$N")

model_sub, get_sample_sub, aggregate_sub = WWFS.make_fns(
    dx,
    surfdem,        zero_uc(),       # surface deterministic
    beddem,         beddem_uc,        # bedrock uncertain
    floatfrac_ice,  zero_uc(),        # flotation deterministic
    source_ice,     source_uc,
    sink_areas,
    rmask_ice
)

aggr_sub = map_mc(model_sub, get_sample_sub, aggregate_sub, N)
outfile_sub = joinpath(output_dir, "stoch_routes_subglacial_bedonly_n$(N)_$(run_name).jls")
serialize(outfile_sub, aggr_sub)
println("✅ Saved ", outfile_sub)

# =============================================================================
# 2) STOCHASTIC ICE-FREE ROUTING (BEDROCK-ONLY “SURFACE”) — bedrock uncertainty only
# =============================================================================
# ice-free trick: use beddem as "surface" too; keep same grid and sink point
surfdem_free   = beddem
beddem_free    = beddem
rmask_free     = trues(size(beddem))              # full domain active
floatfrac_free = 1.0 .* ones(size(beddem))
source_free    = ones(size(beddem))

println("🔄 Running WWFS stochastic ICE-FREE (bedrock uncertainty only) for $run_name | N=$N")

model_free, get_sample_free, aggregate_free = WWFS.make_fns(
    dx,
    surfdem_free,    zero_uc(),     # "surface" deterministic (bedrock)
    beddem_free,     beddem_uc,     # bedrock uncertain (same uncertainty)
    floatfrac_free,  zero_uc(),
    source_free,     source_uc,
    sink_areas,
    rmask_free
)

aggr_free = map_mc(model_free, get_sample_free, aggregate_free, N)
outfile_free = joinpath(output_dir, "stoch_routes_icefree_bedonly_n$(N)_$(run_name).jls")
serialize(outfile_free, aggr_free)
println("✅ Saved ", outfile_free)

# =============================================================================
# Reconstruct routing maps from aggregated outputs and plot both (two colorbars)
# =============================================================================
# Note: aggr*.areas is the aggregated upslope area (same thing you already plot successfully)
area_sub  = Raster(aggr_sub.areas,  dims(thickness))
area_free = Raster(aggr_free.areas, dims(thickness))

function plot_two_routings(
    thickness::Raster,
    area_sub::Raster,
    area_free::Raster,
    savepath::String;
    area_threshold::Float64 = 1e5
)
    glacier_mask = thickness .> 0
    mask_array   = collect(Bool.(glacier_mask))

    x, y, _ = get_axes_and_matrix(thickness)
    _, _, Z_sub_raw  = get_axes_and_matrix(area_sub)
    _, _, Z_free_raw = get_axes_and_matrix(area_free)

    # mask & clip for readability (same threshold for both)
    Z_sub  = copy(Z_sub_raw)
    Z_free = copy(Z_free_raw)

    Z_sub[.!mask_array]  .= NaN
    Z_free[.!mask_array] .= NaN

    Z_sub[Z_sub .<= 0]   .= NaN
    Z_free[Z_free .<= 0] .= NaN

    Z_sub_clip  = clamp.(Z_sub,  0, area_threshold)
    Z_free_clip = clamp.(Z_free, 0, area_threshold)

    # pick common min for nice comparable scaling (keep simple)
    finite_sub  = Z_sub_clip[isfinite.(Z_sub_clip) .& (Z_sub_clip .> 0)]
    finite_free = Z_free_clip[isfinite.(Z_free_clip) .& (Z_free_clip .> 0)]
    if isempty(finite_sub) || isempty(finite_free)
        error("No finite routing values to plot (check masks/threshold).")
    end

    vmin = min(minimum(finite_sub), minimum(finite_free))
    vmax = area_threshold

    fig = Figure(size=(900, 650))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "X (m)",
        ylabel = "Y (m)",
        titlesize = 18,
        xlabelsize = 16,
        ylabelsize = 16,
        xticklabelsize = 16,
        yticklabelsize = 16,
    )

    # --- Plot both as threshold contours OR as semi-transparent heatmaps ---
    # User asked for "stochasticly ... using two colorbar ... blue for subglacial and red for bedrock-only"
    # -> heatmaps with two independent colorbars.

    cmap_sub  = cgrad([:white, "#8EC1FF", "#1F78B4", "#08306B"])      # blue scale
    cmap_free = cgrad([:white, "#FDBBA1", "#FB6A4A", "#CB181D"])      # red scale

    hm_sub = heatmap!(ax, x, y, Z_sub_clip;
        colormap   = cmap_sub,
        colorrange = (vmin, vmax)
    )

    # overlay ice-free routing in red (slightly transparent so blue still visible)
    hm_free = heatmap!(ax, x, y, Z_free_clip;
        colormap   = cmap_free,
        colorrange = (vmin, vmax),
        transparency = true
    )
    hm_free.attributes.alpha[] = 0.55

    # ---- Glacier outline ----
    contour!(ax, x, y, Int.(mask_array); levels=[0.5], color=:black, linewidth=1)

    # ---- Two colorbars (same scaling) ----
    ticks_vals   = range(vmin, vmax; length=3)
    ticks_labels = string.(Int.(round.(ticks_vals)))

    Colorbar(fig[1, 2], hm_sub;
        ticks        = (ticks_vals, ticks_labels),
        label        = "Subglacial routing (upslope area, m²)",
        height       = 250,
        labelsize    = 14,
        ticklabelsize = 13
    )

    Colorbar(fig[1, 3], hm_free;
        ticks        = (ticks_vals, ticks_labels),
        label        = "Ice-free routing (upslope area, m²)",
        height       = 250,
        labelsize    = 14,
        ticklabelsize = 13
    )

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved routes comparison plot to: $savepath")
    return fig
end

plotfile = joinpath(output_dir, "routes_stoch_subglacial_vs_icefree_bedonly_n$(N)_$(run_name).png")
plot_two_routings(thickness, area_sub, area_free, plotfile; area_threshold=area_threshold)
