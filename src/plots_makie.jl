#plot_makie.jl

using CairoMakie
include("LakeAnalysis.jl")
using .LakeAnalysis
using GeoInterface

# --- Your JoG style setup ---
two_column_cm   = 17.8
one_column_cm   = 8.6
font_size_pt    = 10
label_font_size = 12

pt_in_cm = 28.3465

two_column_pt = two_column_cm * pt_in_cm
one_column_pt = one_column_cm * pt_in_cm

dpi        = 600
cm_in_inch = 2.54

px_per_unit = two_column_cm / cm_in_inch * dpi / two_column_pt
pt_per_unit = 1

makie_theme = merge(theme_latexfonts(),
                    Theme(; fontsize=font_size_pt,
                          Axis=(spinewidth=0.5,
                                xtickwidth=0.5,
                                ytickwidth=0.5,
                                xticksize=3,
                                yticksize=3),
                          Colorbar=(spinewidth=0.5, tickwidth=0.5, ticksize=3, size=7),
                          Label=(fontsize=label_font_size, font=:bold),
                          Legend=(rowgap=-8, labelsize=8, framewidth=0.25, padding=(2, 2, 2, 2), margin=(4, 4, 4, 4))))

# Apply the theme globally
#set_theme!(makie_theme)

function get_axes_and_matrix(rt::Raster)
    x, y = collect.(dims(rt))
    x = Float32.(x)
    y = Float32.(y)

    raw = Matrix(rt)
    Z = Array{Float32}(undef, size(raw))

    for j in axes(raw, 2), i in axes(raw, 1)
        val = raw[i, j]
        Z[i, j] = ismissing(val) ? NaN32 : Float32(val)
    end

    return x, y, Z
end

function finite_minmax(Z)
    vals = filter(isfinite, vec(Z))
    if isempty(vals)
        error("No finite values in array.")
    end
    return minimum(vals), maximum(vals)
end


function plot_ice_thickness(rt; gpr_points=nothing, savepath=nothing)
    # Mask zero or negative thickness
    rt_masked = copy(rt)
    rt_masked[rt .<= 0] .= NaN

    x, y, Z = get_axes_and_matrix(rt_masked)
    _, _, Z_full = get_axes_and_matrix(rt)

    vmin, vmax = finite_minmax(Z_full)
    max_idx = argmax(Z_full)
    x_max = x[max_idx[1]]
    y_max = y[max_idx[2]]

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title="Ice thickness (m)")

    # Heatmap and contours
    hm = heatmap!(ax, x, y, Z; colormap=Reverse(:ice), colorrange=(0, vmax))
    # ice thickness countour
    contour!(ax, x, y, Z; levels=0:20:vmax, color=:black)
    glacier_mask = .!isnan.(Z_full) .& (Z_full .> 0)
    contour!(ax, x, y, Int.(glacier_mask); levels=[0.5], color=:black, linewidth=0.8)

    # Max point annotation
    #scatter!(ax, [x_max], [y_max]; color=:red, marker=:xcross, markersize=8)
    #text!(ax, x_max, y_max, text=string(round(vmax, digits=1), " m"), align=(:left, :bottom), fontsize=9,color=:red)

    # GPR points
    if gpr_points !== nothing
        sc = scatter!(ax, gpr_points[:, 1], gpr_points[:, 2]; color=:red, markersize=2, label="GPR measurements")
        #Legend(fig[1, 1], [sc], ["GPR points"], framevisible=false, patchsize=(15,15), labelsize=10)

    end

    # Generate nice intermediate ticks between vmin and vmax, e.g. 5 ticks total
    nticks = 5
    ticks_vals = range(0, vmax, length=nticks)
    ticks_labels = string.(Int.(round.(ticks_vals)))

    # Put colorbar into fig[1, 2], with the same height as ax by linking its height
    cb = Colorbar(fig[1, 2], hm; ticks=(ticks_vals, ticks_labels), label="Ice thickness (m)")
    # Match colorbar height to axis height
    cb.height[] = 350  # or whatever pixel height fits your layout better

    if savepath !== nothing
        #save("$savepath.pdf", fig)
        save("$savepath.png", fig; px_per_unit=4)
    end

    return fig
end

function plot_bedrock(rt; gpr_points=nothing, savepath=nothing, glacier_outline_raster=nothing)
    x, y, Z = get_axes_and_matrix(rt)
    
    vmin, vmax = finite_minmax(Z)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title="Bedrock elevation (m a.s.l.)")

    hm = heatmap!(ax, x, y, Z; colormap=:thermal, colorrange=(vmin, vmax))
    # 20m contour lines
    contour!(ax, x, y, Z; levels=range(vmin, stop=vmax, step=20), linewidth=0.5, color=:black)

    if glacier_outline_raster !== nothing
        x_ice, y_ice, Z_ice = get_axes_and_matrix(glacier_outline_raster)
        #contour!(ax, x_ice, y_ice, Z_ice; levels=0.1:0.1, linewidth=1.0, color=:black)
    end

    if gpr_points !== nothing
    sc = scatter!(ax, gpr_points[:, 1], gpr_points[:, 2]; color=:black, markersize=2, label="GPR measurements")
    #legend_gpr = Legend(fig[1, 1], sc; framevisible=false, labelsize=8)
    # Optionally add a text label for the legend title above the legend:
    #text!(fig[1, 1], "Legend"; position = :topleft, align = (:left, :top), fontsize=8)
    end

    # Generate nice intermediate ticks between vmin and vmax, e.g. 5 ticks total
    nticks = 5
    ticks_vals = range(vmin, vmax, length=nticks)
    ticks_labels = string.(Int.(round.(ticks_vals)))

    # Put colorbar into fig[1, 2], with the same height as ax by linking its height
    cb = Colorbar(fig[1, 2], hm; ticks=(ticks_vals, ticks_labels), label="Elevation (m a.s.l.)")
    # Match colorbar height to axis height
    cb.height[] = 350  # or whatever pixel height fits your layout better


    if savepath !== nothing
        #save("$savepath.pdf", fig)
        save(savepath, fig; px_per_unit=4)
    end

    return fig
end

function plot_uncertainty_bed(r1::Raster, r2::Raster, title::String, subtitle1::String, subtitle2::String, savepath::String)
    x, y, Z1 = get_axes_and_matrix(r1)
    _, _, Z2 = get_axes_and_matrix(r2)

    finite_vals = vcat(Z1[isfinite.(Z1)], Z2[isfinite.(Z2)])
    if isempty(finite_vals)
        @warn "No valid values to plot for $savepath"
        return nothing
    end

    # Define actual data bounds and symmetric color range for white at 0
    vmin_data = minimum(finite_vals)
    vmax_data = maximum(finite_vals)
    vmax_abs = ceil(max(abs(vmin_data), abs(vmax_data)))
    colorrange = (-vmax_abs, vmax_abs)
    cmap = cgrad(:balance, scale=colorrange)

    # Define ticks (integers, always including 0, and both extrema)
    nticks = 5
    ticks_vals = collect(round.(range(vmin_data, vmax_data; length=nticks)))
    if 0 ∉ ticks_vals
        push!(ticks_vals, 0)
        sort!(ticks_vals)
    end
    ticks_labels = string.(Int.(ticks_vals))

    # Begin plotting
    fig = Figure(size=(900, 450), fontsize=10)

    ax1 = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title=subtitle1)
    hm1 = heatmap!(ax1, x, y, Z1; colormap=cmap, colorrange=colorrange)

    ax2 = Axis(fig[1, 2]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title=subtitle2)
    heatmap!(ax2, x, y, Z2; colormap=cmap, colorrange=colorrange)

    # Colorbar with clean integer ticks, 0 centered
    Colorbar(fig[1, 3], hm1;
        label = "Uncertainty (m)",
        height = 300,
        ticks = (ticks_vals, ticks_labels)
    )

    Label(fig[0, :], title; fontsize=12, font=:bold)

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved uncertainty plot to: $savepath")

    return fig
end

# Draw outer rings from a Polygon/MultiPolygon vector file (SHP/GPKG).
# Reprojects to the raster CRS if needed.
function _overlay_depressions!(ax, vec_path::AbstractString, raster_for_crs::Raster)
    ArchGDAL.read(vec_path) do ds
        lyr = ArchGDAL.getlayer(ds, 0)

        # CRS handling
        src_srs = ArchGDAL.getspatialref(lyr)
        dst_srs = ArchGDAL.importWKT(String(crs(raster_for_crs)))
        transf = (src_srs !== nothing && ArchGDAL.toWKT(src_srs) != ArchGDAL.toWKT(dst_srs)) ?
                 ArchGDAL.createcoordinatetransform(src_srs, dst_srs) : nothing

        for feat in lyr
            g = ArchGDAL.getgeom(feat)

            # Uniformly iterate polygons → rings via GeoInterface
            # Wrap to MultiPolygon to handle both Polygon & MultiPolygon
            mp = GeoInterface.MultiPolygon(g)
            for poly in GeoInterface.getgeom(mp)
                # exterior ring is ring index 1; but draw all rings (exterior + holes) anyway
                for ring in GeoInterface.getgeom(poly)
                    coords = GeoInterface.coordinates(ring)  # Vector of (x,y)
                    xs = Vector{Float64}(undef, length(coords))
                    ys = Vector{Float64}(undef, length(coords))
                    @inbounds for i in eachindex(coords)
                        x, y = coords[i]
                        if transf !== nothing
                            x, y, _ = ArchGDAL.transform_point(transf, x, y, 0.0)
                        end
                        xs[i] = x; ys[i] = y
                    end
                    lines!(ax, xs, ys; color = (:black, 0.35), linewidth = 1.0)
                end
            end
        end
    end
end


function plot_lake_depth(
    lakes::Raster,
    thickness::Raster,
    analysis::LakeAnalysisResult,
    phi::Union{Raster, Nothing},
    savepath::String;
    min_depth::Float64 = 2.0,
    show_all_lakes::Bool = false,
    area::Union{Raster, Nothing} = nothing,
    area_threshold::Float64 = 1e5,
    depressions_path::Union{Nothing,String} = nothing
)
    glacier_mask = thickness .> 0
    x, y, Z = get_axes_and_matrix(lakes)
    mask_array = collect(Bool.(glacier_mask))

    # Apply min_depth and glacier mask
    Z_lake = copy(Z)
    Z_lake[Z_lake .< min_depth] .= NaN
    Z_lake[.!mask_array] .= NaN
    Z_lake[Z_lake .== 0] .= NaN

    
    vmin, vmax = finite_minmax(Z_lake)
    vmin = min_depth

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "X (m)",
        ylabel = "Y (m)",
        title = "Water pocket depth (m > $(min_depth))"
    )

    # Plot lake depth
    hm = heatmap!(ax, x, y, Z_lake;
        colormap = :blues,
        colorrange = (vmin, vmax)
    )

    # Plot largest lake outline
    if any(analysis.LargestLake.mask)
        labeled = Int.(analysis.LargestLake.mask)
        #contour!(ax, x, y, labeled; levels=[0.5], color=:red, linewidth=1.5)
    end

    # Optionally plot all lake masks
    if show_all_lakes
        for (_, mask) in analysis.lake_masks
            if any(mask)
                labeled = Int.(mask)
                contour!(ax, x, y, labeled; levels=[0.5], color=:red, linewidth=1)
            end
        end
    end

    # Overlay upslope area outline (if provided)
    if area !== nothing
        _, _, Z_area = get_axes_and_matrix(area)
        area_mask = (Z_area .> area_threshold) .& mask_array  
        area_int = Int.(area_mask)
        contour!(ax, x, y, area_int;
            levels = [0.5],
            color = (:darkblue,0.4),
            linewidth = 1.0
        )
    end

   # Overlay outlines 
   if depressions_path !== nothing
        _overlay_depressions!(ax, depressions_path, lakes) 
   end

    # Overlay hydraulic head contours
    if phi !== nothing
        _, _, Z_phi = get_axes_and_matrix(phi)
        Z_phi[.!mask_array] .= NaN
        vmin_phi = floor(minimum(Z_phi[isfinite.(Z_phi)]), digits=0)
        vmax_phi = ceil(maximum(Z_phi[isfinite.(Z_phi)]), digits=0)
        levels = collect(vmin_phi:10:vmax_phi)
        contour!(ax, x, y, Z_phi; levels=levels, linewidth=0.8, color=:black)
    end

    # Glacier outline
    contour!(ax, x, y, mask_array; levels=[0.5], color=:black, linewidth=1.2)

    # Volume annotations
    total_vol = round(Int, sum(analysis.stats.volume))
    max_vol = round(Int, analysis.LargestLake.volume)

    text!(
        ax, x[1], y[end],
        text = "Total volume: $(total_vol) m³\nLargest water pocket: $(max_vol) m³",
         fontsize = 10, color = :black
    )

    # Generate nice intermediate ticks between vmin and vmax, e.g. 5 ticks total
    nticks = 4
    ticks_vals = range(vmin, vmax, length=nticks)
    ticks_labels = string.(Int.(round.(ticks_vals)))
    # colorbar
    cb = Colorbar(fig[1, 2], hm; ticks=(ticks_vals, ticks_labels), label = "Water pocket depth (m)")
    cb.height[] = 350



    # Legends
    lines!(ax, [NaN], [NaN]; color = :black, linewidth = 0.8, label = "Hydraulic head (10m intervals)")
    lines!(ax, [NaN], [NaN]; color = :red, linewidth = 1.0, label = "Water pocket outlines")
    if area !== nothing
        lines!(ax, [NaN], [NaN]; color = :darkblue, linewidth = 1.0, label = "Upslope area > $(Int(area_threshold)) m²")
    end

    Legend(fig, ax; tellwidth = false, tellheight = false, halign = :left, valign = :top, framevisible = false)

    save(savepath, fig; px_per_unit = 4)
    println("✅ Saved lake depth plot with outlines and annotations to: $savepath")
    return fig
end


function plot_hydraulic_head_and_flux(
    phi::Raster,
    upslope_area::Raster,
    thickness::Raster,
    savepath::String;
    min_threshold::Float64 = 1e4,
    max_threshold::Float64 = 1e6
)
    x, y, Z_phi = get_axes_and_matrix(phi)
    _, _, Z_flux = get_axes_and_matrix(upslope_area)

    # Hydraulic head contours
    vmin_phi, vmax_phi = finite_minmax(Z_phi)

    levels_phi = collect(vmin_phi:20:vmax_phi)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
        aspect = DataAspect(),
        xlabel = "X (m)",
        ylabel = "Y (m)",
        title = "Upslope area and hydraulic head"
    )

    # Base heatmap
    hm = heatmap!(ax, x, y, Z_flux;
        colormap = :blues,
        colorrange = (min_threshold, max_threshold),
        lowclip = :white,
        highclip = :black
    )

    # Build binary mask where upslope area exceeds max_threshold
    highlight_mask = Z_flux .> min_threshold
    highlight_int = Int.(highlight_mask)

    # Add contour outline — adjust linewidth for visual thickness
    contour!(ax, x, y, highlight_int;
        levels = [0.5],       # Contour between 0 and 1
        color = :red,
        linewidth = 1.5        # Try 2.0 or 3.0 for thicker effect
    )

    # Contours of hydraulic head
    contour!(ax, x, y, Z_phi; levels=levels_phi, linewidth=0.6, color=:black)

    # Glacier outline
    glacier_mask = thickness .> 0
    mask_array = collect(Bool.(glacier_mask))
    contour!(ax, x, y, mask_array; levels=[0.5], color=:black, linewidth=1.2)

    # Colorbar synced with main heatmap
    cb = Colorbar(fig[1, 2], hm;
        label = "Upslope area (m²)",
        ticks = [min_threshold, max_threshold]
    )
    cb.height[] = 350

    # Add dummy line for legend entry
    lines!(ax, [NaN], [NaN];
        color = :red,
        linewidth = 1.5,
        label = "Upslope area > $(Int(min_threshold)) m²"
    )

    # Add legend in top-left of axis
    Legend(fig, ax;
        tellwidth = false,
        tellheight = false,
        halign = :left,
        valign = :top,
        framevisible = false
    )


    save(savepath, fig; px_per_unit=4)
    println("✅ Saved hydraulic head and flux plot to: $savepath")
    return fig
end

function boxplot_lake_vol_stoch(
    aggr_main;
    aggr2 = nothing,
    aggr3 = nothing,
    aggr4 = nothing,
    labels::Vector{String} = ["all unc.", "bed. unc.", "surf. unc.", "flot. unc."],
    savepath::String = "lake_fs_volume_boxplot.png"
)
    lake_vols = [aggr_main.lake_fs_vol]

    if aggr2 !== nothing push!(lake_vols, aggr2.lake_fs_vol) end
    if aggr3 !== nothing push!(lake_vols, aggr3.lake_fs_vol) end
    if aggr4 !== nothing push!(lake_vols, aggr4.lake_fs_vol) end

    fig = Figure(size = (100 * length(lake_vols) + 300, 400))
    ax = Axis(fig[1, 1],
        title = "Total lake (>2m) volume distribution (Monte Carlo)",
        ylabel = "Volume (m³)",
        xticks = (1:length(labels), labels)
    )

    x = vcat([fill(i, length(v)) for (i, v) in enumerate(lake_vols)]...)
    y = vcat(lake_vols...)
    boxplot!(ax, x, y)

    save(savepath, fig)
    println("✅ Saved stochastic lake volume boxplot to: $savepath")
    return fig
end

# Helper: turn extract(...) output into (distance, values) vectors
_profile_vals(r::Raster, pts) = begin
    rows = extract(r, pts; geometry=false, skipmissing=false)
    vraw = getproperty.(rows, name(r))              # Vector{Union{Missing, T}}
    vals = Float64.(coalesce.(vraw, NaN))           # keep length, convert missings to NaN
    return vals
end


"""
plot_profiles(bedrock, surface_raw, phi, lakes_free_surf, fname;
              surface_smooth=nothing, surface_fill=nothing)

Plots 1D profiles along the fixed transect:
 A = (962611.20, 6431486.12)  (upstream, right on map)
 B = (962397.32, 6431558.42)  (downstream, left on map)

Saves figure to `fname`.
"""
function plot_profiles(bedrock::Raster, surface_raw::Raster, phi::Raster, lakes_free_surf::Raster, fname::AbstractString;
                       surface_smooth::Union{Raster,Nothing}=nothing,
                       surface_fill::Union{Raster,Nothing}=nothing)

    A = (962632.47,6431521.78) # upstream # picked manually to fit the uper head and the seal (more or less)
    B = (962420.11,6431549.94) # dowstream (stream exit in june 2024)
    dx = 2.0 # sampling, larger values smooth the lines
    x1, y1 = B
    x2, y2 = A
    L = hypot(x2 - x1, y2 - y1)
    n = max(1, floor(Int, L/dx)) + 1
    ts = range(0.0, 1.0; length=n)
    pts = [(x1 + t*(x2 - x1), y1 + t*(y2 - y1)) for t in ts]
    dist = collect(range(0.0, L; length=n))

    # Use the helper for each raster
    z_bed   = _profile_vals(bedrock,         pts)
    z_surf  = _profile_vals(surface_raw,     pts)
    z_phi   = _profile_vals(phi,             pts)
    h_lake  = _profile_vals(lakes_free_surf, pts)
    z_smooth = isnothing(surface_smooth) ? nothing : _profile_vals(surface_smooth, pts)
    z_fill   = isnothing(surface_fill)   ? nothing : _profile_vals(surface_fill,   pts)


    d_smooth, z_smooth = nothing, nothing
    if surface_smooth !== nothing
        z_smooth = _profile_vals(surface_smooth, pts)
    end
    d_fill, z_fill = nothing, nothing
    if surface_fill !== nothing
        z_fill = _profile_vals(surface_fill, pts)
    end

    # --- Bedrock uncertainty shading ---
    err_minus_r = Raster("/scratch-3/cogier/data/BonnePierre_input/WWFS_input/bedrock_err_minus_1m.tif")
    err_plus_r  = Raster("/scratch-3/cogier/data/BonnePierre_input/WWFS_input/bedrock_err_plus_1m.tif")

    z_bed_minus = z_bed + _profile_vals(err_minus_r, pts)
    z_bed_plus  = z_bed + _profile_vals(err_plus_r,  pts)

    # Figure
    fig = Figure(size=(800, 600))
    ax  = Axis(fig[1,1], xlabel="Distance along transect (m)", ylabel="Elevation (m)",
               title="Profiles along A→B")

    # Plot main profiles
    lines!(ax, dist, z_surf,  label="surface (LiDAR)")
    if z_smooth !== nothing
        lines!(ax, dist, z_smooth, label="surface (smoothed)", linestyle=:dash)
    end
    
    if z_fill !== nothing
        lines!(ax, dist, z_fill,   label="supraglacial lake surface (ice equivalent)", linestyle=:dot)
    end
    lines!(ax, dist, z_bed,   label="bedrock (+/- sigma)")

    # Shaded uncertainty band
    poly!(ax,
          vcat(dist, reverse(dist)),
          vcat(z_bed_minus, reverse(z_bed_plus)),
          color=(:gray, 0.3), strokewidth=0)

    lines!(ax, dist, z_phi,   label="hydraulic head φ")

    # Lake free-surface depth: plot as bed = depth
    lines!(ax, dist, z_bed + h_lake,  label="lake_free_surface (depth)")

    axislegend(ax, position=:rb, framevisible=false)

    # Force y-axis limits
    y_min = minimum(z_bed)  - 15 # m
    y_max = maximum(z_surf) + 5  # m
    ylims!(ax, y_min, y_max)

    vlines!(ax, [0, dist[end]]; color=:gray, linestyle=:dash, linewidth=1)
    # DOUBLE CHECK A AND B 
    text!(ax, dist[end]-5, 0, text="A (upstream)", align=(:right, :bottom), space=:data)
    text!(ax, 5, 0, text="B (downstream)", align=(:right, :bottom), space=:data)

    save(fname, fig; px_per_unit=3)
    println("✅ saved transect profile to: $fname")
    return fig
end
