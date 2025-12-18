#plot_makie.jl

using CairoMakie
include("LakeAnalysis.jl")
using .LakeAnalysis
using GeoInterface
using Shapefile

# --- Your JoG style setup ---
two_column_cm   = 17.8
one_column_cm   = 8.6
font_size_pt    = 12
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
                          Legend=(rowgap=-8, labelsize=12, framewidth=0.25, padding=(2, 2, 2, 2), margin=(4, 4, 4, 4))))

# Apply the theme globally

#set_theme!(makie_theme)
set_theme!(Theme(fontsize = 16))   # ← increase this number for bigger text everywhere


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
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)")#, title="Ice thickness (m)")

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
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)")#, title="Bedrock elevation (m a.s.l.)")

    hm = heatmap!(ax, x, y, Z; colormap=:heat, colorrange=(vmin, vmax))
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

# --- Hillshade helper (simple finite differences) ---
function hillshade(Z::AbstractMatrix, x::AbstractVector, y::AbstractVector;
                   azimuth_deg::Float64 = 315.0, altitude_deg::Float64 = 60.0) # reduce altitude_deg for stronger shadows

    # spacing (assumes regular grid)
    dx = abs(x[2] - x[1])
    dy = abs(y[2] - y[1])

    # gradients
    dzdx = similar(Z)
    dzdy = similar(Z)

    # finite differences (central in interior, forward/backward at edges)
    @inbounds for j in axes(Z,2), i in axes(Z,1)
        if !isfinite(Z[i,j])
            dzdx[i,j] = NaN
            dzdy[i,j] = NaN
            continue
        end

        i0 = i == first(axes(Z,1)) ? i : i-1
        i1 = i == last(axes(Z,1))  ? i : i+1
        j0 = j == first(axes(Z,2)) ? j : j-1
        j1 = j == last(axes(Z,2))  ? j : j+1

        z_im1 = Z[i0, j]
        z_ip1 = Z[i1, j]
        z_jm1 = Z[i, j0]
        z_jp1 = Z[i, j1]

        dzdx[i,j] = (z_ip1 - z_im1) / ( (i1 - i0) * dx )
        dzdy[i,j] = (z_jp1 - z_jm1) / ( (j1 - j0) * dy )
    end

    # illumination
    az = deg2rad(azimuth_deg)
    alt = deg2rad(altitude_deg)

    # slope/aspect
    slope = atan.(sqrt.(dzdx.^2 .+ dzdy.^2))
    aspect = atan.(dzdy, .-dzdx)  # aspect measured clockwise from north-ish convention

    hs = sin(alt) .* cos.(slope) .+ cos(alt) .* sin.(slope) .* cos.(az .- aspect)

    # clip and normalize to [0,1]
    hs = clamp.(hs, 0.0, 1.0)
    return hs
end


"""
plot_surface_hillshade(surface; outline_raster=nothing, savepath=nothing)

- `surface` :: Raster (e.g. surface_2024_june.tif)
- `outline_raster` :: Raster where >0 defines glacier (e.g. thickness raster)
Produces a hillshade of the surface DEM + glacier outline.
"""
function plot_surface_hillshade(surface::Raster;
                                outline_raster::Union{Raster,Nothing}=nothing,
                                savepath::Union{String,Nothing}=nothing)

    x, y, Zs = get_axes_and_matrix(surface)

    # --- Hillshade over the full DEM (no clipping) ---
    hs = hillshade(Zs, x, y)

    fig = Figure(size=(800, 600))
    ax  = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)")

    hm = heatmap!(ax, x, y, hs; colormap=:grays, colorrange=(0, 1))

    # --main lake mask overlay (hard-coded shapefile) ---
    lake_mask_path = "/scratch-3/cogier/data/BonnePierre_input/Lac_plein.shp"
    if isfile(lake_mask_path)
        _overlay_depressions!(ax, lake_mask_path;
            linecolor = (:green, 0.9),
            fillcolor = (:green, 0.25),
            lw        = 1.0
        )
    else
        @warn "Lake mask shapefile not found: $lake_mask_path"
    end

    # --- 20 m elevation contours on top of hillshade ---
    vmin, vmax = finite_minmax(Zs)
    step = 20.0
    lmin = floor(vmin/step) * step
    lmax = ceil(vmax/step)  * step
    levels = collect(lmin:step:lmax)

    contour!(ax, x, y, Zs;
        levels    = levels,
        linewidth = 0.6,
        color     = :black
    )

    # --- Glacier outline (from thickness raster) ---
    if outline_raster !== nothing
        _, _, Z_out = get_axes_and_matrix(outline_raster)
        mask_outline = Int.(Z_out .> 0)
        contour!(ax, x, y, mask_outline; levels=[0.5], color=:blue, linewidth=1.5)
    end

    # Optional colorbar (keep/remove as you prefer)
    #cb = Colorbar(fig[1, 2], hm; ticks=([0, 0.5, 1.0], ["0", "0.5", "1"]), label="Hillshade")
    #cb.height[] = 350

    if savepath !== nothing
        save(savepath, fig; px_per_unit=4)
        println("✅ Saved surface hillshade plot to: $savepath")
    end

    return fig
end


function plot_uncertainty_bed(
    r1::Raster;
    r2::Union{Raster,Nothing}=nothing,
    title::String = "",
    subtitle1::String = "",
    subtitle2::String = "",
    savepath::String
)

    # Extract axes + matrices
    x, y, Z1 = get_axes_and_matrix(r1)
    finite_vals1 = Z1[isfinite.(Z1)]

    if isempty(finite_vals1)
        @warn "No valid values for r1."
        return nothing
    end

    # ============ CASE 1: Only ONE raster (single panel) ============
    if isnothing(r2)
        vmin = 0
        vmax = maximum(finite_vals1)
        vmax = ceil(vmax)  # nicer upper bound

        cmap = cgrad(:reds)  # red-toned colormap
        ticks_vals = collect(range(vmin, vmax; length=5))
        ticks_labels = string.(Int.(round.(ticks_vals)))

        fig = Figure(size=(800, 600))
        ax = Axis(fig[1, 1];
            aspect = DataAspect(),
            xlabel = "X (m)",
            ylabel = "Y (m)",
            title = subtitle1
        )

        hm = heatmap!(ax, x, y, Z1;
            colormap = cmap,
            colorrange = (vmin, vmax)
        )

        # --- lines of uncertainty ---
        step   = 5
        lmin   = floor(vmin/step) * step
        lmax   = ceil(vmax/step) * step
        levels = collect(lmin:step:lmax)
        contour!(ax, x, y, Z1; levels=levels, linewidth=0.5, color=:black)


        Colorbar(fig[1, 2], hm;
            label = "Standard deviation (m)",
            height = 350,
            ticks = (ticks_vals, ticks_labels)
        )

        if title != ""
            Label(fig[0, :], title; fontsize=12, font=:bold)
        end

        save(savepath, fig; px_per_unit=4)
        println("✅ Saved single uncertainty plot to: $savepath")

        return fig
    end

    # ============ CASE 2: TWO rasters (your original logic) ============
    _, _, Z2 = get_axes_and_matrix(r2)
    finite_vals = vcat(finite_vals1, Z2[isfinite.(Z2)])

    vmin_data = minimum(finite_vals)
    vmax_data = maximum(finite_vals)
    vmax_abs = ceil(max(abs(vmin_data), abs(vmax_data)))
    colorrange = (-vmax_abs, vmax_abs)
    cmap = cgrad(:balance, scale=colorrange)

    nticks = 5
    ticks_vals = collect(round.(range(vmin_data, vmax_data; length=nticks)))
    if 0 ∉ ticks_vals
        push!(ticks_vals, 0); sort!(ticks_vals)
    end
    ticks_labels = string.(Int.(ticks_vals))

    fig = Figure(size=(900, 450))

    ax1 = Axis(fig[1, 1];
        aspect=DataAspect(),
        xlabel="X (m)",
        ylabel="Y (m)",
        title=subtitle1
    )
    hm1 = heatmap!(ax1, x, y, Z1;
        colormap=cmap, colorrange=colorrange
    )

    ax2 = Axis(fig[1, 2];
        aspect=DataAspect(),
        xlabel="X (m)",
        ylabel="Y (m)",
        title=subtitle2
    )
    heatmap!(ax2, x, y, Z2;
        colormap=cmap, colorrange=colorrange
    )

    Colorbar(fig[1, 3], hm1;
        label="Uncertainty (m)",
        height=300,
        ticks=(ticks_vals, ticks_labels)
    )

    Label(fig[0, :], title; fontsize=12, font=:bold)

    save(savepath, fig; px_per_unit=4)
    println("✅ Saved TWO-raster uncertainty plot to: $savepath")

    return fig
end

function _overlay_depressions!(ax, shp_path::String; linecolor=:black, fillcolor=(:black, 0.2), lw=1.0)
    table = Shapefile.Table(shp_path)

    for geom in table.geometry
        # Each geom can be Polygon or MultiPolygon
        for poly in GeoInterface.getgeom(geom)
            for ring in GeoInterface.getgeom(poly)  # outer + inner rings
                coords = GeoInterface.coordinates(ring)
                xs = first.(coords)
                ys = last.(coords)
                poly!(ax, xs, ys; color=fillcolor, strokecolor=linecolor, strokewidth=lw)
            end
        end
    end
end


function plot_lake_depth(
    lakes::Union{Raster,Nothing},                # ← allow `nothing`
    thickness::Raster,
    analysis::LakeAnalysisResult,
    phi::Union{Raster, Nothing},
    savepath::String;
    min_depth::Float64 = 2.0,
    show_all_lakes::Bool = false,
    area::Union{Raster, Nothing} = nothing,
    area_threshold::Float64 = 1e5,               
    depressions_path::Union{Nothing,String} = nothing,
    stochastic::Bool = true                      # NEW: true = shading (MC), false = contour outline (main)
)
    glacier_mask = thickness .> 0
    mask_array   = collect(Bool.(glacier_mask))

    # pick axes from lakes if present, otherwise from thickness
    base = isnothing(lakes) ? thickness : lakes
    x, y, _ = get_axes_and_matrix(base)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1];
    aspect = DataAspect(),
    xlabel = "X (m)",
    ylabel = "Y (m)",
    #title  = "Water pocket depth (m > $(min_depth))",

    titlesize = 18,
    xlabelsize = 16,
    ylabelsize = 16,
    xticklabelsize = 16,   # size for subplots
    yticklabelsize = 16,
    )

    
    # --- Zoom-out / fixed view window ---
    x_left = minimum(x)         
    xlims!(ax, x_left, 964000.0)
    ylims!(ax, 6431000.0, 6432000.0)

    # ↓↓↓ reduce clutter
    ax.xticks = WilkinsonTicks(4)  # ALSO for subplots
    ax.yticks = WilkinsonTicks(4)

    # --- Optional depressions overlay ---
    if depressions_path !== nothing
        _overlay_depressions!(ax, depressions_path;
            linecolor = (:green, 0.8),
            fillcolor = (:green, 0.2),
            lw        = 1.0
        )
    end

    # --- Upslope area: either heatmap (stochastic) or contour (deterministic) ---
    hm_area = nothing

    if area !== nothing
        _, _, Z_area_raw = get_axes_and_matrix(area)

        if stochastic
            # ---- STOCHASTIC CASE: shaded upslope area heatmap ----
            Z_area = copy(Z_area_raw)

            # Glacier mask and positive-only
            Z_area[.!mask_array] .= NaN
            Z_area[Z_area .<= 0] .= NaN

            # Clip at max threshold for readability
            Z_area_clipped = clamp.(Z_area, 0, area_threshold)

            finite_pos = Z_area_clipped[isfinite.(Z_area_clipped) .& (Z_area_clipped .> 0)]
            if !isempty(finite_pos)
                area_min = minimum(finite_pos)
                area_max = area_threshold

                # Darker blue colormap: white → light blue → mid blue → dark blue
                cmap_area = cgrad([:white, "#8EC1FF", "#1F78B4", "#08306B"])


                hm_area = heatmap!(ax, x, y, Z_area_clipped;
                    colormap   = cmap_area,
                    colorrange = (area_min, area_max)
                )
            end
        else
            # ---- DETERMINISTIC CASE: simple contour outline above threshold ----
            Z_area = copy(Z_area_raw)
            Z_area[.!mask_array] .= 0.0

            area_mask = (Z_area .> area_threshold) .& mask_array
            area_int  = Int.(area_mask)

            contour!(ax, x, y, area_int;
                levels    = [0.5],
                color     = (:darkblue, 0.8),
                linewidth = 1.0
            )
        end
    end

    # --- Lake depth heatmap (viridis) ---
    hm_lake = nothing
    if lakes !== nothing
        _, _, Z = get_axes_and_matrix(lakes)
        Z_lake = copy(Z)
        Z_lake[Z_lake .< min_depth] .= NaN
        Z_lake[.!mask_array]        .= NaN
        Z_lake[Z_lake .== 0]        .= NaN

        vmin, vmax = finite_minmax(Z_lake)
        vmin = min_depth
        vmax = 25.0   # ← FIXED maximum for better comparability in SUBPLOTS (25 correspond to determisitc and also bedrock unc. max)

        # Draw lakes on top of upslope area
        hm_lake = heatmap!(ax, x, y, Z_lake;
            colormap   = :viridis,
            colorrange = (vmin, vmax)
        )
    end

    # --- Largest lake outline (still optional) ---
    if any(analysis.LargestLake.mask)
        labeled = Int.(analysis.LargestLake.mask)
        # contour!(ax, x, y, labeled; levels=[0.5], color=:red, linewidth=1.5)
    end

    # --- All lake masks (optional) ---
    if show_all_lakes
        for (_, mask) in analysis.lake_masks
            if any(mask)
                labeled = Int.(mask)
                contour!(ax, x, y, labeled; levels=[0.5], color=:purple, linewidth=1)
            end
        end
    end

    # --- Outlines of lakes with volume > 1000 m³ ---
    big_volume_threshold = 1000.0

    # Find indices of lakes above threshold
    big_inds = findall(v -> v > big_volume_threshold, analysis.stats.volume)

    if !isempty(big_inds)
        for label in big_inds
            mask = analysis.lake_masks[label]

            if any(mask)
                labeled = Int.(mask)
                contour!(
                    ax, x, y, labeled;
                    levels    = [0.5],
                    color     = (:red, 0.9),
                    linewidth = 2,
                )
            end
        end
    end

    # --- Hydraulic head contours ---
    if phi !== nothing
        _, _, Z_phi = get_axes_and_matrix(phi)
        Z_phi[.!mask_array] .= NaN
        vmin_phi = floor(minimum(Z_phi[isfinite.(Z_phi)]), digits=0)
        vmax_phi = ceil(maximum(Z_phi[isfinite.(Z_phi)]),  digits=0)
        levels   = collect(vmin_phi:10:vmax_phi)
        contour!(ax, x, y, Z_phi;
            levels    = levels,
            linewidth = 0.8,
            color     = :black
        )
    end

    # --- Volume annotations ---
    total_vol = round(Int, sum(analysis.stats.volume))
    max_vol   = round(Int, analysis.LargestLake.volume)
    text!(
        ax, x[1], y[end],
        text     = "Total volume: $(total_vol) m³\nLargest water pocket: $(max_vol) m³",
        fontsize = 10,
        color    = :black
    )

    # --- Colorbar for lake depth (right) ---
    if hm_lake !== nothing
        cr = hm_lake.attributes.colorrange[]
        nticks      = 4
        ticks_vals  = range(cr[1], cr[2], length=nticks)
        ticks_labels = string.(Int.(round.(ticks_vals)))
        Colorbar(fig[1, 2], hm_lake;
            ticks  = (ticks_vals, ticks_labels),
            label  = "Water pockets height (m)",
            height = 200 ## PUT 300 for SINGLE PLOT
        )
    end

    # --- Colorbar for upslope area (only in stochastic mode) ---
    if stochastic && (hm_area !== nothing)
        ar = hm_area.attributes.colorrange[]
        nticks      = 2
        ticks_vals  = range(ar[1], ar[2], length=nticks)
        ticks_labels = string.(Int.(round.(ticks_vals)))
        Colorbar(fig[1, 3], hm_area;
            ticks  = (ticks_vals, ticks_labels),
            label  = "Upslope catchment area (m²)",
            height = 200,  ## PUT 300 for SINGLE PLOT
            labelsize      = 15,
            ticklabelsize  = 14
        )
    end

      # --- GPR water-pick evidence points (only in stochastic plots) ---
    if stochastic
        picks_path = "/scratch-3/cogier/data/BonnePierre_input/Water_picks_gpr_evidence.csv"  # hard coded
        if isfile(picks_path)
            picks_df = CSV.read(picks_path, DataFrame)
            x_picks = Float64.(picks_df.xcoord)
            y_picks = Float64.(picks_df.ycoord)

            scatter!(ax, x_picks, y_picks;
                     color = (:black, 0.4) ,
                     marker = :circle,
                     markersize = 3,
                     label = "GPR water evidence")
        else
            @warn "skipping GPR evidence points."
        end
    end

    # --- Legend entries ---
    lines!(ax, [NaN], [NaN]; color = :black, linewidth = 0.8,
           label = "Hydraulic head (10 m intervals)")
    lines!(ax, [NaN], [NaN]; color = :red, linewidth = 1.0,
           label = "Water pocket outlines")

    if area !== nothing
        if stochastic
            lines!(ax, [NaN], [NaN]; color = :dodgerblue, linewidth = 1.0,
                   label = "Upslope area shading")
        else
            lines!(ax, [NaN], [NaN]; color = :darkblue, linewidth = 1.0,
                   label = "Upslope area > $(Int(area_threshold)) m²")
        end
    end

    Legend(fig, ax;
        tellwidth    = false,
        tellheight   = false,
        halign       = :left,
        valign       = :top,
        framevisible = false
    )

    # --- Glacier outlines (Oct thickness + current mask) ---
    outline_thk = Raster("/scratch-3/cogier/data/BonnePierre_input/WWFS_input/ice_thickness_2024_oct.tif")
    x_out, y_out, Z_out = get_axes_and_matrix(outline_thk)
    mask_outline = Int.(Z_out .> 0)
    contour!(ax, x_out, y_out, mask_outline; levels=[0.5], color=:black, linewidth=1)
    contour!(ax, x_out, y_out, mask_array;   levels=[0.5], color=:black, linewidth=1)

    save(savepath, fig; px_per_unit = 4)
    println("✅ Saved lake depth plot with upslope area to: $savepath")

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
    lines!(ax, dist, z_bed + h_lake,  label="lake_surface (depth)")

    axislegend(ax, position=:rb, framevisible=false)

    # Force y-axis limits
    y_min = minimum(z_bed)  - 15 # m
    y_max = maximum(z_surf) + 5  # m
    #ylims!(ax, y_min, y_max)

    vlines!(ax, [0, dist[end]]; color=:gray, linestyle=:dash, linewidth=1)
    # DOUBLE CHECK A AND B 
    text!(ax, dist[end]-5, 0, text="A (upstream)", align=(:right, :bottom), space=:data)
    text!(ax, 5, 0, text="B (downstream)", align=(:right, :bottom), space=:data)

    save(fname, fig; px_per_unit=3)
    println("✅ saved transect profile to: $fname")
    return fig
end

## plot stochastic ensemble of φ profiles

function plot_transect_spaghetti(
    aggr_spag,                      # e.g. aggr2
    aggr_mean;                      # e.g. aggr5 (reference)
    case_label_spag::AbstractString = "spaghetti case",
    case_label_mean::AbstractString = "aggr5 (mean)",
    plot_phi::Bool = true,
    plot_lake::Bool = false,
    savepath::AbstractString
)
    # --- Basic geometry / static profiles (same for all cases) ---
    @assert hasproperty(aggr_spag, :dist_profile)   "aggr_spag has no dist_profile field"
    @assert hasproperty(aggr_spag, :z_bed_profile)  "aggr_spag has no z_bed_profile field"
    @assert hasproperty(aggr_spag, :z_surf_profile) "aggr_spag has no z_surf_profile field"

    dist           = Float64.(aggr_spag.dist_profile)
    z_bed_profile  = Float64.(aggr_spag.z_bed_profile)
    z_surf_profile = Float64.(aggr_spag.z_surf_profile)
    n_pts          = length(dist)

    # ===================== φ PROFILES =====================
    phi_mat_spag      = nothing   # all realizations for spaghetti case
    phi_spag_mean     = nothing   # mean of aggr_spag ensemble
    phi_ref_mean      = nothing   # mean of aggr_mean ensemble

    if plot_phi
        @assert hasproperty(aggr_spag, :phi_profiles)  "aggr_spag has no phi_profiles field"
        @assert hasproperty(aggr_mean, :phi_profiles)  "aggr_mean has no phi_profiles field"

        phi_profiles_spag = aggr_spag.phi_profiles
        phi_profiles_ref  = aggr_mean.phi_profiles

        # --- spaghetti matrix for aggr_spag ---
        n_real_spag = length(phi_profiles_spag)
        @assert n_real_spag > 0 "phi_profiles in aggr_spag is empty"

        phi_mat_spag = Array{Float64}(undef, n_real_spag, n_pts)
        for i in 1:n_real_spag
            phi_mat_spag[i, :] = Float64.(phi_profiles_spag[i])
        end

        # --- mean φ of aggr_spag (case we compare) ---
        phi_spag_mean = zeros(Float64, n_pts)
        for j in 1:n_pts
            s = 0.0
            for i in 1:n_real_spag
                s += phi_profiles_spag[i][j]
            end
            phi_spag_mean[j] = s / n_real_spag
        end

        # --- mean φ of reference case (aggr_mean, e.g. aggr5) ---
        n_real_ref = length(phi_profiles_ref)
        @assert n_real_ref > 0 "phi_profiles in aggr_mean is empty"

        phi_ref_mean = zeros(Float64, n_pts)
        for j in 1:n_pts
            s = 0.0
            for i in 1:n_real_ref
                s += phi_profiles_ref[i][j]
            end
            phi_ref_mean[j] = s / n_real_ref
        end
    end

    # ================== LAKE FREE-SURFACE ==================
    lake_mat_spag      = nothing
    lake_spag_mean     = nothing
    lake_ref_mean      = nothing

    if plot_lake
        @assert hasproperty(aggr_spag, :lake_profiles) "aggr_spag has no lake_profiles field"
        @assert hasproperty(aggr_mean, :lake_profiles) "aggr_mean has no lake_profiles field"

        lake_profiles_spag = aggr_spag.lake_profiles
        lake_profiles_ref  = aggr_mean.lake_profiles

        # --- spaghetti matrix for aggr_spag ---
        n_real_spag = length(lake_profiles_spag)
        @assert n_real_spag > 0 "lake_profiles in aggr_spag is empty"

        lake_mat_spag = Array{Float64}(undef, n_real_spag, n_pts)
        for i in 1:n_real_spag
            lake_mat_spag[i, :] = Float64.(lake_profiles_spag[i])
        end

        # --- mean lake depth of aggr_spag (case we compare) ---
        lake_spag_mean = zeros(Float64, n_pts)
        for j in 1:n_pts
            s = 0.0
            for i in 1:n_real_spag
                s += lake_profiles_spag[i][j]
            end
            lake_spag_mean[j] = s / n_real_spag
        end

        # --- mean lake depth of reference case (aggr_mean) ---
        n_real_ref = length(lake_profiles_ref)
        @assert n_real_ref > 0 "lake_profiles in aggr_mean is empty"

        lake_ref_mean = zeros(Float64, n_pts)
        for j in 1:n_pts
            s = 0.0
            for i in 1:n_real_ref
                s += lake_profiles_ref[i][j]
            end
            lake_ref_mean[j] = s / n_real_ref
        end
    end

    # ======================= FIGURE ========================
    fig = Figure(size = (800, 500))
    ax  = Axis(fig[1, 1];
        xlabel = "Distance along transect (m)",
        ylabel = "Elevation / hydraulic head (m a.s.l.)",
        title  = "Transect A→B – $(case_label_spag) vs $(case_label_mean)"
    )

    # --- φ spaghetti for aggr_spag (grey) ---
    if plot_phi && (phi_mat_spag !== nothing)
        n_real_spag = size(phi_mat_spag, 1)
        for i in 1:n_real_spag
            lines!(ax, dist, phi_mat_spag[i, :];
                   color = (:gray, 0.25), linewidth = 0.7)
        end
    end

    # --- lake spaghetti for aggr_spag (semi-transparent blue, plotted as bed + depth) ---
    if plot_lake && (lake_mat_spag !== nothing)
        n_real_spag = size(lake_mat_spag, 1)
        for i in 1:n_real_spag
            lines!(ax, dist, z_bed_profile .+ lake_mat_spag[i, :];
                   color = (:dodgerblue, 0.15), linewidth = 0.7)
        end
    end

    # --- bed & surface for reference ---
    lines!(ax, dist, z_bed_profile;
           color = :black, linestyle = :dot,  linewidth = 0.8, label = "bed")
    lines!(ax, dist, z_surf_profile;
           color = :black, linestyle = :dash, linewidth = 0.8, label = "surface")

    # --- means of the CASE we compare (aggr_spag) in color ---
    if plot_phi && (phi_spag_mean !== nothing)
        lines!(ax, dist, phi_spag_mean;
               color = :dodgerblue, linewidth = 2.0,
               label = "φ mean ($(case_label_spag))")
    end

    if plot_lake && (lake_spag_mean !== nothing)
        lines!(ax, dist, z_bed_profile .+ lake_spag_mean;
               color = :dodgerblue, linestyle = :solid, linewidth = 2.0,
               label = "lake_free_surface mean ($(case_label_spag))")
    end

    # --- means of the REFERENCE case (aggr_mean) in black ---
    if plot_phi && (phi_ref_mean !== nothing)
        lines!(ax, dist, phi_ref_mean;
               color = :black, linewidth = 2.5,
               label = "φ mean ($(case_label_mean))")
    end

    if plot_lake && (lake_ref_mean !== nothing)
        lines!(ax, dist, z_bed_profile .+ lake_ref_mean;
               color = :black, linestyle = :solid, linewidth = 2.5,
               label = "lake_free_surface mean ($(case_label_mean))")
    end

    # --- y-limits ---
    y_min = minimum(z_bed_profile)
    y_max = maximum(z_surf_profile)

    if plot_phi && (phi_mat_spag !== nothing)
        y_min = min(y_min, minimum(phi_mat_spag))
        y_max = max(y_max, maximum(phi_mat_spag))
    end
    if plot_phi && (phi_spag_mean !== nothing)
        y_min = min(y_min, minimum(phi_spag_mean))
        y_max = max(y_max, maximum(phi_spag_mean))
    end
    if plot_phi && (phi_ref_mean !== nothing)
        y_min = min(y_min, minimum(phi_ref_mean))
        y_max = max(y_max, maximum(phi_ref_mean))
    end

    if plot_lake && (lake_mat_spag !== nothing)
        lake_elev_min = minimum(z_bed_profile .+ minimum(lake_mat_spag, dims=1))
        lake_elev_max = maximum(z_bed_profile .+ maximum(lake_mat_spag, dims=1))
        y_min = min(y_min, lake_elev_min)
        y_max = max(y_max, lake_elev_max)
    end
    if plot_lake && (lake_spag_mean !== nothing)
        lake_elev_mean_spag = z_bed_profile .+ lake_spag_mean
        y_min = min(y_min, minimum(lake_elev_mean_spag))
        y_max = max(y_max, maximum(lake_elev_mean_spag))
    end
    if plot_lake && (lake_ref_mean !== nothing)
        lake_elev_mean_ref = z_bed_profile .+ lake_ref_mean
        y_min = min(y_min, minimum(lake_elev_mean_ref))
        y_max = max(y_max, maximum(lake_elev_mean_ref))
    end

    ylims!(ax, y_min - 10, y_max + 5)

    axislegend(ax, position = :rb, framevisible = false)

    save(savepath, fig; px_per_unit = 3)
    println("✅ Saved transect spaghetti plot to: ", savepath)

    return fig
end

"""
    plot_lake_volume_boxplot(
        lake_vols, largest_lake_vols, labels, savepath;
        plot_largest=true
    )

Plots a boxplot for stochastic lake volumes.

Arguments
---------
- `lake_vols` :: Vector of vectors  
      Each element contains the N realizations of total lake volume for one case.
- `largest_lake_vols` :: Vector of vectors  
      Same but for largest-lake volume.
- `labels` :: Vector of String  
      Labels below each category.
- `savepath` :: String  
      Where to save the figure.

Keyword
-------
- `plot_largest=true`  
      If `false`, only total lake volume is plotted.
"""

function plot_lake_volume_boxplot(
    lake_vols::Vector{<:Vector},
    largest_lake_vols::Vector{<:Vector},
    labels::Vector{String},
    savepath::String;
    plot_largest::Bool = true,
    logscale::Bool = false
)

    n = length(labels)
    centers = 1:n

    fig = Figure(size = (420, 420))

    ax = Axis(fig[1, 1];
        ylabel = "Total water pockets volume (m³)",
        xlabel = "Field perturbed",
        xticks  = (centers, labels),
        yscale  = logscale ? log10 : identity
    )

    if plot_largest
        # --- Two-category plot: total & largest ---
        offset = 0.12
        w = 0.3

        pos_tot = centers .- offset
        pos_lrg = centers .+ offset

        x_tot = vcat([fill(pos_tot[i], length(v)) for (i,v) in enumerate(lake_vols)]...)
        y_tot = vcat(lake_vols...)

        x_lrg = vcat([fill(pos_lrg[i], length(v)) for (i,v) in enumerate(largest_lake_vols)]...)
        y_lrg = vcat(largest_lake_vols...)

        boxplot!(ax, x_tot, y_tot; color=:dodgerblue, width=w)
        boxplot!(ax, x_lrg, y_lrg; color=:orange, width=w)

        scatter!(ax, pos_tot, [mean(v) for v in lake_vols];
                 color=:black, marker=:cross, markersize=9)
        scatter!(ax, pos_lrg, [mean(v) for v in largest_lake_vols];
                 color=:black, marker=:cross, markersize=9)

        lines!(ax, [NaN], [NaN]; color=:dodgerblue, label="Total pockets volume")
        lines!(ax, [NaN], [NaN]; color=:orange,     label="Largest pocket volume")

    else
        # --- Single-category plot: total only ---
        w = 0.5
        pos_tot = centers

        x_tot = vcat([fill(pos_tot[i], length(v)) for (i,v) in enumerate(lake_vols)]...)
        y_tot = vcat(lake_vols...)

        boxplot!(ax, x_tot, y_tot; color=:dodgerblue, width=w)

        scatter!(ax, pos_tot, [mean(v) for v in lake_vols];
                 color=:black, marker=:cross, markersize=9)

        lines!(ax, [NaN], [NaN]; color=:dodgerblue, label="Total pockets volume")
    end

    save(savepath, fig)
    println("✅ Saved boxplot to: $savepath")
    return fig
end
