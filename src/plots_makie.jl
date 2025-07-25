#plot_makie.jl

using CairoMakie
include("LakeAnalysis.jl")
using .LakeAnalysis

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
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title="Bedrock elevation")

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



function plot_lake_depth(
    lakes::Raster,
    thickness::Raster,
    analysis::LakeAnalysisResult,
    phi::Union{Raster, Nothing},
    savepath::String;
    min_depth::Float64 = 2.0,
    show_all_lakes::Bool = false,
    area::Union{Raster, Nothing} = nothing,
    area_threshold::Float64 = 1e5
)
    glacier_mask = thickness .> 0
    x, y, Z = get_axes_and_matrix(lakes)
    mask_array = collect(Bool.(glacier_mask))

    # Apply min_depth and glacier mask
    Z_lake = copy(Z)
    Z_lake[Z_lake .< min_depth] .= NaN
    Z_lake[.!mask_array] .= NaN
    Z_lake[Z_lake .== 0] .= NaN

    
    vmin, vmax = finite_minmax(Z_full)
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
            color = (:darkblue,0.8),
            linewidth = 1.0
        )
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
        halign = :left, valign = :top, fontsize = 10, color = :black
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
        label = "Upslope area > $(Int(max_threshold)) m²"
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


