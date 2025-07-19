using CairoMakie


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
    Z = Float32.(Matrix(rt))  # Ensure Z[y, x]

    # Flip y and Z if y is descending (top to bottom)
    #if y[2] < y[1]
        #y = reverse(y)
        #Z = reverse(Z, dims=1)
    #end

    return x, y, Z
end



function plot_ice_thickness(rt; gpr_points=nothing, savepath=nothing)
    # Mask zero or negative thickness
    rt_masked = copy(rt)
    rt_masked[rt .<= 0] .= NaN

    x, y, Z = get_axes_and_matrix(rt_masked)
    _, _, Z_full = get_axes_and_matrix(rt)

    vmax = ceil(maximum(Z_full))
    max_idx = argmax(Z_full)
    x_max = x[max_idx[1]]
    y_max = y[max_idx[2]]

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="X (m)", ylabel="Y (m)", title="Ice Thickness (m)")

    # Heatmap and contours
    hm = heatmap!(ax, x, y, Z; colormap=Reverse(:ice), colorrange=(0, vmax))
    # ice thickness countour
    contour!(ax, x, y, Z; levels=10:20:vmax, color=:black)
    contour!(ax, x, y, Z_full; levels=0.1:0.1, linewidth=0.5, color=:black)

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
    vmin = floor(minimum(filter(x -> !isnan(x),skipmissing(Z))), digits=0)
    vmax = ceil(maximum(filter(x -> !isnan(x),skipmissing(Z))), digits=0)

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
    cb = Colorbar(fig[1, 2], hm; ticks=(ticks_vals, ticks_labels), label="Elevation (m)")
    # Match colorbar height to axis height
    cb.height[] = 350  # or whatever pixel height fits your layout better


    if savepath !== nothing
        #save("$savepath.pdf", fig)
        save(savepath, fig; px_per_unit=4)
    end

    return fig
end
