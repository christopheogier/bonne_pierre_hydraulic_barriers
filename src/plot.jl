"""
plot.jl

Functions to plot ice thickness and bedrock elevation using Raster.jl and Plots.jl.

Optional GPR points can be added.
"""

ENV["GKSwstype"] = "100"  # Enable headless plotting (no GUI required)

using Plots
using Rasters
using DelimitedFiles  # for loading GPR points if needed

"""
    plot_ice_thickness(rt::Raster; gpr_points=nothing, savepath=nothing)

Plot ice thickness raster with contour lines and optional GPR points.
Highlights the location of maximum ice thickness.
Saves the figure if `savepath` is provided.

TODO: Do not plot ice thickness of 0m (basically mask the data outside the glacier)
"""
function plot_ice_thickness(rt::Raster; gpr_points=nothing, savepath=nothing, year = "unknown")
    # Set default image size and resolution
    default(size=(1000, 800), dpi=300)

    # Compute min/max for contours and colorbar
    vmin = floor(minimum(rt), digits=0)
    vmax = ceil(maximum(rt), digits=0)

    # Find max thickness location
    max_val = maximum(rt)
    max_idx = argmax(rt)
    xdim, ydim = dims(rt)
    x_max = xdim[max_idx[1]]
    y_max = ydim[max_idx[2]]

    # Plot heatmap
    heatmap(
        rt;
       
        xlabel = "X (m)",
        ylabel = "Y (m)",
        color = :viridis,
        aspect_ratio = :equal,
        legend = :topleft

    )

    # Overlay contour lines
    contour!(
        rt;
        levels = range(vmin, stop=vmax, step=10),
        linewidth = 1.0,
        linecolor = :black,
        label = true,
    )

    # Mark max thickness location
    scatter!(
        [x_max], [y_max];
        color = :red,
        markersize = 4,
        label = "Max: $(round(max_val, digits=2)) m"
    )

    # Overlay GPR points if provided
    if gpr_points !== nothing
        scatter!(
            gpr_points[:, 1],
            gpr_points[:, 2];
            color = :black,
            marker = :circle,
            markersize = 1.5,
            markerstrokecolor = :black,
            label = "GPR points", 
            title = "Ice thickness - $year",
            colorbar_title = "Ice thickness (m)"
        )
    end

    # Save figure
    if savepath !== nothing
        savefig(savepath)
    end
end


"""
    plot_bedrock(rt::Raster; savepath=nothing)

Plot bedrock elevation raster with contour lines.
Save the figure if savepath is provided.
"""
function plot_bedrock(rt::Raster; savepath=nothing)
    heatmap(
        rt;
        title = "Bedrock Elevation - Bonne Pierre",
        xlabel = "X (m)",
        ylabel = "Y (m)",
        color = :thermal,
        aspect_ratio = :equal,
    )
    contour!(
        rt;
        levels = range(2400, stop=3200, step=20),
        linewidth = 1,
        linecolor = :black,
        label = false,
    )
    if savepath !== nothing
        savefig(savepath)
    end
end

