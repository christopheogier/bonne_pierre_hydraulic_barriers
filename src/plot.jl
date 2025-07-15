"""
plot.jl

Functions to plot ice thickness and bedrock elevation using Raster.jl and Plots.jl.

Optional GPR points can be added.

"""

using Plots
using Rasters
using DelimitedFiles  # for loading GPR points if needed
using Plots, Statistics
using Contour  # for contour lines
using StatsPlots



ENV["GKSwstype"] = "100"  # Enable headless plotting (no GUI required)



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

"""
    plot_lake_depth(lakes::Raster, savepath::String; largest_mask::Union{Nothing, BitMatrix}=nothing)

Plot a heatmap of lake depths. If `largest_mask` is provided, the outline of the largest lake is overlaid.
"""
function plot_lake_depth(lakes::Raster, savepath::String)
    plt = heatmap(
        lakes;
        title = "Lake depth (m)",
        colorbar_title = "Lake depth (m)",
        color = :blues,
        axis = false,
        ticks = false,
        size = (800, 700),
        aspect_ratio = :equal
    )

    savefig(plt, savepath)
    return plt
end



function plot_hydraulic_head(phi::Raster, savepath::String)
    phi_clean = replace(phi, NaN => missing)
    vals = vec(collect(phi_clean))
    vals = filter(!isnan, skipmissing(vals))

    plt = heatmap(
        phi_clean;
        title="Hydraulic head (m)",
        colorbar_title="Hydraulic head (m)",
        color=:thermal,
        axis=false,
        ticks=false,
        size=(800, 700),
        aspect_ratio=:equal
    )

    if !isempty(vals)
        vmin = floor(minimum(vals), digits=0)
        vmax = ceil(maximum(vals), digits=0)

        if vmin < vmax
            contour!(
                plt,
                phi_clean;
                levels=range(vmin, vmax; step=20),
                linewidth=1.0,
                linecolor=:black,
                label=false
            )
        else
            @warn "Skipping contour: vmin == vmax ($(vmin))"
        end
    else
        @warn "Skipping contour: no valid phi values found"
    end

    savefig(plt, savepath)
    return plt
end

function heatmap_mindepth_vs_smoothing(summaries::DataFrame, run_name::AbstractString, zvar::Symbol, savepath::String)
    # Filter DataFrame for the specified run
    run_df = filter(:run => r -> occursin(run_name, r), summaries)

    # Sort to ensure reshaping works
    sort!(run_df, [:smooth_surface_ice_fraction, :min_depth_m])

    # Extract unique x and y axis values
    smoothing_vals = unique(run_df.smooth_surface_ice_fraction)
    min_depth_vals = unique(run_df.min_depth_m)

    # Check reshaping consistency
    if length(smoothing_vals) * length(min_depth_vals) != nrow(run_df)
        error("Inconsistent data: missing combinations of smoothing and min_depth.")
    end

    # Reshape the Z matrix (zvar values)
    Z = reshape(run_df[!, zvar], (length(min_depth_vals), length(smoothing_vals)))'

    # Plot heatmap
    plt = heatmap(
        min_depth_vals,
        smoothing_vals,
        Z,
        xlabel = "Minimum Lake Depth [m]",
        ylabel = "Surface Smoothing Fraction",
        title = string(zvar) * " – " * run_name,
        colorbar_title = string(zvar),
        c = :blues,
        size = (700, 500),
        dpi = 300
    )

    savefig(plt, savepath)
end

function plot_run_summary(summaries::DataFrame, run_name::String, output_dir::String)
    # Filter data for the selected run
    df = filter(row -> row.run == run_name, summaries)

    # Convert numeric parameters to string for plotting
    df.smooth_str = string.(df.smooth_surface_ice_fraction)
    df.depth_str = string.(df.min_depth_m)
    df.fill_str = string.(df.fill_frac)

    # Set categorical x-axis labels
    df.label = "d=" .* df.depth_str .* ", s=" .* df.smooth_str .* ", f=" .* df.fill_str

    # Sort labels for consistent plotting
    sortperm = sortperm(df.label)
    df = df[sortperm, :]

    # Plot 1: Number of lakes
    bar1 = bar(
        df.label,
        df.n_lakes,
        legend = false,
        title = "Number of Lakes - $(run_name)",
        xlabel = "Parameters (min_depth, smoothing, fill)",
        ylabel = "Number of Lakes",
        xticks = :auto,
        rotation = 45,
        bar_width = 0.6,
        color = :steelblue,
        size = (900, 400),
        dpi = 200
    )
    savefig(bar1, joinpath(output_dir, "$(run_name)_barplot_n_lakes.png"))

    # Plot 2: Total Volume
    bar2 = bar(
        df.label,
        df.total_volume_m3 ./ 1e6,
        legend = false,
        title = "Total Lake Volume - $(run_name)",
        xlabel = "Parameters (min_depth, smoothing, fill)",
        ylabel = "Volume [10⁶ m³]",
        xticks = :auto,
        rotation = 45,
        bar_width = 0.6,
        color = :darkgreen,
        size = (900, 400),
        dpi = 200
    )
    savefig(bar2, joinpath(output_dir, "$(run_name)_barplot_total_volume.png"))

    println("  📊 Saved summary plots for run: ", run_name)
end

