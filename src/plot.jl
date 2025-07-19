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
using Images, ImageComponentAnalysis
using DataFrames



ENV["GKSwstype"] = "100"  # Enable headless plotting (no GUI required)



"""
    plot_ice_thickness(rt::Raster; gpr_points=nothing, savepath=nothing)

Plot ice thickness raster with contour lines and optional GPR points.
Highlights the location of maximum ice thickness.
Saves the figure if `savepath` is provided.

TODO: Do not plot ice thickness of 0m (basically mask the data outside the glacier)
"""
function plot_ice_thickness_former(rt::Raster; gpr_points=nothing, savepath=nothing, year = "unknown")
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


function plot_ice_thickness(
    rt::Raster; 
    gpr_points::Union{Matrix{Float64}, Nothing}=nothing,
    savepath::Union{String, Nothing}=nothing,
    year::String="")

    # Extract raster data and coordinates
    data = Matrix(rt)
    xs = coordinates(rt, dims=1)
    ys = coordinates(rt, dims=2)

    # Create figure and axis
    fig = Figure(resolution=(1200, 1000))
    ax = Axis(fig[1, 1]; title="Ice thickness $(year)", aspect=DataAspect())

    # Plot ice thickness as heatmap
    heatmap!(ax, xs, ys, data; colormap=:ice, colorrange=(minimum(data), maximum(data)))

    # Add GPR points if provided
    if gpr_points !== nothing
        scatter!(ax, gpr_points[:, 1], gpr_points[:, 2]; color=:black, markersize=3)
    end

    Colorbar(fig[1, 2], ax; label="Ice thickness (m)")

    if savepath !== nothing
        save(savepath, fig)
    else
        display(fig)
    end

    return fig
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


function plot_lake_depth_with_outlines(lakes::Raster, surface::Raster, thickness::Raster, savepath::String; min_depth::Float64=2.0)
    # Get coordinates and data
    x = coordinates(lakes, 1)
    y = coordinates(lakes, 2)
    z = permutedims(parent(lakes))  # Makie expects (y,x) layout

    # Compute lake outlines where depth > min_depth
    lake_mask = (lakes .> min_depth) .& (thickness .> 0) .& .!ismissing.(surface)
    lake_labels = label_components(collect(Bool.(lake_mask)))

    # Glacier outline
    glacier_mask = (thickness .> 0) .& .!ismissing.(surface)

    fig = Figure(resolution=(800, 700))
    ax = Axis(fig[1,1], aspect=DataAspect(), title="Lake depth (>{min_depth} m)", xlabel="x", ylabel="y")

    # Heatmap of lake depth
    heatmap!(ax, x, y, z; colormap=:blues, colorrange=(0, maximum(z)), interpolate=false)

    # Plot lake outlines
    labeled_array = parent(lake_labels)
    for label in 1:maximum(labeled_array)
        mask = labeled_array .== label
        if count(mask) == 0
            continue
        end
        C = contours(mask; levels=[0.5])
        for c in C
            for level in c
                lines!(ax, x[level[:,1]], y[level[:,2]], color=:red, linewidth=1.5)
            end
        end
    end

    # Plot glacier outline
    glacier_mask_array = collect(Bool.(glacier_mask))
    Cg = contours(glacier_mask_array; levels=[0.5])
    for c in Cg
        for level in c
            lines!(ax, x[level[:,1]], y[level[:,2]], color=:black, linewidth=1.2)
        end
    end

    Colorbar(fig[1,2], ax, label="Lake depth (m)")

    savefig(fig, savepath)
    println("✅ Saved lake depth plot with outlines to: ", savepath)
    return fig
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

function plot_selected_scenarios(summaries::DataFrame, output_path::String)

    # Filter scenarios
    case1 = filter(row -> row.min_depth_m == 0.0 && row.smooth_surface_ice_fraction == 0.0, summaries)
    case2 = filter(row -> row.min_depth_m == 2.0 && row.smooth_surface_ice_fraction == 0.1, summaries)

    case1.label .= "min=0.0 / smooth=0.0"
    case2.label .= "min=2.0 / smooth=0.1"

    df = vcat(case1, case2)
    sort!(df, [:run, :label])

    labels = unique(df.label)
    runs = unique(df.run)

    # Create long-format DataFrame: one row per (run, label, metric)
    long_df = DataFrame(run = String[], label = String[], metric = String[], value = Float64[], n_lakes = Int[])
    metrics = [:total_volume_m3, :largest_single_volume_m3]

    for metric in metrics
        append!(long_df, DataFrame(
            run = df.run,
            label = df.label,
            metric = fill(string(metric), nrow(df)),
            value = df[!, metric],
            n_lakes = df.n_lakes,
        ))
    end

    # Create x-axis label that groups by run and metric
    long_df.group = string.(long_df.run, " / ", long_df.metric)

    # Plot: group by run/metric, color by scenario
    p = groupedbar(
        long_df.group,
        long_df.value,
        group = long_df.label,
        bar_position = :dodge,
        xlabel = "Run / Metric",
        ylabel = "Volume [m³]",
        legend = :topright,
        title = "Lake Volume Summary (Selected Scenarios)",
        size = (1100, 500),
        rotation = 45,
        color = [:dodgerblue :orangered]
    )

    # Add number of lakes as annotations on top of bars
    xpos = 1
    bar_count_per_group = length(labels)
    for row in eachrow(long_df)
        annotate!(xpos, row.value + 0.03 * row.value, text("$(row.n_lakes)", :black, 8, :center))
        xpos += 1
    end

    # Save plot
    outfile = joinpath(output_path, "_lakes_summary.png")
    savefig(p, outfile)
    println("✅ Saved plot: ", outfile)
end




