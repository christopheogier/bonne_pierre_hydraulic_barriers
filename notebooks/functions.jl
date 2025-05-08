"""
    compute_ice_thickness(surface::Raster, bed::Raster) -> Raster

Computes ice thickness as the difference between a glacier surface raster and a bedrock raster.
The bed raster is resampled to match the surface raster grid using bilinear interpolation.
Negative thickness values are set to 0, and invalid surface elevations (e.g., < 0) are masked as NaN.

# Arguments
- `surface::Raster`: Glacier surface elevation raster (e.g. from LiDAR or photogrammetry).
- `bed::Raster`: Bedrock elevation raster (may have a different resolution or grid).

# Returns
- `Raster`: Ice thickness raster (same extent and resolution as `surface`).
"""
function compute_ice_thickness(surface::Raster, bed::Raster)
    # Clean, crop surface: set invalid elevation values to NaN
    surface = crop(surface; to=bed)

    surface[surface .< 0] .= NaN

    # Resample bed to match surface resolution and grid
    bed_resamp = resample(bed; to=surface, method=:bilinear)

    # Compute ice thickness
    ice_thickness = surface .- bed_resamp
    ice_thickness[ice_thickness .< 0] .= 0

    return ice_thickness
end


"""
    plot_ice_thickness(ice_thickness::Raster; plot_title::String="Ice Thickness", step::Real=10)

Plot a heatmap of glacier ice thickness with contours and highlight the location of maximum thickness.

# Arguments
- `ice_thickness::Raster`: A 2D raster of ice thickness values (in meters).
- `plot_title::String="Ice Thickness"`: Title for the plot.
- `step::Real=10`: Contour interval in meters.

# Behavior
- Negative thickness values are set to zero.
- The maximum ice thickness is identified and marked on the plot with a red dot and label.
- Contours are drawn from the minimum to the maximum ice thickness with a given step.
using Rasters, Plots
"""

function plot_ice_thickness(
    ice_thickness::Raster,
    plot_title::String = "Ice Thickness",
    step::Real = 10)

    # Max value and coordinates
    #max_value = maximum(ice_thickness[.!isnan.(ice_thickness)]) # filter
    #max_index = argmax(ice_thickness[.!isnan.(ice_thickness)])  # (row, col) tuple
    #x, y = dims(ice_thickness)
    #x_max, y_max = x[max_index[2]], y[max_index[1]]  # Note: column, row order

    # Contour limits
    vmin = floor(minimum(ice_thickness), digits=0)
    vmax = ceil(maximum(ice_thickness), digits=0)

    # Plot
    Plots.heatmap(
        ice_thickness,
        title = plot_title,
        colorbar_title = "Ice thickness [m]",
        xlabel = "Easting",
        ylabel = "Northing",
        c = :blues,
        aspect_ratio = :equal
    )

    # Add max point
    #scatter!([x_max], [y_max],color = :red, markersize = 4, label = "Max: $(round(max_value, digits=2)) m")

    # Add contours
    contour!(
        ice_thickness;
        levels = range(start=vmin, stop=vmax, step=step),
        linewidth = 1.0,
        linecolor = :black,
        label = false
    )
end


