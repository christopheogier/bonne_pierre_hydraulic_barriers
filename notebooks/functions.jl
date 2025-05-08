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

# structure at the top level to store attributes
"""
`LakeAnalysis` is a structure to store the results of a lake analysis. It holds the labeled image of detected lakes and 
the statistical measurements for each lake, including their volume.

# Fields:
- `labels::Matrix{Int}`: A matrix where each pixel is assigned a label corresponding to a detected lake.
- `stats::DataFrame`: A data frame containing combined measurements, such as area, for each labeled lake.
"""
struct LakeAnalysis
    labels::Matrix{Int}    # Labeled image of lakes
    stats::DataFrame       # Combined measurements and volumes
end

"""
`analyze_lakes` takes two input rasters (lakes and ice mask), processes them to detect individual lakes, 
calculates their volumes, and returns a `LakeAnalysis` structure.

# Arguments:
- `lakes::Matrix{Float64}`: A raster representing the lake areas where non-zero values indicate the presence of lakes.
- `ice_mask::Matrix{Bool}`: A binary mask representing the ice-covered areas, where `true` indicates ice presence.

# Returns:
- `LakeAnalysis`: A structure containing two fields:
    - `labels`: A matrix of labeled connected components, where each lake has a unique integer label.
    - `stats`: A data frame containing various measurements for each lake, including the computed volume.

# Example:
    result = analyze_lakes(lakes, ice_mask)
    println(result.labels)
    println(result.stats)
"""
function analyze_lakes(lakes, ice_mask)
    # Build binary lake mask
    mask_lakes = (lakes .> 0) .& (ice_mask .> 0)

    # Spatial resolution (assuming square pixels)
    dx = step(dims(lakes)[1])

    # Label connected components (individual lakes)
    labeled_image = Images.label_components(mask_lakes)

    # Analyze the labeled components to get stats (area in pixels, etc.)
    measurements = analyze_components(labeled_image, BasicMeasurement())

    # Compute volumes for each labeled lake
    lake_volumes = Float64[]
    for label in 1:maximum(labeled_image)
        mask_current = labeled_image .== label
        volume = sum(lakes[mask_current]) * dx^2  # Volume in m³
        push!(lake_volumes, volume)
    end

    # Add volumes as a new column to the measurements DataFrame
    measurements.volume = lake_volumes

    return LakeAnalysis(labeled_image, measurements)
end



