module LakeAnalysis

using Images, DataFrames, Statistics, Rasters, ImageComponentAnalysis

export LakeAnalysisResult, analyze_lakes, boxplot_lakes

struct LakeAnalysisResult
    labels::Matrix{Int}
    stats::DataFrame
    min_depth::Float64
end

function analyze_lakes(
    lakes::Raster,
    surface::Raster,
    thickness::Raster;
    min_depth::Float64 = 0.0
)
    # Use raster geometry to compute area
    dx = step(dims(lakes)[1])
    pixel_area = dx^2

    # Create ice mask from thickness
    ice_mask = (thickness .> 0) .& .!ismissing.(surface)

    # Mask lake pixels to be only under ice and deeper than min_depth
    mask_lakes = (lakes .> min_depth) .& ice_mask

    # Label connected components
    labeled_image = Images.label_components(collect(Bool.(mask_lakes)))

    # Analyze components
    measurements = analyze_components(labeled_image, BasicMeasurement())

    lake_volumes, lake_areas_m2 = Float64[], Float64[]
    for label in 1:maximum(labeled_image)
        mask_current = labeled_image .== label
        volume = sum(collect(lakes)[mask_current]) * pixel_area
        area = sum(mask_current) * pixel_area
        push!(lake_volumes, volume)
        push!(lake_areas_m2, area)
    end

    measurements.volume = lake_volumes
    measurements.area_m2 = lake_areas_m2

    return LakeAnalysisResult(labeled_image, measurements, min_depth)
end

function boxplot_lakes(lakes::LakeAnalysisResult, col_str::AbstractString)
    col = Symbol(col_str)
    values = lakes.stats[!, col]
    return (
        mean = mean(values),
        median = median(values),
        q1 = quantile(values, 0.25),
        q3 = quantile(values, 0.75),
        iqr = quantile(values, 0.75) - quantile(values, 0.25),
        max = maximum(values),
        count = length(values)
    )
end

end # module
