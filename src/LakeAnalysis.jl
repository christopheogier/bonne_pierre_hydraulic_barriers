module LakeAnalysis

using Images, DataFrames, Statistics, Rasters, ImageComponentAnalysis

export LakeAnalysisResult, LargestLake, analyze_lakes, boxplot_lakes

# Struct to hold data about the largest lake
struct LargestLake
    volume::Float64
    mask::BitMatrix
end

# Main analysis result
struct LakeAnalysisResult
    labels::Matrix{Int}
    stats::DataFrame
    min_depth::Float64
    LargestLake::LargestLake
end

function analyze_lakes(
    lakes::Raster,
    surface::Raster,
    thickness::Raster;
    min_depth::Float64 = 0.0
)
    dx = step(dims(lakes)[1])
    pixel_area = dx^2
    ice_mask = (thickness .> 0) .& .!ismissing.(surface)
    mask_lakes = (lakes .> min_depth) .& ice_mask

    labeled_image = Images.label_components(collect(Bool.(mask_lakes)))
    measurements = analyze_components(labeled_image, BasicMeasurement())

    lake_volumes, lake_areas_m2 = Float64[], Float64[]
    max_volume, max_label = -Inf, -1

    lakes_array = collect(lakes)  # avoid repeated conversion
    for label in 1:maximum(labeled_image)
        mask_current = labeled_image .== label
        volume = sum(lakes_array[mask_current]) * pixel_area
        area = sum(mask_current) * pixel_area
        push!(lake_volumes, volume)
        push!(lake_areas_m2, area)
        if volume > max_volume
            max_volume = volume
            max_label = label
        end
    end

    measurements.volume = lake_volumes
    measurements.area_m2 = lake_areas_m2

    largest = LargestLake(max_volume, labeled_image .== max_label)

    return LakeAnalysisResult(labeled_image, measurements, min_depth, largest)
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
