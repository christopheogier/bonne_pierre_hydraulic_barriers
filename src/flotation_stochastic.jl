### stochastic_flotation.jl ###
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using Rasters, Serialization
using WhereTheWaterFlowsSubglacially
const WWFS = WhereTheWaterFlowsSubglacially

include("functions.jl")

datadir = "/scratch-3/cogier/data/BonnePierre_input"
inputdir = joinpath(datadir, "WWFS_input")
outputdir = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis_flot"

surface = clean_raster(Raster(joinpath(inputdir, "surface_2024_oct_resamp_1m.tif")))
bedrock = clean_raster(Raster(joinpath(inputdir, "bedrock_resamp_1m.tif")))
thickness = clean_raster(Raster(joinpath(inputdir, "ice_thickness_2024_oct.tif")))

surfdem = surface
beddem = bedrock
rmask = thickness .> 0
floatfrac = 1 .* ones(size(surfdem))
source = ones(size(surfdem))

cov_fn = WWFS.GRF.gaussian_kernel
rel_unc = 0.1
corr_lengths = [10.0, 100.0, 1000.0]

for (i, L) in enumerate(corr_lengths)
    println("🔄 Running flotation-only WWFS with L = $L m")

    surfdem_uc = Uncertainty(absuc=0.0, reluc=0.0)
    beddem_uc = Uncertainty(absuc=0.0, reluc=0.0)
    floatfrac_uc = Uncertainty(absuc=0.0, reluc=rel_unc, correlation_length=L, covariance_fn=cov_fn)
    source_uc = Uncertainty()

    sink_areas = (
        outlet = [CartesianIndices((1:10, 1:length(dims(surface)[2])))[:]]
    )

    model, get_sample, aggregate = WWFS.make_fns(
        step(dims(surface)[1]),
        surfdem, surfdem_uc,
        beddem, beddem_uc,
        floatfrac, floatfrac_uc,
        source, source_uc,
        sink_areas,
        rmask
    )

    aggr = map_mc(model, get_sample, aggregate, 20)
    serialize(joinpath(outputdir, "aggr_flot_L$(Int(L)).jls"), aggr)
    println("✅ Saved aggr_flot_L$(Int(L)).jls")
end
