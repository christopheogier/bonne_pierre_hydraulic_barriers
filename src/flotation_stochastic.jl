### stochastic_flotation_june2024.jl ###
using Pkg
Pkg.activate("/scratch-3/cogier/hydraulic_barriers/")

using Rasters, Serialization
using WhereTheWaterFlowsSubglacially, WhereTheWaterFlows
const WWFS = WhereTheWaterFlowsSubglacially
const WWF  = WhereTheWaterFlows

include("functions.jl")
include("LakeAnalysis.jl")
using .LakeAnalysis

# --- Paths ---
datadir_root   = "/scratch-3/cogier/data/BonnePierre_input"
inputdir       = joinpath(datadir_root, "WWFS_input")
outputdir      = "/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis_flot"
mkpath(outputdir)

# --- June 2024 data ---
surface   = clean_raster(Raster("/scratch-3/cogier/data/BonnePierre_output/WWFS_analysis/2024_June_smoothed_surface_with_lakefilled.tif"))
bedrock   = clean_raster(Raster(joinpath(inputdir, "bedrock_resamp_1m.tif")))
thickness = clean_raster(Raster(joinpath(inputdir, "ice_thickness_2024_june.tif")))

# --- Base fields ---
surfdem   = surface
beddem    = bedrock
rmask     = thickness .> 0
floatfrac = ones(size(surfdem))           # baseline f = 1 everywhere
source    = ones(size(surfdem))           # uniform water input (not used for flotation-only sensitivity)

# --- MC setup ---
N        = 20                              # realizations
cov_fn   = WWFS.GRF.gaussian_kernel
rel_unc  = 0.1                             # 10% relative uncertainty on f
corr_Ls  = [50.0, 100.0, 500.0]            # correlation lengths to test (m)

# convenient grid spacing & sink definition
x, y = dims(surface)
dx   = step(x)
# Sink definitions
    sink_areas = (
    outlet = [CartesianIndices((1:10, 1:length(y)))[:],
              CartesianIndices((1:10, 1:(length(y)÷2)))[:]])  # arbitrary

for L in corr_Ls
    println("🔄 Flotation-only WWFS • June 2024 • corr length L = $(Int(L)) m")

    # Uncertainties: only flotation varies
    surfdem_uc   = Uncertainty(absuc=0.0, reluc=0.0)
    beddem_uc    = Uncertainty(absuc=0.0, reluc=0.0)
    floatfrac_uc = Uncertainty(absuc=0.0, reluc=rel_unc, correlation_length=L, covariance_fn=cov_fn)
    source_uc    = Uncertainty()

    model, get_sample, aggregate = WWFS.make_fns(
        dx,
        surfdem,   surfdem_uc,
        beddem,    beddem_uc,
        floatfrac, floatfrac_uc,
        source,    source_uc,
        sink_areas,
        rmask;
        # keep defaults for gamma/rho etc.
        min_lake_depth = 2.0
    )

    # Aggregate maps/statistics over N runs
    aggr = map_mc(model, get_sample, aggregate, N)


    # Save
    outfile = joinpath(outputdir, "aggr_flot_L$(Int(L))_2024_June.jls")
    serialize(outfile, aggr)
    println("✅ Saved $(outfile)")
end
