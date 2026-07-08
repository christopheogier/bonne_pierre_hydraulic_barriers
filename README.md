# Bonne Pierre hydraulic barriers

This repository contains the workflow used to reproduce the hydraulic-barrier analysis for the Bonne Pierre case study discussed in the paper:

Ogier et al. (2026), *Potential glacier contributions to the 2024 La Bérarde flood*.
Preprint: https://egusphere.copernicus.org/preprints/2026/egusphere-2026-466/

## Scripts

- prepare.jl: prepares the DEMs, thickness fields, and uncertainty rasters needed as input for the routing workflow.
- main.jl: runs the deterministic subglacial routing analysis and writes water pockets and hydraulic-head outputs.
- main_stochastic.jl: runs the stochastic subglacial workflow with uncertainty propagation and saves Monte Carlo summaries.
- analyse_stoch.jl: analyzes the stochastic outputs, builds summary tables, and produces the main figures.
- synthetic_phi.jl: generates synthetic hydraulic-potential examples.
- lbocv_gpr.py: builds leave-ball-out GPR cross-validation subsets for estimating interpolation uncertainty.
- analyze_lbocv_res_gpr.py: analyzes leave-ball-out residuals to estimate mass-conservation interpolation error and spatial correlation.

The scripts write outputs into the Bonne Pierre data/output folders. The underlying input and output data for the workflow are available in the research collection: http://hdl.handle.net/20.500.11850/800814


