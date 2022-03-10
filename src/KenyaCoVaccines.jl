__precompile__(true)

module KenyaCoVaccines

using DynamicHMC,Dates,OrdinaryDiffEq,DiffEqCallbacks,JLD2,Parameters,Distributions
using NamedArrays,TransformVariables,LogDensityProblems,LinearAlgebra,BlackBoxOptim,Random
using DynamicHMC.Diagnostics,MCMCDiagnostics,MCMCChains, Interpolations, FileIO, CSV, DataFrames
using RecursiveArrayTools
using DataInterpolations
using StatsPlots
using Plots.PlotMeasures
using ForwardDiff

relative_testing_rate = vcat([0.0 for d in Date(2020,2,21):Day(1):Date(2020,3,14)],
        [t/120 for (t,d) in enumerate(Date(2020,3,15):Day(1):Date(2020,7,12))],
        [1.0 for k = 1:500])

include("transmissionmodel.jl")
include("types.jl")
include("likelihood.jl")
include("priors.jl")
include("inference.jl")
include("dataflow.jl")
include("vacc_scenarios.jl")

#Load defaults for the PCR and serological sensitivity after infection
@load("data/rel_sero_detection_after_infection.jld2")
@load("data/smoothedPCR_detection_after_infection.jld2")

#Load the prior for the p_test and chi parameters
@load("data/non_urban_detection_fits.jld2")

#@load("data/semi_urban_rural_detection_fits.jld2")
#@load("data/rural_fits.jld2")


end # module
