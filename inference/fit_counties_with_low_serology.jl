## Fit counties with significant sero data

## Do inference in parallel

using Distributed

## add number of processors desired
# addprocs(2)
params = (exename = `nice -19 /Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia`,
    dir = "./Github/KenyaCoVaccinesPrivate")
addprocs([("samandfi@10.0.0.6", 6)]; params...)
addprocs(6)
# params = (exename = `nice -19 /Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia`,
#     dir = "./Github/KenyaCoVaccinesPrivate")
# addprocs([("samandfi@10.0.0.6", :auto)]; params...)
# addprocs(2)

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

## Load relevant code into scope for each worker

@everywhere using JLD2, LinearAlgebra, TransformVariables, LogDensityProblems, MCMCChains, Dates
@everywhere import KenyaCoVaccines

## Load data
# @everywhere begin
using DataFrames
@load("data/linelist_data_with_pos_neg_20feb_to_27sept_c_age__pcr.jld2")
@load("data/cleaned_serology_data_c_age__combined.jld2")
@load("data/deaths_20210919.jld2")
@load("data/p_IH.jld2")
@load("data/p_ID.jld2")

#Load contact structure and population size data
M_Kenya_ho, M_Kenya_other, M_Kenya_school, M_Kenya_work = KenyaCoVaccines.get_rescaled_contact_matrices("data/agemixingmatrix_Kenya_all_types.jld2")
N_kenya = KenyaCoVaccines.get_population_size_matrix("data/2019_census_age_pyramids_counties.csv")


# Parameter transformations and initial KE guess for HMC

trans = as((
    β₀ = as(Real, 0.0, 5.0),  #β₀= as(Real, 0.0, 7.0),
    β_home = as(Real, 0.0, 5.0),
    β_school = as(Real, 0.0, 5.0),
    β_other = as(Real, 0.0, 5.0),
    β_work = as(Real, 0.0, 5.0),
    ϵ = as(Real, 0.0, 0.999),
    χ = as(Real, 0.1, 30.0),
    p_test = as(Real, 0.01, 5),
    E₀ = as(Real, 0.0, 10e3),
    inc_R_αvar = as(Real, 0.0, 1.0),
    time_scale_αvar = as(Real, 0.01, 0.5),
    mid_point_αvar = as(Real, 0.0, 4.0),
    inc_R_δvar = as(Real, 0.001, 1.5),
    time_scale_δvar = as(Real, 0.01, 0.5),
    mid_point_δvar = as(Real, 0.0, 4.0),
    init_scale = as(Real, 0.0, 2.0)
))

D = Diagonal(0.1 * ones(TransformVariables.dimension(trans)))

# end


## Name for each county with substantial serological data, and matching priors

semiurbancounties_lowserodata = ["Kajiado", "Kiambu", "Machakos", "Isiolo", "Kirinyaga", "Murang'a",
    "Trans Nzoia", "Nyandarua", "Meru", "Bungoma", "Nyamira", "Kakamega",
    "Laikipia", "Taita Taveta", "Makueni", "Elgeyo-Marakwet",
    "Tharaka-Nithi", "Kericho", "Bomet", "Homa Bay", "Busia", "Migori",
    "Nandi", "Vihiga"]

ruralcounties_lowserodata = ["Garissa", "Kitui",
    "Lamu", "Marsabit", "Narok", "Tana River", "Turkana", "Wajir",
    "Baringo", "Mandera", "Samburu", "West Pokot"]
counties_with_informed_priors = [semiurbancounties_lowserodata; ruralcounties_lowserodata]

# counties_with_informed_priors = ["Kiambu"
# # , "Machakos", "Murang'a", "Meru",
#     # "Kakamega", "Makueni", "Kericho", "Busia", "Vihiga"]

#Baseline seropos for low density coverage counties at December 1st
#Based on the unweighted seroprevalence from Ngere et al in Nairobi mid-Sept scaled down 33% -> 25% unweighted seroprevalence
init_sero_baseline = [(34 + 74) / (179 + 244), (100 + 95 + 51) / (265 + 241 + 134), 21 / 61, 9 / 40, 9 / 40, 9 / 40] .* (25 / 33)


## Parallelised loop for fitting each county

@distributed for i = 1:36
    name = counties_with_informed_priors[i]
    model = KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, KenyaCoVaccines.priors_onegroup_alpha_delta_variant_noncities_fitted;
        linelist_data_with_pos_neg = linelist_data_with_pos_neg,
        serology_data = serology_data,
        average_sero_init = init_sero_baseline,
        deaths = deaths_data,
        pop_data = N_kenya,
        M_county_ho = M_Kenya_ho,
        M_county_other = M_Kenya_other,
        M_county_school = M_Kenya_school,
        M_county_work = M_Kenya_work,
        rapid = false,
        scenario = 0,
        startdate = Date(2020, 12, 1),
        enddate = Date(2021, 9, 24))

    model.init_sero = init_sero_baseline #Enforcing this because maybe more reliable with low data

    KenyaCoVaccines.inferparameters!(model, 2000, trans, 0.05, D; serowaningrate = 0.0)
    @save("modelfits/$(name)_model.jld2", model)
    println("Finished inference for county $(name)")

end
##
rmprocs(workers())




