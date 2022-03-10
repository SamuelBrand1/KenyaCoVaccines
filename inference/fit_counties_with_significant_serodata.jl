## Fit counties with significant sero data

## Do inference in parallel

using Distributed

## add number of processors desired
# addprocs(2)
params = (exename = `nice -19 /Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia`,
      dir = "./Github/KenyaCoVaccinesPrivate")
addprocs([("samandfi@10.0.0.6", 7)]; params...)
addprocs(4)

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

##Load contact structure and population size data
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

counties_with_serodata = ["Nairobi", "Mombasa",
      "Nakuru", "Uasin Gishu", "Embu", "Kisumu", "Siaya", "Kisii", "Nyeri",
      "Kilifi", "Kwale"]

priors = [fill(KenyaCoVaccines.priors_onegroup_alpha_delta_variant_cities, 2)
      fill(KenyaCoVaccines.priors_onegroup_alpha_delta_variant_outside_cities, 9)]

#Baseline seropos for low density coverage counties at December 1st
#Based on the unweighted seroprevalence from Ngere et al in Nairobi mid-Sept

init_sero_baseline = [(34 + 74) / (179 + 244), (100 + 95 + 51) / (265 + 241 + 134), 21 / 61, 9 / 40, 9 / 40, 9 / 40]

## parallelised loop for fitting each county
@distributed for i = 1:11
      name = counties_with_serodata[i]
      prior = priors[i]
      model = KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, prior;
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

      KenyaCoVaccines.inferparameters!(model, 2000, trans, 0.05, D; serowaningrate = 0.0)
      @save("modelfits/$(name)_model.jld2", model)

end
##
rmprocs(workers())


##Tests
@everywhere counties_with_serodata = ["Nairobi", "Mombasa",
      "Nakuru", "Uasin Gishu", "Embu", "Kisumu", "Siaya", "Kisii", "Nyeri",
      "Kilifi", "Kwale"]



@everywhere priors = [fill(KenyaCoVaccines.priors_onegroup_alpha_delta_variant_cities, 2)
      fill(KenyaCoVaccines.priors_onegroup_alpha_delta_variant_outside_cities, 9)]


r = @spawnat 2 begin
      name = counties_with_serodata[1]
      prior = priors[1]
      model = KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, prior;
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

      # KenyaCoVaccines.inferparameters!(model, 2000, trans, 0.05, D; serowaningrate = 0.0)
      # @save("modelfits/$(name)_model.jld2", model)
end

model = fetch(r)

name = counties_with_serodata[1]
prior = priors[1]
model = KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, prior;
      PCR_cases = linelist_data_with_pos_neg,
      sero_cases = serology_data,
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
KenyaCoVaccines.inferparameters!(model, 2000, trans, 0.05, D; serowaningrate = 0.0)
