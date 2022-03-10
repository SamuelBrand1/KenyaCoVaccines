using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra, LogExpFunctions
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks
using RecursiveArrayTools, Interpolations, MCMCChains, Suppressor, Interpolations
using RCall
import KenyaCoVaccines

##


@load("data/p_ID.jld2")
@load("data/N_kenya.jld2")
model = @suppress_err load("modelfits/Nairobi_model.jld2")["model"]
age_symptoms = model.p_symp
veff1 = 1 - model.VE_severe_disease_risk[2]
veff2 = 1 - model.VE_severe_disease_risk[3]
veff_death1 = 0.2
veff_death2 = 0.1

include("generation_methods.jl");
include("inference_methods.jl");

@load("data/fitted_mort_by_county_age_variant.jld2")
@load("data/fitted_rel_risk_ICU.jld2")
@load("data/fitted_rel_risk_hosp.jld2")

## Get vaccine scenarios

t_0 = (Date(2021, 9, 27) - Date(2020, 2, 20)).value #Vaccines start Sept 2021
t_end = (Date(2021, 9, 27) + Month(18) - Date(2020, 2, 20)).value #Option where vaccines finish in 18 months
t_end_rapid = (Date(2021, 9, 27) + Month(6) - Date(2020, 2, 20)).value #Option where vaccines finish in 6 months
ts = [0; t_0; t_end; t_end + 100] .+ 14 # Time points for 18 month deployment lagged 14 days to account for delay between jab and maximum effectiveness
ts_rapid = [0; t_0; t_end_rapid; t_end_rapid + 100] .+ 14 # Time points for 6 month deployment

# Number of people DOUBLE jabbed by time points under three coverages
fully_vaccinated_30perc = [0; 0; 15.9e6 / 2; 15.9e6 / 2]
fully_vaccinated_50perc = [0; 0; 25e6 / 2; 25e6 / 2]
fully_vaccinated_70perc = [0; 0; 35e6 / 2; 35e6 / 2]

# Projection to daily rate of people becoming double jabbed

fully_vaccinated_interp_30perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_30perc, ts)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_50perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_50perc, ts)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_70perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_70perc, ts)(t) for t in 1:(t_end+100)]
daily_new_fully_vaccinated_30perc = diff(fully_vaccinated_interp_30perc)[:]
daily_new_fully_vaccinated_50perc = diff(fully_vaccinated_interp_50perc)[:]
daily_new_fully_vaccinated_70perc = diff(fully_vaccinated_interp_70perc)[:]

fully_vaccinated_interp_30perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_30perc, ts_rapid)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_50perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_50perc, ts_rapid)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_70perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_70perc, ts_rapid)(t) for t in 1:(t_end+100)]
daily_new_fully_vaccinated_30perc_rapid = diff(fully_vaccinated_interp_30perc_rapid)[:]
daily_new_fully_vaccinated_50perc_rapid = diff(fully_vaccinated_interp_50perc_rapid)[:]
daily_new_fully_vaccinated_70perc_rapid = diff(fully_vaccinated_interp_70perc_rapid)[:]

scenarios_vac_rates = [zeros(length(daily_new_fully_vaccinated_30perc)),
    daily_new_fully_vaccinated_30perc, daily_new_fully_vaccinated_50perc, daily_new_fully_vaccinated_70perc,
    daily_new_fully_vaccinated_30perc_rapid, daily_new_fully_vaccinated_50perc_rapid, daily_new_fully_vaccinated_70perc_rapid]

scenarios_names = ["no_vac", "30_perc", "50_perc", "70_perc", "30_perc_rapid", "50_perc_rapid", "70_perc_rapid"]

## Save the number of vaccines by each day 

dates = [Date(2020, 2, 20) + Day(t) for t = 1:(t_end+100)]
n = length(dates)
number_first_dose_30_perc = diff([fully_vaccinated_interp_30perc[(43):end]; fill(fully_vaccinated_interp_30perc[end], 56 - 14)])
number_second_dose_30_perc = diff([fully_vaccinated_interp_30perc[(15):end]; fill(fully_vaccinated_interp_30perc[end], 14)])
number_first_dose_50_perc = diff([fully_vaccinated_interp_50perc[(43):end]; fill(fully_vaccinated_interp_50perc[end], 56 - 14)])
number_second_dose_50_perc = diff([fully_vaccinated_interp_50perc[(15):end]; fill(fully_vaccinated_interp_50perc[end], 14)])
number_first_dose_70_perc = diff([fully_vaccinated_interp_70perc[(43):end]; fill(fully_vaccinated_interp_70perc[end], 56 - 14)])
number_second_dose_70_perc = diff([fully_vaccinated_interp_70perc[(15):end]; fill(fully_vaccinated_interp_70perc[end], 14)])

number_first_dose_30_perc_rapid = diff([fully_vaccinated_interp_30perc_rapid[(43):end]; fill(fully_vaccinated_interp_30perc_rapid[end], 56 - 14)])
number_second_dose_30_perc_rapid = diff([fully_vaccinated_interp_30perc_rapid[(15):end]; fill(fully_vaccinated_interp_30perc_rapid[end], 14)])
number_first_dose_50_perc_rapid = diff([fully_vaccinated_interp_50perc_rapid[(43):end]; fill(fully_vaccinated_interp_50perc_rapid[end], 56 - 14)])
number_second_dose_50_perc_rapid = diff([fully_vaccinated_interp_50perc_rapid[(15):end]; fill(fully_vaccinated_interp_50perc_rapid[end], 14)])
number_first_dose_70_perc_rapid = diff([fully_vaccinated_interp_70perc_rapid[(43):end]; fill(fully_vaccinated_interp_70perc_rapid[end], 56 - 14)])
number_second_dose_70_perc_rapid = diff([fully_vaccinated_interp_70perc_rapid[(15):end]; fill(fully_vaccinated_interp_70perc_rapid[end], 14)])


vaccination_rates_kenya = DataFrame(dates = dates[2:end],
    number_first_dose_30_perc = number_first_dose_30_perc,
    number_second_dose_30_perc = number_second_dose_30_perc,
    number_first_dose_50_perc = number_first_dose_50_perc,
    number_second_dose_50_perc = number_second_dose_50_perc,
    number_first_dose_70_perc = number_first_dose_70_perc,
    number_second_dose_70_perc = number_second_dose_70_perc,
    number_first_dose_30_perc_rapid = number_first_dose_30_perc_rapid,
    number_second_dose_30_perc_rapid = number_second_dose_30_perc_rapid,
    number_first_dose_50_perc_rapid = number_first_dose_50_perc_rapid,
    number_second_dose_50_perc_rapid = number_second_dose_50_perc_rapid,
    number_first_dose_70_perc_rapid = number_first_dose_70_perc_rapid,
    number_second_dose_70_perc_rapid = number_second_dose_70_perc_rapid)

plot(vaccination_rates_kenya.dates, vaccination_rates_kenya.number_second_dose_70_perc_rapid)

@rput vaccination_rates_kenya

R"""
    save(vaccination_rates_kenya,file = "generation_of_data_for_econ_analysis/vaccinations_per_day.rda")
"""

# vaccine_projections = 


## Get basic models for each county

fitfiles = readdir("modelfits", join = true)
model_vect = [@suppress_err load(filename)["model"] for filename in fitfiles]

##
for scen_num = 1:7
    println("Beginning Scenario $(scen_num)")
    for k = 1:length(model_vect)
        perc_over_18 = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])
        KenyaCoVaccines.change_prob!(model_vect[k], scenarios_vac_rates[scen_num].*perc_over_18; startdate = Date(2020, 12, 1), enddate = Date(2021, 12, 1))
    end


    infs_types = mapreduce((model, mort_scale) -> gather_infection_type_for_econ_analysis(model, Date(2023, 6, 1);
            mort_scale = mort_scale, mort_age = mort_age, mort_var = mort_var,
            H_ICU_to_death_fit = H_ICU_to_death_fit, H_hosp = H_hosp, H_var_ab = H_var_ab, H_var_delta = H_var_delta),
        add_infs, model_vect, mort_scales)

    kenya_inc_A = infs_types.asymptomatic_infs
    kenya_inc_M = infs_types.mild_infections
    kenya_inc_sev = infs_types.severe_infs
    kenya_inc_crit = infs_types.critical_infs
    kenya_inc_deaths = infs_types.deadly_infs

    savefilename_inc_A = "generation_of_data_for_econ_analysis/scenario_$(scen_num)/kenya_inc_A_" * scenarios_names[scen_num] * ".rda"
    savefilename_inc_M = "generation_of_data_for_econ_analysis/scenario_$(scen_num)/kenya_inc_M_" * scenarios_names[scen_num] * ".rda"
    savefilename_inc_sev = "generation_of_data_for_econ_analysis/scenario_$(scen_num)/kenya_inc_sev_" * scenarios_names[scen_num] * ".rda"
    savefilename_inc_crit = "generation_of_data_for_econ_analysis/scenario_$(scen_num)/kenya_inc_crit_" * scenarios_names[scen_num] * ".rda"
    savefilename_inc_deaths = "generation_of_data_for_econ_analysis/scenario_$(scen_num)/kenya_inc_deaths_" * scenarios_names[scen_num] * ".rda"


    @rput kenya_inc_A kenya_inc_M kenya_inc_sev kenya_inc_crit kenya_inc_deaths
    @rput savefilename_inc_A savefilename_inc_M savefilename_inc_sev savefilename_inc_crit savefilename_inc_deaths

    R"""
        save(kenya_inc_A,file = savefilename_inc_A)
        save(kenya_inc_M,file = savefilename_inc_M)
        save(kenya_inc_sev,file = savefilename_inc_sev)
        save(kenya_inc_crit,file = savefilename_inc_crit)
        save(kenya_inc_deaths,file = savefilename_inc_deaths)
    """
end


