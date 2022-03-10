using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra, LogExpFunctions
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks, Sundials
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

## Create vaccine callbacks

include("create_vac_callbacks.jl");

## Projection to daily rate of people becoming double jabbed

fully_vaccinated_interp_30perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_30perc, ts)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_50perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_50perc, ts)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_70perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_70perc, ts)(t) for t in 1:(t_end+100)]
daily_new_fully_vaccinated_30perc = clamp!(diff(fully_vaccinated_interp_30perc)[:], 0.0, Inf)
daily_new_fully_vaccinated_50perc = clamp!(diff(fully_vaccinated_interp_50perc)[:], 0.0, Inf)
daily_new_fully_vaccinated_70perc = clamp!(diff(fully_vaccinated_interp_70perc)[:], 0.0, Inf)

fully_vaccinated_interp_30perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_30perc, ts_rapid)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_50perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_50perc, ts_rapid)(t) for t in 1:(t_end+100)]
fully_vaccinated_interp_70perc_rapid = [DataInterpolations.LinearInterpolation(fully_vaccinated_70perc, ts_rapid)(t) for t in 1:(t_end+100)]
daily_new_fully_vaccinated_30perc_rapid = clamp!(diff(fully_vaccinated_interp_30perc_rapid)[:], 0.0, Inf)
daily_new_fully_vaccinated_50perc_rapid = clamp!(diff(fully_vaccinated_interp_50perc_rapid)[:], 0.0, Inf)
daily_new_fully_vaccinated_70perc_rapid = clamp!(diff(fully_vaccinated_interp_70perc_rapid)[:], 0.0, Inf)

scenarios_vac_rates = [zeros(length(daily_new_fully_vaccinated_30perc)),
    daily_new_fully_vaccinated_30perc, daily_new_fully_vaccinated_50perc, daily_new_fully_vaccinated_70perc,
    daily_new_fully_vaccinated_30perc_rapid, daily_new_fully_vaccinated_50perc_rapid, daily_new_fully_vaccinated_70perc_rapid]

scenarios_names = ["no_vac", "30_perc", "50_perc", "70_perc", "30_perc_rapid", "50_perc_rapid", "70_perc_rapid"]

## Save the number of vaccines by each day 

# dates = [Date(2020, 2, 20) + Day(t) for t = 1:(t_end+100)]
# n = length(dates)
# number_first_dose_30_perc = diff([fully_vaccinated_interp_30perc[(43):end]; fill(fully_vaccinated_interp_30perc[end], 56 - 14)])
# number_second_dose_30_perc = diff([fully_vaccinated_interp_30perc[(15):end]; fill(fully_vaccinated_interp_30perc[end], 14)])
# number_first_dose_50_perc = diff([fully_vaccinated_interp_50perc[(43):end]; fill(fully_vaccinated_interp_50perc[end], 56 - 14)])
# number_second_dose_50_perc = diff([fully_vaccinated_interp_50perc[(15):end]; fill(fully_vaccinated_interp_50perc[end], 14)])
# number_first_dose_70_perc = diff([fully_vaccinated_interp_70perc[(43):end]; fill(fully_vaccinated_interp_70perc[end], 56 - 14)])
# number_second_dose_70_perc = diff([fully_vaccinated_interp_70perc[(15):end]; fill(fully_vaccinated_interp_70perc[end], 14)])

# number_first_dose_30_perc_rapid = diff([fully_vaccinated_interp_30perc_rapid[(43):end]; fill(fully_vaccinated_interp_30perc_rapid[end], 56 - 14)])
# number_second_dose_30_perc_rapid = diff([fully_vaccinated_interp_30perc_rapid[(15):end]; fill(fully_vaccinated_interp_30perc_rapid[end], 14)])
# number_first_dose_50_perc_rapid = diff([fully_vaccinated_interp_50perc_rapid[(43):end]; fill(fully_vaccinated_interp_50perc_rapid[end], 56 - 14)])
# number_second_dose_50_perc_rapid = diff([fully_vaccinated_interp_50perc_rapid[(15):end]; fill(fully_vaccinated_interp_50perc_rapid[end], 14)])
# number_first_dose_70_perc_rapid = diff([fully_vaccinated_interp_70perc_rapid[(43):end]; fill(fully_vaccinated_interp_70perc_rapid[end], 56 - 14)])
# number_second_dose_70_perc_rapid = diff([fully_vaccinated_interp_70perc_rapid[(15):end]; fill(fully_vaccinated_interp_70perc_rapid[end], 14)])


# vaccination_rates_kenya = DataFrame(dates = dates[2:end],
#     number_first_dose_30_perc = number_first_dose_30_perc,
#     number_second_dose_30_perc = number_second_dose_30_perc,
#     number_first_dose_50_perc = number_first_dose_50_perc,
#     number_second_dose_50_perc = number_second_dose_50_perc,
#     number_first_dose_70_perc = number_first_dose_70_perc,
#     number_second_dose_70_perc = number_second_dose_70_perc,
#     number_first_dose_30_perc_rapid = number_first_dose_30_perc_rapid,
#     number_second_dose_30_perc_rapid = number_second_dose_30_perc_rapid,
#     number_first_dose_50_perc_rapid = number_first_dose_50_perc_rapid,
#     number_second_dose_50_perc_rapid = number_second_dose_50_perc_rapid,
#     number_first_dose_70_perc_rapid = number_first_dose_70_perc_rapid,
#     number_second_dose_70_perc_rapid = number_second_dose_70_perc_rapid)

# plot(vaccination_rates_kenya.dates, vaccination_rates_kenya.number_second_dose_70_perc_rapid)

# @rput vaccination_rates_kenya

# R"""
#     save(vaccination_rates_kenya,file = "generation_of_data_for_econ_analysis/vaccinations_per_day.rda")
# """


##make the immune escape callback


## Get basic models for each county

fitfiles = readdir("modelfits", join = true)
model_vect = [@suppress_err load(filename)["model"] for filename in fitfiles]

##Test model

model = model_vect[30]
perc_over_18 = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])
KenyaCoVaccines.change_prob_immune_escape!(model, perc_over_18; startdate = Date(2020, 12, 1), enddate = Date(2021, 12, 1))

##
# scenarios_vac_rates[7]
# first_day_1vac = findfirst(scenarios_vac_rates[7] .> 0) - (Date(2020, 12, 1) - Date(2020, 2, 20)).value - 56
# first_day_2vac = findfirst(scenarios_vac_rates[7] .> 0) - (Date(2020, 12, 1) - Date(2020, 2, 20)).value
# last_day_1vac = findlast(scenarios_vac_rates[7] .> 0) - (Date(2020, 12, 1) - Date(2020, 2, 20)).value - 56
# last_day_2vac = findlast(scenarios_vac_rates[7] .> 0) - (Date(2020, 12, 1) - Date(2020, 2, 20)).value



θs = model |> flatten_chain
function immune_escape!(integrator)
    println(size(integrator.u))
    R = @view integrator.u[:, 7, :]
    #Reduce generation time by 30%
    gen_time_scale = 0.7
    integrator.p[1] = 1 * (1 / gen_time_scale) * integrator.p[1]# β₀
    integrator.p[6] = (1 / gen_time_scale) * integrator.p[6] # α
    integrator.p[8:10] = (1 / gen_time_scale) * integrator.p[8:10]# αP, γA, γM

    #Reduce vaccine protection against acquiring and transmitting infection by 50%
    integrator.p[19:24] = 0.5 * integrator.p[19:24] #ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3
    #Reduce protection from reinfection by 50%
    reinf_red = 0.5
    integrator.p[28] = (1.0 * reinf_red + 0.16 * (1 - reinf_red)) / 0.16
    #50% of R -> W at emergence of immune escape variant but set loss of complete immunity to 0
    integrator.u[:, 7, :] .-= 0.5 .* R
    integrator.u[:, 8, :] .+= 0.5 .* R
    integrator.p[11] = 0.0 # ω

end

fs = falses(6, 10, 3)
fs[:, 7, :] .= true
idxs_R_cvode = fs[:]
fs = falses(6, 10, 3)
fs[:, 8, :] .= true
idxs_W_cvode = fs[:]

function immune_escape_CVODE!(integrator)

    #Reduce generation time by 30%
    gen_time_scale = 0.7
    integrator.p[1] = 1 * (1 / gen_time_scale) * integrator.p[1]# β₀
    integrator.p[6] = (1 / gen_time_scale) * integrator.p[6] # α
    integrator.p[8:10] = (1 / gen_time_scale) * integrator.p[8:10]# αP, γA, γM

    #Reduce vaccine protection against acquiring and transmitting infection by 50%
    integrator.p[19:24] = 0.5 * integrator.p[19:24] #ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3
    #Reduce protection from reinfection by 50%
    reinf_red = 0.5
    integrator.p[28] = (1.0 * reinf_red + 0.16 * (1 - reinf_red)) / 0.16
    #50% of R -> W at emergence of immune escape variant but set loss of complete immunity to 0
    integrator.u[idxs_W_cvode] .+= 0.5 .* integrator.u[idxs_R_cvode]
    integrator.u[idxs_R_cvode] .-= 0.5 .* integrator.u[idxs_R_cvode]
    integrator.p[11] = 0.0 # ω

end

# cb_imm_esc = PresetTimeCallback([349.0], immune_escape!)
cb_imm_esc_cvode = PresetTimeCallback([349.0], immune_escape_CVODE!, save_positions = (false, false))

# cb_imm_esc = CallbackSet()
solver = CVODE_BDF(linear_solver = :GMRES)
@time sol_ie = solve_model_immune_escape(θs[2], model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[2], solver)
sol_ie.u[1]
# E = [sum(u[:, 2, :]) for u in sol_ie.u]
# R = [sum(u[:, 7, :]) for u in sol_ie.u]
# W = [sum(u[:, 8, 1]) for u in sol_ie.u]
C = [sum(u[:, 10, :]) for u in sol_ie.u]
V1 = [sum(u[:, 1:8, 2]) for u in sol_ie.u]
V2 = [sum(u[:, 1:8, 3]) for u in sol_ie.u]
# plot(R)
# plot(W)
# plot(E)

plot(diff(C))
plot(V1)
plot(V2)
plot((V2 .+ V1) ./ sum(model.N[2:6]))
##
# findall(sol_ie.u[end] .< 0)
##

plot(diff(C))
vline!([349], lab = "")
vline!([427], lab = "")
vline!([278.62], lab = "")
##
@time infs = project_daily_infections_by_age_and_vaccine_immune_escape(θs[6], model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[7])
total_infections = sum(infs[2], dims = [2, 3])[:]
plot(total_infections)

##
infs_gathered = gather_infection_type_for_econ_analysis_immune_escape(model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[1];
    mort_scale = 1.0, mort_age = mort_age, mort_var = mort_var,
    H_ICU_to_death_fit = H_ICU_to_death_fit, H_hosp = H_hosp, H_var_ab = H_var_ab, H_var_delta = H_var_delta)

##
A_pred = get_credible_intervals(sum(infs_gathered.asymptomatic_infs, dims = 2)[:, 1, :])
M_pred = get_credible_intervals(sum(infs_gathered.mild_infections, dims = 2)[:, 1, :])
sev_pred = get_credible_intervals(sum(infs_gathered.severe_infs, dims = 2)[:, 1, :])
crit_pred = get_credible_intervals(sum(infs_gathered.critical_infs, dims = 2)[:, 1, :])
deaths_pred = get_credible_intervals(sum(infs_gathered.deadly_infs, dims = 2)[:, 1, :])

## plot(A_pred.pred)
plot!(deaths_pred.pred)

##
for scen_num = 1:2
    println("Beginning immune escape scenario $(scen_num)")
    for k = 1:length(model_vect)
        perc_over_18 = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])
        KenyaCoVaccines.change_prob_immune_escape!(model_vect[k],  perc_over_18; startdate = Date(2020, 12, 1), enddate = Date(2023, 12, 1))
    end

    infs_types = mapreduce((model, mort_scale) -> gather_infection_type_for_econ_analysis_immune_escape(model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[scen_num];
            mort_scale = mort_scale, mort_age = mort_age, mort_var = mort_var,
            H_ICU_to_death_fit = H_ICU_to_death_fit, H_hosp = H_hosp, H_var_ab = H_var_ab, H_var_delta = H_var_delta),
        add_infs, model_vect, mort_scales)

    kenya_inc_A = infs_types.asymptomatic_infs
    kenya_inc_M = infs_types.mild_infections
    kenya_inc_sev = infs_types.severe_infs
    kenya_inc_crit = infs_types.critical_infs
    kenya_inc_deaths = infs_types.deadly_infs

    savefilename_inc_A = "generation_of_data_for_econ_analysis/immune_escape_scenario_$(scen_num)/kenya_inc_A_" * scenarios_names[scen_num] * "_imm_esc.rda"
    savefilename_inc_M = "generation_of_data_for_econ_analysis/immune_escape_scenario_$(scen_num)/kenya_inc_M_" * scenarios_names[scen_num] * "_imm_esc.rda"
    savefilename_inc_sev = "generation_of_data_for_econ_analysis/immune_escape_scenario_$(scen_num)/kenya_inc_sev_" * scenarios_names[scen_num] * "_imm_esc.rda"
    savefilename_inc_crit = "generation_of_data_for_econ_analysis/immune_escape_scenario_$(scen_num)/kenya_inc_crit_" * scenarios_names[scen_num] * "_imm_esc.rda"
    savefilename_inc_deaths = "generation_of_data_for_econ_analysis/immune_escape_scenario_$(scen_num)/kenya_inc_deaths_" * scenarios_names[scen_num] * "_imm_esc.rda"


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

##Test randomized draws



