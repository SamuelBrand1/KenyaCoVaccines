## As for plots, these are the ones I can think of:
## * Averted cases by vaccine scenarios over the simulation period(so maybe the time series may be needed after all as an outcome – how fine-grained it depends on you)
## * Averted deaths by vaccine scenarios over the simulation period.
## * Assumed vaccine dose coverage by scenario and age

## * cumulative severe, critical and deaths plots

using Dates, Statistics
using StatsPlots, JLD2, Plots.PlotMeasures, CSV, DataFrames
using Parameters, Distributions, LinearAlgebra, Tullio
using RecursiveArrayTools, MCMCChains, Interpolations, Suppressor
using RCall
import KenyaCoVaccines

## Load background data

@load("data/p_ID.jld2")
@load("data/p_IH.jld2")
@load("data/N_kenya.jld2")
model = @suppress_err load("modelfits/Nairobi_model.jld2")["model"]
age_symptoms = model.p_symp
veff1 = 1 - model.VE_severe_disease_risk[2]
veff2 = 1 - model.VE_severe_disease_risk[3]
veff_death1 = 0.2
veff_death2 = 0.1

include("generation_methods.jl");
include("inference_methods.jl");
include("projection_methods.jl");
@load("data/fitted_sol.jld2")
η_death = sol.u[end]
@load("data/fitted_mort_by_county_age_variant.jld2")
@load("data/fitted_rel_risk_ICU.jld2")
@load("data/fitted_rel_risk_hosp.jld2")

##

p_ID_fitted = normalize!([p * exp(η_death * t) for (t, p) in enumerate(p_ID)], 1)
p_IH_sev = normalize([zeros(5); [exp(η_hosp * t) * p_IH[t] for t = 1:45]], 1)
p_IH_crit = normalize([zeros(5); [exp(η_ICU * t) * p_IH[t] for t = 1:45]], 1)
# p_IH_sev = normalize([exp(0.0 * t) * p_IH[t] for t = 1:50], 1)
# p_IH_crit = normalize([exp(0.0 * t) * p_IH[t] for t = 1:50], 1)


##Gather hospitalisation data

kenya_MoH = CSV.File("data/2022-01-10_MOH_Kenya_user_tweets__EXTRACTED_cleaned_dates.csv") |> DataFrame
kenya_MoH.dates

kenya_MoH.rel_time = [(d - Date(2020, 2, 19)).value for d in kenya_MoH.dates]
idxs_hosp_valid = kenya_MoH.in_hospital .> 0
idxs_ICU_valid = (kenya_MoH.intensive_care_unit .> 0) .& (kenya_MoH.intensive_care_unit .< 4000)
idxs_suppoxy_valid = kenya_MoH.supplementary_oxygen .> 0
idxs_march1 = kenya_MoH.rel_time .>= (Date(2021, 3, 1) - Date(2020, 2, 19)).value

###
### CONSTRUCT DELAY DISTRIBUTIONS
###

#Data distributions
d_incubation = LogNormal(1.644, 0.363)#Lauer estimate
d_duration_in_hosp = Gamma(8.0472378052287681, 0.0017022772867586302 * 365)
mean(d_duration_in_hosp)
std(d_duration_in_hosp)
d_ICUstay = Gamma(53.273947475158401, 0.0003599898843749761 * 365) #From Cameline
mean(d_ICUstay)
std(d_ICUstay)
d_recovery = Exponential(2.4)
#lag distributions
p_IS = [cdf(d_incubation, t) - cdf(d_incubation, t - 1) for t = 1:100] #infection to symptoms
p_ICU = [cdf(d_ICUstay, t) - cdf(d_ICUstay, t - 1) for t = 1:100] #ICU to leave ICU
p_HR = [cdf(d_duration_in_hosp, t) - cdf(d_duration_in_hosp, t - 1) for t = 1:100] #Hospital to discharge
p_ICUR = KenyaCoVaccines.simple_conv(p_ICU, p_HR) #ICU to discharge assuming its the sum of the two
p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
p_R = [cdf(d_recovery, t) - cdf(d_recovery, t - 1) for t = 1:1000]
#Upper tail functions
Q_HR = vcat([1.0], [1 - cumsum(p_HR)[t] for t = 1:100])
Q_ICUH = vcat([1.0], [1 - cumsum(p_ICU)[t] for t = 1:100])
Q_ICUR = vcat([1.0], [1 - cumsum(p_ICUR)[t] for t = 1:100])
Q_R = vcat([1.0], [1 - cumsum(p_R)[t] for t = 1:1000])
F_R = 1 .- Q_R


## Load death data

@load("data/cleaned_linelist20220207_deaths_c_age__date_of_lab_confirmation.jld2")

t0 = findfirst(deaths_data.dates .== Date(2021, 1, 1))
t1 = findfirst(deaths_data.dates .== Date(2021, 10, 31))


deathsdata = cat(sum(deaths_data.deaths[t0:t1, :, 1:4], dims=3),
    sum(deaths_data.deaths[t0:t1, :, 5:10], dims=3),
    sum(deaths_data.deaths[t0:t1, :, 11:12], dims=3),
    sum(deaths_data.deaths[t0:t1, :, 13:14], dims=3),
    sum(deaths_data.deaths[t0:t1, :, 15:16], dims=3),
    sum(deaths_data.deaths[t0:t1, :, 17], dims=3), dims=3)

kenya_deaths = sum(deathsdata, dims=2)[:, 1, :]

N_kenya_agegrp = [sum(N_kenya[1:4, :]),
    sum(N_kenya[5:10, :]),
    sum(N_kenya[11:12, :]),
    sum(N_kenya[13:14, :]),
    sum(N_kenya[15:16, :]),
    sum(N_kenya[17, :])]


## Gather the inf data for severe, critical and deaths from the RDA files

include("gather_inf_data_severe_vac_waning.jl");
include("gather_inf_data_crit_vac_waning.jl");
include("gather_inf_data_deaths_vac_waning.jl");

## Gather the inferred mean daily incidences
# μ_deaths = sum(kenya_inc_deaths_sc_1[:, 4:6, :], dims = 2)[:, 1, :]
μ_deaths_vect = [Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_1[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_1, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_2[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_2, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_3[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_3, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_4[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_4, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_5[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_5, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_6[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_6, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_deaths_sc_7[:, j], p_ID_fitted) for j = 1:size(kenya_inc_deaths_sc_7, 2)])')]

α_deaths = fill(0.05, size(μ_deaths_vect[4]))

t_dec1 = (Date(2020, 12, 1) - Date(2020, 2, 19)).value
t_sept11_model = (Date(2021, 9, 1) - Date(2020, 12, 1)).value
total_deaths_before_dec2020 = sum(kenya_deaths[1:t_dec1, :])
# total_deaths_before_jan2021 = sum(kenya_deaths[1:(t_dec1+30), :])

nb_deaths = [map((μ, α) -> create_NB_dist(μ, α), μ_deaths, fill(0.05, size(μ_deaths))) .|> d -> rand(d) for μ_deaths in μ_deaths_vect] .|> deaths -> get_credible_intervals(1e5 .* deaths[31:end, :] ./ sum(N_kenya_agegrp[1:end]))
nb_cumdeaths = [map((μ, α) -> create_NB_dist(μ, α), μ_deaths, fill(0.05, size(μ_deaths))) .|> d -> rand(d) for μ_deaths in μ_deaths_vect] .|> deaths -> get_credible_intervals(1e5 .* cumsum(deaths[31:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))
nb_cumdeaths_from_sept1 = [map((μ, α) -> create_NB_dist(μ, α), μ_deaths, fill(0.05, size(μ_deaths))) .|> d -> rand(d) for μ_deaths in μ_deaths_vect] .|> deaths -> get_credible_intervals(1e5 .* cumsum(deaths[274:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))

## Gather the inferred mean critical daily incidence 

crit_inc_vect = [Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_1[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_1, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_2[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_2, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_3[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_3, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_4[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_4, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_5[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_5, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_6[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_6, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_crit_sc_7[:, j], p_IH_crit) for j = 1:size(kenya_inc_crit_sc_7, 2)])')]

α_ICU = fill(alpha_ICU, size(crit_inc_vect[1]))

nb_crit = [map((μ, α) -> create_NB_dist(μ, α), ICU_infs, α_ICU) .|> d -> rand(d) for ICU_infs in crit_inc_vect] .|> ICUs -> get_credible_intervals(1e5 .* ICUs[31:end, :] ./ sum(N_kenya_agegrp[1:end]))
nb_cumcrit = [map((μ, α) -> create_NB_dist(μ, α), ICU_infs, α_ICU) .|> d -> rand(d) for ICU_infs in crit_inc_vect] .|> ICUs -> get_credible_intervals(1e5 .* cumsum(ICUs[31:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))

nb_cumcrit_from_sept1 = [map((μ, α) -> create_NB_dist(μ, α), ICU_infs, α_ICU) .|> d -> rand(d) for ICU_infs in crit_inc_vect] .|> ICUs -> get_credible_intervals(1e5 .* cumsum(ICUs[274:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))

## Gather the inferred mean severe daily incidence

sev_inc_vect = [Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_1[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_1, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_2[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_2, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_3[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_3, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_4[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_4, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_5[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_5, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_6[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_6, 2)])'),
    Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(kenya_inc_sev_sc_7[:, j], p_IH_sev) for j = 1:size(kenya_inc_sev_sc_7, 2)])')]

α_hosp = fill(alpha_hosp, size(sev_inc_vect[1]))

nb_sev = [map((μ, α) -> create_NB_dist(μ, α), sev_infs, α_hosp) .|> d -> rand(d) for sev_infs in sev_inc_vect] .|> hosps -> get_credible_intervals(1e5 .* hosps[31:end, :] ./ sum(N_kenya_agegrp[1:end]))
nb_cumsev = [map((μ, α) -> create_NB_dist(μ, α), sev_infs, α_hosp) .|> d -> rand(d) for sev_infs in sev_inc_vect] .|> hosps -> get_credible_intervals(1e5 .* cumsum(hosps[31:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))

nb_cumsev_from_sept1 = [map((μ, α) -> create_NB_dist(μ, α), sev_infs, α_hosp) .|> d -> rand(d) for sev_infs in sev_inc_vect] .|> hosps -> get_credible_intervals(1e5 .* cumsum(hosps[274:end, :], dims=1) ./ sum(N_kenya_agegrp[1:end]))

## Calculate ICU occupancy


ICU_occup_vect = map(mat -> Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(mat[:, j], Q_ICUH) for j = 1:size(mat, 2)])'), crit_inc_vect)
α_ICU = fill(alpha_ICU, size(ICU_occup_vect[1]))

nb_ICU_occup = [map((μ, α) -> create_NB_dist(μ, α), ICU_occups, α_ICU) .|> d -> rand(d) for ICU_occups in ICU_occup_vect] .|> ICUs -> get_credible_intervals(1e5 .* ICUs[31:end, :] ./ sum(N_kenya_agegrp[1:end]))

## Calculate hosp occupancy

ICU_or_hosp_occup_vect = map(mat -> Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(mat[:, j], Q_ICUR) for j = 1:size(mat, 2)])'), crit_inc_vect)
hosp_occup_from_ICU_vect = ICU_or_hosp_occup_vect .- ICU_occup_vect
hosp_occup_from_sev_inc_vect = map(mat -> Matrix(vecvec_to_mat([KenyaCoVaccines.simple_conv(mat[:, j], Q_ICUH) for j = 1:size(mat, 2)])'), sev_inc_vect)
hosp_occup_vect = hosp_occup_from_ICU_vect .+ hosp_occup_from_sev_inc_vect

α_hosp = fill(alpha_hosp, size(ICU_occup_vect[1]))

nb_hosp_occup = [map((μ, α) -> create_NB_dist(μ, α), hosp_occups, α_hosp) .|> d -> rand(d) for hosp_occups in hosp_occup_vect] .|> hosps -> get_credible_intervals(1e5 .* hosps[31:end, :] ./ sum(N_kenya_agegrp[1:end]))


## X ticks and dates

xtickdates = [Date(2020, 2, 1) + Month(k) for k = 19:29]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]
t_sept_model = (Date(2021, 9, 1) - Date(2021, 1, 1)).value
t_sept_data = (Date(2021, 9, 1) - Date(2020, 2, 19)).value
t_end_plot = (Date(2022, 7, 1) - Date(2020, 2, 19)).value

t_start_model = (Date(2021, 1, 1) - Date(2020, 2, 19)).value
t_end_model = (Date(2022, 7, 1) - Date(2020, 2, 19)).value
t_end_model_fit = (Date(2021, 9, 1) - Date(2021, 1, 1)).value
t_end_fit = (Date(2021, 9, 1) - Date(2020, 2, 19)).value

xs_model_proj = (t_end_fit+1):t_end_model
xs_model_fitted_to_data = t_start_model:t_end_fit

n_model_proj = length(xs_model_proj)
n_model_fit = length(xs_model_fitted_to_data)

## Cumulative deaths

plot(xticks=(xtickpos, xticklab),
    size=(800, 400),
    left_margin=5mm,
    xlims=(t_end_fit, t_end_plot + 5),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Cum. reported deaths per 100,000",
    title="Projected cumulative reported deaths",
    tickfont=11, guidefont=13, titlefont=20)

plot!(xs_model_proj, nb_cumdeaths[1].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[1].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[1].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=:grey,
    lab="No vaccination")

plot!(xs_model_proj, nb_cumdeaths[2].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[2].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[2].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1,
    lab="30% coverage")

plot!(xs_model_proj, nb_cumdeaths[3].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[3].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[3].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2,
    lab="50% coverage")

plot!(xs_model_proj, nb_cumdeaths[4].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[4].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[4].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3,
    lab="70% coverage")

plot!(xs_model_proj, nb_cumdeaths[5].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[5].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[5].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(xs_model_proj, nb_cumdeaths[6].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[6].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[6].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(xs_model_proj, nb_cumdeaths[7].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumdeaths[7].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumdeaths[7].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

vline!([(Date(2021, 11, 15) - Date(2020, 2, 19)).value], color=:black, ls=:dash, lab="Scenario: variant introduction date")
# savefig("generation_of_data_for_econ_analysis/main_plots/cum_deaths.pdf")



## Plot cumulative ICU admissions

plt_cum_ICU = plot(xticks=(xtickpos, xticklab),
    size=(800, 400),
    left_margin=5mm,
    xlims=(t_end_fit, t_end_plot + 5),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Cum. ICU admissions per 100,000",
    title="Projected cumulative reported ICU admissions",
    tickfont=11, guidefont=13, titlefont=20)

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[1].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[1].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[1].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=:grey,
    lab="No vaccination")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[2].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[2].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[2].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1,
    lab="30% coverage")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[3].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[3].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[3].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2,
    lab="50% coverage")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[4].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[4].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[4].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3,
    lab="70% coverage")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[5].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[5].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[5].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[6].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[6].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[6].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(plt_cum_ICU, xs_model_proj, nb_cumcrit[7].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumcrit[7].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumcrit[7].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

vline!(plt_cum_ICU, [(Date(2021, 11, 15) - Date(2020, 2, 19)).value], color=:black, ls=:dash, lab="Scenario: variant introduction date")

# savefig(plt_cum_ICU,"generation_of_data_for_econ_analysis/main_plots/cum_ICUadmission.pdf")


## Plot cumulative hosp admissions

plt_cum_hosp = plot(xticks=(xtickpos, xticklab),
    size=(800, 400),
    left_margin=5mm,
    xlims=(t_end_fit, t_end_plot + 5),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Cum. hosp. admissions per 100,000",
    title="Projected cumulative reported health facility admissions",
    tickfont=11, guidefont=13, titlefont=16)

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[1].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[1].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[1].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=:grey,
    lab="No vaccination")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[2].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[2].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[2].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1,
    lab="30% coverage")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[3].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[3].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[3].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2,
    lab="50% coverage")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[4].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[4].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[4].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3,
    lab="70% coverage")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[5].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[5].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[5].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[6].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[6].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[6].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(plt_cum_hosp, xs_model_proj, nb_cumsev[7].pred[(n_model_fit+1):(n_model_fit+n_model_proj)],
    ribbon=(nb_cumsev[7].lb[(n_model_fit+1):(n_model_fit+n_model_proj)], nb_cumsev[7].ub[(n_model_fit+1):(n_model_fit+n_model_proj)]),
    fillalpha=0.1, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

vline!(plt_cum_hosp, [(Date(2021, 11, 15) - Date(2020, 2, 19)).value], color=:black, ls=:dash, lab="Scenario: variant introduction date")

# savefig(plt_cum_hosp,"generation_of_data_for_econ_analysis/main_plots/cum_hospadmission.pdf")

## Reported deaths
# nb_deaths
xtickdates = [Date(2020, 2, 1) + Month(k) for k = 11:28]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]

plt_deaths = plot(xticks=(xtickpos, xticklab),
    size=(1000, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Daily reported deaths per 100,000",
    title="Reported daily deaths vs model",
    tickfont=10, guidefont=13, titlefont=20)
# (n_model_fit+1):(n_model_fit+n_model_proj)
plot!(plt_deaths, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+n_model_proj))], nb_deaths[1].pred[1:(n_model_fit+n_model_proj)],
    ribbon=(nb_deaths[1].lb[1:(n_model_fit+n_model_proj)], nb_deaths[1].ub[1:(n_model_fit+n_model_proj)]),
    fillalpha=0.3, color=:grey, lw=3,
    lab="Model projection")


xs_data = xs_model_fitted_to_data[1]:(xs_model_fitted_to_data[1]-1+size(kenya_deaths, 1))
scatter!(plt_deaths, xs_data, 1e5 .* sum(kenya_deaths, dims=2) ./ sum(N_kenya_agegrp),
    lab="Data used in fitting")


## Alternate daily deaths plot

xs_data = xs_model_fitted_to_data[1]:(xs_model_fitted_to_data[1]-1+size(kenya_deaths, 1))

xtickdates = [Date(2020, 2, 1) + Month(k) for k = 11:29]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]

t_proj = (Date(2021, 9, 2) - Date(2021, 1, 1)).value
t_proj_axis = (Date(2021, 9, 2) - Date(2020, 2, 19)).value


plt_deaths2 = plot(xticks=(xtickpos, xticklab),
    size=(1000, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj + 10),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Daily reported deaths per 100,000",
    title="Projected reported deaths by scenario",
    xtickfont=9, guidefont=13, titlefont=20, lw=1)


plot!(plt_deaths2, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+n_model_proj))], nb_deaths[1].pred[1:(n_model_fit+n_model_proj)],
    ribbon=(nb_deaths[1].lb[1:(n_model_fit+n_model_proj)], nb_deaths[1].ub[1:(n_model_fit+n_model_proj)]),
    fillalpha=0.2, color=:grey,
    lab="No vaccination")

plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[2].pred[t_proj:end])),
    nb_deaths[2].pred[t_proj:end],
    ribbon=(nb_deaths[2].lb[t_proj:end], nb_deaths[2].ub[t_proj:end]),
    fillalpha=0.2, color=1,
    lab="30% coverage")
plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[3].pred[t_proj:end])),
    nb_deaths[3].pred[t_proj:end],
    ribbon=(nb_deaths[3].lb[t_proj:end], nb_deaths[3].ub[t_proj:end]),
    fillalpha=0.2, color=2,
    lab="50% coverage")

plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[4].pred[t_proj:end])),
    nb_deaths[4].pred[t_proj:end],
    ribbon=(nb_deaths[4].lb[t_proj:end], nb_deaths[4].ub[t_proj:end]),
    fillalpha=0.2, color=3,
    lab="70% coverage")

plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[5].pred[t_proj:end])),
    nb_deaths[5].pred[t_proj:end],
    ribbon=(nb_deaths[5].lb[t_proj:end], nb_deaths[5].ub[t_proj:end]),
    fillalpha=0.2, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[6].pred[t_proj:end])),
    nb_deaths[6].pred[t_proj:end],
    ribbon=(nb_deaths[6].lb[t_proj:end], nb_deaths[6].ub[t_proj:end]),
    fillalpha=0.2, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(plt_deaths2,
    (t_proj_axis+1):(t_proj_axis+length(nb_deaths[7].pred[t_proj:end])),
    nb_deaths[7].pred[t_proj:end],
    ribbon=(nb_deaths[7].lb[t_proj:end], nb_deaths[7].ub[t_proj:end]),
    fillalpha=0.2, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

scatter!(plt_deaths2, xs_data, 1e5 .* sum(kenya_deaths, dims=2) ./ sum(N_kenya_agegrp),
    lab="Data used in fitting", color=:black)



bar_lab_strs = reverse(["0%" "30%" "50%" "70%" "30% (rapid)" "50% (rapid)" "70% (rapid)"])
# bar_cum_deaths = bar(orientation=:horizontal,
#                     yticks = (1:7,bar_lab_strs))
bar!(7:-1:1, [nb_cumdeaths_from_sept1[k].pred[end] for k = 1:7],
    yticks=(1:7, bar_lab_strs),
    inset=(1, bbox(0.05, 0.65, 0.2, 0.25, :bottom, :right)),
    subplot=2, orientation=:horizontal,
    lab="", fillcolor=[:grey, 1, 2, 3, 1, 2, 3], fillalpha=[1, 1, 1, 1, 0.5, 0.5, 0.5],
    title="Cum. deaths from 01/09/2021 by cov.",
    xlabel="Deaths per 100,000",
    titlefont=9, guidefont=8,
    bg_inside=nothing,
    grid=nothing)
scatter!([nb_cumdeaths_from_sept1[k].pred[end] for k = 1:7], 7:-1:1, ms=0,
    xerr=([nb_cumdeaths_from_sept1[k].lb[end] for k = 1:7], [nb_cumdeaths_from_sept1[k].ub[end] for k = 1:7]),
    lw=3, lab="", subplot=2)

# savefig(plt_deaths2, "generation_of_data_for_econ_analysis/main_plots/reported_deaths_vs2.png")

## Plot ICU occupancy

xtickdates = [Date(2020, 2, 1) + Month(k) for k = 11:26]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]

t_proj = (Date(2021, 9, 2) - Date(2021, 1, 1)).value
t_proj_axis = (Date(2021, 9, 2) - Date(2020, 2, 19)).value


plt_ICU_occup = plot(xticks=(xtickpos, xticklab),
    size=(800, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj + 10),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="ICU occupancy per 100,000",
    title="Reported ICU occupancy vs model",
    tickfont=11, guidefont=13, titlefont=20)

plot!(plt_ICU_occup, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+100))], nb_ICU_occup[1].pred[1:(n_model_fit+100)],
    ribbon=(nb_ICU_occup[1].lb[1:(n_model_fit+100)], nb_ICU_occup[1].ub[1:(n_model_fit+100)]),
    fillalpha=0.3, color=:grey, lw=3,
    lab="Model projection")

scatter!(plt_ICU_occup, kenya_MoH.rel_time[idxs_ICU_valid.*idxs_march1], 1e5 .* kenya_MoH.intensive_care_unit[idxs_ICU_valid.*idxs_march1] ./ sum(N_kenya_agegrp),
    lab="Data used in fitting")

# savefig(plt_ICU_occup, "generation_of_data_for_econ_analysis/main_plots/ICU_occupancy.pdf")
## Alternate ICU occupancy plot

plt_ICU_occup2 = plot(xticks=(xtickpos, xticklab),
    size=(1000, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj + 10),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="ICU occupancy per 100,000",
    title="Projected reported ICU occupancy by scenario",
    xtickfont=9, guidefont=13, titlefont=20, lw=1)

plot!(plt_ICU_occup2, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+n_model_proj))], nb_ICU_occup[1].pred[1:(n_model_fit+n_model_proj)],
    ribbon=(nb_ICU_occup[1].lb[1:(n_model_fit+n_model_proj)], nb_ICU_occup[1].ub[1:(n_model_fit+n_model_proj)]),
    fillalpha=0.2, color=:grey,
    lab="No vaccination")

plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[2].pred[t_proj:end])),
    nb_ICU_occup[2].pred[t_proj:end],
    ribbon=(nb_ICU_occup[2].lb[t_proj:end], nb_ICU_occup[2].ub[t_proj:end]),
    fillalpha=0.2, color=1,
    lab="30% coverage")
plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[3].pred[t_proj:end])),
    nb_ICU_occup[3].pred[t_proj:end],
    ribbon=(nb_ICU_occup[3].lb[t_proj:end], nb_ICU_occup[3].ub[t_proj:end]),
    fillalpha=0.2, color=2,
    lab="50% coverage")

plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[4].pred[t_proj:end])),
    nb_ICU_occup[4].pred[t_proj:end],
    ribbon=(nb_ICU_occup[4].lb[t_proj:end], nb_ICU_occup[4].ub[t_proj:end]),
    fillalpha=0.2, color=3,
    lab="70% coverage")

plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[5].pred[t_proj:end])),
    nb_ICU_occup[5].pred[t_proj:end],
    ribbon=(nb_ICU_occup[5].lb[t_proj:end], nb_ICU_occup[5].ub[t_proj:end]),
    fillalpha=0.2, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[6].pred[t_proj:end])),
    nb_ICU_occup[6].pred[t_proj:end],
    ribbon=(nb_ICU_occup[6].lb[t_proj:end], nb_ICU_occup[6].ub[t_proj:end]),
    fillalpha=0.2, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(plt_ICU_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_ICU_occup[7].pred[t_proj:end])),
    nb_ICU_occup[7].pred[t_proj:end],
    ribbon=(nb_ICU_occup[7].lb[t_proj:end], nb_ICU_occup[7].ub[t_proj:end]),
    fillalpha=0.2, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

scatter!(plt_ICU_occup2, kenya_MoH.rel_time[idxs_ICU_valid.*idxs_march1][1:118], 1e5 .* kenya_MoH.intensive_care_unit[idxs_ICU_valid.*idxs_march1][1:118] ./ sum(N_kenya_agegrp),
    lab="Data used in fitting", color=:black)


bar_lab_strs = reverse(["0%" "30%" "50%" "70%" "30% (rapid)" "50% (rapid)" "70% (rapid)"])
# bar_cum_deaths = bar(orientation=:horizontal,
#                     yticks = (1:7,bar_lab_strs))
bar!(7:-1:1, [nb_cumcrit_from_sept1[k].pred[end] for k = 1:7],
    yticks=(1:7, bar_lab_strs),
    inset=(1, bbox(0.05, 0.65, 0.2, 0.25, :bottom, :right)),
    subplot=2, orientation=:horizontal,
    lab="", fillcolor=[:grey, 1, 2, 3, 1, 2, 3], fillalpha=[1, 1, 1, 1, 0.5, 0.5, 0.5],
    title="Cum. critical from 01/09/2021 by cov.",
    xlabel="Crit. cases per 100,000",
    titlefont=9, guidefont=8,
    bg_inside=nothing,
    grid=nothing)
scatter!([nb_cumcrit_from_sept1[k].pred[end] for k = 1:7], 7:-1:1, ms=0,
    xerr=([nb_cumcrit_from_sept1[k].lb[end] for k = 1:7], [nb_cumcrit_from_sept1[k].ub[end] for k = 1:7]),
    lw=3, lab="", subplot=2)

# scatter!(plt_ICU_occup2, kenya_MoH.rel_time[idxs_ICU_valid.*idxs_march1], 1e5 .* kenya_MoH.intensive_care_unit[idxs_ICU_valid.*idxs_march1] ./ sum(N_kenya_agegrp),
#     lab="Data used in fitting", color=:black)

# savefig(plt_ICU_occup2, "generation_of_data_for_econ_analysis/main_plots/ICU_occupancy_v2.png")



## Plot hospital occupancy

xtickdates = [Date(2020, 2, 1) + Month(k) for k = 11:21]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]

plt_hosp = plot(xticks=(xtickpos, xticklab),
    size=(800, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj + 10),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Hospital occupancy per 100,000",
    title="Reported hosp. occupancy vs model",
    tickfont=11, guidefont=13, titlefont=20)

plot!(plt_hosp, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+100))], nb_hosp_occup[1].pred[1:(n_model_fit+100)],
    ribbon=(nb_hosp_occup[1].lb[1:(n_model_fit+100)], nb_hosp_occup[1].ub[1:(n_model_fit+100)]),
    fillalpha=0.3, color=:grey, lw=3,
    lab="Model projection")

scatter!(plt_hosp, kenya_MoH.rel_time[idxs_hosp_valid.*idxs_march1], 1e5 .* kenya_MoH.in_hospital[idxs_hosp_valid.*idxs_march1] ./ sum(N_kenya_agegrp),
    lab="Data used in fitting")

# savefig(plt_hosp, "generation_of_data_for_econ_analysis/main_plots/hosp_occupancy.pdf")

## Alternate hosp plot

# xtickdates = [Date(2020, 2, 1) + Month(k) for k = 11:26]
# xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
# xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]

# t_proj = (Date(2021, 9, 2) - Date(2021, 1, 1)).value
# t_proj_axis = (Date(2021, 9, 2) - Date(2020, 2, 19)).value


plt_hosp_occup2 = plot(xticks=(xtickpos, xticklab),
    size=(1000, 400),
    left_margin=5mm,
    xlims=(xs_model_fitted_to_data[1] - 3, xs_model_fitted_to_data[end] + n_model_proj + 10),
    # ylims = (10, 17),
    legend=:topleft,
    ylabel="Hospital occupancy per 100,000",
    title="Projected reported hospital occupancy by scenario",
    xtickfont=9, guidefont=13, titlefont=20, lw=1)

plot!(plt_hosp_occup2, [xs_model_fitted_to_data; collect((xs_model_fitted_to_data[end]+1):(xs_model_fitted_to_data[end]+n_model_proj))], nb_hosp_occup[1].pred[1:(n_model_fit+n_model_proj)],
    ribbon=(nb_hosp_occup[1].lb[1:(n_model_fit+n_model_proj)], nb_hosp_occup[1].ub[1:(n_model_fit+n_model_proj)]),
    fillalpha=0.2, color=:grey,
    lab="No vaccination")

plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[2].pred[t_proj:end])),
    nb_hosp_occup[2].pred[t_proj:end],
    ribbon=(nb_hosp_occup[2].lb[t_proj:end], nb_hosp_occup[2].ub[t_proj:end]),
    fillalpha=0.2, color=1,
    lab="30% coverage")
plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[3].pred[t_proj:end])),
    nb_hosp_occup[3].pred[t_proj:end],
    ribbon=(nb_hosp_occup[3].lb[t_proj:end], nb_hosp_occup[3].ub[t_proj:end]),
    fillalpha=0.2, color=2,
    lab="50% coverage")

plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[4].pred[t_proj:end])),
    nb_hosp_occup[4].pred[t_proj:end],
    ribbon=(nb_hosp_occup[4].lb[t_proj:end], nb_hosp_occup[4].ub[t_proj:end]),
    fillalpha=0.2, color=3,
    lab="70% coverage")

plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[5].pred[t_proj:end])),
    nb_hosp_occup[5].pred[t_proj:end],
    ribbon=(nb_hosp_occup[5].lb[t_proj:end], nb_hosp_occup[5].ub[t_proj:end]),
    fillalpha=0.2, color=1, ls=:dash,
    lab="30% coverage, rapid rollout")

plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[6].pred[t_proj:end])),
    nb_hosp_occup[6].pred[t_proj:end],
    ribbon=(nb_hosp_occup[6].lb[t_proj:end], nb_hosp_occup[6].ub[t_proj:end]),
    fillalpha=0.2, color=2, ls=:dash,
    lab="50% coverage, rapid rollout")

plot!(plt_hosp_occup2,
    (t_proj_axis+1):(t_proj_axis+length(nb_hosp_occup[7].pred[t_proj:end])),
    nb_hosp_occup[7].pred[t_proj:end],
    ribbon=(nb_hosp_occup[7].lb[t_proj:end], nb_hosp_occup[7].ub[t_proj:end]),
    fillalpha=0.2, color=3, ls=:dash,
    lab="70% coverage, rapid rollout")

scatter!(plt_hosp_occup2, kenya_MoH.rel_time[idxs_hosp_valid.*idxs_march1][1:118], 1e5 .* kenya_MoH.in_hospital[idxs_hosp_valid.*idxs_march1][1:118] ./ sum(N_kenya_agegrp),
    lab="Data used in fitting",color = :black)
# scatter!(plt_hosp_occup2, kenya_MoH.rel_time[idxs_hosp_valid.*idxs_march1], 1e5 .* kenya_MoH.in_hospital[idxs_hosp_valid.*idxs_march1] ./ sum(N_kenya_agegrp),
#     lab="Data used in fitting")
bar_lab_strs = reverse(["0%" "30%" "50%" "70%" "30% (rapid)" "50% (rapid)" "70% (rapid)"])
# bar_cum_deaths = bar(orientation=:horizontal,
#                     yticks = (1:7,bar_lab_strs))
bar!(7:-1:1, [nb_cumsev_from_sept1[k].pred[end] for k = 1:7],
    yticks=(1:7, bar_lab_strs),
    inset=(1, bbox(0.05, 0.65, 0.2, 0.25, :bottom, :right)),
    subplot=2, orientation=:horizontal,
    lab="", fillcolor=[:grey, 1, 2, 3, 1, 2, 3], fillalpha=[1, 1, 1, 1, 0.5, 0.5, 0.5],
    title="Cum. severe from 01/09/2021 by cov.",
    xlabel="Sev. cases per 100,000",
    titlefont=9, guidefont=8,
    bg_inside=nothing,
    grid=nothing)
scatter!([nb_cumsev_from_sept1[k].pred[end] for k = 1:7], 7:-1:1, ms=0,
    xerr=([nb_cumsev_from_sept1[k].lb[end] for k = 1:7], [nb_cumsev_from_sept1[k].ub[end] for k = 1:7]),
    lw=3, lab="", subplot=2)
# savefig(plt_hosp_occup2, "generation_of_data_for_econ_analysis/main_plots/hosp_occupancy_vs2.png")

## Checks

