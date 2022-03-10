using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra, LogExpFunctions
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks, Turing
using RecursiveArrayTools, Interpolations, MCMCChains, Suppressor, Interpolations
import KenyaCoVaccines

##Load methods

include("generation_methods.jl");
include("inference_methods.jl");
@load("data/p_ID.jld2")
@load("data/N_kenya.jld2")

## Load the first and fully vaccinated data from OWID, then lag by 14 days to account for delay to full effectiveness plus use interpolation to fill in gaps in reporting
## Considering three strategies: 
# Scenario 2: Current government vaccination strategy: 30% of population over 18, 15.9 million vaccine doses
# Scenario 3: Higher population coverage 50% of the population over 18 years, 25 million doses
# Scenario 4: Higher population coverage 70% of the population over 18 years, 35 million doses
# owid_kenya = CSV.File("data/owid_20_jan_2022.csv") |> DataFrame
# t_end = (Date(2022, 6, 1) - Date(2020, 2, 20)).value

# _ts = [(Date(d, DateFormat("dd/mm/yyyy")) - Date(2020, 2, 20)).value for d in owid_kenya.date]
# _fully_vaccinated = owid_kenya.people_fully_vaccinated
# idxs2 = (.~ismissing.(_fully_vaccinated))
# ts2 = [0; _ts[idxs2][1] - 1; _ts[idxs2]; t_end; t_end + 100] .+ 14
# fully_vaccinated_30perc = [0; 0; _fully_vaccinated[idxs2]; 15.9e6 / 2; 15.9e6 / 2]
# fully_vaccinated_50perc = [0; 0; _fully_vaccinated[idxs2]; 25e6 / 2; 25e6 / 2]
# fully_vaccinated_70perc = [0; 0; _fully_vaccinated[idxs2]; 35e6 / 2; 35e6 / 2]


# fully_vaccinated_interp_30perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_30perc, ts2)(t) for t in 1:(t_end+100)]
# fully_vaccinated_interp_50perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_50perc, ts2)(t) for t in 1:(t_end+100)]
# fully_vaccinated_interp_70perc = [DataInterpolations.LinearInterpolation(fully_vaccinated_70perc, ts2)(t) for t in 1:(t_end+100)]

# daily_new_fully_vaccinated = diff(fully_vaccinated_interp_30perc)[:]
# plot(daily_new_fully_vaccinated)
# plot!(clamp!(smooth_vector(daily_new_fully_vaccinated), 0.0, Inf), lw = 2)

t_0 = (Date(2021, 9, 27) - Date(2020, 2, 20)).value
t_end = (Date(2021, 9, 27) + Month(18) - Date(2020, 2, 20)).value
t_end_rapid = (Date(2021, 9, 27) + Month(6) - Date(2020, 2, 20)).value
ts = [0; t_0; t_end; t_end + 100] .+ 14
ts_rapid = [0; t_0; t_end_rapid; t_end_rapid + 100] .+ 14

fully_vaccinated_30perc = [0; 0; 15.9e6 / 2; 15.9e6 / 2]
fully_vaccinated_50perc = [0; 0; 25e6 / 2; 25e6 / 2]
fully_vaccinated_70perc = [0; 0; 35e6 / 2; 35e6 / 2]

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

plot(daily_new_fully_vaccinated_30perc)
plot!(daily_new_fully_vaccinated_50perc)
plot!(daily_new_fully_vaccinated_30perc_rapid)


##
sol_mort_fit = load("fitted_sol.jld2")["sol"]
x = sol_mort_fit.u
mort_scales = [logistic.(x[1:29]); exp(x[30]); logistic.(x[31:47])]
mort_age = logistic.(x[(47+1):(47+6)])
bar(mort_age)

mort_var = exp.(x[54:55])
η = x[56]

##Make projections 
prop_vac = sum(N_kenya[5:end, "Baringo"]) / sum(N_kenya[5:end, :])

model = @suppress_err load("modelfits/Nairobi_model.jld2")["model"]
model_novac = deepcopy(model)
KenyaCoVaccines.change_prob!(model_novac, zeros(length(daily_new_fully_vaccinated_30perc)) .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_30perc = deepcopy(model)
KenyaCoVaccines.change_prob!(model_30perc, daily_new_fully_vaccinated_30perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_50perc = deepcopy(model)
KenyaCoVaccines.change_prob!(model_50perc, daily_new_fully_vaccinated_50perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_70perc = deepcopy(model)
KenyaCoVaccines.change_prob!(model_70perc, daily_new_fully_vaccinated_70perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_30perc_rapid = deepcopy(model)
KenyaCoVaccines.change_prob!(model_30perc_rapid, daily_new_fully_vaccinated_30perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_50perc_rapid = deepcopy(model)
KenyaCoVaccines.change_prob!(model_50perc_rapid, daily_new_fully_vaccinated_50perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
model_70perc_rapid = deepcopy(model)
KenyaCoVaccines.change_prob!(model_70perc_rapid, daily_new_fully_vaccinated_70perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))


inf_predictions_novac = gather_post_mean_infection_rates_by_type(model_novac, Date(2023, 6, 1))
pred_novac = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions_novac, model.VE_severe_disease_risk, η, p_ID)
inf_predictions = gather_post_mean_infection_rates_by_type(model_30perc, Date(2023, 6, 1))
pred_30perc = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
inf_predictions = gather_post_mean_infection_rates_by_type(model_50perc, Date(2023, 6, 1))
pred_50perc = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
inf_predictions = gather_post_mean_infection_rates_by_type(model_70perc, Date(2023, 6, 1))
pred_70perc = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)


inf_predictions = gather_post_mean_infection_rates_by_type(model_30perc_rapid, Date(2023, 6, 1))
pred_30perc_rapid = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
inf_predictions = gather_post_mean_infection_rates_by_type(model_50perc_rapid, Date(2023, 6, 1))
pred_50perc_rapid = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
inf_predictions = gather_post_mean_infection_rates_by_type(model_70perc_rapid, Date(2023, 6, 1))
pred_70perc_rapid = pred_death_rate_county(mort_scales[1], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)

fitfiles = readdir("modelfits", join = true)

for k = 2:47
    model = @suppress_err load(fitfiles[k])["model"]
    println("Making vaccine projections for $(model.areaname)")
    prop_vac = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])

    model_novac = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_novac, zeros(length(daily_new_fully_vaccinated_30perc)) .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_30perc = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_30perc, daily_new_fully_vaccinated_30perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_50perc = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_50perc, daily_new_fully_vaccinated_50perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_70perc = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_70perc, daily_new_fully_vaccinated_70perc .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_30perc_rapid = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_30perc_rapid, daily_new_fully_vaccinated_30perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_50perc_rapid = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_50perc_rapid, daily_new_fully_vaccinated_50perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))
    model_70perc_rapid = deepcopy(model)
    KenyaCoVaccines.change_prob!(model_70perc_rapid, daily_new_fully_vaccinated_70perc_rapid .* prop_vac; startdate = Date(2020, 12, 1), enddate = Date(2023, 6, 1))

    inf_predictions_novac = gather_post_mean_infection_rates_by_type(model_novac, Date(2023, 6, 1))
    pred_novac .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions_novac, model.VE_severe_disease_risk, η, p_ID)
    inf_predictions = gather_post_mean_infection_rates_by_type(model_30perc, Date(2023, 6, 1))
    pred_30perc .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
    inf_predictions = gather_post_mean_infection_rates_by_type(model_50perc, Date(2023, 6, 1))
    pred_50perc .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
    inf_predictions = gather_post_mean_infection_rates_by_type(model_70perc, Date(2023, 6, 1))
    pred_70perc = pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)


    inf_predictions = gather_post_mean_infection_rates_by_type(model_30perc_rapid, Date(2023, 6, 1))
    pred_30perc_rapid .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
    inf_predictions = gather_post_mean_infection_rates_by_type(model_50perc_rapid, Date(2023, 6, 1))
    pred_50perc_rapid .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)
    inf_predictions = gather_post_mean_infection_rates_by_type(model_70perc_rapid, Date(2023, 6, 1))
    pred_70perc_rapid .+= pred_death_rate_county(mort_scales[k], mort_age, mort_var, 0.15, inf_predictions, model.VE_severe_disease_risk, η, p_ID)


end


##
plot(sum(pred_novac, dims = 1)[:])
pred_70perc_rapid
plot!(sum(pred_70perc_rapid, dims = 1)[:])

@save("death_projs.jld2", pred_novac, pred_30perc, pred_50perc, pred_70perc, pred_30perc_rapid, pred_50perc_rapid, pred_70perc_rapid)