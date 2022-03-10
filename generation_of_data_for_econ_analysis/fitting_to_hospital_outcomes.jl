using StatsPlots, Dates, JLD2, Statistics, Optim, Parameters, Distributions, DataFrames, CSV, LinearAlgebra
using Plots.PlotMeasures, GalacticOptim, Test, Flux, NamedArrays, OrdinaryDiffEq, DiffEqCallbacks
using TransformVariables, LogDensityProblems, MCMCChains, DynamicHMC, ForwardDiff, LogExpFunctions, Random, StatsBase
import KenyaCoVaccines

## Load relevant data

M_Kenya_ho, M_Kenya_other, M_Kenya_school, M_Kenya_work = KenyaCoVaccines.get_rescaled_contact_matrices("data/agemixingmatrix_Kenya_all_types.jld2")
N_kenya = KenyaCoVaccines.get_population_size_matrix("data/2019_census_age_pyramids_counties.csv")
@load("data/p_ID.jld2")#NB mean 19.06 days, var 119.78 days^2
# @load("data/cleaned_linelist20211109_deaths_c_age__date_of_lab_confirmation.jld2")
@load("data/deaths_20210919.jld2")
@load("data/N_kenya.jld2")
@load("data/linelist_data_with_pos_neg_20feb_to_27sept_c_age__pcr.jld2")
@load("data/serological_data_with_20feb_to_10Mar2021_age.jld2")

## Load relevant methods

include("generation_methods.jl")
include("fitting_methods.jl")

##



fitfiles = readdir("modelfits", join = true)
filename = fitfiles[30]
fit_dict = load(filename)
model = fit_dict[first(keys(fit_dict))]
name = model.areaname

#Baseline seropos for low density coverage counties at December 1st
#Based on the unweighted seroprevalence from Ngere et al in Nairobi mid-Sept
init_sero_baseline = [(34 + 74) / (179 + 244), (100 + 95 + 51) / (265 + 241 + 134), 21 / 61, 9 / 40, 9 / 40, 9 / 40]



model_vac, vacc_rate_1, vacc_rate_2 = KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, KenyaCoVaccines.priors_onegroup_alpha_delta_variant_cities;
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
    scenario = 3,#<---- activate vaccines
    startdate = Date(2020, 12, 1),
    enddate = Date(2021, 9, 24))




model_no_vac, vacc_rate_1_nv, vacc_rate_2_nv= KenyaCoVaccines.CoVAreaModel(name, KenyaCoVaccines.ll_onegroup_newvariant_infboost, KenyaCoVaccines.priors_onegroup_alpha_delta_variant_cities;
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
    scenario = 0,#<---- activate vaccines
    startdate = Date(2020, 12, 1),
    enddate = Date(2021, 9, 24))

model_vac.MCMC_results = deepcopy(model.MCMC_results)
model_no_vac.MCMC_results = deepcopy(model.MCMC_results)

## TEST CODE



θ_sample = get(model_vac.MCMC_results.chain[1], [:β₀, :ϵ, :β_home, :β_school, :β_other, :β_work, :E₀, :init_scale]; flatten = true)
u0 = KenyaCoVaccines.create_initial_conditions(θ_sample, model_vac)

schoolsclosed_time = Float64((Date(2021, 1, 1) - Date(2020, 12, 1)).value)
function affect_schoolsclosed!(integrator)
    integrator.p[3] = 0.0 #Schools shut
end
schoolsclosed_cb = PresetTimeCallback([schoolsclosed_time], affect_schoolsclosed!; save_positions = (false, false))

#Set up callback for reopening Schools
schoolsreopening_time = Float64((Date(2021, 1, 6) - Date(2020, 12, 1)).value)
function affect_schoolsopening!(integrator)
    integrator.p[3] = β_school # Schools open #
end
schoolsreopening_cb = PresetTimeCallback([schoolsreopening_time], affect_schoolsopening!; save_positions = (false, false))

cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb)
chn = model_vac.MCMC_results.chain
@unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other, ω = model_vac

k = 1
β₀ = chn[:β₀][:][k]
β_home = chn[:β_home][:][k]
β_school = chn[:β_school][:][k]
β_other = chn[:β_other][:][k]
β_work = chn[:β_work][:][k]
ϵ = chn[:ϵ][:][k]
χ = chn[:χ][:][k]
E₀ = chn[:E₀][:][k]
inc_R_αvar = chn[:inc_R_αvar][:][k]
time_scale_αvar = chn[:time_scale_αvar][:][k]
mid_point_αvar = chn[:mid_point_αvar][:][k]
inc_R_δvar = chn[:inc_R_δvar][:][k]
time_scale_δvar = chn[:time_scale_δvar][:][k]
mid_point_δvar = chn[:mid_point_δvar][:][k]
init_scale = chn[:init_scale][:][k]
ve_ac1 = VE_acquisition[1]
ve_ac2 = VE_acquisition[2]
ve_ac3 = VE_acquisition[3]
ve_inf1 = VE_infectiousness[1]
ve_inf2 = VE_infectiousness[2]
ve_inf3 = VE_infectiousness[3]
ve_dis1 = VE_severe_disease_risk[1]
ve_dis2 = VE_severe_disease_risk[2]
ve_dis3 = VE_severe_disease_risk[3]
#Set up parameter vector for the model
p = [β₀, β_home, β_school, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]



sol = solve(model_vac.prob, BS3(), tspan = (0.0, 365.0),
    callback = cb,
    u0 = u0, p = p,
    saveat = 1,
    # isoutofdomain = (u, p, t) -> any(x -> x < -1e-3, u),
    verbose = true)
du = similar(u0)
model_vac.prob.f(du, u0, p, 0.0)
sol.t
sol.u[end]
plot(sol.t, [sum(u[6, 1:8, 3]) / sum(u[6, 1:8, :]) for u in sol.u])

u0 .- u02

model.init_sero
model_no_vac.init_sero


inc_all_with_vac = get_daily_incidence_all_infections(model_vac, Date(2021, 12, 1))
inc_pred_with_vac = get_credible_intervals(sum(inc_all_with_vac.incidence, dims = 2)[:, 1, :])

inc_all_without_vac = get_daily_incidence_all_infections(model_no_vac, Date(2021, 12, 1))
inc_pred_without_vac = get_credible_intervals(sum(inc_all_without_vac.incidence, dims = 2)[:, 1, :])

inc_all_model = get_daily_incidence_all_infections(model, Date(2021, 12, 1))
inc_pred_model = get_credible_intervals(sum(inc_all_model.incidence, dims = 2)[:, 1, :])


plot(inc_pred_with_vac.pred)
plot!(inc_pred_without_vac.pred)
plot!(inc_pred_model.pred)

inc_sev = get_daily_critical_incidence_by_variant_scaled_by_risk(model_vac, Date(2021, 12, 1))
inc_sev_nv = get_daily_critical_incidence_by_variant_scaled_by_risk(model_no_vac, Date(2021, 12, 1))

crit_delta = get_credible_intervals(inc_sev.critical_infections_delta)
crit_delta_nv = get_credible_intervals(inc_sev_nv.critical_infections_delta)
plot(crit_delta.pred)
plot!(crit_delta_nv.pred)


##
x = rand(100, 3, 3)
y = [1, 2, 3]
z = [4, 5, 6]

x[1, :, :] .* y' .* z
## Load methods

include("generation_methods.jl");

## Weight incidence by proportion test positive
I = sortperm(linelist_data_with_pos_neg.areas)
# pos_test_weights = normalize([sum(linelist_data_with_pos_neg.cases[:, I[k], :, 1]) for k = 1:47], 1)
pos_test_weights = normalize([sum(linelist_data_with_pos_neg.cases[:, I[k], :, 1]) for k = 1:47], 1)

# bar(pos_test_weights,
#     orientation = :horizontal,
#     yticks = (1:47, linelist_data_with_pos_neg.areas[I]),
#     size = (500,800))
##
fitfiles = readdir("modelfits_delta", join = true)
filename = fitfiles[1]
fit_dict = load(filename)
model = fit_dict[first(keys(fit_dict))]
name = model.areaname

incidence_pred = get_daily_incidence(model, Date(2022, 1, 1))
variant_pred = get_variant_prediction(incidence_pred, model)
# PCR_pred = get_PCR_prediction
model.relative_testing_rate = [model.relative_testing_rate; ones(1000)]
N̄ = normalize(N_kenya[:, name], 1)
baseline_IHR = pos_test_weights[1] * baseline_IFR_estimate(N̄, risk_inf, IHR_age)
rel_risk_reinf = 0.15

weighted_inc_wt = baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_wt
weighted_inc_wt .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_wt
weighted_inc_alpha_beta = baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_alpha_beta
weighted_inc_alpha_beta .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_alpha_beta
weighted_inc_delta = baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_delta
weighted_inc_delta .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_delta


for (k, filename) in enumerate(fitfiles[2:end])
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    name = model.areaname

    incidence_pred = get_daily_incidence(model, Date(2022, 1, 1))
    variant_pred = get_variant_prediction(incidence_pred, model)
    model.relative_testing_rate = [model.relative_testing_rate; ones(1000)]
    N̄ = normalize(N_kenya[:, name], 1)
    baseline_IHR = pos_test_weights[k+1] * baseline_IFR_estimate(N̄, risk_inf, IHR_age)
    rel_risk_reinf = 0.15

    weighted_inc_wt .+= baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_wt
    weighted_inc_wt .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_wt
    weighted_inc_alpha_beta .+= baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_alpha_beta
    weighted_inc_alpha_beta .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_alpha_beta
    weighted_inc_delta .+= baseline_IHR .* (incidence_pred.firstincidence₁ .+ rel_risk_reinf .* (incidence_pred.incidence₁ .- incidence_pred.firstincidence₁)) .* variant_pred.prop_inc_delta
    weighted_inc_delta .+= baseline_IHR .* (incidence_pred.firstincidence₂ .+ rel_risk_reinf .* (incidence_pred.incidence₂ .- incidence_pred.firstincidence₂)) .* variant_pred.prop_inc_delta
end

## Gather Hospital data
kenya_MoH = CSV.File("data/2021-12-07_MOH_Kenya_user_tweets__EXTRACTED_cleaned_dates.csv") |> DataFrame
kenya_MoH.rel_time = [(d - Date(2020, 2, 20)).value for d in kenya_MoH.dates]
idxs_hosp_valid = kenya_MoH.in_hospital .> 0
idxs_ICU_valid = (kenya_MoH.intensive_care_unit .> 0) .& (kenya_MoH.intensive_care_unit .< 4000)
idxs_suppoxy_valid = kenya_MoH.supplementary_oxygen .> 0
scatter(kenya_MoH.rel_time[idxs_suppoxy_valid], kenya_MoH.supplementary_oxygen[idxs_suppoxy_valid],
    title = "Number on supp. oxy.", lab = "")
scatter(kenya_MoH.rel_time[kenya_MoH.in_hospital.>0], kenya_MoH.in_hospital[kenya_MoH.in_hospital.>0],
    lab = "", title = "Number in hospital")
scatter(kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid],
    lab = "", title = "Number in ICU")
plot(mean(weighted_inc_wt, dims = 2), lab = "wt")
plot!(mean(weighted_inc_alpha_beta, dims = 2), lab = "alpha and beta")
plot!(mean(weighted_inc_delta, dims = 2), lab = "Delta")
vline!([(Date(2021, 2, 1) - Date(2020, 2, 20)).value])

###
### CONSTRUCT DELAY DISTRIBUTIONS
###
##
d_incubation = LogNormal(1.644, 0.363)#Lauer estimate
# d_duration_in_hosp = Weibull(1.572, 17.819)#Haw et al
d_duration_in_hosp = Gamma(8.0472378052287681, 0.0017022772867586302 * 365)
mean(d_duration_in_hosp)
std(d_duration_in_hosp)
# shape = 53.273947475158401, scale = 0.0003599898843749761

# d_ICUstay = Gamma(2.66, 3.42)#Fit to IQR
d_ICUstay = Gamma(53.273947475158401, 0.0003599898843749761 * 365) #From Cameline
mean(d_ICUstay)
std(d_ICUstay)
d_recovery = Exponential(2.4)
#lag distributions
p_IS = [cdf(d_incubation, t) - cdf(d_incubation, t - 1) for t = 1:100]#infection to symptoms
p_ICU = [cdf(d_ICUstay, t) - cdf(d_ICUstay, t - 1) for t = 1:100]#ICU to leave ICU
p_HR = [cdf(d_duration_in_hosp, t) - cdf(d_duration_in_hosp, t - 1) for t = 1:100]#Hospital to discharge
p_ICUR = KenyaCoVSD.simple_conv(p_ICU, p_HR)#ICU to discharge assuming its the sum of the two
p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
p_R = [cdf(d_recovery, t) - cdf(d_recovery, t - 1) for t = 1:1000]
#Upper tail functions
Q_HR = vcat([1.0], [1 - cumsum(p_HR)[t] for t = 1:100])
Q_ICUH = vcat([1.0], [1 - cumsum(p_ICU)[t] for t = 1:100])
Q_ICUR = vcat([1.0], [1 - cumsum(p_ICUR)[t] for t = 1:100])
Q_R = vcat([1.0], [1 - cumsum(p_R)[t] for t = 1:1000])
F_R = 1 .- Q_R

###
### FIT ICU
###
##


## Construct optimisation functions with delays
t_max = (Date(2021, 11, 4) - Date(2020, 2, 20)).value
function sample_hosp_dev(H, H_var,
    weighted_inc_wt,
    weighted_inc_alpha_beta,
    weighted_inc_delta,
    valid_times,
    hosp_data, d, clustering_factor)

    T = eltype(H)
    k = rand(1:size(weighted_inc_wt, 2))
    pred = convert.(T, H .* (H_var[1] .* weighted_inc_wt[:, k] .+ H_var[2] .* weighted_inc_alpha_beta[:, k] .+ H_var[3] .* weighted_inc_delta[:, k]))
    hosp_deviance = T(0)
    n = length(d)
    _pred = [sum(pred[max(1, t - n + 1):t] .* d[min(n, t):-1:1]) for t = 1:length(pred)]
    critical_prev = clamp!(KenyaCoVSD.simple_conv(_pred, Q_ICUH), 0.0, Inf)

    for (i, t) in enumerate(valid_times)
        if t <= t_max
            μ = critical_prev[t]
            σ² = μ + clustering_factor * μ^2
            p_negbin = 1 - (clustering_factor * (μ^2) / σ²)
            r_negbin = 1 / clustering_factor

            hosp_deviance -= logpdf(NegativeBinomial(r_negbin, p_negbin), hosp_data[i])
        end
    end

    return hosp_deviance, critical_prev
end


H = 0.5
H_var = [1.0, 1.5, 3.0]
p_start = 0.117
r_start = 1.870
d = delay(logit(p_start), log(r_start))


dev, pred = sample_hosp_dev(H, H_var, weighted_inc_wt, weighted_inc_alpha_beta, weighted_inc_delta,
    kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid], d, 0.25)
ts = kenya_MoH.rel_time[idxs_ICU_valid]
scatter(kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid])
scatter!(ts, pred[ts])

## 
p = (weighted_inc_wt, weighted_inc_alpha_beta, weighted_inc_delta, kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid])
scatter(p[4], p[5])

function reported_ICU_occupany(x, p)
    weighted_inc_wt, weighted_inc_alpha_beta, weighted_inc_delta, valid_times, hosp_data = p
    T = eltype(x)
    ICU_rel_risk = exp(x[1])
    ICU_rel_risk_var = [T(1), 1 + exp(x[2]), 1 + exp(x[3])]
    logit_p = x[4]
    log_r = x[5]
    cf = exp(x[6])
    # cf = 0.05
    d = delay(logit_p, log_r)
    dev, pred = sample_hosp_dev(ICU_rel_risk, ICU_rel_risk_var, weighted_inc_wt, weighted_inc_alpha_beta, weighted_inc_delta,
        valid_times, hosp_data, d, cf)

    return dev, pred, valid_times
end

x₀ = [H; log.([0.5, 0.5]); logit(p_start); log(r_start); log(0.25)]

dev, pred, ts = reported_ICU_occupany(x₀, p)

scatter(ts, kenya_MoH.intensive_care_unit[idxs_ICU_valid])
scatter!(ts, pred[ts])

dev_vect = fill(-99.0, 3000)
iter = 0
function cb_ICU_fit(x, dev, pred, ts)
    global iter
    iter += 1
    if iter <= length(dev_vect)
        dev_vect[iter] = dev
    end
    println(dev)
    plt = scatter(ts, kenya_MoH.intensive_care_unit[idxs_ICU_valid], lab = "Data", color = :black, ms = 3, fillalpha = 0.3)
    scatter!(ts, pred[ts], lab = "pred", color = :red)
    display(plt)
    return false
end

##
opt_func_ICU = OptimizationFunction(reported_ICU_occupany, GalacticOptim.AutoForwardDiff())
opt_func_ICU = OptimizationProblem(opt_func_ICU, x₀, p)
sol = solve(opt_func_ICU, ADAM(0.01); cb = cb_ICU_fit, maxiters = 3000)

##opt_res_IFR
saved_dev = copy(dev_vect)
x̂ = [opt_res_IFR.u; log(0.25)]
##
opt_prob_IFR = OptimizationProblem(opt_func_IFR, opt_res_IFR.u, p)
opt_res_IFR = solve(opt_prob_IFR, ADAM(0.001); cb = cb_ICU_fit, maxiters = 3000)
##
saved_dev = [saved_dev; copy(dev_vect)]
scatter(saved_dev)
plot!(mvav_cols(saved_dev), lw = 5)
hline!([530], lw = 5)

##
ICU_rel_risk = exp.(sol[1])
ICU_rel_var = vcat([1.0], 1 .+ exp.(sol[2:3]))
p_hat_ICU = logistic(sol[4])
r_hat_ICU = exp(sol[5])
cf = exp(sol[6])
d_fitted = delay(sol[4], sol[5])

norm = sum(d_fitted)
@save("outcome_parameter_fits/ICU_fits.jld2", ICU_rel_risk, ICU_rel_var, p_hat_ICU, r_hat_ICU, cf, d_fitted)

###
### FIT HOSPITALISATIONS
###

## Get critical incidence best estimate, this gives adjustment to hospitalisation prevalence due to people leaving ICU


n = length(d_fitted)

pred_array = ICU_rel_risk .* (ICU_rel_var[1] .* weighted_inc_wt .+ ICU_rel_var[2] .* weighted_inc_alpha_beta .+ ICU_rel_var[3] .* weighted_inc_delta)
for k = 1:2000
    pred_array[:, k] .= [sum(pred_array[max(1, t - n + 1):t, k] .* d_fitted[min(n, t):-1:1]) for t = 1:size(pred_array, 1)]
end
fitted_critical_incidence = pred_array
fitted_critical_prev = copy(pred_array)
for k = 1:2000
    fitted_critical_prev[:, k] .= KenyaCoVSD.simple_conv(fitted_critical_incidence[:, k], Q_ICUH)
end
fitted_critical_in_hospital_prev = copy(pred_array)
for k = 1:2000
    fitted_critical_in_hospital_prev[:, k] .= KenyaCoVSD.simple_conv(fitted_critical_in_hospital_prev[:, k], Q_ICUR)
end

fitted_in_hosp_ward_from_ICU = fitted_critical_in_hospital_prev .- fitted_critical_prev

## Get delayed predictions based on the delay fitted from criticals

delayed_hosp_inc_wt = weighted_inc_wt
delayed_hosp_inc_alpha_beta = weighted_inc_alpha_beta
delayed_hosp_inc_delta = weighted_inc_delta

for k = 1:2000
    delayed_hosp_inc_wt[:, k] .= [sum(delayed_hosp_inc_wt[max(1, t - n + 1):t, k] .* d_fitted[min(n, t):-1:1]) for t = 1:size(delayed_hosp_inc_wt, 1)]
    delayed_hosp_inc_alpha_beta[:, k] .= [sum(delayed_hosp_inc_alpha_beta[max(1, t - n + 1):t, k] .* d_fitted[min(n, t):-1:1]) for t = 1:size(delayed_hosp_inc_alpha_beta, 1)]
    delayed_hosp_inc_delta[:, k] .= [sum(delayed_hosp_inc_delta[max(1, t - n + 1):t, k] .* d_fitted[min(n, t):-1:1]) for t = 1:size(delayed_hosp_inc_delta, 1)]
end

##
p2 = (delayed_hosp_inc_wt, delayed_hosp_inc_alpha_beta, delayed_hosp_inc_delta, kenya_MoH.rel_time[kenya_MoH.in_hospital.>0], kenya_MoH.in_hospital[kenya_MoH.in_hospital.>0])
# kenya_MoH.rel_time[kenya_MoH.in_hospital.>0], kenya_MoH.in_hospital[kenya_MoH.in_hospital.>0]
scatter(p2[4], p2[5])
scatter!(p2[4], delayed_hosp_inc_delta[p[4], 1])


function sample_ward_dev(W, W_var,
    delayed_hosp_inc_wt,
    delayed_hosp_inc_alpha_beta,
    delayed_hosp_inc_delta,
    valid_times,
    ward_data, clustering_factor)

    T = eltype(W_var)
    k = rand(1:size(delayed_hosp_inc_wt, 2))
    pred = convert.(T, W .* (W_var[1] .* delayed_hosp_inc_wt[:, k] .+ W_var[2] .* delayed_hosp_inc_alpha_beta[:, k] .+ W_var[3] .* delayed_hosp_inc_delta[:, k]))
    hosp_deviance = T(0)
    ward_prev = clamp!(KenyaCoVSD.simple_conv(pred, Q_HR), 0.0, Inf) .+ fitted_in_hosp_ward_from_ICU[:, k]

    for (i, t) in enumerate(valid_times)
        if t <= t_max
            μ = ward_prev[t]
            σ² = μ + clustering_factor * μ^2
            p_negbin = 1 - (clustering_factor * (μ^2) / σ²)
            r_negbin = 1 / clustering_factor

            hosp_deviance -= logpdf(NegativeBinomial(r_negbin, p_negbin), ward_data[i])
        end
    end

    return hosp_deviance, ward_prev
end

dev, pred = sample_ward_dev(1.0, [1.0, 1.5, 1.5], p2[1], p2[2], p2[3], p2[4], p2[5], 0.25)
scatter(p2[4], p2[5])
scatter!(p2[4], pred[p2[4]])

function reported_hosp_occupany(x, p)
    delayed_hosp_inc_wt, delayed_hosp_inc_alpha_beta, delayed_hosp_inc_delta, valid_times, ward_data = p
    T = eltype(x)
    hosp_rel_risk = exp(x[1])
    hosp_rel_risk_var = [T(1), 1 + exp(x[2]), 1 + exp(x[3])]
    cf = exp(x[4])

    dev, pred = sample_ward_dev(hosp_rel_risk, hosp_rel_risk_var, delayed_hosp_inc_wt, delayed_hosp_inc_alpha_beta, delayed_hosp_inc_delta,
        valid_times, ward_data, cf)

    return dev, pred, valid_times
end

dev_vect_H = fill(-99.0, 3000)
iter = 0
function cb_hosp_fit(x, dev, pred, ts)
    global iter
    iter += 1
    if iter <= length(dev_vect)
        dev_vect_H[iter] = dev
    end
    println(dev)
    plt = scatter(ts, kenya_MoH.in_hospital[kenya_MoH.in_hospital.>0], lab = "Data", color = :black, ms = 3, fillalpha = 0.3)
    scatter!(ts, pred[ts], lab = "pred", color = :red)
    display(plt)
    return false
end


x₀_hosp = [0.0; log.([0.5, 0.5]); log(0.2)]
opt_func_hosp = OptimizationFunction(reported_hosp_occupany, GalacticOptim.AutoForwardDiff())
opt_prob_hosp = OptimizationProblem(opt_func_hosp, x₀_hosp, p2)
sol_hosp_5 = solve(opt_prob_hosp, ADAM(0.01); cb = cb_hosp_fit, maxiters = 3000)

# saved_dev = [saved_dev; copy(dev_vect)]
saved_dev = copy(dev_vect)
scatter(saved_dev)
plot!(mvav_cols(saved_dev), lw = 5)
hline!([530], lw = 5)

##
hosp_rel_risk_5days = exp(sol_hosp_5[1])
hosp_rel_var_5days = [1.0, 1 + exp(sol_hosp_5[2]), 1 + exp(sol_hosp_5[3])]
cf_hosp_5days = exp(sol_hosp_5[4])

@save("outcome_parameter_fits/hosp_fits_5days.jld2", hosp_rel_risk_5days, hosp_rel_var_5days, cf_hosp_5days)
# hosp_rel_risk_16days = exp(sol_hosp_16[1])
# hosp_rel_var_16days = [1.0, 1 + exp(sol_hosp_16[2]), 1 + exp(sol_hosp_16[3])]
# cf_hosp_16days = exp(sol_hosp_16[4])
# @save("outcome_parameter_fits/hosp_fits_16days.jld2", hosp_rel_risk_16days, hosp_rel_var_16days, cf_hosp_16days)

