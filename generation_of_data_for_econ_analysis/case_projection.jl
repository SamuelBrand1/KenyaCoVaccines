
fit = perc_pos_to_num_tests_scale(model; lookback = 150)
θs = model |> flatten_chain

## Project PCR test results from the infection projections

sampled_PCRs = map((inf, θ) -> generate_PCR_predictions(inf, θ, model, fit; M_PCR = 150), infections, θs)
sampled_pos_tests = [sample[:, 1] for sample in sampled_PCRs]
sampled_prop_pos_tests = [[sample[i, 2] == 0 ? 0.0 : sample[i, 1] / sample[i, 2] for i = 1:size(sample, 1)] for sample in sampled_PCRs]

post_pred_pos_tests = get_credible_intervals(sampled_pos_tests)
post_pred_prop_pos_tests = get_credible_intervals(sampled_prop_pos_tests)

scatter(model.PCR_cases[model.model_start_time:end, 1])
plot!(post_pred_pos_tests.pred, ribbon = (post_pred_pos_tests.lb, post_pred_pos_tests.ub))

scatter(model.PCR_cases[model.model_start_time:end, 1] ./ model.PCR_cases[model.model_start_time:end, 2])
plot!(post_pred_prop_pos_tests.pred, ribbon = (post_pred_prop_pos_tests.lb, post_pred_prop_pos_tests.ub))

## Project serological test results from the infection projections

n = length(sero_prevalence_by_age_preds[1].pred)
uncertaintyfirstinfections_by_age = [[infs[:, a] for infs in firstinfections_by_age] for a = 1:size(firstinfections_by_age[1], 2)]
sero_prevalence_by_age = [map((ι_sero, θ) -> generate_sero_predictions(ι_sero, θ, age_grp, model), uncertaintyfirstinfections_by_age[age_grp], θs)
                          for age_grp = 1:6]

sero_data_by_age = [model.sero_cases[:, age_grp, :] for age_grp = 1:6]
dates = [Date(2020, 2, 19) + Day(k) for k = 1:size(model.sero_cases, 1)]
sero_prop_by_age = [monthly_proportions_with_errors(sero, dates, Date(2020, 2, 19)) for sero in sero_data_by_age]

##
k = 1
plot(model.model_start_time:(model.model_start_time-1+n), sero_prevalence_by_age_preds[k].pred,
    lab = "",
    ribbon = (sero_prevalence_by_age_preds[k].lb, sero_prevalence_by_age_preds[k].ub))

scatter!(sero_prop_by_age[k].ts, sero_prop_by_age[k].prop,
    yerrors = (sero_prop_by_age[k].lb, sero_prop_by_age[k].ub),
    lab = "")

##Project variant frequency

n = size(post_pred_pos_tests.pred, 1)
var_frqs = θs .|> θ -> get_variant_prediction(θ, model, n)

wt_frqs = [data.prop_inc_wt for data in var_frqs]
alpha_beta_frqs = [data.prop_inc_alpha_beta for data in var_frqs]
delta_frqs = [data.prop_inc_delta for data in var_frqs]

wt_frq_preds = get_credible_intervals(wt_frqs)
alpha_beta_frq_preds = get_credible_intervals(alpha_beta_frqs)
delta_frq_preds = get_credible_intervals(delta_frqs)

plot(wt_frq_preds.pred, ribbon = (wt_frq_preds.lb, wt_frq_preds.ub), lab = "WT")
plot!(alpha_beta_frq_preds.pred, ribbon = (alpha_beta_frq_preds.lb, alpha_beta_frq_preds.ub), lab = "Alpha/Beta")
plot!(delta_frq_preds.pred, ribbon = (delta_frq_preds.lb, delta_frq_preds.ub), lab = "Delta")

##

X = Matrix(RecursiveArrayTools.vecvec_to_mat(infections)')
pred = get_credible_intervals(X)

plot(pred.pred, ribbon = (pred.lb, pred.ub))

param_names = names(model.MCMC_results.chain)

get(model.MCMC_results.chain[1, :, :], param_names; flatten = true)

##
include("generation_methods.jl");

## Combine incidence data
fitfiles = readdir("modelfits", join = true)
filename = fitfiles[1]
fit_dict = load(filename)
model = fit_dict[first(keys(fit_dict))]
name = model.areaname
inc_by_type = get_daily_incidence_by_type(model, Date(2022, 6, 1),
    1.2,
    1.2,
    0.1,
    0.1,
    0.1)

kenya_inc_asymp = inc_by_type.asymp_infections
kenya_inc_M = inc_by_type.mild_infections
kenya_inc_severe = inc_by_type.severe_infections
kenya_inc_critical = inc_by_type.critical_infections


for filename in fitfiles[2:end]
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    name = model.areaname
    inc_by_type = get_daily_incidence_by_type(model, Date(2022, 6, 1), 1.2, 1.2, 0.1, 0.1, 0.1)
    kenya_inc_asymp .+= inc_by_type.asymp_infections
    kenya_inc_M .+= inc_by_type.mild_infections
    kenya_inc_severe .+= inc_by_type.severe_infections
    kenya_inc_critical .+= inc_by_type.critical_infections
end
