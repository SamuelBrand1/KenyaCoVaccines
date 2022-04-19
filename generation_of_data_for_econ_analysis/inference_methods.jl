"""
    function gather_post_mean_infection_rates_by_type(model::KenyaCoVaccines.CoVAreaModel)

Gather the posterior mean predictions for infections by: age, dose, primary/reinfection episode, variant.        
"""
function gather_post_mean_infection_rates_by_type(model::KenyaCoVaccines.CoVAreaModel, projection_date::Date)
    infections_by_age_dose = model |> flatten_chain .|> θ -> project_daily_infections_by_age_and_vaccine(θ, model, projection_date)

    uncertaintyinfections_by_age_unvac_first = [[infs[1][:, a, 1] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]
    uncertaintyinfections_by_age_unvac_reinf = [[infs[2][:, a, 1] .- infs[1][:, a, 1] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]
    uncertaintyinfections_by_age_1dose_first = [[infs[1][:, a, 2] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]
    uncertaintyinfections_by_age_1dose_reinf = [[infs[2][:, a, 2] .- infs[1][:, a, 2] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]
    uncertaintyinfections_by_age_2dose_first = [[infs[1][:, a, 3] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]
    uncertaintyinfections_by_age_2dose_reinf = [[infs[2][:, a, 3] .- infs[1][:, a, 3] for infs in infections_by_age_dose] for a = 1:size(infections_by_age_dose[1][1], 2)]

    n = size(uncertaintyinfections_by_age_unvac_first[1][1], 1)
    θs = model |> flatten_chain
    var_frqs = θs .|> θ -> get_variant_prediction(θ, model, n)
    wt_frqs = [data.prop_inc_wt for data in var_frqs]
    alpha_beta_frqs = [data.prop_inc_alpha_beta for data in var_frqs]
    delta_frqs = [data.prop_inc_delta for data in var_frqs]



    return (infections_by_age_unvac_first_pred = [get_credible_intervals(uncertaintyinfections_by_age_unvac_first[age]).pred for age = 1:6],
        infections_by_age_1dose_first_pred = [get_credible_intervals(uncertaintyinfections_by_age_1dose_first[age]).pred for age = 1:6],
        infections_by_age_2dose_first_pred = [get_credible_intervals(uncertaintyinfections_by_age_2dose_first[age]).pred for age = 1:6],
        infections_by_age_unvac_reinf_pred = [get_credible_intervals(uncertaintyinfections_by_age_unvac_reinf[age]).pred for age = 1:6],
        infections_by_age_1dose_reinf_pred = [get_credible_intervals(uncertaintyinfections_by_age_1dose_reinf[age]).pred for age = 1:6],
        infections_by_age_2dose_reinf_pred = [get_credible_intervals(uncertaintyinfections_by_age_2dose_reinf[age]).pred for age = 1:6],
        wt_frq_pred = get_credible_intervals(wt_frqs).pred,
        alpha_beta_frq_pred = get_credible_intervals(alpha_beta_frqs).pred,
        delta_frq_pred = get_credible_intervals(delta_frqs).pred)

end

"""
    function gather_post_mean_infection_rates_by_type_novac(model::KenyaCoVaccines.CoVAreaModel)

Gather the posterior mean predictions for infections by: age, primary/reinfection episode, variant. Only to be used with zero vaccination rates.
"""
function gather_post_mean_infection_rates_by_type_novac(model::KenyaCoVaccines.CoVAreaModel, projection_date::Date)
    firstinfections_by_age_dose = model |> flatten_chain .|> θ -> project_daily_infections_by_age_and_vaccine(θ, model, projection_date)

    uncertaintyfirstinfections_by_age_first = [[infs[1][:, a, 1] for infs in firstinfections_by_age_dose] for a = 1:size(firstinfections_by_age_dose[1][1], 2)]
    uncertaintyfirstinfections_by_age_reinf = [[infs[2][:, a, 1] .- infs[1][:, a, 1] for infs in firstinfections_by_age_dose] for a = 1:size(firstinfections_by_age_dose[1][1], 2)]

    n = size(uncertaintyfirstinfections_by_age_first[1][1], 1)
    θs = model |> flatten_chain
    var_frqs = θs .|> θ -> get_variant_prediction(θ, model, n)
    wt_frqs = [data.prop_inc_wt for data in var_frqs]
    alpha_beta_frqs = [data.prop_inc_alpha_beta for data in var_frqs]
    delta_frqs = [data.prop_inc_delta for data in var_frqs]

    return (infections_by_age_first_pred = [get_credible_intervals(uncertaintyfirstinfections_by_age_first[age]).pred for age = 1:6],
        infections_by_age_reinf_pred = [get_credible_intervals(uncertaintyfirstinfections_by_age_reinf[age]).pred for age = 1:6],
        wt_frq_pred = get_credible_intervals(wt_frqs).pred,
        alpha_beta_frq_pred = get_credible_intervals(alpha_beta_frqs).pred,
        delta_frq_pred = get_credible_intervals(delta_frqs).pred)

end

"""
    function pred_death_rate_county(mort_scale, mort_age, mort_var, mort_reinf, predictions, VE_severe_disease_risk, η, p_ID)

Given various parameters, convert the infection predictions into a death rate predictions.       
"""
function pred_death_rate_county(mort_scale, mort_age, mort_var, mort_reinf, predictions, VE_severe_disease_risk, η, p_ID)
    pred_mean = vecvec_to_mat((1 - VE_severe_disease_risk[1]) .* mort_age .* predictions.infections_by_age_unvac_first_pred .+ (1 - VE_severe_disease_risk[2]) .* mort_age .* predictions.infections_by_age_1dose_first_pred .+ (1 - VE_severe_disease_risk[3]) .* mort_age .* predictions.infections_by_age_2dose_first_pred)
    pred_mean .+= mort_reinf * vecvec_to_mat((1 - VE_severe_disease_risk[1]) .* mort_age .* predictions.infections_by_age_unvac_reinf_pred .+ (1 - VE_severe_disease_risk[2]) .* mort_age .* predictions.infections_by_age_1dose_reinf_pred .+ (1 - VE_severe_disease_risk[3]) .* mort_age .* predictions.infections_by_age_2dose_reinf_pred)
    pred_mean .*= mort_scale .* (predictions.wt_frq_pred .+ mort_var[1] .* predictions.alpha_beta_frq_pred .+ mort_var[2] .* predictions.delta_frq_pred)'
    _p_ID = normalize(p_ID .* [exp(η * t) for t = 1:length(p_ID)], 1)
    for a = 1:6
        pred_mean[a, :] = [sum(pred_mean[a, max(1, t - 99):t] .* reverse(_p_ID[1:min(t, 100)])) for t = 1:size(pred_mean, 2)]
    end
    return pred_mean
end


"""
    function pred_death_rate_county_novac(mort_scale, mort_age, mort_var, mort_reinf, predictions, VE_severe_disease_risk, η, p_ID)

Given various parameters, convert the infection predictions into a death rate predictions. Only to be used with zero vaccination rates.
"""
function pred_death_rate_county_novac(mort_scale, mort_age, mort_var, mort_reinf, predictions, η, p_ID)
    pred_mean = vecvec_to_mat(mort_age .* predictions.infections_by_age_first_pred)
    pred_mean .+= mort_reinf * vecvec_to_mat(mort_age .* predictions.infections_by_age_reinf_pred)
    pred_mean .*= mort_scale .* (predictions.wt_frq_pred .+ mort_var[1] .* predictions.alpha_beta_frq_pred .+ mort_var[2] .* predictions.delta_frq_pred)'
    _p_ID = normalize(p_ID .* [exp(η * t) for t = 1:length(p_ID)], 1)
    for a = 1:6
        pred_mean[a, :] = [sum(pred_mean[a, max(1, t - 99):t] .* reverse(_p_ID[1:min(t, 100)])) for t = 1:size(pred_mean, 2)]
    end
    return pred_mean
end


"""
    function poisson_loss_death_pred_county(county_deaths, mort_scale, mort_age, mort_var, mort_reinf, predictions, η, p_ID)

Return the negative log-likelihood of the death data given the daily predictions.        
"""
function poisson_loss_death_pred_county(county_deaths, mort_scale, mort_age, mort_var, mort_reinf, predictions, η, p_ID)
    preds = pred_death_rate_county_novac(mort_scale, mort_age, mort_var, mort_reinf, inf_predictions[1], η, p_ID)
    T = eltype(mort_scale)
    loss = T(0.0)
    for a = 1:size(county_deaths, 1), t = 31:size(county_deaths, 2)
        loss -= logpdf(Poisson(preds[a, t]), county_deaths[a, t])
    end
    return loss, sum(preds, dims = 1)[:]
end
