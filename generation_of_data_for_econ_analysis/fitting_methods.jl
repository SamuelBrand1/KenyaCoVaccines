"""
    function get_daily_incidence(model::KenyaCoVaccines.CoVAreaModel)

Look into the `model` struct and for each MCMC draw of the underlying transmission parameters generate 
    a projection upto `projection_date` of the daily incidence and first incidence in each age group.                
"""
function get_daily_incidence_all_infections(model::KenyaCoVaccines.CoVAreaModel, projection_date::Date)
    println("Generating incidence predictions for $(model.areaname)")
    @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other, ω = model
    chn = model.MCMC_results.chain
    sim_duration = (projection_date - Date(2020, 12, 1)).value

    # Vaccine efficacies
    ve_ac1 = VE_acquisition[1]
    ve_ac2 = VE_acquisition[2]
    ve_ac3 = VE_acquisition[3]
    ve_inf1 = VE_infectiousness[1]
    ve_inf2 = VE_infectiousness[2]
    ve_inf3 = VE_infectiousness[3]
    ve_dis1 = VE_severe_disease_risk[1]
    ve_dis2 = VE_severe_disease_risk[2]
    ve_dis3 = VE_severe_disease_risk[3]

    #Declare arrays to put sampled incidences in
    incidence = zeros(sim_duration, size(model.prob.u0, 1), size(chn, 1))


    for k = 1:size(chn, 1)
        #For each draw from MCMC set up parameters for solving the underlying transmission model
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

        #Set up parameter vector for the model
        p = [β₀, β_home, β_school, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]

        #Set up the callback for lockdown
        schoolsclosure_times = [Float64((Date(2021, 3, 19) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 7, 16) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 10, 1) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 12, 23) - Date(2020, 12, 1)).value)]
        function affect_schoolsclosed!(integrator)
            integrator.p[3] = 0.0 #Schools shut
        end
        schoolsclosed_cb = PresetTimeCallback(schoolsclosure_times, affect_schoolsclosed!, save_positions = (false, false))

        #Set up callback for reopening Schools
        schoolsreopening_times = [Float64((Date(2021, 1, 4) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 5, 10) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 7, 26) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 10, 11) - Date(2020, 12, 1)).value)]
        function affect_schoolsopening!(integrator)
            integrator.p[3] = β_school # Schools open #
        end
        schoolsreopening_cb = PresetTimeCallback(schoolsreopening_times, affect_schoolsopening!, save_positions = (false, false))

        cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb)

        #Set up initial conditions
        θ_sample = get(model.MCMC_results.chain[k], [:β₀, :ϵ, :β_home, :β_school, :β_other, :β_work, :E₀, :init_scale]; flatten = true)
        u0 = KenyaCoVaccines.create_initial_conditions(θ_sample, model)

        #Solve model
        sol = solve(prob, BS3();
            tspan = (0, sim_duration),
            callback = cb,
            u0 = u0, p = p,
            saveat = 1,
            isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
            verbose = true)

        #Get all infections 
        incidence[:, :, k] = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age(sol, 10), 0.0, Inf)

    end

    return (; incidence)

end

# """
#     function perc_pos_to_num_tests_scale(model;lookback=21)

# Fit a linear realationship between the % PCR positive (7 day running average) and the number of PCR tests performed on each day (7-day running average),
# over the last `lookback` days of available data. Returns a `NamedTuple` with the regression coefficients and the std over error.
# """
# function perc_pos_to_num_tests_scale(model; lookback = 21)
#     prop_pos = model.PCR_cases[:, 1] ./ model.PCR_cases[:, 2]
#     prop_pos_mv_av = [mean(prop_pos[(t-3):t+3]) for t = 4:(length(prop_pos)-3)]
#     number_pos_mv_av = [mean(model.PCR_cases[:, 1][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 1])-3)]
#     number_tests_mv_av = [mean(model.PCR_cases[:, 2][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 2])-3)]

#     idxs = (.~isnan.(prop_pos_mv_av)) .& (prop_pos_mv_av .> 0) .& (1:length(prop_pos_mv_av) .>= length(prop_pos_mv_av) - lookback) .& (number_tests_mv_av .> 0)
#     xs = prop_pos_mv_av[idxs]
#     ys = number_tests_mv_av[idxs]
#     β = [ones(length(xs)) xs] \ ys
#     pred_ys = [ones(length(xs)) xs] * β
#     std_ys = std(ys .- pred_ys)
#     return (β = β, σ = std_ys)
# end

# """
#     function get_PCR_prediction(incidence_pred,model::KenyaCoVSD.CoVAreaModel;p_test_scale=5e-4,lookback = 21)

# Given a set of daily incidence predictions in `incidence_pred`, then generate predictions of observed PCR test results.
# NB: This assumes that the incidence predictions are each generated from a MCMC sample in the same order as the observation parameters.        
# """
# function get_PCR_prediction(daily_incidences_across_samples, model::KenyaCoVaccines.CoVAreaModel; p_test_scale = 5e-4, lookback = 21, M_PCR = 30.0, clustering_factor_PCR = 0.5)
#     println("Generating PCR test predictions for $(model.areaname)")
#     @unpack PCR_cases, N, PCR_array, relative_testing_rate = model
#     chn = model.MCMC_results.chain

#     #Set up arrays in memory
#     PCR_pred = zeros(size(daily_incidences_across_samples.incidence, 1), size(chn, 1))
#     prop_PCR_pred = zeros(size(daily_incidences_across_samples.incidence, 1), size(chn, 1))
#     prop_PCR_pred_uncorrected = zeros(size(daily_incidences_across_samples.incidence, 1), size(chn, 1))
#     forecast_testing_rate = zeros(size(daily_incidences_across_samples.incidence, 1), size(chn, 1))
#     test_fit = perc_pos_to_num_tests_scale(model; lookback = lookback)

#     start_idx = (Date(2020, 12, 1) - Date(2020, 2, 20)).value

#     for k = 1:size(chn, 1)
#         #Get detection rate parameter draws
#         p_test = chn[:p_test][:][k]
#         χ = chn[:χ][:][k]
#         #Get the overall incidence rate and project a PCR positive numbers time series
#         incidence = sum(daily_incidences_across_samples.incidence[:, :, k], dims = 2)[:]
#         PCR = KenyaCoVaccines.simple_conv(incidence, PCR_array)

#         #PCR test predictions
#         p_PCR_pred_num = p_test_scale .* p_test .* PCR
#         p_PCR_pred_denom = p_PCR_pred_num .+ ((p_test_scale * p_test / χ) * (sum(N) .- PCR))

#         PCR_pred[:, k] = p_PCR_pred_num
#         prop_PCR_pred[:, k] = p_PCR_pred_num ./ p_PCR_pred_denom
#         prop_PCR_pred_uncorrected[:, k] = PCR ./ sum(N)
#         forecast_testing_rate[:, k] .= [ones(length(prop_PCR_pred[:, k])) prop_PCR_pred[:, k]] * test_fit.β
#     end

#     #Construct PCR prediction = predicted %PCR pos x number of tests
#     PCR_numbers_forecast = similar(prop_PCR_pred)
#     #Differentiate between when neg PCR tests are not available and forecasting periods
#     num_tests = PCR_cases[start_idx:end, 2]
#     n = length(num_tests)
#     idxs_neg_tests = num_tests .>= 0

#     for j = 1:size(PCR_numbers_forecast, 2)
#         past = PCR_pred[1:n, j] .* (.~idxs_neg_tests) .+ prop_PCR_pred[1:n, j] .* num_tests .* (idxs_neg_tests)
#         PCR_numbers_forecast[:, j] .= vcat(past[1:(n-8)],
#             prop_PCR_pred[(n-7):end, j] .* forecast_testing_rate[(n-7):end, j])
#     end

#     ##Sample from the posterior distribution of actual cases
#     sampled_PCR_pos = zeros(Int64, size(PCR_numbers_forecast))

#     for i = 1:(n-8), j = 1:size(PCR_numbers_forecast, 2)
#         if idxs_neg_tests[i]
#             p = prop_PCR_pred[i, j] + 0.001
#             sampled_PCR_pos[i, j] = rand(BetaBinomial(num_tests[i], p * M_PCR, (1 - p) * M_PCR))
#         else
#             μ = PCR_pred[i, j] + 0.001
#             σ² = μ + clustering_factor_PCR * μ^2
#             p_negbin = 1 - (clustering_factor_PCR * μ^2 / σ²)
#             r_negbin = 1 / clustering_factor_PCR
#             sampled_PCR_pos[i, j] = rand(NegativeBinomial(r_negbin, p_negbin))
#         end
#     end

#     for i = (n-7):size(sampled_PCR_pos, 1), j = 1:size(PCR_numbers_forecast, 2)
#         projected_tests = round(Int64, mean(num_tests[idxs_neg_tests]))
#         if !isnan(test_fit.σ)
#             projected_tests = max(round(Int64, test_fit.β[1] + test_fit.β[2] * prop_PCR_pred[i, j] + test_fit.σ * randn()), 0)
#         end
#         p = prop_PCR_pred[i, j] + 0.001
#         sampled_PCR_pos[i, j] = rand(BetaBinomial(projected_tests, p * M_PCR, (1 - p) * M_PCR))
#     end

#     return (PCR_numbers_forecast = PCR_numbers_forecast,
#         sampled_PCR_pos = sampled_PCR_pos,
#         PCR_pred = PCR_pred,
#         prop_PCR_pred_uncorrected = prop_PCR_pred_uncorrected,
#         prop_PCR_pred = prop_PCR_pred,
#         forecast_testing_rate = forecast_testing_rate,
#         test_fit = test_fit)

# end


# function group_prop_pos(data; removenans = true,window = 7, last_window = true, lwr = 0.025, upr = 0.975)
#     @assert any(data .>= 0) "Negative numbers indicates that some days don't have number of tests."

#     n = size(data, 1)
#     num_pos = [sum(data[idxs, 1]) for idxs in Iterators.partition(1:n, window)]
#     num_total = [sum(data[idxs, 2]) for idxs in Iterators.partition(1:n, window)]

#     if !last_window
#         num_pos = num_pos[1:(end-1)]
#         num_total = num_total[1:(end-1)]
#     end
#     prop = num_pos ./ num_total
#     if removenans
#         prop[isnan.(prop)] .= 0.0 #Deals with division by 0 if no data available
#     end
#     lb = [p - invlogcdf(Beta(num_pos[t] + 0.5, num_total[t] - num_pos[t] + 0.5), log(lwr)) for (t, p) in enumerate(prop)]
#     ub = [invlogcdf(Beta(num_pos[t] + 0.5, num_total[t] - num_pos[t] + 0.5), log(upr)) - p for (t, p) in enumerate(prop)]

#     return (; prop, lb, ub)
# end



# ## Methods below are C&P from non-age structured version of the model
# # """
# #     function perc_pos_to_num_tests_scale(model;lookback=21)

# # Fit a linear realationship between the % PCR positive (7 day running average) and the number of PCR tests performed on each day (7-day running average),
# # over the last `lookback` days of available data. Returns a `NamedTuple` with the regression coefficients and the std over error.
# # """
# # function perc_pos_to_num_tests_scale(model; lookback = 21)
# #     prop_pos = model.PCR_cases[:, 1] ./ model.PCR_cases[:, 2]
# #     prop_pos_mv_av = [mean(prop_pos[(t-3):t+3]) for t = 4:(length(prop_pos)-3)]
# #     number_pos_mv_av = [mean(model.PCR_cases[:, 1][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 1])-3)]
# #     number_tests_mv_av = [mean(model.PCR_cases[:, 2][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 2])-3)]

# #     idxs = (.~isnan.(prop_pos_mv_av)) .& (prop_pos_mv_av .> 0) .& (1:length(prop_pos_mv_av) .>= length(prop_pos_mv_av) - lookback) .& (number_tests_mv_av .> 0)
# #     xs = prop_pos_mv_av[idxs]
# #     ys = number_tests_mv_av[idxs]
# #     β = [ones(length(xs)) xs] \ ys
# #     pred_ys = [ones(length(xs)) xs] * β
# #     std_ys = std(ys .- pred_ys)
# #     return (β = β, σ = std_ys)
# # end

# # """
# #     function get_unscaled_predictions(ι₁,ι₂,p_ID,rel_test_rate;fitlength=430)

# # Useful method for collating a matrix of of unnormalised (i.e. IFR not applied) estimates of the daily death rate from the model incidence estimates.
# # """
# # function get_unscaled_predictions(ι₁, ι₂, p_ID, rel_test_rate; fitlength = 430)
# #     unscaled_deaths1 = similar(ι₁)[1:fitlength, :]
# #     unscaled_deaths2 = similar(ι₂)[1:fitlength, :]
# #     for j = 1:size(unscaled_deaths1, 2)
# #         unscaled_deaths1[:, j] .= KenyaCoVSD.simple_conv(ι₁[:, j], p_ID)[1:fitlength] .* rel_test_rate[1:fitlength]
# #         unscaled_deaths2[:, j] .= KenyaCoVSD.simple_conv(ι₂[:, j], p_ID)[1:fitlength] .* rel_test_rate[1:fitlength]
# #     end
# #     return unscaled_deaths1, unscaled_deaths2
# # end


# # """
# #     function get_daily_incidence(model::KenyaCoVSD.CoVAreaModel,
# #                         projection_date::Date;
# #                         σ = 0.16,
# #                         ω = 1/180,
# #                         ct_min2 = 0.555)

# # Look into the `model` struct and for each MCMC draw of the underlying transmission parameters generate 
# #     a projection upto `projection_date` of the daily incidence and first incidence in each SES group.                
# # """
# # function get_daily_incidence(model::KenyaCoVSD.CoVAreaModel,
# #     projection_date::Date;
# #     σ = 0.16,
# #     ω = 1 / 180,
# #     ct_min2 = 0.555)
# #     println("Generating incidence predictions for $(model.areaname)")
# #     @unpack N, M_BB, prob, α, γ = model
# #     MCMCdraws = get(model.MCMC_results.chain, [:ct_min1, :R₀, :ϵ, :χ₁, :χ₂, :p_test₁, :p_test₂, :P_eff, :schooleffect, :Rα, :time_scale_αvar, :mid_point_αvar, :Rδ, :time_scale_δvar, :mid_point_δvar, :E₀])
# #     incidence₁ = zeros((projection_date - Date(2020, 2, 20)).value, length(first(MCMCdraws)))
# #     incidence₂ = zeros((projection_date - Date(2020, 2, 20)).value, length(first(MCMCdraws)))
# #     firstincidence₁ = zeros((projection_date - Date(2020, 2, 20)).value, length(first(MCMCdraws)))
# #     firstincidence₂ = zeros((projection_date - Date(2020, 2, 20)).value, length(first(MCMCdraws)))

# #     for k = 1:length(first(MCMCdraws))
# #         #For each draw from MCMC set up parameters for solving the underlying transmission model
# #         N₁ = MCMCdraws.P_eff[k] * N
# #         N₂ = N - N₁
# #         p = [MCMCdraws.R₀[k], model.α, model.γ, MCMCdraws.ϵ[k], σ, N₁, N₂, ω, MCMCdraws.schooleffect[k], MCMCdraws.ct_min1[k], ct_min2, MCMCdraws.Rα[k], MCMCdraws.time_scale_αvar[k], 347.0 + MCMCdraws.mid_point_αvar[k] * 30.0, MCMCdraws.Rδ[k], MCMCdraws.time_scale_δvar[k], 406.0 + MCMCdraws.mid_point_δvar[k] * 30.0]
# #         u0 = [N₁, MCMCdraws.E₀[k], 0.0, 0.0, 0.0, 0.0, N₂, MCMCdraws.E₀[k], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# #         #Solve model
# #         sol = solve(prob, BS3(); tspan = (0, (projection_date - Date(2020, 2, 20)).value),
# #             u0 = u0,
# #             p = p,
# #             saveat = 1,
# #             isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
# #             force_dtmin = true)

# #         ι₁ = clamp!(diff(sol[:C₁]), 0.0, Inf)
# #         ι₂ = clamp!(diff(sol[:C₂]), 0.0, Inf)
# #         ι_sero₁ = clamp!(diff(sol[:C_sero₁]), 0.0, Inf)
# #         ι_sero₂ = clamp!(diff(sol[:C_sero₂]), 0.0, Inf)

# #         incidence₁[:, k] = ι₁
# #         incidence₂[:, k] = ι₂
# #         firstincidence₁[:, k] = ι_sero₁
# #         firstincidence₂[:, k] = ι_sero₂

# #     end

# #     return (incidence₁ = incidence₁, incidence₂ = incidence₂,
# #         firstincidence₁ = firstincidence₁, firstincidence₂ = firstincidence₂)

# # end

# # """
# #     function get_PCR_prediction(incidence_pred,model::KenyaCoVSD.CoVAreaModel;p_test_scale=5e-4,lookback = 21)

# # Given a set of daily incidence predictions in `incidence_pred`, then generate predictions of observed PCR test results.
# # NB: This assumes that the incidence predictions are each generated from a MCMC sample in the same order as the observation parameters.        
# # """
# # function get_PCR_prediction(incidence_pred, model::KenyaCoVSD.CoVAreaModel; p_test_scale = 5e-4, lookback = 21, M_PCR = 30.0, clustering_factor_PCR = 0.5)
# #     println("Generating PCR test predictions for $(model.areaname)")
# #     @unpack PCR_cases, N, PCR_array, relative_testing_rate = model
# #     MCMCdraws = get(model.MCMC_results.chain, [:ct_min1, :R₀, :ϵ, :χ₁, :χ₂, :p_test₁, :p_test₂, :P_eff, :schooleffect, :Rα, :time_scale_αvar, :mid_point_αvar, :Rδ, :time_scale_δvar, :mid_point_δvar, :E₀])

# #     PCR_pred = zeros(size(incidence_pred.incidence₁, 1), length(first(MCMCdraws)))
# #     prop_PCR_pred = zeros(size(incidence_pred.incidence₁, 1), length(first(MCMCdraws)))
# #     forecast_testing_rate = zeros(size(incidence_pred.incidence₁, 1), length(first(MCMCdraws)))
# #     test_fit = perc_pos_to_num_tests_scale(model; lookback = lookback)

# #     for k = 1:length(first(MCMCdraws))
# #         N₁ = MCMCdraws.P_eff[k] * N
# #         N₂ = N - N₁
# #         ι₁ = incidence_pred.incidence₁[:, k]
# #         ι₂ = incidence_pred.incidence₂[:, k]
# #         ι_sero₁ = incidence_pred.firstincidence₁[:, k]
# #         ι_sero₂ = incidence_pred.firstincidence₂[:, k]
# #         #Generate implied number of PCR-positive in each group
# #         PCR₁ = KenyaCoVSD.simple_conv(ι₁, PCR_array)
# #         PCR₂ = KenyaCoVSD.simple_conv(ι₂, PCR_array)
# #         #PCR test predictions
# #         PCR_pred[:, k] = p_test_scale * relative_testing_rate[1:length(ι₁)] .* (MCMCdraws.p_test₁[k] .* PCR₁ .+ MCMCdraws.p_test₂[k] .* PCR₂)
# #         p_PCR_pred_num = p_test_scale .* (MCMCdraws.p_test₁[k] .* PCR₁ .+ MCMCdraws.p_test₂[k] .* PCR₂)
# #         p_PCR_pred_denom = p_PCR_pred_num .+ p_test_scale .* ((MCMCdraws.p_test₁[k] / MCMCdraws.χ₁[k]) .* (N₁ .- PCR₁) .+ (MCMCdraws.p_test₂[k] / MCMCdraws.χ₂[k]) .* (N₂ .- PCR₂))
# #         prop_PCR_pred[:, k] .= p_PCR_pred_num ./ p_PCR_pred_denom
# #         #PCR forecast accounting for test rates in last lookback days of period
# #         forecast_testing_rate[:, k] .= [ones(length(prop_PCR_pred[:, k])) prop_PCR_pred[:, k]] * test_fit.β
# #     end
# #     #Construct PCR prediction = predicted %PCR pos x number of tests
# #     PCR_forecast = similar(prop_PCR_pred)
# #     #Differentiate between when neg PCR tests are not available and forecasting periods
# #     num_tests = PCR_cases[:, 2]
# #     n = length(num_tests)
# #     idxs_neg_tests = num_tests .>= 0

# #     for j = 1:size(PCR_forecast, 2)
# #         past = PCR_pred[1:n, j] .* (.~idxs_neg_tests) .+ prop_PCR_pred[1:n, j] .* num_tests .* (idxs_neg_tests)
# #         PCR_forecast[:, j] .= vcat(past[1:(n-8)],
# #             prop_PCR_pred[(n-7):end, j] .* forecast_testing_rate[(n-7):end, j])
# #     end

# #     ##Sample from the posterior distribution of actual cases
# #     sampled_PCR_pos = zeros(Int64, size(PCR_forecast))

# #     for i = 1:(n-8), j = 1:size(PCR_forecast, 2)
# #         if idxs_neg_tests[i]
# #             p = prop_PCR_pred[i, j] + 0.01
# #             sampled_PCR_pos[i, j] = rand(BetaBinomial(num_tests[i], p * M_PCR, (1 - p) * M_PCR))
# #         else
# #             μ = PCR_pred[i, j] + 0.001
# #             σ² = μ + clustering_factor_PCR * μ^2
# #             p_negbin = 1 - (clustering_factor_PCR * μ^2 / σ²)
# #             r_negbin = 1 / clustering_factor_PCR
# #             sampled_PCR_pos[i, j] = rand(NegativeBinomial(r_negbin, p_negbin))
# #         end
# #     end

# #     for i = (n-7):size(sampled_PCR_pos, 1), j = 1:size(PCR_forecast, 2)
# #         projected_tests = round(Int64, mean(num_tests[idxs_neg_tests]))
# #         if !isnan(test_fit.σ)
# #             projected_tests = max(round(Int64, test_fit.β[1] + test_fit.β[2] * prop_PCR_pred[i, j] + test_fit.σ * randn()), 0)
# #         end
# #         p = prop_PCR_pred[i, j] + 0.01
# #         sampled_PCR_pos[i, j] = rand(BetaBinomial(projected_tests, p * M_PCR, (1 - p) * M_PCR))
# #     end

# #     return (PCR_forecast = PCR_forecast,
# #         sampled_PCR_pos = sampled_PCR_pos,
# #         PCR_pred = PCR_pred,
# #         prop_PCR_pred = prop_PCR_pred,
# #         forecast_testing_rate = forecast_testing_rate,
# #         test_fit = test_fit)

# # end


# # """
# #     function get_variant_prediction(incidence_pred,model)

# # Given a set of daily incidence predictions in `incidence_pred`, then generate predictions of proportion infected by each variant.
# # NB: This assumes that the incidence predictions are each generated from a MCMC sample in the same order as the observation parameters.   
# # """
# # function get_variant_prediction(incidence_pred, model)
# #     println("Generating variant frequency predictions for $(model.areaname)")
# #     MCMCdraws = get(model.MCMC_results.chain, [:ct_min1, :R₀, :ϵ, :χ₁, :χ₂, :p_test₁, :p_test₂, :P_eff, :schooleffect, :Rα, :time_scale_αvar, :mid_point_αvar, :Rδ, :time_scale_δvar, :mid_point_δvar, :E₀])

# #     prop_inc_delta = [KenyaCoVSD.gen_logistic(t; κ = MCMCdraws.time_scale_δvar[k], x₀ = 406.0 + MCMCdraws.mid_point_δvar[k] * 30.0) for t = 1:size(incidence_pred.incidence₁, 1), k = 1:size(incidence_pred.incidence₁, 2)]
# #     prop_inc_alpha_beta = (1 .- prop_inc_delta) .* [KenyaCoVSD.gen_logistic(t; κ = MCMCdraws.time_scale_αvar[k], x₀ = 406.0 + MCMCdraws.mid_point_αvar[k] * 30.0) for t = 1:size(incidence_pred.incidence₁, 1), k = 1:size(incidence_pred.incidence₁, 2)]
# #     prop_inc_wt = (1 .- prop_inc_delta) .* (1 .- [KenyaCoVSD.gen_logistic(t; κ = MCMCdraws.time_scale_αvar[k], x₀ = 406.0 + MCMCdraws.mid_point_αvar[k] * 30.0) for t = 1:size(incidence_pred.incidence₁, 1), k = 1:size(incidence_pred.incidence₁, 2)])

# #     return (prop_inc_wt = prop_inc_wt, prop_inc_alpha_beta = prop_inc_alpha_beta, prop_inc_delta = prop_inc_delta)
# # end

# # """
# #     function get_unscaled_deaths_predictions(incidence_pred,variant_pred,p_ID,rel_test_rate)

# #     Given a set of daily incidence predictions in `incidence_pred`, then generate predictions of unscaled trend of a delay with delay disttibution `p_ID` by each variant.
# # NB: This assumes that the incidence predictions are each generated from a MCMC sample in the same order as the observation parameters.   
# # """
# # function get_unscaled_delayed_predictions(incidence_pred, variant_pred, p_ID, rel_test_rate)
# #     println("Generating unscaled delayed trend predictions for $(model.areaname)")

# #     unscaled_trend1_wt = incidence_pred.incidence₁ .* variant_pred.prop_inc_wt
# #     unscaled_trend1_alpha_beta = incidence_pred.incidence₁ .* variant_pred.prop_inc_alpha_beta
# #     unscaled_trend1_delta = incidence_pred.incidence₁ .* variant_pred.prop_inc_delta
# #     unscaled_trend2_wt = incidence_pred.incidence₂ .* variant_pred.prop_inc_wt
# #     unscaled_trend2_alpha_beta = incidence_pred.incidence₂ .* variant_pred.prop_inc_alpha_beta
# #     unscaled_trend2_delta = incidence_pred.incidence₂ .* variant_pred.prop_inc_delta

# #     rel_test_rate_trunc = copy(rel_test_rate)[1:size(unscaled_trend1_wt, 1)]

# #     for j = 1:size(unscaled_trend1_wt, 2)
# #         unscaled_trend1_wt[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend1_wt[:, j], p_ID) .* rel_test_rate_trunc
# #         unscaled_trend1_alpha_beta[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend1_alpha_beta[:, j], p_ID) .* rel_test_rate_trunc
# #         unscaled_trend1_delta[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend1_delta[:, j], p_ID) .* rel_test_rate_trunc
# #         unscaled_trend2_wt[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend2_wt[:, j], p_ID) .* rel_test_rate_trunc
# #         unscaled_trend2_alpha_beta[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend2_alpha_beta[:, j], p_ID) .* rel_test_rate_trunc
# #         unscaled_trend2_delta[:, j] .= KenyaCoVSD.simple_conv(unscaled_trend2_delta[:, j], p_ID) .* rel_test_rate_trunc
# #     end

# #     return (unscaled_trend1_wt = unscaled_trend1_wt,
# #         unscaled_trend1_alpha_beta = unscaled_trend1_alpha_beta,
# #         unscaled_trend1_delta = unscaled_trend1_delta,
# #         unscaled_trend2_wt = unscaled_trend2_wt,
# #         unscaled_trend2_alpha_beta = unscaled_trend2_alpha_beta,
# #         unscaled_trend2_delta = unscaled_trend2_delta)
# # end

# # function mvav_cols(X)
# #     _X = zeros(size(X))
# #     for j = 1:size(X, 2)
# #         _X[:, j] .= [X[1:3, j]
# #             [mean(X[(t-3):(t+3), j]) for t = 4:(size(X, 1)-3)]
# #             X[(end-2):end, j]]
# #     end
# #     return _X
# # end

# # function get_credible_intervals(X)
# #     pred = mean(X, dims = 2)
# #     lb = pred .- [quantile(X[t, :], 0.025) for t = 1:size(X, 1)]
# #     ub = [quantile(X[t, :], 0.975) for t = 1:size(X, 1)] .- pred
# #     return (pred = pred, lb = lb, ub = ub)
# # end

# # """
# #     function sampled_dev_deaths_data_for_county(μ,μ_var,unscaled_pred_wt,unscaled_pred_alpha_beta,unscaled_pred_delta,reported_deaths)

# # Return Poisson deviance (neg. log-likelihood) of county with IFR `μ` and variant risk scalars `μ_var`, and precalculated unscaled trends. Deviance is sampled from posterior distribution of unscaled trends.        
# # """
# # function sampled_dev_deaths_data_for_county(μ, μ_var, unscaled_pred_wt, unscaled_pred_alpha_beta, unscaled_pred_delta, reported_deaths)
# #     sample = rand(1:size(unscaled_pred_wt, 2))
# #     pred = μ .* μ_var[1] .* unscaled_pred_wt[:, sample]
# #     pred .+= μ .* μ_var[2] .* unscaled_pred_alpha_beta[:, sample]
# #     pred .+= μ .* μ_var[3] .* unscaled_pred_delta[:, sample]
# #     T = eltype(μ_var)
# #     deaths_deviance = T(0)
# #     for t = 1:length(reported_deaths)
# #         deaths_deviance -= logpdf(Poisson(pred[t] .+ 0.001), reported_deaths[t])
# #     end
# #     return deaths_deviance, pred
# # end

# # """
# #     function sampled_dev_deaths_data_for_county(μ,μ_var,unscaled_pred_wt,unscaled_pred_alpha_beta,unscaled_pred_delta,reported_deaths,clustering_factor)

# # Return Neg. Bin. deviance (neg. log-likelihood) of county with IFR `μ` and variant risk scalars `μ_var`, and precalculated unscaled trends. Deviance is sampled from posterior distribution of unscaled trends.        
# # """
# # function sampled_dev_deaths_data_for_county(μ, μ_var, unscaled_pred_wt, unscaled_pred_alpha_beta, unscaled_pred_delta, reported_deaths, clustering_factor)
# #     sample = rand(1:size(unscaled_pred_wt, 2))
# #     pred = μ .* μ_var[1] .* unscaled_pred_wt[:, sample]
# #     pred .+= μ .* μ_var[2] .* unscaled_pred_alpha_beta[:, sample]
# #     pred .+= μ .* μ_var[3] .* unscaled_pred_delta[:, sample]
# #     T = eltype(μ_var)
# #     deaths_deviance = T(0)
# #     for t = 1:length(reported_deaths)
# #         neg_bin_mean = pred[t] +0.01
# #         σ² = neg_bin_mean + clustering_factor * neg_bin_mean^2
# #         p_negbin = 1 - (clustering_factor * neg_bin_mean^2 / σ²)
# #         r_negbin = 1 / clustering_factor
# #         deaths_deviance -= logpdf(NegativeBinomial(r_negbin, p_negbin), reported_deaths[t])
# #     end
# #     return deaths_deviance, pred
# # end


# # """
# #     function sampled_dev_deaths_data(μ_vect,μ_var,unscaled_pred_wt_vect,unscaled_pred_alpha_beta_vect,unscaled_pred_delta_vect,reported_deaths_vect)

# # Return Poisson deviance (neg. log-likelihood) of all Kenya with IFRs `μ_vect` and variant risk scalars `μ_var`, and precalculated unscaled trends. Deviance is sampled from posterior distribution of unscaled trends.        
# # """
# # function sampled_dev_deaths_data(μ_vect, μ_var, unscaled_pred_wt_vect, unscaled_pred_alpha_beta_vect, unscaled_pred_delta_vect, reported_deaths_vect)
# #     deaths_deviance, pred = sampled_dev_deaths_data_for_county(μ_vect[1], μ_var, unscaled_pred_wt_vect[1], unscaled_pred_alpha_beta_vect[1], unscaled_pred_delta_vect[1], reported_deaths_vect[1])
# #     for k = 2:length(unscaled_pred_wt_vect)
# #         dev, _pred = sampled_dev_deaths_data_for_county(μ_vect[k], μ_var, unscaled_pred_wt_vect[k], unscaled_pred_alpha_beta_vect[k], unscaled_pred_delta_vect[k], reported_deaths_vect[k])
# #         deaths_deviance += dev
# #         pred .+= _pred
# #     end
# #     return deaths_deviance, pred
# # end

# # """
# #     function sampled_dev_deaths_data(μ_vect,μ_var,unscaled_pred_wt_vect,unscaled_pred_alpha_beta_vect,unscaled_pred_delta_vect,reported_deaths_vect,clustering_factor)

# # Return Neg bin. deviance (neg. log-likelihood) of all Kenya with IFRs `μ_vect` and variant risk scalars `μ_var`, and precalculated unscaled trends. Deviance is sampled from posterior distribution of unscaled trends.        
# # """
# # function sampled_dev_deaths_data(μ_vect, μ_var, unscaled_pred_wt_vect, unscaled_pred_alpha_beta_vect, unscaled_pred_delta_vect, reported_deaths_vect, clustering_factor)
# #     deaths_deviance, pred = sampled_dev_deaths_data_for_county(μ_vect[1], μ_var, unscaled_pred_wt_vect[1], unscaled_pred_alpha_beta_vect[1], unscaled_pred_delta_vect[1], reported_deaths_vect[1], clustering_factor)
# #     for k = 2:length(unscaled_pred_wt_vect)
# #         dev, _pred = sampled_dev_deaths_data_for_county(μ_vect[k], μ_var, unscaled_pred_wt_vect[k], unscaled_pred_alpha_beta_vect[k], unscaled_pred_delta_vect[k], reported_deaths_vect[k], clustering_factor)
# #         deaths_deviance += dev
# #         pred .+= _pred
# #     end
# #     return deaths_deviance, pred
# # end