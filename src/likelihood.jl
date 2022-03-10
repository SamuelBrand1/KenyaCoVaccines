"""
    function simple_conv(incidence,kern)

Direct implementation of discrete convolution [incidence ⋆ kern](t).
"""
function simple_conv(incidence, kern)
    z = similar(incidence)
    n = length(incidence)
    m = length(kern)
    for t = 1:n
        z[t] = 0.0
        for s = 1:(t-1)
            if t - s <= m && t - s >= 1
                z[t] += kern[t-s] * incidence[s]
            end
        end
    end
    return z
end


"""
    function gen_logistic(x;L = 1,κ = 1, x₀ = 0)
Generalised logistic function (defaults to 'StatsFuns.logistic').
"""
function gen_logistic(x; L = 1, κ = 1, x₀ = 0)
    L / (1 + exp(-κ * (x - x₀)))
end


"""

"""
function pred_exposure_distribution_using_K(θ, model)
    @unpack init_sero, p_symp, p_sus, γA, M_county_ho, M_county_work, M_county_other = model
    @unpack β₀, β_home, β_other, β_work, E₀, ϵ, init_scale = θ

    γ = γA  # recovery from asymptomatic and symptomatic infectious states are assumed equal

    d1, = size(M_county_ho)
    K = zeros(eltype(β₀), d1, d1)
    for a = 1:d1, b = 1:d1
        K[a, b] = (1 - init_sero[a] * init_scale) * p_sus[a] * (p_symp[b] + ϵ * (1 - p_symp[b])) * β₀ * (1 / γ) * (β_home * M_county_ho[a, b] + β_other * M_county_other[a, b] + β_work * M_county_work[a, b]) # school matrix not included: schools had been closed on 1st Jan 2021
    end
    v = (K^10) * ones(eltype(β₀), d1) # 10 generations of infections used

    return v / sum(v) # normalise
end

"""
        function create_initial_conditions(initnum,N)

Create an initial conditions (on 1st Jan) to mimic conditions around time of vaccina introduction 6th march.
Fixes the initial susceptible based on round three seroprevalence by age, then
init_scale : scale that seroprevalence to a optimal initial attack rate based on data, this reduces the number of parameters
to be estiated for the intial conditions
"""
function create_initial_conditions(θ, model)
    @unpack init_sero, p_symp, p_sus, hosp_rate_by_age, α, αP, γA, γM, M_county_ho, M_county_work, M_county_other, N = model
    @unpack β₀, β_home, β_school, β_other, β_work, E₀, init_scale = θ

    γV = γM
    T = eltype(init_scale)
    u0 = zeros(T, 6, 10, 3)

    δ = p_symp
    υ = hosp_rate_by_age

    u0[:, 2, 1] .= E₀ .* KenyaCoVaccines.pred_exposure_distribution_using_K(θ, model) #β₀,β_home,β_other,β_work,E₀,ϵ,init_scale,init_sero,p_symp,p_sus,γA, M_county_ho, M_county_work, M_county_other;
    u0[:, 3, 1] .= α .* u0[:, 2, 1] .* (1 .- δ) ./ (1 + γA)
    u0[:, 4, 1] .= α .* u0[:, 2, 1] .* δ ./ (1 + αP)
    u0[:, 5, 1] .= αP .* u0[:, 4, 1] .* (1 .- υ) ./ (1 + γM)
    u0[:, 6, 1] .= αP .* u0[:, 4, 1] .* υ ./ (1 + γV)
    u0[:, 7, 1] .= N .* init_sero .* init_scale .* 0.95
    u0[:, 8, 1] .= N .* init_sero .* init_scale .* 0.05 # assume litte wanning of immunity by Jan, only 5% of those recovered (R)

    u0[:, 1, 1] .= N .- sum(u0[:, 2:8, 1], dims = [2, 3])[:]

    return u0
end

"""
    function ll_twogroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180,ct_min2 = 0.445)
Log-likelihood of the case and serology data in `model` assuming that new variant is more transmissible rather than immune evading.
"""
function ll_onegroup_newvariant_infboost(θ, model::KenyaCoVaccines.CoVAreaModel, seroreversionrate; ω = 1 / 180, p_test_scale = 5e-4) ##Scale is set so that p_test 1 is about 1% overall chance of detection
    @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other = model
    @unpack β₀, β_home, β_school, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ  #

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

    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.0
    T = eltype(β₀)

    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = convert.(T, [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3])

    u0 = convert.(T, KenyaCoVaccines.create_initial_conditions(θ, model))

    #Set up the callback for lockdown
    schoolsclosure_times = [Float64((Date(2021, 3, 19) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 7, 16) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 10, 1) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 12, 23) - Date(2020, 12, 1)).value)]
    function affect_schoolsclosed!(integrator)
        integrator.p[3] = 0.0 #Schools shut
    end
    schoolsclosed_cb = PresetTimeCallback(schoolsclosure_times, affect_schoolsclosed!)

    #Set up callback for reopening Schools
    schoolsreopening_times = [Float64((Date(2021, 1, 4) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 5, 10) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 7, 26) - Date(2020, 12, 1)).value),
        Float64((Date(2021, 10, 11) - Date(2020, 12, 1)).value)]
    function affect_schoolsopening!(integrator)
        integrator.p[3] = β_school # Schools open #
    end
    schoolsreopening_cb = PresetTimeCallback(schoolsreopening_times, affect_schoolsopening!)

    cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb)

    #Sero-waning
    sero_array = vcat(baseline_sero_array[1:30], [(1 - seroreversionrate)^k for k = 1:length(baseline_sero_array[31:end])])

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    LL = T(0.0)

    sol = solve(prob, BS3();
        tspan = (0, model_end_time),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
        verbose = false) #

    if sol.retcode == :Success
        ι = KenyaCoVaccines.get_incidence_time_array(sol, 10)  # Sum over age and vaccine doses
        ι_sero = KenyaCoVaccines.get_incidence_time_array_by_age(sol, 9)
        PCR = KenyaCoVaccines.simple_conv(ι, PCR_array)

        sero = zeros(T, size(ι_sero))
        for a = 1:size(sero, 2)
            sero[:, a] = sero_sensitivity .* KenyaCoVaccines.simple_conv(ι_sero[:, a], sero_array)
        end

        #Calculate log-likelihood for PCR testing
        for t = janstarttime:model_end_time
            if PCR_cases[model_start_time+t, 2] < 0 #Negative tests unavailable
                #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
                μ = p_test_scale * relative_testing_rate[model_start_time+t] * (p_test * PCR[t]) + 0.001
                σ² = μ + clustering_factor_PCR * μ^2
                p_negbin = 1 - (clustering_factor_PCR * μ^2 / σ²)
                r_negbin = 1 / clustering_factor_PCR
                LL += logpdf(NegativeBinomial(r_negbin, p_negbin), PCR_cases[model_start_time+t, 1])#likelihood contribution from PCR testing --- positive case detection
            else  #Negative tests available
                #p_PCR_pred = (p_test*PCR[t])/(p_test*(((χ-1)*PCR[t] + sum(N))/χ)) + 0.001
                p_PCR_pred_num = p_test_scale * p_test * PCR[t]
                p_PCR_pred_denom = p_PCR_pred_num + p_test_scale * ((p_test / χ) * (sum(N) - PCR[t]))
                p_PCR_pred = p_PCR_pred_num / p_PCR_pred_denom + 0.001
                #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
                LL += logpdf(BetaBinomial(PCR_cases[model_start_time+t, 2], p_PCR_pred * M_PCR, (1 - p_PCR_pred) * M_PCR), PCR_cases[model_start_time+t, 1])#likelihood contribution from PCR testing --- proportion postive
            end
        end

        #Calculate log-likelihood contribution from serological testing
        for t = janstarttime:(min(size(sero_cases[(model_start_time:end), :, :], 1), size(sero, 1))-1), a = 1:size(sero_cases[(model_start_time:end), :, :], 2)
            #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
            #Overall serological testing
            if sero_cases[model_start_time+t, a, 1] + sero_cases[model_start_time+t, a, 2] > 0
                p_sero_pred = (init_sero[a] * init_scale) + (sero[t, a] ./ N[a])#CAPTURE THAT INITIALLY THERE WAS SEROPOSITIVITY
                p_hat = p_sero_pred + (1 - p_sero_pred) * (1 - sero_specificity)
                LL += logpdf(BetaBinomial(sero_cases[model_start_time+t, a, 1] + sero_cases[model_start_time+t, a, 2], M_BB * p_hat, M_BB * (1 - p_hat)), sero_cases[model_start_time+t, a, 1])#Likelihood contribution from sero testing
            end
        end

    else
        LL = T(-Inf)
    end

    return LL::T
end

# """
#     function ll_twogroup_newvariant_infboost(θ,model::KenyaCoVSD.CoVAreaModel,seroreversionrate;σ = 0.16,ω = 1/180,ct_min2 = 0.445)
# Log-likelihood of the case and serology data in `model` assuming that new variant is more transmissible rather than immune evading.
# """
# function ll_onegroup_newvariant_infboost(θ, model::KenyaCoVaccines.CoVAreaModel, seroreversionrate; ω = 1 / 180, p_test_scale = 5e-4) ##Scale is set so that p_test 1 is about 1% overall chance of detection
#     @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other = model
#     @unpack β₀, β_home, β_school, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ  #

#     # Vaccine efficacies
#     ve_ac1 = VE_acquisition[1]
#     ve_ac2 = VE_acquisition[2]
#     ve_ac3 = VE_acquisition[3]
#     ve_inf1 = VE_infectiousness[1]
#     ve_inf2 = VE_infectiousness[2]
#     ve_inf3 = VE_infectiousness[3]
#     ve_dis1 = VE_severe_disease_risk[1]
#     ve_dis2 = VE_severe_disease_risk[2]
#     ve_dis3 = VE_severe_disease_risk[3]

#     #Set variance scalers
#     clustering_factor_PCR = 0.5
#     M_PCR = 30.0
#     T = eltype(β₀)

#     #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
#     p = convert.(T, [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3])

#     u0 = convert.(T, KenyaCoVaccines.create_initial_conditions(θ, model))

#     #Set up the callback for lockdown
#     schoolsclosure_times = [Float64((Date(2021, 3, 19) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 7, 16) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 10, 1) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 12, 23) - Date(2020, 12, 1)).value)]
#     function affect_schoolsclosed!(integrator)
#         integrator.p[3] = 0.0 #Schools shut
#     end
#     schoolsclosed_cb = PresetTimeCallback(schoolsclosure_times, affect_schoolsclosed!)

#     #Set up callback for reopening Schools
#     schoolsreopening_times = [Float64((Date(2021, 1, 4) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 5, 10) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 7, 26) - Date(2020, 12, 1)).value),
#         Float64((Date(2021, 10, 11) - Date(2020, 12, 1)).value)]
#     function affect_schoolsopening!(integrator)
#         integrator.p[3] = β_school # Schools open #
#     end
#     schoolsreopening_cb = PresetTimeCallback(schoolsreopening_times, affect_schoolsopening!)

#     cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb)

#     #Sero-waning
#     sero_array = vcat(baseline_sero_array[1:30], [(1 - seroreversionrate)^k for k = 1:length(baseline_sero_array[31:end])])

#     #Solve for daily incidence by age (not currently broken down by severity of disease)
#     LL = T(0.0)

#     sol = solve(prob, BS3();
#         tspan = (0, model_end_time),
#         callback = cb,
#         u0 = u0, p = p,
#         saveat = 1,
#         isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
#         verbose = false) #

#     if sol.retcode == :Success
#         ι = KenyaCoVaccines.get_incidence_time_array(sol, 10)  # Sum over age and vaccine doses
#         ι_sero = KenyaCoVaccines.get_incidence_time_array_by_age(sol, 9)
#         PCR = KenyaCoVaccines.simple_conv(ι, PCR_array)

#         sero = zeros(T, size(ι_sero))
#         for a = 1:size(sero, 2)
#             sero[:, a] = sero_sensitivity .* KenyaCoVaccines.simple_conv(ι_sero[:, a], sero_array)
#         end

#         #Calculate log-likelihood for PCR testing
#         for t = janstarttime:model_end_time
#             if PCR_cases[model_start_time+t, 2] < 0 #Negative tests unavailable
#                 #Convert from μ,α parameterisation to p,r parameterisation for using Neg. binomial model
#                 μ = p_test_scale * relative_testing_rate[model_start_time+t] * (p_test * PCR[t]) + 0.001
#                 σ² = μ + clustering_factor_PCR * μ^2
#                 p_negbin = 1 - (clustering_factor_PCR * μ^2 / σ²)
#                 r_negbin = 1 / clustering_factor_PCR
#                 LL += logpdf(NegativeBinomial(r_negbin, p_negbin), PCR_cases[model_start_time+t, 1])#likelihood contribution from PCR testing --- positive case detection
#             else  #Negative tests available
#                 #p_PCR_pred = (p_test*PCR[t])/(p_test*(((χ-1)*PCR[t] + sum(N))/χ)) + 0.001
#                 p_PCR_pred_num = p_test_scale * p_test * PCR[t]
#                 p_PCR_pred_denom = p_PCR_pred_num + p_test_scale * ((p_test / χ) * (sum(N) - PCR[t]))
#                 p_PCR_pred = p_PCR_pred_num / p_PCR_pred_denom + 0.001
#                 #Convert from p_PCR,inv_M_PCR parameterisation of Beta-Binomial to standard α,β parameterisation
#                 LL += logpdf(BetaBinomial(PCR_cases[model_start_time+t, 2], p_PCR_pred * M_PCR, (1 - p_PCR_pred) * M_PCR), PCR_cases[model_start_time+t, 1])#likelihood contribution from PCR testing --- proportion postive
#             end
#         end

#         #Calculate log-likelihood contribution from serological testing
#         for t = janstarttime:(min(size(sero_cases[(model_start_time:end), :, :], 1), size(sero, 1))-1), a = 1:size(sero_cases[(model_start_time:end), :, :], 2)
#             #Convert from the p_hat,M_BB parameterisation of the Betabinomial distribution to the α_BB,β_BB parameterisation
#             #Overall serological testing
#             if sero_cases[model_start_time+t, a, 1] + sero_cases[model_start_time+t, a, 2] > 0
#                 p_sero_pred = (init_sero[a] * init_scale) + (sero[t, a] ./ N[a])#CAPTURE THAT INITIALLY THERE WAS SEROPOSITIVITY
#                 p_hat = p_sero_pred + (1 - p_sero_pred) * (1 - sero_specificity)
#                 LL += logpdf(BetaBinomial(sero_cases[model_start_time+t, a, 1] + sero_cases[model_start_time+t, a, 2], M_BB * p_hat, M_BB * (1 - p_hat)), sero_cases[model_start_time+t, a, 1])#Likelihood contribution from sero testing
#             end
#         end

#     else
#         LL = T(-Inf)
#     end

#     return LL::T
# end

