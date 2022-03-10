"""
    function flatten_chain(chn::Chains)

Provide a vector of all MCMC draws as `NamedTuple` objects.
"""
function flatten_chain(chn::Chains)
    param_names = names(chn)
    [get(chn[i, :, k], param_names; flatten = true) for i = 1:size(chn, 1), k = 1:size(chn, 3)][:]
end

"""
    function flatten_chain(chn::Chains,n::Integer)

Provide a vector of `n` sampled MCMC draws as `NamedTuple` objects.
"""
function flatten_chain(chn::Chains, n::Integer)
    sampled_chn = sample(chn, n)
    param_names = names(chn)
    [get(sampled_chn[i, :, k], param_names; flatten = true) for i = 1:size(sampled_chn, 1), k = 1:size(sampled_chn, 3)][:]
end

"""
    function flatten_chain(model::KenyaCoVaccines.CoVAreaModel)

Provide a vector of all MCMC draws as `NamedTuple` objects from a `KenyaCoVaccines.CoVAreaModel` object
"""
function flatten_chain(model::KenyaCoVaccines.CoVAreaModel)
    flatten_chain(model.MCMC_results.chain)
end

"""
    function flatten_chain(model::KenyaCoVaccines.CoVAreaModel,n::Integer)

Provide a vector of `n` sampled MCMC draws as `NamedTuple` objects from a `KenyaCoVaccines.CoVAreaModel` object
"""
function flatten_chain(model::KenyaCoVaccines.CoVAreaModel, n::Integer)
    flatten_chain(model.MCMC_results.chain, n)
end

"""
    project_daily_infections_all(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)

Return total infections as a timeseries, this sums over age and vaccine status. Simulation use the parameters in `θ` with waning immunity rate `ω`.    
"""
function project_daily_infections_all(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, BS3(); abstol = 1e-3,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        # isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
        verbose = true)
    ι = clamp!(KenyaCoVaccines.get_incidence_time_array(sol, 10), 0.0, Inf)  # Sum over age and vaccine doses
    return ι
end

"""
    project_daily_infections_by_age_and_vaccine(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)

Return total infections as a timeseries, this sums over age and vaccine status. Simulation use the parameters in `θ` with waning immunity rate `ω`.    
"""
function project_daily_infections_by_age_and_vaccine(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, AutoTsit5(Rosenbrock23()); abstol = 1e-3,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        # isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u),
        verbose = true)

    ι_first = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol, 9), 0.0, Inf)  # Sum over age and vaccine doses
    ι_all = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol, 10), 0.0, Inf)
    return ι_first, ι_all
end

"""
    function project_daily_infections_by_age_and_vaccine_immune_escape(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date,immune_esc_cb; ω = 1 / 180)

Project daily infections, broken down by age and vaccination status, in a structure that can incoporate immune escape variants.

The specific effect of the immune escape variant is encoded in the 
"""
function project_daily_infections_by_age_and_vaccine_immune_escape(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date, immune_esc_cb, cb_vac; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3, 1.0, 0.0, 0.0]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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
    cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb, immune_esc_cb, cb_vac)

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, CVODE_BDF(linear_solver = :GMRES);# abstol = 1e-15,reltol = 1e-15,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        # isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u),
        verbose = true,
        force_dtmin = true)
    # global solves
    # solves += 1
    # println(solves)        
    # sol = solve(prob, Tsit5();tspan = (0, sim_duration),AutoTsit5(Rodas5())
    #     callback = cb,
    #     u0 = u0, p = p,
    #     saveat = 1,
    #     # isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u),
    #     verbose = true)

    ι_first = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol, 9), 0.0, Inf)  # Sum over age and vaccine doses
    ι_all = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol, 10), 0.0, Inf)
    if size(ι_first, 1) < 912
        ι_first = vcat(ι_first, zeros(912 - size(ι_first, 1), 6, 3))
    end
    if size(ι_all, 1) < 912
        ι_all = vcat(ι_first, zeros(912 - size(ι_first, 1), 6, 3))
    end


    return ι_first, ι_all
end



"""
    function solve_model(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)

Return `ODESolution` object. Simulation use the parameters in `θ` with waning immunity rate `ω`.
NB: this saves every aspect of the simulation, and, therefore, costs a large amount of memory is not used temporarily.
"""
function solve_model(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, AutoTsit5(Rosenbrock23()); abstol = 1e-3,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u),
        verbose = true)

    return sol

end

function solve_model_immune_escape(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date, immune_esc_cb, vac_cb, solver; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3, 1.0, 0.0, 0.0]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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
    # vaccine_rollout_days = PresetTimeCallback([258, 314, 889, 945], init -> 0, save_positions = (false, false))
    cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb, immune_esc_cb, vac_cb)

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, solver; abstol = 1e-3, reltol = 1e-3,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        # isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u),
        verbose = true)

    return sol

end


"""
    function project_daily_infections_first_by_age(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)

Return total first (i.e. no reinfections) infections as a timeseries, this sums over vaccine status, but seperates by age. Simulation use the parameters in `θ` with waning immunity rate `ω`.    
"""
function project_daily_infections_first_by_age(θ, model::KenyaCoVaccines.CoVAreaModel, projection_date::Date; ω = 1 / 180)
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
    #Set transmission parameters --- SCHOOLS START SHUT IN THIS MODEL
    p = [β₀, β_home, 0.0, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, alpha_variant_time + mid_point_αvar * 30, inc_R_δvar, time_scale_δvar, delta_variant_time + mid_point_δvar * 30, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3]
    u0 = KenyaCoVaccines.create_initial_conditions(θ, model)
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

    #Solve for daily incidence by age (not currently broken down by severity of disease)
    sim_duration = (projection_date - Date(2020, 12, 1)).value
    sol = solve(prob, BS3(); abstol = 1e-3,
        tspan = (0, sim_duration),
        callback = cb,
        u0 = u0, p = p,
        saveat = 1,
        isoutofdomain = (u, p, t) -> any(x -> x < 0, u),
        verbose = false)
    ι = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age(sol, 9), 0.0, Inf)  # Sum over vaccine doses
    return ι
end


"""
    function perc_pos_to_num_tests_scale(model;lookback=21)

Fit a linear realationship between the % PCR positive (7 day running average) and the number of PCR tests performed on each day (7-day running average),
over the last `lookback` days of available data. Returns a `NamedTuple` with the regression coefficients and the std over error.
"""
function perc_pos_to_num_tests_scale(model::KenyaCoVaccines.CoVAreaModel; lookback = 21)
    prop_pos = model.PCR_cases[:, 1] ./ model.PCR_cases[:, 2]
    prop_pos_mv_av = [mean(prop_pos[(t-3):t+3]) for t = 4:(length(prop_pos)-3)]
    number_pos_mv_av = [mean(model.PCR_cases[:, 1][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 1])-3)]
    number_tests_mv_av = [mean(model.PCR_cases[:, 2][(t-3):t+3]) for t = 4:(length(model.PCR_cases[:, 2])-3)]

    idxs = (.~isnan.(prop_pos_mv_av)) .& (prop_pos_mv_av .> 0) .& (1:length(prop_pos_mv_av) .>= length(prop_pos_mv_av) - lookback) .& (number_tests_mv_av .> 0)
    xs = prop_pos_mv_av[idxs]
    ys = number_tests_mv_av[idxs]
    β = [ones(length(xs)) xs] \ ys
    pred_ys = [ones(length(xs)) xs] * β
    std_ys = std(ys .- pred_ys)
    return (β = β, σ = std_ys)
end


function generate_PCR_predictions(ι::AbstractVector, θ, model::KenyaCoVaccines.CoVAreaModel, fit; M_PCR = 30.0, p_test_scale = 5e-4)
    @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other = model
    @unpack β₀, β_home, β_school, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ  #

    PCR = KenyaCoVaccines.simple_conv(ι, PCR_array)
    p_PCR_pred_num = p_test_scale .* p_test .* PCR
    p_PCR_pred_denom = p_PCR_pred_num .+ p_test_scale .* ((p_test / χ) .* (sum(N) .- PCR))
    p_PCR_pred = (p_PCR_pred_num ./ p_PCR_pred_denom) .+ 0.001
    _PCR_cases = [PCR_cases[model_start_time:end, :]; fill(-1, length(ι) - size(PCR_cases[model_start_time:end, :], 1), 2)]
    n_tests = _PCR_cases[:, 2]
    for i = 1:length(n_tests)
        if n_tests[i] < 0
            n_tests[i] = round(Int64, fit.β[1] + fit.β[2] * p_PCR_pred[i] + fit.σ * randn())
        end
    end
    clamp!(n_tests, 0, Inf)
    sampled_PCR_pos_tests = [rand(BetaBinomial(n_tests[i], p_PCR_pred[i] * M_PCR, (1.0 - p_PCR_pred[i]) * M_PCR)) for i = 1:length(n_tests)]
    return [sampled_PCR_pos_tests n_tests]
end

function generate_sero_predictions(ι_sero::AbstractVector, θ, age_grp, model::KenyaCoVaccines.CoVAreaModel; seroreversionrate = 0.0)
    @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other = model
    @unpack β₀, β_home, β_school, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ  #

    sero_array = vcat(baseline_sero_array[1:30], [(1 - seroreversionrate)^k for k = 1:600])
    rand_sensitivity = rand(Beta(M_BB * sero_sensitivity, M_BB * (1 - sero_sensitivity)))
    unscaled_sero = (init_sero[age_grp] * init_scale) .+ KenyaCoVaccines.simple_conv(ι_sero, sero_array) ./ N[age_grp]
    sero = rand_sensitivity .* unscaled_sero .+ (1 - sero_specificity) .* (1 .- unscaled_sero)

    return sero
end

"""
    function get_variant_prediction(θ, model, n::Integer)

For times `1:n` give the model prediction of frequency of different variants conditional on parameter `θ`.
"""
function get_variant_prediction(θ, model, n::Integer)
    @unpack PCR_cases, sero_cases, init_sero, baseline_sero_array, PCR_array, sero_sensitivity, sero_specificity, p_symp, p_sus, hosp_rate_by_age, N, M_BB, prob, α, αP, γA, γM, relative_testing_rate, alpha_variant_time, delta_variant_time, model_start_time, model_end_time, janstarttime, VE_acquisition, VE_infectiousness, VE_severe_disease_risk, M_county_ho, M_county_school, M_county_work, M_county_other = model
    @unpack β₀, β_home, β_school, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ  #

    prop_inc_delta = [KenyaCoVaccines.gen_logistic(t; κ = time_scale_δvar, x₀ = delta_variant_time + mid_point_δvar * 30) for t = 1:n]
    prop_inc_alpha_beta = (1 .- prop_inc_delta) .* [KenyaCoVaccines.gen_logistic(t; κ = time_scale_αvar, x₀ = alpha_variant_time + mid_point_αvar * 30) for t = 1:n]
    prop_inc_wt = (1 .- prop_inc_delta) .* (1 .- [KenyaCoVaccines.gen_logistic(t; κ = time_scale_αvar, x₀ = alpha_variant_time + mid_point_αvar * 30) for t = 1:n])

    return (prop_inc_wt = prop_inc_wt, prop_inc_alpha_beta = prop_inc_alpha_beta, prop_inc_delta = prop_inc_delta)
end

"""
    function weekly_proportions_with_errors(X::AbstractMatrix; lwr = 0.025, upr = 0.975)

Generate a time series of weekly proportions of pos/neg test results.        
"""
function weekly_proportions_with_errors(X::AbstractMatrix; lwr = 0.025, upr = 0.975)
    @assert size(X, 2) == 2 "This method only works for two column matrices: first column being positive tests and second column being negative tests"

    weekly_pos_tests = [sum(X[grp, 1]) for grp in Iterators.partition(1:size(X, 1), 7)]
    weekly_neg_tests = [sum(X[grp, 2]) for grp in Iterators.partition(1:size(X, 1), 7)]
    weekly_total_tests = weekly_pos_tests + weekly_neg_tests

    ts = collect(0:7:(7*(length(weekly_total_tests)-1)))
    prop = [weekly_total_tests[i] > 0 ? weekly_pos_tests[i] / weekly_total_tests[i] : 0.0 for i = 1:length(weekly_total_tests)]
    lb = [weekly_total_tests[i] > 0 ? prop[i] - invlogcdf(Beta(0.5 + weekly_pos_tests[i], 0.5 + weekly_neg_tests[i]), log(lwr)) : prop[i] for i = 1:length(weekly_total_tests)]
    ub = [weekly_total_tests[i] > 0 ? invlogcdf(Beta(0.5 + weekly_pos_tests[i], 0.5 + weekly_neg_tests[i]), log(upr)) - prop[i] : 1 - prop[i] for i = 1:length(weekly_total_tests)]
    idxs = weekly_total_tests .> 0
    return (ts = ts[idxs], prop = prop[idxs], lb = lb[idxs], ub = ub[idxs])
end

"""
    function monthly_proportions_with_errors(X::AbstractMatrix, dates::Vector{Date}, date0::Date; lwr = 0.025, upr = 0.975)

Generate a time series of monthly proportions of pos/neg test results.          
"""
function monthly_proportions_with_errors(X::AbstractMatrix, dates::Vector{Date}, date0::Date; lwr = 0.025, upr = 0.975)
    @assert size(X, 2) == 2 "This method only works for two column matrices: first column being positive tests and second column being negative tests"

    samplemonths = month.(dates)
    samplyears = year.(dates)

    ts = Float64[]
    monthly_pos_tests = Int64[]
    monthly_neg_tests = Int64[]

    for y = minimum(samplyears):maximum(samplyears)
        for m = minimum(samplemonths):maximum(samplemonths)
            idxs = (samplemonths .== m) .& (samplyears .== y)
            if any(idxs)
                day1 = Date(y, m, 1)
                dayend = lastdayofmonth(day1)
                push!(ts, mean([(day1 - date0).value, (dayend - date0).value]))
                push!(monthly_pos_tests, sum(X[idxs, 1]))
                push!(monthly_neg_tests, sum(X[idxs, 2]))
            end
        end
    end

    monthly_total_tests = monthly_pos_tests + monthly_neg_tests
    idxs = monthly_total_tests .> 0
    prop = [monthly_total_tests[i] > 0 ? monthly_pos_tests[i] / monthly_total_tests[i] : 0.0 for i = 1:length(monthly_total_tests)]
    lb = [monthly_total_tests[i] > 0 ? prop[i] - invlogcdf(Beta(0.5 + monthly_pos_tests[i], 0.5 + monthly_neg_tests[i]), log(lwr)) : prop[i] for i = 1:length(monthly_total_tests)]
    ub = [monthly_total_tests[i] > 0 ? invlogcdf(Beta(0.5 + monthly_pos_tests[i], 0.5 + monthly_neg_tests[i]), log(upr)) - prop[i] : 1 - prop[i] for i = 1:length(monthly_total_tests)]

    return (ts = ts[idxs], prop = prop[idxs], lb = lb[idxs], ub = ub[idxs])
end



"""
    function get_credible_intervals(X::AbstractMatrix; lwr = 0.025, upr = 0.975)

Get credible intervals of a matrix `X` over time where dim 1 is assumed to be time, and dim 2 is the variation over MCMC samples.
"""
function get_credible_intervals(X::AbstractMatrix; lwr = 0.025, upr = 0.975)
    pred = mean(X, dims = 2)
    lb = pred .- [quantile(X[t, :], lwr) for t = 1:size(X, 1)]
    ub = [quantile(X[t, :], upr) for t = 1:size(X, 1)] .- pred
    return (pred = pred, lb = lb, ub = ub)
end

"""
    function get_credible_intervals(X::Vector{AbstractArray}; lwr = 0.025, upr = 0.975)

Converts the vector of arrays into a matrix then applies `get_credible_intervals`.
"""
function get_credible_intervals(X::AbstractVector; lwr = 0.025, upr = 0.975)
    X = Matrix(RecursiveArrayTools.vecvec_to_mat(X)')
    return get_credible_intervals(X; lwr = lwr, upr = upr)
end


"""
    function mvav_cols(X::AbstractMatrix)

Return 7-day central moving average for each column in Matrix `X`.
"""
function mvav_cols(X::AbstractMatrix)
    _X = zeros(size(X))
    for j = 1:size(X, 2)
        _X[:, j] .= [X[1:3, j]
            [mean(X[(t-3):(t+3), j]) for t = 4:(size(X, 1)-3)]
            X[(end-2):end, j]]
    end
    return _X
end

"""
    function smooth_vector(X)

Return 7-day central moving average for vector `X`.
"""
function smooth_vector(X::AbstractVector)
    return [mean(X[1:3])
        [mean(X[(t-3):(t+3)]) for t = 4:(length(X)-3)]
        mean(X[(end-2):end])]
end



function vecmat_to_array(X)
    _X = zeros(size(X[1], 1), size(X[1], 2), length(X))
    for k = 1:length(X)
        _X[:, :, k] .= X[k]
    end
    return permutedims(_X, [1, 3, 2])
end

function gather_infection_type_for_econ_analysis(model, projection_date::Date;
    mort_scale::Float64, mort_age, mort_var,
    H_ICU_to_death_fit, H_hosp, H_var_ab, H_var_delta,
    veff1 = 0.24, veff2 = 0.187, veff_death1 = 0.2, veff_death2 = 0.1)

    infections_by_age_dose = model |> flatten_chain .|> θ -> project_daily_infections_by_age_and_vaccine(θ, model, Date(2023, 6, 1))
    age_symptoms = model.p_symp

    uncertaintyinfections_by_age_unvac_first = Matrix.([vecvec_to_mat([infs[1][:, a, 1] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_unvac_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 1] .- infs[1][:, a, 1] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_1dose_first = Matrix.([vecvec_to_mat([infs[1][:, a, 2] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_1dose_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 2] .- infs[1][:, a, 2] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_2dose_first = Matrix.([vecvec_to_mat([infs[1][:, a, 3] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_2dose_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 3] .- infs[1][:, a, 3] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])

    θs = model |> flatten_chain
    n = size(uncertaintyinfections_by_age_unvac_first[1], 1)
    var_frqs = θs .|> θ -> get_variant_prediction(θ, model, n)
    wt_frqs = Matrix(vecvec_to_mat([data.prop_inc_wt for data in var_frqs])')
    alpha_beta_frqs = Matrix(vecvec_to_mat([data.prop_inc_alpha_beta for data in var_frqs])')
    delta_frqs = Matrix(vecvec_to_mat([data.prop_inc_delta for data in var_frqs])')

    var_deadliness = wt_frqs .+ mort_var[1] .* alpha_beta_frqs .+ mort_var[2] .* delta_frqs
    var_severe_risk = H_hosp .* wt_frqs .+ H_var_ab .* mort_var[1] .* alpha_beta_frqs .+ H_var_delta .* mort_var[2] .* delta_frqs
    var_critical_risk = H_ICU_to_death_fit .* var_deadliness

    all_infs = [uncertaintyinfections_by_age_unvac_first[a] + uncertaintyinfections_by_age_unvac_reinf[a] +
                uncertaintyinfections_by_age_1dose_first[a] + uncertaintyinfections_by_age_1dose_reinf[a] +
                uncertaintyinfections_by_age_2dose_first[a] + uncertaintyinfections_by_age_2dose_reinf[a]
                for a = 1:6]

    temp_deadly_infs = [uncertaintyinfections_by_age_unvac_first[a] + 0.15 * uncertaintyinfections_by_age_unvac_reinf[a] +
                        veff_death1 * uncertaintyinfections_by_age_1dose_first[a] + veff_death1 * 0.15 * uncertaintyinfections_by_age_1dose_reinf[a] +
                        veff_death2 * uncertaintyinfections_by_age_2dose_first[a] + veff_death2 * 0.15 * uncertaintyinfections_by_age_2dose_reinf[a]
                        for a = 1:6] .* mort_scale .* mort_age

    temp_sev_crit_infs = [uncertaintyinfections_by_age_unvac_first[a] + 0.15 * uncertaintyinfections_by_age_unvac_reinf[a] +
                          veff1 * uncertaintyinfections_by_age_1dose_first[a] + veff1 * 0.15 * uncertaintyinfections_by_age_1dose_reinf[a] +
                          veff2 * uncertaintyinfections_by_age_2dose_first[a] + veff2 * 0.15 * uncertaintyinfections_by_age_2dose_reinf[a]
                          for a = 1:6] .* mort_scale .* mort_age

    asymptomatic_infs = vecmat_to_array(all_infs .* (1.0 .- age_symptoms))
    deadly_infs = vecmat_to_array(temp_deadly_infs .|> x -> x .* var_deadliness)
    severe_infs = vecmat_to_array(temp_sev_crit_infs .|> x -> x .* var_severe_risk)
    critical_infs = vecmat_to_array(temp_sev_crit_infs .|> x -> x .* var_critical_risk)
    mild_infections = vecmat_to_array((all_infs .* age_symptoms)) .- deadly_infs .- severe_infs .- critical_infs

    return (; asymptomatic_infs, mild_infections, severe_infs, critical_infs, deadly_infs)

end

function gather_infection_type_for_econ_analysis_immune_escape(model, projection_date::Date, cb_imm_esc, cb_vac;
    mort_scale::Float64, mort_age, mort_var,
    H_ICU_to_death_fit, H_hosp, H_var_ab, H_var_delta,
    veff1 = 0.24, veff2 = 0.187, veff_death1 = 0.2, veff_death2 = 0.1)

    infections_by_age_dose = model |> flatten_chain .|> θ -> project_daily_infections_by_age_and_vaccine_immune_escape(θ, model, projection_date, cb_imm_esc, cb_vac)
    age_symptoms = model.p_symp

    uncertaintyinfections_by_age_unvac_first = Matrix.([vecvec_to_mat([infs[1][:, a, 1] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_unvac_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 1] .- infs[1][:, a, 1] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_1dose_first = Matrix.([vecvec_to_mat([infs[1][:, a, 2] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_1dose_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 2] .- infs[1][:, a, 2] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_2dose_first = Matrix.([vecvec_to_mat([infs[1][:, a, 3] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])
    uncertaintyinfections_by_age_2dose_reinf = Matrix.([vecvec_to_mat([infs[2][:, a, 3] .- infs[1][:, a, 3] for infs in infections_by_age_dose])' for a = 1:size(infections_by_age_dose[1][1], 2)])

    θs = model |> flatten_chain
    n = size(uncertaintyinfections_by_age_unvac_first[1], 1)
    var_frqs = θs .|> θ -> get_variant_prediction(θ, model, n)
    wt_frqs = Matrix(vecvec_to_mat([data.prop_inc_wt for data in var_frqs])')
    alpha_beta_frqs = Matrix(vecvec_to_mat([data.prop_inc_alpha_beta for data in var_frqs])')
    delta_frqs = Matrix(vecvec_to_mat([data.prop_inc_delta for data in var_frqs])')

    var_deadliness = wt_frqs .+ mort_var[1] .* alpha_beta_frqs .+ mort_var[2] .* delta_frqs
    var_severe_risk = H_hosp .* wt_frqs .+ H_var_ab .* mort_var[1] .* alpha_beta_frqs .+ H_var_delta .* mort_var[2] .* delta_frqs
    var_critical_risk = H_ICU_to_death_fit .* var_deadliness

    all_infs = [uncertaintyinfections_by_age_unvac_first[a] + uncertaintyinfections_by_age_unvac_reinf[a] +
                uncertaintyinfections_by_age_1dose_first[a] + uncertaintyinfections_by_age_1dose_reinf[a] +
                uncertaintyinfections_by_age_2dose_first[a] + uncertaintyinfections_by_age_2dose_reinf[a]
                for a = 1:6]

    temp_deadly_infs = [uncertaintyinfections_by_age_unvac_first[a] + 0.15 * uncertaintyinfections_by_age_unvac_reinf[a] +
                        veff_death1 * uncertaintyinfections_by_age_1dose_first[a] + veff_death1 * 0.15 * uncertaintyinfections_by_age_1dose_reinf[a] +
                        veff_death2 * uncertaintyinfections_by_age_2dose_first[a] + veff_death2 * 0.15 * uncertaintyinfections_by_age_2dose_reinf[a]
                        for a = 1:6] .* mort_scale .* mort_age

    temp_sev_crit_infs = [uncertaintyinfections_by_age_unvac_first[a] + 0.15 * uncertaintyinfections_by_age_unvac_reinf[a] +
                          veff1 * uncertaintyinfections_by_age_1dose_first[a] + veff1 * 0.15 * uncertaintyinfections_by_age_1dose_reinf[a] +
                          veff2 * uncertaintyinfections_by_age_2dose_first[a] + veff2 * 0.15 * uncertaintyinfections_by_age_2dose_reinf[a]
                          for a = 1:6] .* mort_scale .* mort_age

    asymptomatic_infs = vecmat_to_array(all_infs .* (1.0 .- age_symptoms))
    deadly_infs = vecmat_to_array(temp_deadly_infs .|> x -> x .* var_deadliness)
    severe_infs = vecmat_to_array(temp_sev_crit_infs .|> x -> x .* var_severe_risk)
    critical_infs = vecmat_to_array(temp_sev_crit_infs .|> x -> x .* var_critical_risk)
    mild_infections = vecmat_to_array((all_infs .* age_symptoms)) .- deadly_infs .- severe_infs .- critical_infs

    return (; asymptomatic_infs, mild_infections, severe_infs, critical_infs, deadly_infs)

end


function add_infs(infs1, infs2)
    return (asymptomatic_infs = infs1.asymptomatic_infs + infs2.asymptomatic_infs,
        mild_infections = infs1.mild_infections + infs2.mild_infections,
        severe_infs = infs1.severe_infs + infs2.severe_infs,
        critical_infs = infs1.critical_infs + infs2.critical_infs,
        deadly_infs = infs1.deadly_infs + infs2.deadly_infs)
end
