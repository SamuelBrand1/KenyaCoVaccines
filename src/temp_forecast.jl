"""

Predict and plot out vaccination impact in a county given initial condition and vaccination parameter inputs
"""

function predict_vacc_county_new(name::String,linelist_data_with_pos_neg,serology_data,
                    M_Kenya_ho,M_Kenya_other,M_Kenya_school,M_Kenya_work,N_kenya;
                    PCR_array,relative_testing_rate,
                    rapid, scenario,
                    prob_sus_w,
                    prob_sus,
                    prob_symp,
                    H_rate_by_age,
                    p_ID,p_IH, all_IFR_fits,ICU_fit,hosp_fit,
                    deaths_data,
                    priors,
                    T_to_full_efficacy,
                    vacc_end_date = Date(2022,6,30), # date vaccination program ends
                    startdate = Date(2020,12,31),
                    enddate = Date(2022,6,30),
                    obs_factor)



    fit_dict = load("modelfits/$(name)_model.jld2")
    model = fit_dict[first(keys(fit_dict))]
    up_county = uppercase(model.areaname)


    #Set times, all times after model_start_time are in reference to the model_start_time
    model_start_time = model.model_start_time # start in 1st Dec 2020
    model_end_time = model.model_end_time   # end of model run for fitting
    marchendpoint = model.marchendpoint # 10th March when roudn 3 serosurvey data ends
    alpha_variant_time = model.alpha_variant_time  # earliest time of alpha variant introduction
    delta_variant_time = model.delta_variant_time # earliest time of delta variant introduction
    janstarttime = model.janstarttime  # time t start including observations in likelihood
    vacc_end_point =  (vacc_end_date - model.startdate).value # End of model run
    pre_vacc_period = (Date(2021,3,6) - model.startdate).value # vaccination starts on 6th March 2021
    model_end_point =  (enddate - startdate).value # End of model run

    # Variants data
    kenya_variants = CSV.File("data/kenya_gisaid_variants.csv") |> DataFrame
    kenya_variants.date = [Date(d,DateFormat("dd/mm/yyyy")) for d in kenya_variants.date]
    kenya_variants.rel_time = [(d - Date(2020,2,20)).value for d in kenya_variants.date]

    idxs_alpha = kenya_variants.variant .== "Alpha"
    idxs_beta = kenya_variants.variant .== "Beta"
    idxs_delta = kenya_variants.variant .== "Delta"

    #Daily prediction of variant prevalence using linear interpolation
    interp_alpha = DataInterpolations.LinearInterpolation(kenya_variants.perc_sequences[idxs_alpha]./100,kenya_variants.rel_time[idxs_alpha])
    alpha_perc = [interp_alpha(t) for t in 1:maximum(kenya_variants.rel_time[idxs_alpha])]
    interp_beta = DataInterpolations.LinearInterpolation(kenya_variants.perc_sequences[idxs_beta]./100,kenya_variants.rel_time[idxs_beta])
    beta_perc = [interp_beta(t) for t in 1:maximum(kenya_variants.rel_time[idxs_beta])]
    interp_delta = DataInterpolations.LinearInterpolation(kenya_variants.perc_sequences[idxs_delta]./100,kenya_variants.rel_time[idxs_delta])
    delta_perc = [interp_delta(t) for t in 1:maximum(kenya_variants.rel_time[idxs_delta])]

    delta_perc = delta_perc[model_start_time:end]  # pick data from 1st Dec 2020
    alpha_perc = alpha_perc[model_start_time:end]
    beta_perc  = beta_perc[model_start_time:end]

    n = model_end_point

    delta_perc = [delta_perc;ones(n-length(delta_perc))]
    alpha_beta_perc = alpha_perc .+ beta_perc
    alpha_beta_perc = [alpha_beta_perc;zeros(n-length(alpha_beta_perc))]
    other_perc = 1.0 .- alpha_beta_perc .- delta_perc


    # fitted parameters for county
    MCMCdraws = get(model.MCMC_results.chain,[:β₀,:β_home,:β_school,:β_other,:β_work,:ϵ,:χ,:p_test,:E₀,:inc_R_αvar,:time_scale_αvar,:mid_point_αvar,:inc_R_δvar,:time_scale_δvar,:mid_point_δvar,:init_scale])
    prior = priors[countynames.==model.areaname][1]

    # Data sets to populate and save after loop
    crit_dis_inc_samples  = zeros(model_end_point,6,2000) # 6 age groups, 2000 samples from HMC
    severe_dis_inc_samples = zeros(model_end_point,6,2000)
    asymp_dis_inc_samples = zeros(model_end_point,6,2000)
    mild_dis_inc_samples = zeros(model_end_point,6,2000)
    death_inc_samples = zeros(model_end_point,6,2000)

    for k in 1:2000
        α= model.α  # 1/mean latent period
        ϵ= MCMCdraws.ϵ[k] #0.5    # Relative infectiousness of asymptomatic cases
        αP= model.αP # rate of movement from Pre symptomatic state
        γA = model.γA #Recovery rate from Asymptomatic state
        γM= model.γM  #Recovery rate from Mild disease

        relative_testing_rate = KenyaCoVaccines.relative_testing_rate
        PCR_array= KenyaCoVaccines.PCR_array #Relative Sensitivity of PCR test on days post-infection
        PCR_sensitivity = model.PCR_sensitivity #Base sensitivity for RT-PCR test
        PCR_specificity = model.PCR_specificity #Base specificity for RT_PCR test
        baseline_sero_array = KenyaCoVaccines.rel_sero_array_26days #Relative sensitivity of serology test on days post-infection
        sero_sensitivity = model.sero_sensitivity#Base sensitivity for serology test
        sero_specificity= model.sero_specificity #Base specificity for serology test
        M_BB= model.M_BB


        ω= model.ω # Waning of natural protection, rate per day
        σω= model.p_sus_w # relative suceptibility after first episode: 0.16 for all ages for now.
        σ = model.p_sus
        δ = model.p_symp
        υ = model.hosp_rate_by_age


        β₀ = MCMCdraws.β₀[k]
        β_home = MCMCdraws.β_home[k]
        β_school = MCMCdraws.β_school[k]
        β_other = MCMCdraws.β_other[k]
        β_work = MCMCdraws.β_work[k]

        #Set variance scalers
        clustering_factor_PCR = 0.5
        M_PCR = 30.
        T = eltype(β₀)

        inc_R_αvar= MCMCdraws.inc_R_αvar[k]
        time_scale_αvar= MCMCdraws.time_scale_αvar[k]
        mid_point_αvar= MCMCdraws.mid_point_αvar[k]
        inc_R_δvar= MCMCdraws.inc_R_δvar[k]
        time_scale_δvar= MCMCdraws.time_scale_δvar[k]
        mid_point_δvar= MCMCdraws.mid_point_δvar[k]
        init_scale = MCMCdraws.init_scale[k]

        p = convert.(T,[β₀,β_home,0.0,β_other,β_work,α,ϵ,αP,γA,γM,ω,inc_R_αvar,time_scale_αvar,alpha_variant_time + mid_point_αvar*30,inc_R_δvar,time_scale_δvar, delta_variant_time + mid_point_δvar*30,
         init_scale])  # start with schools closed
        p = vcat(p,model.VE_acquisition,model.VE_infectiousness,model.VE_severe_disease_risk) # add vaccine efficacies to the paramter vector


        init_sero = model.init_sero
        γV = γM

        N= model.N
        E₀ = MCMCdraws.E₀[k]

        θ = (β₀=β₀,β_home=β_home,β_school=β_school,β_other=β_other,β_work=β_work,E₀=E₀,ϵ=ϵ,init_scale=init_scale)
        u0 = KenyaCoVaccines.create_initial_conditions(θ,model)


        #Set up the callback for lockdown
        schoolsclosure_times = [Float64((Date(2021, 3, 19) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 7, 16) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 10, 1) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 12, 23) - Date(2020, 12, 1)).value)]
        function affect_schoolsclosed!(integrator)
            integrator.p[3] = 0.0 #Schools shut
        end
        schoolsclosed_cb = PresetTimeCallback(schoolsclosure_times, affect_schoolsclosed!,save_positions=(false,false))

        #Set up callback for reopening Schools
        schoolsreopening_times = [Float64((Date(2021, 1, 4) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 5, 10) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 7, 26) - Date(2020, 12, 1)).value),
            Float64((Date(2021, 10, 11) - Date(2020, 12, 1)).value)]
        function affect_schoolsopening!(integrator)
            integrator.p[3] = β_school # Schools open #
        end
        schoolsreopening_cb = PresetTimeCallback(schoolsreopening_times, affect_schoolsopening!,save_positions=(false,false))

        cb = CallbackSet(schoolsclosed_cb, schoolsreopening_cb)

        immune_escape_variant_time = Float64((Date(2022,1,1) - model.startdate).value) #  December 1 2021, time to implement variant immunity escape

        function immune_escape_variant!(integrator)
            # Reduce vaccine efficacies by 50%
            integrator.p[19:27] *= 0.5

            # Start with 'variant_E' exposed number of individuals
            variant_E = 1
            prop_vacc = sum(integrator.u[:,1:8,:],dims=2)[:,1,:] ./ sum(integrator.u[:,1:8,:],dims=[2,3])[:]
            integrator.u[:,2,:] .= variant_E .* prop_vacc # number exposed in each age group and vaccinatin status
            integrator.u[:,1,:] .= integrator.u[:,1,:] .- integrator.u[:,2,:] #substract E from S

            # move 50% from R and W -> S
            loss_susceptible_R = 0.5 .* integrator.u[:,7,:]
            loss_susceptible_W = 0.5 .* integrator.u[:,8,:]
            integrator.u[:,1,:] .+= loss_susceptible_R
            integrator.u[:,7,:] .-= loss_susceptible_R
            integrator.u[:,1,:] .+= loss_susceptible_W
            integrator.u[:,8,:] .-= loss_susceptible_W
        end
        immune_escape_variant_cb = PresetTimeCallback([immune_escape_variant_time],immune_escape_variant!,save_positions=(false,false))

        # Update callbacks based on scenario
        if  scenario==1 || scenario>4  # For Scenario 4-6
         cb = CallbackSet(cb, immune_escape_variant_cb)
        end

        county_model = KenyaCoVaccines.CoVAreaModel(name,KenyaCoVaccines.ll_onegroup_newvariant_infboost,prior;
                                PCR_cases = linelist_data_with_pos_neg,
                                sero_cases = serology_data,
                                average_sero_init = zeros(6),
                                deaths = deaths_data,
                                pop_data= N_kenya,
                                M_county_ho=M_Kenya_ho,
                                M_county_other=M_Kenya_other,
                                M_county_school=M_Kenya_school,
                                M_county_work=M_Kenya_work,
                                rapid = rapid,
                                scenario= scenario,
                                startdate=startdate,
                                enddate=enddate)

        #Solve dynamical system
        sol = solve(county_model.prob, BS3();tspan = (0,model_end_point),callback = cb,u0=u0, p=p,saveat = 1) #isoutofdomain=(u,p,t)->any(x->x<0,u)

        # Number of unscalled (for testing) cases over time by severity
        ι = KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol,10)

        # Number of tests over time by age
        PCR_cases = sum(linelist_data_with_pos_neg.cases[:,linelist_data_with_pos_neg.areas .== up_county,:,:],dims = [2])[:,1,:,:]
        p_vec = [[1,2,3,4],[5,6,7,8,9,10],[11,12],[13,14],[15,16],[17]]
        PCR_cases_collaped =  zeros(Int64,size(PCR_cases,1),6,2)
        for a in 1:6
            PCR_cases_collaped[:,a:a,:] .= sum(PCR_cases[:,p_vec[a][1]:p_vec[a][end],:], dims=2)
        end
        PCR_cases = PCR_cases_collaped
        num_tests =  PCR_cases[model_start_time:end,:,1] .+ PCR_cases[model_start_time:end,:,2]
        num_tests_age = zeros(size(ι,1),6)
        for a in 1:6
            num_tests_age[:,a] = vcat(num_tests[:,a],fill(mean(num_tests[:,a]),(length(ι)-length(num_tests[:,a]))))  # use mean number of teste for forcast periods, could be changed
        end

        # scalling parameters
        p_test = MCMCdraws.p_test[k]
        χ = MCMCdraws.χ[k]
        p_test_scale = 5e-4


        # Apply vaccine efficacy against disease
        v = (1 .- model.VE_severe_disease_risk') .* υ  # efficacy

        l_inc = [inc .* v for inc in ι]
        cases_wt = similar(ι)
        cases_αβ = similar(ι)
        cases_δ  = similar(ι)
        crit_incidence = similar(ι)
        sev_incidence = similar(ι)


        for i in 1:size(ι,1)
            cases_wt[i] = sum(l_inc[i] .* other_perc[i],dims=2) # sum over dose
            cases_αβ[i] = sum(l_inc[i] .* alpha_beta_perc[i], dims=2)
            cases_δ[i]  = sum(l_inc[i] .* delta_perc[i], dims=2)

            PCR_num_wt = cases_wt[i] .* p_test_scale .* p_test
            PCR_denom_wt = PCR_num_wt .+ p_test_scale*((p_test/χ)*(N .- cases_wt[i]))
            p_PCR_pred_wt = PCR_num_wt./PCR_denom_wt
            cases_wt[i]  =  p_PCR_pred_wt .* num_tests_age[i,:]

            PCR_num_αβ = cases_αβ[i] .* p_test_scale .* p_test
            PCR_denom_αβ = PCR_num_αβ .+ p_test_scale*((p_test/χ)*(N .- cases_αβ[i]))
            p_PCR_pred_αβ = PCR_num_αβ./PCR_denom_αβ
            cases_αβ[i]  =  p_PCR_pred_αβ .* num_tests_age[i,:]

            PCR_num_δ = cases_δ[i] .* p_test_scale .* p_test
            PCR_denom_δ = PCR_num_wt .+ p_test_scale*((p_test/χ)*(N .- cases_δ[i]))
            p_PCR_pred_δ = PCR_num_δ./PCR_denom_δ
            cases_δ[i]  =  p_PCR_pred_δ .* num_tests_age[i,:]

            # need to edit structure on ICU and hosp fit to include kth sample fit rathert than point est (single value)
            crit_incidence[i] = ICU_fit[1]*(cases_wt[i] .+ ICU_fit[2].*cases_αβ[i] .+ ICU_fit[3].*cases_δ[i])
            sev_incidence[i] = hosp_fit[1]*(cases_wt[i] .+ hosp_fit[2].*cases_αβ[i] .+ hosp_fit[3].*cases_δ[i])
        end

        crit_incidence =  Matrix(VectorOfArray([inc[:] for inc in crit_incidence])')  # convert to a matrix
        sev_incidence =  Matrix(VectorOfArray([inc[:] for inc in sev_incidence])')

        for a in 1:6 # convolve with PCR array to get right timing of PCR detectable incidence
            crit_incidence[:,a] = KenyaCoVaccines.simple_conv(crit_incidence[:,a],KenyaCoVaccines.PCR_array)
            sev_incidence[:,a] = KenyaCoVaccines.simple_conv(sev_incidence[:,a],KenyaCoVaccines.PCR_array)
        end

        # Saving severe and critical disease incidence
        crit_dis_inc_samples[:,:,k]   .= crit_incidence
        severe_dis_inc_samples[:,:,k] .= sev_incidence

        # Asymptomatic cases
        asymp_dis_inc = [inc .* (ones(3)'.* (1 .- δ )) for inc in ι]
        asymp_dis_inc = Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in asymp_dis_inc])')  #sum over dose
        asymp_dis_inc_samples[:,:,k] .= asymp_dis_inc

        # Mild disease
        symp_dis_inc  =  [inc .* (ones(3)'.* δ) for inc in ι]
        symp_dis_inc =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in symp_dis_inc])')  #sum over dose
        for a in 1:6
            symp_dis_inc[:,a] =  KenyaCoVaccines.simple_conv(symp_dis_inc[:,a],p_IH)  # convolve with time to hsop function to get right timing of Hospitalisation
        end
        mild_dis_inc = symp_dis_inc .- sev_incidence .- crit_incidence
        mild_dis_inc_samples[:,:,k] .= mild_dis_inc

        # Number of unscalled deaths by age group over time
        IFR = all_IFR_fits[[nm[1] .== name for nm in all_IFR_fits]][1].IFR # EDIT to include HMC samples, so there shoudl be a 'k' index on the IFR
        d = (1 .- VE_severe_disease_risk') .* IFR # efficacy multiplied by risk of death by age, assumes only those with severe disease die

        death_inc =  [inc .* d for inc in ι]
        death_inc =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in death_inc])')  #sum over dose
        for a in 1:6
            death_inc[:,a] =  KenyaCoVaccines.simple_conv(death_inc[:,a],p_ID)  # convolve with time to death function to get right timing of deaths
        end
        death_inc = death_inc ./ obs_factor

    end # end loop over HMC chains

    # pick data from Jan for CEA analysis
    asymp_dis_inc_samples = asymp_dis_inc_samples[janstarttime:end,:,:]
    crit_dis_inc_samples = crit_dis_inc_samples[janstarttime:end,:,:]
    severe_dis_inc_samples = severe_dis_inc_samples[janstarttime:end,:,:]
    mild_dis_inc_samples = mild_dis_inc_samples[janstarttime:end,:,:]
    death_inc_samples = death_inc_samples[janstarttime:end,:,:]

    if rapid  # becasue we now have rapid scenario we save them in different directory: I have created this filing structure in the repository
        rapid_text = "rapid"
    else
        rapid_text = "non_rapid"
    end

    @save("scenario_$(scenario)/$(rapid_text)/mild_cases/severe_cases_$(name).jld2",mild_dis_inc_samples)
    @save("scenario_$(scenario)/$(rapid_text)/severe_cases/severe_cases_$(name).jld2",severe_dis_inc_samples)
    @save("scenario_$(scenario)/$(rapid_text)/critical_cases/critical_cases_$(name).jld2",crit_dis_inc_samples)
    @save("scenario_$(scenario)/$(rapid_text)/asymp_cases/asymp_cases_$(name).jld2",asymp_dis_inc_samples)
    @save("scenario_$(scenario)/$(rapid_text)/deaths/deaths_$(name).jld2",death_inc_samples)

    # Number of vaccines administered over time, this is the same for all 2000 chains
    doses_Kenya_baseline = 0.
    if scenario==1
        doses_Kenya = 0.
    else
        doses_Kenya = (1.25 + 9.76 +4.9) .* 1000000  #Number of doses planned for Kenya over the ynext one year under a 30% cover plan
    end

    age_dist_vacc, doses_daily = KenyaCoVaccines.age_dist([name],doses_Kenya,N_kenya,vacc_end_point,pre_vacc_period,rapid)
    vacc_rate_1,vacc_rate_2 = KenyaCoVaccines.vacc_rate_data(doses_daily,vacc_end_point,pre_vacc_period,age_dist_vacc,N_kenya,name,scenario,rapid)
    num_doses = sum(vacc_rate_1, dims=2)[:] .+ sum(vacc_rate_2, dims=2)[:]
    num_doses = num_doses[janstarttime:end]
    complete_doses_age = vacc_rate_2[janstarttime:end,:]  # number with complete doses


    age_dist_vacc, doses_daily = KenyaCoVaccines.age_dist([name],doses_Kenya_baseline,N_kenya,vacc_end_point,pre_vacc_period,rapid)
    vacc_rate_1,vacc_rate_2 = KenyaCoVaccines.vacc_rate_data(doses_daily,vacc_end_point,pre_vacc_period,age_dist_vacc,N_kenya,name,baseline_scenario,rapid)
    num_doses_baseline = sum(vacc_rate_1, dims=2)[:] .+ sum(vacc_rate_2, dims=2)[:]
    num_doses_baseline = num_doses_baseline[janstarttime:end]
    complete_doses_age_baseline = vacc_rate_2[janstarttime:end,:]

    @save("scenario_$(scenario)/$(rapid_text)/doses/doses_$(name).jld2",num_doses)
    @save("scenario_$(scenario)/$(rapid_text)/doses_age/doses_$(name).jld2",complete_doses_age)

    return nothing
end
