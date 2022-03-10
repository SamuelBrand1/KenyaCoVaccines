

function extra_doses(vacc_data,N_vacc)
    dose_cum = zeros(size(vacc_data))
    for age in 1:6
        dose_cum[:,age] = cumsum(vacc_data[:,age])
    end
    above_N = trues(size(dose_cum))
    for time in 1:size(dose_cum)[1]
        above_N[time,:] .= (dose_cum[time,:] .> N_vacc)
    end
    vacc_extra = vacc_data .- (vacc_data .* (1 .- above_N))  # set doses past full coverage time in agegroup to zero
    return(vacc_extra)
end

"""
Function to generate vaccination rates: The GoK plan taken more generally, not in Phases
"""

function vacc_rate_data(doses_daily,model_days,pre_vacc_period,age_dist_vacc,N_kenya,county,scenario,rapid)
    N = Array(N_kenya[:,county])
    N = vcat(sum(N[1:4]),sum(N[5:10]),sum(N[11:12]),sum(N[13:14]),sum(N[15:16]),N[17])  # reduce age categories to six
    N_vacc = vcat(N[1]*2/20,N[2:end])  # only 2/20 in the first age group need to be vaccinated

    if scenario == 3 || scenario == 6
        doses_daily = doses_daily .* (25/15.9) # Increase factor: Use 25 million doses instead of the 15.9 million baseline doses
    end

    if scenario == 4 || scenario == 7
        doses_daily = doses_daily .* (35/15.9) # Increase factor: Use 35 million doses instead of the 15.9 million baseline doses
    end

    dose_interval = 8*7 # number of days between 1st and second doses
    doses_daily = doses_daily/2  # not all available doses are given have to reserve one doese for each dose given

    vacc_1_data = zeros(model_days,6)  # To hold number of dose 1 vaccinations per day in each age group ## VACC DATA SHOULD BE TIME (ROWS) x AGE GROUP (COLS)
    vacc_2_data = zeros(model_days,6)  # To hold number of dose 2 vaccinations per day in each age group ## VACC DATA SHOULD BE TIME (ROWS) x AGE GROUP (COLS)

    ds = sum(doses_daily[(pre_vacc_period+1): (pre_vacc_period+dose_interval)]) #/2  # doses that wont be issued at the tail end of vaccination period, recover them by giving earlier to keep total doses

    if rapid
       ds = ds/(180 - 56)  # per day addional dose 1 doses 56 days before end of vaccination period
    else
       ds = ds/(model_days - pre_vacc_period - 56)  # per day addional dose 1 doses 56 days before end of vaccination period
    end

    for time in (pre_vacc_period+1): (pre_vacc_period+dose_interval)
      vacc_1_data[time,:] .= (doses_daily[time] + ds).*age_dist_vacc
    end

    for time in (pre_vacc_period+dose_interval+1): model_days
       if rapid
          if  time <= (pre_vacc_period+180 - 56)
              vacc_1_data[time,:] .= (doses_daily[time] + ds).*age_dist_vacc
          end
       else
          if  time <= (model_days - 56)
              vacc_1_data[time,:] .= (doses_daily[time] + ds).*age_dist_vacc
          end
       end
       vacc_2_data[time,:] .= sum(vacc_1_data[time .- dose_interval,:]).*age_dist_vacc
    end

    # limit the number of doses to not exceed population numbers

    vacc_1_extra = KenyaCoVaccines.extra_doses(vacc_1_data,N_vacc)
    vacc_1_data_ed = vacc_1_data

    age = 6
    while age > 1  # push extra doses down the age groups
        vacc_1_data_ed[:,(age-1)] .+= vacc_1_extra[:,age]
        vacc_1_data_ed[:,age] .-= vacc_1_extra[:,age]
        vacc_1_extra = KenyaCoVaccines.extra_doses(vacc_1_data_ed,N_vacc)
        age = age -1
    end

    for age in 1:5
        vacc_1_data_ed[:,(age+1)] .+= vacc_1_extra[:,age]
        vacc_1_data_ed[:,age] .-= vacc_1_extra[:,age]
        vacc_1_extra = KenyaCoVaccines.extra_doses(vacc_1_data_ed,N_vacc)
    end

    # now check if any extra doses and set to zero (overalocation in the county) assume not transfer between counties ATM
    vacc_1_extra = KenyaCoVaccines.extra_doses(vacc_1_data_ed,N_vacc)
    vacc_1_data_ed .-= vacc_1_extra
    vacc_1_data = vacc_1_data_ed

    vacc_2_extra = KenyaCoVaccines.extra_doses(vacc_2_data,N_vacc)
    vacc_2_data_ed = vacc_2_data

    age = 6
    while age > 1  # push extra doses down the age groups
        vacc_2_data_ed[:,(age-1)] .+= vacc_2_extra[:,age]
        vacc_2_data_ed[:,age] .-= vacc_2_extra[:,age]
        vacc_2_extra = KenyaCoVaccines.extra_doses(vacc_2_data_ed,N_vacc)
        age = age -1
    end

    for age in 1:5
        vacc_2_data_ed[:,(age+1)] .+= vacc_2_extra[:,age]
        vacc_2_data_ed[:,age] .-= vacc_2_extra[:,age]
        vacc_2_extra = KenyaCoVaccines.extra_doses(vacc_2_data_ed,N_vacc)
    end

    # now check if any extra doses and set to zero (overalocation in the county) assume not transfer between counties ATM
    vacc_2_extra = KenyaCoVaccines.extra_doses(vacc_2_data_ed,N_vacc)
    vacc_2_data_ed .-= vacc_2_extra
    vacc_2_data = vacc_2_data_ed

    return vacc_1_data,vacc_2_data
end

"""
  Function to generate target age distribution of vaccinees in a given county
"""

function age_dist(county,total_doses,N_kenya,model_days,pre_vacc_period,rapid)

    pop_county =  N_kenya[:,county]
    pop_kenya =  sum(N_kenya, dims=2)

    pop_kenya_over_18 = sum(vcat((2/5)*pop_kenya[4],pop_kenya[5:end]))  # 2/5 is age 18 and 19 among 15-19 age category
    pop_county_over_18 = sum(vcat((2/5)*pop_county[4],pop_county[5:end]))

    vacc_dist_county = pop_county_over_18/pop_kenya_over_18 # proportion od vaccine to a county is proportional to its population above 18 years: assumption

    doses_per_day = zeros(model_days)
    doses_per_day[(pre_vacc_period+1):model_days]  .= total_doses/(model_days-pre_vacc_period) # use vaccination rate (flat rate, using 15.91 million doses planned for Phase I,II and III
    doses_per_day_county = doses_per_day.*vacc_dist_county #

    if rapid
        doses_per_day = zeros(model_days) # to simulate some rapid vaccination over two months
        doses_per_day[(pre_vacc_period+1):(pre_vacc_period+180)] .= total_doses/180 # compress all doses assigned to over six months
        doses_per_day_county = doses_per_day.*vacc_dist_county
    end

    N50  = sum(Array(N_kenya[:,county])[11:17])  # number >50years

    number_doses_over_50 = N50/(model_days - pre_vacc_period - 56)  # The number of doses (dose 1 ) per day to cover 100% of this age group
    prop_over_50 = min(1,number_doses_over_50/(doses_per_day_county[pre_vacc_period+1]/2)) # prioritise vaccination on 50+, this is the fraction of doses that will go to 50+

    if rapid
        number_doses_over_50 = N50/(180 - 56)  # The number of doses (dose 1 ) per day to cover 100% of this age group
        prop_over_50 = min(1,number_doses_over_50/(doses_per_day_county[pre_vacc_period+1]/2)) # prioritise vaccination on 50+, this is the fraction of doses that will go to 50+
    end

    age_dist_vacc = zeros(6)  # age distrbution of vaccinees
    age_dist_vacc[1] =  (2/32)*(1-prop_over_50)  # make 18-50 age group represent remaining doses. min(1,#) above is used because there are some counties wehre doses assigned is not able to cover beyond the pop aged >50: two age years 18 and 19 out of 32 vaccinated
    age_dist_vacc[2] = (30/32)*(1-prop_over_50)
    age_dist_vacc[3:6] .= prop_over_50/4
    sum(age_dist_vacc)

    return age_dist_vacc,doses_per_day_county
end


# Vacc predictions functions

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
                    p_ID,p_IH,
                    deaths_data,
                    priors,
                    T_to_full_efficacy,
                    vacc_end_date = Date(2022,6,30), # date vaccination program ends
                    startdate = Date(2020,12,31),
                    enddate = Date(2022,6,30),
                    obs_factor)

    # load all_IFR_fits,ICU_fit,hosp_fit here, or pass in to the function as arguments

    # Load the model object for the county
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


"""
OLDER FUNTION THAT WAS USING POINT ESTIMATES
Predict and plot out vaccination impact in a county given initial condition and vaccination parameter inputs
"""

function predict_vacc_county(name::String,linelist_data_with_pos_neg,serology_data,
                    M_Kenya_ho,M_Kenya_other,M_Kenya_school,M_Kenya_work,N_kenya;
                    PCR_array,relative_testing_rate,
                    all_param_fits_newvariant,
                    rapid, scenario,baseline_scenario,countynames,
                    prob_sus_w,
                    prob_sus,
                    prob_symp,
                    H_rate_by_age,
                    av_sero_init,p_ID,p_IH, all_IFR_fits,ICU_fit,hosp_fit,
                    deaths_data,
                    priors,
                    T_to_full_efficacy,
                    VE_acquisition,  #Vaccine efficacy against acquisition by dose
                    VE_infectiousness, #Vaccine efficacy against transmission by dose
                    VE_severe_disease_risk,  #Vaccine efficacy against disease by dose
                    vacc_end_date = Date(2022,6,30), # date vaccination program ends
                    startdate = Date(2020,12,31),
                    enddate = Date(2022,6,30),
                    obs_factor)

    up_county = uppercase(name)

    # fitted parameters for county
    fit_p = all_param_fits_newvariant[[fit[1] .== name for fit in all_param_fits_newvariant]][1][2]
    prior = priors[countynames.==name][1]

    α= 1/3.1  # 1/mean latent period
    ϵ= fit_p[6] #0.5    # Relative infectiousness of asymptomatic cases
    αP= 1/1.9  # rate of movement from Pre symptomatic state
    γA = 1/2.4 #Recovery rate from Asymptomatic state
    γM= 1/2.4 #Recovery rate from Mild disease

    relative_testing_rate = KenyaCoVaccines.relative_testing_rate
    PCR_array= KenyaCoVaccines.PCR_array #Relative Sensitivity of PCR test on days post-infection
    PCR_sensitivity = 1. #Base sensitivity for RT-PCR test
    PCR_specificity = 0.995 #Base specificity for RT_PCR test
    baseline_sero_array = KenyaCoVaccines.rel_sero_array_26days #Relative sensitivity of serology test on days post-infection
    sero_sensitivity = 0.825 #Base sensitivity for serology test
    sero_specificity= 0.992 #Base specificity for serology test
    M_BB= 40.74


    ω= 1/180 # Waning of natural protection, rate per day
    σω= prob_sus_w  # relative suceptibility after first episode: 0.16 for all ages for now.
    σ = prob_sus
    δ = prob_symp
    υ = H_rate_by_age


    β₀ = fit_p[1]
    β_home = fit_p[2]
    β_school = fit_p[3]
    β_other = fit_p[4]
    β_work = fit_p[5]

    #Set variance scalers
    clustering_factor_PCR = 0.5
    M_PCR = 30.
    T = eltype(β₀)

    inc_R_αvar= fit_p[15]
    time_scale_αvar= fit_p[16]
    mid_point_αvar= fit_p[17]
    inc_R_δvar= fit_p[18]
    time_scale_δvar= fit_p[19]
    mid_point_δvar= fit_p[20]
    init_scale = fit_p[21]

    #Set times, all times after model_start_time are in reference to the model_start_time
    model_start_time = (startdate - Date(2020,2,20)).value # start in 1st Dec 2020
    model_end_time = (Date(2021,9,24) - startdate).value  # end of model run for fitting
    marchendpoint = (Date(2021,3,10) - startdate).value    # 10th March when roudn 3 serosurvey data ends
    alpha_variant_time = (Date(2021,2,1) - startdate).value  # earliest time of alpha variant introduction
    delta_variant_time = (Date(2021,4,1) - startdate).value  # earliest time of delta variant introduction
    janstarttime = (Date(2021,1,1) - startdate).value # time t start including observations in likelihood
    vacc_end_point =  (vacc_end_date - startdate).value # End of model run
    pre_vacc_period = (Date(2021,3,6) - startdate).value # vaccination starts on 6th March 2021
    model_end_point =  (enddate - startdate).value # End of model run

    p = convert.(T,[β₀,β_home,β_school,β_other,β_work,α,ϵ,αP,γA,γM,ω,inc_R_αvar,time_scale_αvar,alpha_variant_time + mid_point_αvar*30,inc_R_δvar,time_scale_δvar, delta_variant_time + mid_point_δvar*30,
     init_scale])
    p = vcat(p,VE_acquisition,VE_infectiousness,VE_severe_disease_risk) # add vaccine efficacies to the paramter vector

    sero_cases = sum(serology_data.sero[:,serology_data.areas .== up_county,:,:],dims = [2])[:,1,:,:]
    if sum(sero_cases[(model_start_time:(model_start_time+marchendpoint)),:,:]) == 0
       init_sero = av_sero_init  # replace with the average
    else
        p_vec = [[1,2,3,4],[5,6,7,8,9,10],[11,12],[13,14],[15,16],[17]]  # reduce sero cases to 6 age groups
        sero_cases_collaped =  zeros(Int64,size(sero_cases,1),6,2)
        for a in 1:6
            sero_cases_collaped[:,a:a,:] .= sum(sero_cases[:,p_vec[a][1]:p_vec[a][end],:], dims=2)
        end
        sero_cases = sero_cases_collaped
        u_sero = sero_cases[(model_start_time:(model_start_time+marchendpoint)),:,:]
        u_sero = sum(u_sero, dims=1)[1,:,:]
        init_sero = u_sero[:,1] ./ sum(u_sero,dims=2)[:]
    end

    if isnan(init_sero[1])  # action if first value is isnan: replace with next succeding available value
        a = 1
        while isnan(init_sero[1])
            init_sero[1] = init_sero[a+1]
            a += 1
        end
    end


    for a in 2:6
        if isnan(init_sero[a])
            init_sero[a] = init_sero[a-1] # use the seroprevalence in the preceding age groups
        end
    end


    γV = γM

    N= Array(N_kenya[:,name])  #Area's population size
    N = vcat(sum(N[1:4]),sum(N[5:10]),sum(N[11:12]),sum(N[13:14]),sum(N[15:16]),N[17]) # reduce to six age categories

    E_vec = fit_p[9:14]
    u0 = KenyaCoVaccines.create_initial_conditions(init_sero,init_scale,E_vec,N,α,ω,δ,αP,υ,γA,γM,γV) #NOTE_THAT γV = γM

    #Set up the callback for lockdown
    schoolsclosed_time = Float64((Date(2021,1,1) - startdate).value)
    function affect_schoolsclosed!(integrator)
        integrator.p[3] = 0. #Schools shut
    end
    schoolsclosed_cb = PresetTimeCallback([schoolsclosed_time],affect_schoolsclosed!,save_positions=(false,false))

    #Set up callback for reopening Schools
    schoolsreopening_time = Float64((Date(2021,1,6) - startdate).value)
    function affect_schoolsopening!(integrator)
        integrator.p[3] = β_school # Schools open #
    end
    schoolsreopening_cb = PresetTimeCallback([schoolsreopening_time],affect_schoolsopening!,save_positions=(false,false))

    cb = CallbackSet(schoolsclosed_cb,schoolsreopening_cb)
    cb_baseline  = CallbackSet(schoolsclosed_cb,schoolsreopening_cb)

    immune_escape_variant_time = Float64((Date(2021,12,1) - startdate).value) #  December 1 2021, time to implement variant immunity escape

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

    if  baseline_scenario==1 || baseline_scenario>4  # For Scenario 4-6
     cb_baseline = CallbackSet(cb_baseline, immune_escape_variant_cb)
    end

    # Vaccination rates and vaccine uptake
    county_model_baseline = KenyaCoVaccines.CoVAreaModel_proj(name,KenyaCoVaccines.ll_onegroup_newvariant_infboost,prior;
                            PCR_cases = linelist_data_with_pos_neg,
                            sero_cases = serology_data,
                            average_sero_init = av_sero_init,
                            deaths = deaths_data,
                            pop_data= N_kenya,
                            M_county_ho=M_Kenya_ho,
                            M_county_other=M_Kenya_other,
                            M_county_school=M_Kenya_school,
                            M_county_work=M_Kenya_work,
                            rapid = rapid,
                            scenario= baseline_scenario,
                            time_to_full_efficacy = T_to_full_efficacy,
                            hosp_rate_by_age = H_rate_by_age,
                            p_symp = prob_symp,
                            p_sus = prob_sus,
                            p_sus_w=prob_sus_w,
                            startdate=startdate,
                            enddate=enddate)

    county_model = KenyaCoVaccines.CoVAreaModel_proj(name,KenyaCoVaccines.ll_onegroup_newvariant_infboost,prior;
                            PCR_cases = linelist_data_with_pos_neg,
                            sero_cases = serology_data,
                            average_sero_init = av_sero_init,
                            deaths = deaths_data,
                            pop_data= N_kenya,
                            M_county_ho=M_Kenya_ho,
                            M_county_other=M_Kenya_other,
                            M_county_school=M_Kenya_school,
                            M_county_work=M_Kenya_work,
                            rapid = rapid,
                            scenario= scenario,
                            time_to_full_efficacy = T_to_full_efficacy,
                            hosp_rate_by_age = H_rate_by_age,
                            p_symp = prob_symp,
                            p_sus = prob_sus,
                            p_sus_w=prob_sus_w,
                            startdate=startdate,
                            enddate=enddate)

    #Solve dynamical system
    sol = solve(county_model.prob, BS3();tspan = (0,model_end_point),callback = cb,u0=u0, p=p,saveat = 1) #isoutofdomain=(u,p,t)->any(x->x<0,u)
    sol_baseline = solve(county_model_baseline.prob, BS3();tspan = (0,model_end_point),callback = cb_baseline,u0=u0, p=p, saveat = 1) # isoutofdomain=(u,p,t)->any(x->x<0,u)

    # Number of unscalled (for testing) cases over time by severity
    ι = KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol,10)
    ι_baseline = KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol_baseline,10)

    # Variants data
    kenya_variants = CSV.File("data/kenya_gisaid_variants.csv") |> DataFrame
    kenya_variants.date = [Date(d,DateFormat("dd/mm/yyyy")) for d in kenya_variants.date]
    kenya_variants.rel_time = [(d - Date(2020,2,20)).value for d in kenya_variants.date]

    idxs_alpha = kenya_variants.variant .== "Alpha"
    idxs_beta = kenya_variants.variant .== "Beta"
    idxs_delta = kenya_variants.variant .== "Delta"

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
        num_tests_age[:,a] = vcat(num_tests[:,a],fill(mean(num_tests[:,a]),(length(ι)-length(num_tests[:,a]))))  # use mean number of teste for forcast periods
    end

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

    n = size(ι ,1)

    delta_perc = [delta_perc;ones(n-length(delta_perc))]
    alpha_beta_perc = alpha_perc .+ beta_perc
    alpha_beta_perc = [alpha_beta_perc;zeros(n-length(alpha_beta_perc))]
    other_perc = 1.0 .- alpha_beta_perc .- delta_perc

    # scalling parameters
    p_test = fit_p[8]
    χ = fit_p[7]
    p_test_scale = 5e-4

    # Apply vaccine efficacy against disease
    υ = (1 .- VE_severe_disease_risk') .* ones(6)  # efficacy
    if baseline_scenario==0
        υ_baseline = ones(3)' .* ones(6) # No vaccine efficacy
    else
        υ_baseline = (1 .- VE_severe_disease_risk') .* ones(6)  # in case of a baselise scenarion with vaccination
    end

    l_inc = [inc .* υ for inc in ι]
    cases_wt = similar(l_inc)
    cases_αβ = similar(l_inc)
    cases_δ  = similar(l_inc)
    crit_incidence = similar(l_inc)
    sev_incidence = similar(l_inc)

    l_baseline_inc = [inc .* υ_baseline for inc in ι_baseline]
    cases_baseline_wt = similar(l_baseline_inc)
    cases_baseline_αβ = similar(l_baseline_inc)
    cases_baseline_δ  = similar(l_baseline_inc)
    crit_incidence_baseline = similar(l_baseline_inc)
    sev_incidence_baseline = similar(l_baseline_inc)

    for i in 1:size(l_inc,1)
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

        crit_incidence[i] = ICU_fit[1]*(cases_wt[i] .+ ICU_fit[2].*cases_αβ[i] .+ ICU_fit[3].*cases_δ[i])
        sev_incidence[i] = hosp_fit[1]*(cases_wt[i] .+ hosp_fit[2].*cases_αβ[i] .+ hosp_fit[3].*cases_δ[i])

        # baseline  iteration
        cases_baseline_wt[i] = sum(l_baseline_inc[i] .* other_perc[i],dims=2) # sum over dose
        cases_baseline_αβ[i] = sum(l_baseline_inc[i] .* alpha_beta_perc[i], dims=2)
        cases_baseline_δ[i]  = sum(l_baseline_inc[i] .* delta_perc[i], dims=2)

        PCR_num_wt = cases_baseline_wt[i] .* p_test_scale .* p_test
        PCR_denom_wt = PCR_num_wt .+ p_test_scale*((p_test/χ)*(N .- cases_baseline_wt[i]))
        p_PCR_pred_wt = PCR_num_wt./PCR_denom_wt
        cases_baseline_wt[i]  =  p_PCR_pred_wt .* num_tests_age[i,:]

        PCR_num_αβ = cases_baseline_αβ[i] .* p_test_scale .* p_test
        PCR_denom_αβ = PCR_num_αβ .+ p_test_scale*((p_test/χ)*(N .- cases_baseline_αβ[i]))
        p_PCR_pred_αβ = PCR_num_αβ./PCR_denom_αβ
        cases_baseline_αβ[i]  =  p_PCR_pred_αβ .* num_tests_age[i,:]

        PCR_num_δ = cases_baseline_δ[i] .* p_test_scale .* p_test
        PCR_denom_δ = PCR_num_wt .+ p_test_scale*((p_test/χ)*(N .- cases_baseline_δ[i]))
        p_PCR_pred_δ = PCR_num_δ./PCR_denom_δ
        cases_baseline_δ[i]  =  p_PCR_pred_δ .* num_tests_age[i,:]

        crit_incidence_baseline[i] = ICU_fit[1]*(cases_baseline_wt[i] .+ ICU_fit[2].*cases_baseline_αβ[i] .+ ICU_fit[3].*cases_baseline_δ[i])
        sev_incidence_baseline[i] = hosp_fit[1]*(cases_baseline_wt[i] .+ hosp_fit[2].*cases_baseline_αβ[i] .+ hosp_fit[3].*cases_baseline_δ[i])
    end

    crit_incidence =  Matrix(VectorOfArray([inc[:] for inc in crit_incidence])')  # convert to a matrix
    sev_incidence =  Matrix(VectorOfArray([inc[:] for inc in sev_incidence])')
    crit_incidence_baseline=  Matrix(VectorOfArray([inc[:] for inc in crit_incidence_baseline])')  # convert to a matrix
    sev_incidence_baseline =  Matrix(VectorOfArray([inc[:] for inc in sev_incidence_baseline])')

    for a in 1:6 # convolve with PCR array to get right timing of PCR detectable incidence
        crit_incidence[:,a] = KenyaCoVaccines.simple_conv(crit_incidence[:,a],KenyaCoVaccines.PCR_array)
        sev_incidence[:,a] = KenyaCoVaccines.simple_conv(sev_incidence[:,a],KenyaCoVaccines.PCR_array)
        crit_incidence_baseline[:,a] = KenyaCoVaccines.simple_conv(crit_incidence_baseline[:,a],KenyaCoVaccines.PCR_array)
        sev_incidence_baseline[:,a] = KenyaCoVaccines.simple_conv(sev_incidence_baseline[:,a],KenyaCoVaccines.PCR_array)
    end

    # Saving severe and critical disease incidence
    crit_dis_inc   = crit_incidence[janstarttime:end,:]
    severe_dis_inc = sev_incidence[janstarttime:end,:]
    crit_dis_inc_baseline   = crit_incidence_baseline[janstarttime:end,:]
    severe_dis_inc_baseline = sev_incidence_baseline[janstarttime:end,:]

    @save("scenario_$(scenario)/severe_cases/severe_cases_$(name).jld2",severe_dis_inc)
    @save("scenario_$(baseline_scenario)/severe_cases/severe_cases_$(name).jld2",severe_dis_inc_baseline)
    @save("scenario_$(scenario)/critical_cases/critical_cases_$(name).jld2",crit_dis_inc)
    @save("scenario_$(baseline_scenario)/critical_cases/critial_cases_$(name).jld2",crit_dis_inc_baseline)

    # Asymptomatic cases
    asymp_dis_inc = [inc .* (ones(3)'.* (1 .- prob_symp)) for inc in ι]
    asymp_dis_inc = Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in asymp_dis_inc])')  #sum over dose
    asymp_dis_inc = asymp_dis_inc[janstarttime:end,:]

    asymp_dis_inc_baseline  =  [inc .* (ones(3)'.* (1 .- prob_symp)) for inc in ι_baseline]
    asymp_dis_inc_baseline = Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in asymp_dis_inc_baseline])')  #sum over dose
    asymp_dis_inc_baseline = asymp_dis_inc_baseline[janstarttime:end,:]

    @save("scenario_$(scenario)/asymp_cases/asymp_cases_$(name).jld2",asymp_dis_inc)
    @save("scenario_$(baseline_scenario)/asymp_cases/asymp_cases_$(name).jld2",asymp_dis_inc_baseline)

    # Mild disease
    symp_dis_inc  =  [inc .* (ones(3)'.* prob_symp) for inc in ι]
    symp_dis_inc =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in symp_dis_inc])')  #sum over dose
    for a in 1:6
        symp_dis_inc[:,a] =  KenyaCoVaccines.simple_conv(symp_dis_inc[:,a],p_IH)  # convolve with time to hsop function to get right timing of Hospitalisation
    end
    symp_dis_inc = symp_dis_inc[janstarttime:end,:]
    mild_dis_inc = symp_dis_inc .- severe_dis_inc .- crit_dis_inc


    symp_dis_inc_baseline =  [inc .* (ones(3)'.* prob_symp) for inc in ι_baseline]
    symp_dis_inc_baseline =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in symp_dis_inc_baseline])')  #sum over dose
    for a in 1:6
        symp_dis_inc_baseline[:,a] =  KenyaCoVaccines.simple_conv(symp_dis_inc_baseline[:,a],p_IH)  # convolve with time to hsop function to get right timing of Hospitalisation
    end
    symp_dis_inc_baseline = symp_dis_inc_baseline[janstarttime:end,:]
    mild_dis_inc_baseline = symp_dis_inc_baseline .- severe_dis_inc_baseline .- crit_dis_inc_baseline

    @save("scenario_$(scenario)/mild_cases/mild_cases_$(name).jld2",mild_dis_inc)
    @save("scenario_$(baseline_scenario)/mild_cases/mild_cases_$(name).jld2",mild_dis_inc_baseline)

    # Number of unscalled deaths by age group over time
    IFR = all_IFR_fits[[nm[1] .== name for nm in all_IFR_fits]][1].IFR
    υ = (1 .- VE_severe_disease_risk') .* IFR # efficacy multiplied by risk of death by age, assumes only those with severe disease die
    υ_baseline = ones(3)' .* IFR # No vaccine efficacy

    death_inc =  [inc .* υ for inc in ι]
    death_inc =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in death_inc])')  #sum over dose
    for a in 1:6
        death_inc[:,a] =  KenyaCoVaccines.simple_conv(death_inc[:,a],p_ID)  # convolve with time to death function to get right timing of deaths
    end
    death_inc = death_inc[janstarttime:end,:] ./ obs_factor

    death_inc_baseline =  [inc .* υ for inc in ι_baseline]
    death_inc_baseline =  Matrix(VectorOfArray([sum(inc,dims=2)[:] for inc in death_inc_baseline])')  #sum over dose
    for a in 1:6
        death_inc_baseline[:,a] =  KenyaCoVaccines.simple_conv(death_inc_baseline[:,a],p_ID)  # convolve with time to death function to get right timing of deaths
    end
    death_inc_baseline = death_inc_baseline[janstarttime:end,:] ./ obs_factor

    @save("scenario_$(scenario)/deaths/deaths_$(name).jld2",death_inc)
    @save("scenario_$(baseline_scenario)/deaths/deaths_$(name).jld2",death_inc_baseline)

    # Number of vaccines administered over time
    # number of doses planned to cover about 30% of pop>18y

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

    @save("scenario_$(scenario)/doses/doses_$(name).jld2",num_doses)
    @save("scenario_$(baseline_scenario)/doses/doses_$(name).jld2",num_doses_baseline)

    @save("scenario_$(scenario)/doses_age/doses_$(name).jld2",complete_doses_age)
    @save("scenario_$(baseline_scenario)/doses_age/doses_$(name).jld2",complete_doses_age_baseline)

    # Plot vaccination impact for county

    ι = KenyaCoVaccines.get_incidence_time_array(sol,10)
    ι_baseline = KenyaCoVaccines.get_incidence_time_array(sol_baseline,10)
    p_test = fit_p[8]
    χ = fit_p[7]
    p_test_scale = 5e-4

    PCR = KenyaCoVaccines.simple_conv(ι,PCR_array)
    PCR_cases = sum(linelist_data_with_pos_neg.cases[:,linelist_data_with_pos_neg.areas .== up_county,:,:],dims = [2,3])[:,1,1,:]
    PCR_cases_pos = PCR_cases[model_start_time:end,1]

    p_PCR_pred_num = p_test_scale*p_test*PCR
    p_PCR_pred_denom = p_PCR_pred_num .+ p_test_scale*((p_test/χ)*(sum(N) .- PCR))
    p_PCR_pred = p_PCR_pred_num./p_PCR_pred_denom
    num_tests =  PCR_cases[model_start_time:end,1] + PCR_cases[model_start_time:end,2]

    # Might later change number of tests function from mean
    num_tests = vcat(num_tests,fill(mean(num_tests),(length(p_PCR_pred)-length(num_tests))))
    pred = num_tests .* p_PCR_pred

    PCR = KenyaCoVaccines.simple_conv(ι_baseline,PCR_array)
    p_PCR_pred_num = p_test_scale*p_test*PCR
    p_PCR_pred_denom = p_PCR_pred_num .+ p_test_scale*((p_test/χ)*(sum(N) .- PCR))
    p_PCR_pred = p_PCR_pred_num./p_PCR_pred_denom
    pred_baseline = num_tests .* p_PCR_pred

    pcr_pos_ma = KenyaCoVaccines.sma(PCR_cases_pos,7)  # add 7-day movinf average to the plot

    cover = 0
    if scenario == 2 || scenario == 5
        cover =30
    end
    if scenario == 3 || scenario == 6
        cover = 50
    end
    if scenario == 4 || scenario == 7
        cover = 70
    end

    baseline_cover = 0
    if baseline_scenario == 2 || baseline_scenario == 5
        baseline_cover =30
    end
    if baseline_scenario == 3 || baseline_scenario == 6
        baseline_cover = 50
    end
    if baseline_scenario == 4 || baseline_scenario == 7
        baseline_cover = 70
    end

    xticks_monthly = [(Date(y,m,1) - Date(2021,1,1)).value for m = 1:12,y = 2021:2022][1:end]
    xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2021:2022][1:end]


    plt = plot(pred[janstarttime:end],
            lab = "$(cover)% coverage",
            lw=1,
            legend = :topright,
            title = "PCR prevalence $(name)",
            xticks = (xticks_monthly,xticklabs_monthly),
            ylabel = "Daily cases",
            color=:blue,
            size = (800,600),dpi = 300)
    plot!(plt,pred_baseline[janstarttime:end],color=:red,lab = "Baseline:$(baseline_cover)% coverage",lw=1)
    plot!(plt,pcr_pos_ma[janstarttime:end],lw=1,color=:green,lab="7-day MA PCR prevalence")
    scatter!(plt,PCR_cases_pos[janstarttime:end],color=:green,lab="Observed PCR prevalence")

    savefig(plt,"scenario_$(scenario)/plots/vaccination_impact/predicted_PCR_incidence_$(name).png")

    # Save scalled PCR incidence for later plotting Kenyawide model fit
    if  baseline_scenario==0
        pred_baseline_fit =  pred_baseline[janstarttime:model_end_time]
        @save("scenario_$(baseline_scenario)/scaled_PCR_pred/scaled_PCR_pred_$(name).jld2",pred_baseline_fit)
    end

    return nothing
end


"""

Pedict vaccination impact of multiple counties using the `predict_vacc_county` method.
"""
function predict_vacc_for_multiple_counties(namelist::Vector{String},
                    linelist_data_with_pos_neg,serology_data,
                    M_Kenya_ho,M_Kenya_other,M_Kenya_school,M_Kenya_work,N_kenya;
                    PCR_array,
                    relative_testing_rate,
                    all_param_fits_newvariant,
                    rapid,
                    scenario,
                    baseline_scenario,
                    countynames,
                    prob_sus_w,
                    prob_sus,
                    prob_symp,
                    H_rate_by_age,
                    av_sero_init,
                    p_ID,
                    p_IH,
                    all_IFR_fits,
                    ICU_fit,
                    hosp_fit,
                    deaths_data,
                    priors,
                    T_to_full_efficacy,
                    VE_acquisition,  #Vaccine efficacy against acquisition by dose
                    VE_infectiousness, #Vaccine efficacy against transmission by dose
                    VE_severe_disease_risk,  #Vaccine efficacy against disease by dose
                    vacc_end_date, # date vaccination program ends
                    startdate,
                    enddate,
                    obs_factor)

    for name in namelist
        println("Starting projection for $(name)")
        pred = KenyaCoVaccines.predict_vacc_county(name,linelist_data_with_pos_neg,serology_data,
                            M_Kenya_ho,M_Kenya_other,M_Kenya_school,M_Kenya_work,N_kenya;
                            PCR_array=PCR_array,
                            relative_testing_rate=relative_testing_rate,
                            all_param_fits_newvariant = all_param_fits_newvariant,
                            rapid=rapid,
                            scenario=scenario,
                            baseline_scenario=baseline_scenario,
                            countynames=countynames,
                            prob_sus_w=prob_sus_w,
                            prob_sus=prob_sus,
                            prob_symp=prob_symp,
                            H_rate_by_age=H_rate_by_age,
                            av_sero_init=av_sero_init,
                            p_ID=p_ID,
                            p_IH=p_IH,
                            all_IFR_fits=all_IFR_fits,
                            ICU_fit=ICU_fit,
                            hosp_fit=hosp_fit,
                            deaths_data = deaths_data,
                            priors=priors,
                            T_to_full_efficacy=T_to_full_efficacy,
                            VE_acquisition=VE_acquisition,  #Vaccine efficacy against acquisition by dose
                            VE_infectiousness=VE_infectiousness, #Vaccine efficacy against transmission by dose
                            VE_severe_disease_risk=VE_severe_disease_risk,  #Vaccine efficacy against disease by dose
                            vacc_end_date = vacc_end_date, # date vaccination program ends
                            startdate = startdate,
                            enddate = enddate,
                            obs_factor=obs_factor)
    end
    return nothing
end


function sma(a::Vector, n::Int)  # moving average function
    vals = zeros(size(a,1) - (n-1))

    for i in 1:size(a,1) - (n-1)
      vals[i] = mean(a[i:i+(n-1)])
    end
    return vcat(missings(n-1),vals)
end
