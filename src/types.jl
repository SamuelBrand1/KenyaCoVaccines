"""
struct MCMCResults
Struct for holding the MCMC chain (a vector of NamedTuples which represent the accepted parameter sets), and the MCMC tree statistics
"""
struct MCMCResults
    chain
    logL::Vector{Float64}
    treestatistics::Vector{DynamicHMC.TreeStatisticsNUTS}
end


"""
    mutable struct CoVArea

Type for containing the information required to calculate the log-likelihood contribution of one sub-area.
"""

Base.@kwdef mutable struct CoVAreaModel
    areaname::String #Name of area
    PCR_cases::Matrix{Int64} = zeros(Int64,365,2) #Daily PCR cases: dim 1: day, dim 2: num pos and total
    sero_cases::Array{Int64,3} = zeros(Int64,365,17,2) #Daily serology tests: dim 1: day, dim2 : age group, dim 3: num pos and neg
    average_sero_init::Vector{Float64} # seroprevalence to replace roond three seroprevalence used for initial starting conditions for counties without any sero data in round 3
    init_sero ::Vector{Float64} # initial (round 3) seroprevelance by age
    deaths::Matrix{Int64} = zeros(Int64,365,17) #Daily deaths
    pop_data::NamedMatrix{Float64} # Population by age and county
    N  #Area's population size
    α::Float64 = 1/3.1  # 1/mean latent period
    ϵ::Float64 = 0.5    # Relative infectiousness of asymptomatic cases
    αP::Float64 = 1/1.9  # rate of movement from Pre symptomatic state
    γA::Float64 = 1/2.4 #Recovery rate from Asymptomatic state
    γM::Float64 = 1/2.4 #Recovery rate from Mild disease

    relative_testing_rate::Vector{Float64} = relative_testing_rate
    PCR_array::Vector{Float64} = PCR_array #Relative Sensitivity of PCR test on days post-infection
    PCR_sensitivity::Float64 = 1. #Base sensitivity for RT-PCR test
    PCR_specificity::Float64 = 0.995 #Base specificity for RT_PCR test
    baseline_sero_array::Vector{Float64} = rel_sero_array_26days #Relative sensitivity of serology test on days post-infection
    sero_sensitivity::Float64 = 0.927 #Base sensitivity for serology test
    sero_specificity::Float64 = 0.992 #Base specificity for serology test
    M_BB::Float64 = 194.0 #Shrinkage factor for BetaBinomial distribution that incorporates the uncertainty in the sensitivity (fitted to serology group)
    log_likelihood::Function = θ -> 0.
    log_priors::Function = θ -> 0.

    M_county_ho::Matrix{Float64}     # "M_county_ho: Normalised age mixing matrix, home setting"
    M_county_other::Matrix{Float64}  # "M_county_other: Normalised age mixing matrix, social/other settings"
    M_county_school::Matrix{Float64} # "M_county_school: Normalised age mixing matrix, schools setting"
    M_county_work::Matrix{Float64}   # "M_county_work: Normalised age mixing matrix, work setting"

    rapid::Bool = false  # introduce vaccination in a rapid manner: finish target doses in 180 days
    scenario::Int64      # Scenario to run: Choices from 1-4
    time_to_full_efficacy::Float64  = 14. # Time to full vaccine vacc_efficacy
    hosp_rate_by_age:: Vector{Float64} = zeros(17) # Hospitalisation rate by age
    p_symp:: Vector{Float64}  # risk of symptoms by age
    p_sus:: Vector{Float64} # susceptibility by age
    p_sus_w:: Vector{Float64} # relative susceptibility following primary infection

    ω::Float64 = 1/180 # Waning of natural protection, rate per day

    # Include vaccination parameters in the vector of parameters
    # Sources: vacce effiecay against infections and severe VE_severe_disease_risk
        # Voysey M, Costa Clemens SA, Madhi SA, Weckx LY, Folegatti PM, et al. Single-dose administration
        #and the influence of the timing of the booster dose on immunogenicity and efficacy of
        #ChAdOx1 nCoV-19 (AZD1222) vaccine: a pooled analysis of four randomised trials. Lancet
        #(2021). doi:10.1016/S0140-6736(21)00432-3.

        # https://www.medrxiv.org/content/10.1101/2021.09.28.21264260v1

    VE_acquisition = [0,0.639,0.599] #Vaccine efficacy against acquisition by dose
    VE_infectiousness = [0,0.18,0.63] #Vaccine efficacy against transmission by dose
    VE_severe_disease_risk = [0,0.76,0.813] #Vaccine efficacy against disease by dose

    model_start_time # start in 1st Dec 2020
    model_end_time  # end of model run for fitting
    janstarttime  # time t start including observations in likelihood
    marchendpoint  # 10th March when roudn 3 serosurvey data ends
    alpha_variant_time  # earliest time of alpha variant introduction
    delta_variant_time  # earliest time of delta variant introduction

    startdate # date to start model run for projections/fitting
    enddate  # date to end model run for projections/fitting

    prob::ODEProblem

    MCMC_results::Union{MCMCResults,Nothing} = nothing #This field gets added after the MCMC run
end


"""
    CoVAreaModel(name::String,loglk_func,log_prior_func;
                        case_data,sero_data,death_data,pop_data::NamedArray)
Constructor for a `CoVAreaModel` struct using country data objects.
"""
function CoVAreaModel(name::String, log_likelihood, log_priors;
    linelist_data_with_pos_neg,
    serology_data,
    average_sero_init,
    deaths,
    pop_data::NamedMatrix,
    M_county_ho::Matrix{Float64},
    M_county_other::Matrix{Float64},
    M_county_school::Matrix{Float64},
    M_county_work::Matrix{Float64},
    rapid::Bool,
    scenario::Int64,
    startdate,
    enddate) # relative susceptibility by age following first episeode

    # Data input area
    #--------------------------------------------------------------------------------------------------------------------------
    time_to_full_efficacy = 14.0  # assumed 14 days post vaccination. This is the duration of stay in the partial vaccine protection states D1p and D2p

    #Set times, all times after model_start_time are in reference to the model_start_time
    model_start_time = (Date(2020, 12, 1) - Date(2020, 2, 20)).value # start in 1st Dec 2020
    model_end_time = (Date(2021, 9, 24) - Date(2020, 12, 1)).value  # end of model run for fitting
    janstarttime = (Date(2021, 1, 1) - Date(2020, 12, 1)).value # time t start including observations in likelihood
    marchendpoint = (Date(2021, 3, 10) - Date(2020, 12, 1)).value # 10th March when roudn 3 serosurvey data ends
    alpha_variant_time = (Date(2021, 2, 1) - Date(2020, 12, 1)).value  # earliest time of alpha variant introduction
    delta_variant_time = (Date(2021, 4, 1) - Date(2020, 12, 1)).value  # earliest time of delta variant introduction

    model_days = (Date(2022, 6, 30) - Date(2020, 12, 1)).value   # Spread vaccination till 2022 June
    pre_vacc_period = (Date(2021, 3, 6) - Date(2020, 12, 1)).value # vaccination starts on 6th March 2021


    # relative susceptibility by age
    prob_sus_w_17 = [0.636022, 0.37745, 0.424016, 0.456355, 0.633717, 0.725298, 0.739908, 0.714375, 0.7324, 0.790052, 0.917996, 1.02082, 1.06435, 1.31636, 1.39356, 1.31638, 1.531594]
    # probability of developing symptoms
    prop_symp_17 = [0.0638201, 0.0139376, 0.0195662, 0.0242427, 0.0631479, 0.093603, 0.0992078, 0.0895518, 0.0963006, 0.120111, 0.18606, 0.253572, 0.286403, 0.532203, 0.628416, 0.532232, 0.8299886]  # From Matt Keeling
    # relative susceptibility forllowinf first episode
    p_sus_w = ones(6) .* 0.16

    #reduce susceptibility and symptom risk vectors to six age groups
    p_sus = vcat(mean(prob_sus_w_17[1:4]), mean(prob_sus_w_17[5:10]), mean(prob_sus_w_17[11:12]), mean(prob_sus_w_17[13:14]), mean(prob_sus_w_17[15:16]), prob_sus_w_17[17])
    p_symp = vcat(mean(prop_symp_17[1:4]), mean(prop_symp_17[5:10]), mean(prop_symp_17[11:12]), mean(prop_symp_17[13:14]), mean(prop_symp_17[15:16]), prop_symp_17[17])

    # Age specific risk of severe disease
    #Estimate from CDC COVID-19 Response Team. Severe Outcomes Among Patients with Coronavirus Disease 2019 (COVID-19) - United States, February 12-March 16, 2020. MMWR Morb. Mortal. Wkly. Rep. 69, 343–346 (2020).
    #Model age groups: [0-19], [20-49], [50-59], [60-69], [70-79], 80+
    hosp_rate_CDC = [mean([1.6, 2.5]),#0-19 year olds
        mean([14.3, 20.8]),#20-44 yos
        mean([21.2, 28.3]),#45-54
        mean([20.5, 30.1]),#55-64
        mean([28.6, 43.5]),#65-74
        mean([30.5, 58.7]),#75-84
        mean([31.3, 70.3])] ./ 100 #85+

    # number of doses planned to cover about 30% of pop>18y
    if scenario == 0 || scenario == 1
        doses_Kenya = 0.0
    else
        doses_Kenya = (1.25 + 9.76 + 4.9) .* 1000000  #Number of doses planned for Kenya over the ynext one year
    end

    # END Data input area
    #--------------------------------------------------------------------------------------------------------------------------

    #Hsopitalisation rate by age. These rates dont actually affect the hospitalisations incidence reported from this work. We fit an IHR and ICU rates based on modelled incidence and observed hospital and ICU admissions
    H_rate_by_age = [hosp_rate_CDC[1], hosp_rate_CDC[2], hosp_rate_CDC[3], hosp_rate_CDC[4], hosp_rate_CDC[5], hosp_rate_CDC[6]]
    hosp_rate_by_age = H_rate_by_age .* 0.1 # Taking 10% of Verity rates. The IFR estimate from out model and death data was 10% of the expected from Verity data

    uprname = uppercase(name)
    PCR_cases = sum(linelist_data_with_pos_neg.cases[:, linelist_data_with_pos_neg.areas.==uprname, :, :], dims = [2, 3])[:, 1, 1, :]
    sero_cases = sum(serology_data.sero[:, serology_data.areas.==uprname, :, :], dims = [2])[:, 1, :, :]

    # reduce sero cases to 6 age groups
    p_vec = [[1, 2, 3, 4], [5, 6, 7, 8, 9, 10], [11, 12], [13, 14], [15, 16], [17]]
    sero_cases_collaped = zeros(Int64, size(sero_cases, 1), 6, 2)
    for a in 1:6
        sero_cases_collaped[:, a:a, :] .= sum(sero_cases[:, p_vec[a][1]:p_vec[a][end], :], dims = 2)
    end
    sero_cases = sero_cases_collaped

    # Put the initial serology as the default value --- this is based on Ngere et al and therefore more recent that rd 2 serology

    init_sero = average_sero_init

    # #Change from default if more than 10 tests
    # if sum(sero_cases[(model_start_time:(model_start_time+marchendpoint)), :, :]) == 0  # no round 3 data for some counties in the group of 36 counties without significant serodata
    #     init_sero = average_sero_init  # replace with the average
    # else
    #     u_sero = sero_cases[(model_start_time:(model_start_time+marchendpoint)), :, :]
    #     u_sero = sum(u_sero, dims = 1)[1, :, :]
    #     init_sero = u_sero[:, 1] ./ sum(u_sero, dims = 2)[:]
    # end

    # if isnan(init_sero[1])  # action if first value is isnan: replace with next succeding available value
    #     a = 1
    #     while isnan(init_sero[1])
    #         init_sero[1] = init_sero[a+1]
    #         a += 1
    #     end
    # end

    # for a in 2:6
    #     if isnan(init_sero[a])
    #         init_sero[a] = init_sero[a-1] # use the seroprevalence in the preceding age groups
    #     end
    # end

    #

    deaths = deaths.deaths[:, deaths.areas.==uprname, :][:, 1, :][:, :]

    #dates = [Date(2020,12,1) + Day(k) for k = 1:size(PCR_cases[(model_start_time:(model_start_time + model_end_time)),:],1)]
    N = pop_data[:, name]
    N = vcat(sum(N[1:4]), sum(N[5:10]), sum(N[11:12]), sum(N[13:14]), sum(N[15:16]), N[17]) # reduce to six age categories

    age_dist_vacc, doses_daily = KenyaCoVaccines.age_dist([name], doses_Kenya, pop_data, model_days, pre_vacc_period, rapid)
    vacc_rate_1, vacc_rate_2 = KenyaCoVaccines.vacc_rate_data(doses_daily, model_days, pre_vacc_period, age_dist_vacc, pop_data, name, scenario, rapid)

    # Reduce contact matrix from 17 to 6 age groups
    M_county_ho = KenyaCoVaccines.reduce_age_categories(M_county_ho, pop_data, name)
    M_county_other = KenyaCoVaccines.reduce_age_categories(M_county_other, pop_data, name)
    M_county_school = KenyaCoVaccines.reduce_age_categories(M_county_school, pop_data, name)
    M_county_work = KenyaCoVaccines.reduce_age_categories(M_county_work, pop_data, name)

    prob = build_prob_vaccine_modelling(N, vacc_rate_1, vacc_rate_2,
        M_county_ho, M_county_other, M_county_school, M_county_work;
        time_to_full_efficacy = time_to_full_efficacy,
        hosp_rate_by_age = hosp_rate_by_age,
        p_symp = p_symp,
        p_sus = p_sus,
        p_sus_w = p_sus_w,
        startdate = startdate,
        enddate = enddate)

    return KenyaCoVaccines.CoVAreaModel(areaname = name,
        log_likelihood = log_likelihood,
        log_priors = log_priors;
        PCR_cases = PCR_cases,
        sero_cases = sero_cases,
        init_sero = init_sero,
        average_sero_init = average_sero_init,
        deaths = deaths,
        startdate = startdate,
        enddate = enddate,
        model_start_time = model_start_time,
        model_end_time = model_end_time,
        janstarttime = janstarttime,
        alpha_variant_time = alpha_variant_time,
        delta_variant_time = delta_variant_time,
        marchendpoint = marchendpoint,
        prob = prob,
        pop_data = pop_data,
        N = N,
        M_county_ho = M_county_ho,
        M_county_other = M_county_other,
        M_county_school = M_county_school,
        M_county_work = M_county_work,
        rapid = rapid,
        scenario = scenario,
        time_to_full_efficacy = time_to_full_efficacy,
        hosp_rate_by_age = hosp_rate_by_age,
        p_symp = p_symp,
        p_sus = p_sus,
        p_sus_w = p_sus_w)

end

function (model::KenyaCoVaccines.CoVAreaModel)(θ)
    model.log_likelihood(θ,model,1/180) + model.log_priors(θ)
end

function change_prob!(model::CoVAreaModel, fully_vacc_rate::AbstractVector;
    startdate::Date,
    enddate::Date,
    maximum_uptake = 0.8)
    _prob = build_prob_vaccine_modelling_simplified(model.N, fully_vacc_rate, model.M_county_ho, model.M_county_other, model.M_county_school, model.M_county_work;
        maximum_uptake = maximum_uptake,
        startdate = startdate,
        enddate = enddate)
    model.prob = _prob
    return nothing
end

function change_prob_immune_escape!(model::CoVAreaModel, vaccine_prop;
    startdate::Date,
    enddate::Date,
    maximum_uptake = 0.8)
    _prob = build_prob_vaccine_modelling_vaccine_escape(model.N, vaccine_prop, model.M_county_ho, model.M_county_other, model.M_county_school, model.M_county_work;
        maximum_uptake = maximum_uptake,
        startdate = startdate,
        enddate = enddate)
    model.prob = _prob
    return nothing
end
