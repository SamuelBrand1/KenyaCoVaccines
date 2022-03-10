"""
    predict_incidence_and_prev(;pred_incidence::Vector{Float64},
                                        pred_hospitalisations::Vector{Float64},
                                        pred_deaths::Vector{Float64},
										symptomaticrate::Float64,
										prop_critical::Float64,
										p_IS::Vector{Float64},
										Q_R::Vector{Float64},
                                        Q_ICUH::Vector{Float64},
                                        Q_HR::Vector{Float64},
                                        Q_ICUR::Vector{Float64},
										F_R::Vector{Float64})

Augment a prediction of incidence of infection, hospitalisation and death in an area with predictions of symptomatic/asymptomatic incidence, critical care (ICU) incidence

                                        # Arguments
- `pred_incidence::Vector{Float64}`: Daily predictions of new infections
- `pred_hospitalisations::Vector{Float64}`: Daily predictions of new hospitalisations.
- `pred_deaths::Vector{Float64}`: Daily predictions of new deaths.
- `symptomaticrate::Float64`: Broad population symptomatic rate.
- `prop_critical::Float64`: Proportion of hospitalised cases that require ICU treatment
- `p_IS::Vector{Float64}`: Distribution of days between infection and symptoms (for symptomatic cases).
- `Q_R::Vector{Float64}`: Upper tail function for days between infection and recovery
- `Q_ICUH::Vector{Float64}`: Upper tail function for days between arriving in ICU, surviving, and returning to less specialist wards
- `Q_HR::Vector{Float64}`: Upper tail function for days between arriving at hospital, remaining severe but not critical and then getting discharged
- `Q_ICUR::Vector{Float64}`: Upper tail function for days between arriving in ICU, surviving, returning to less specialist wards and then discharge
- `F_R::Vector{Float64}`: Distribution function for days until recovery from infection
"""
function predict_incidence_and_prev(;pred_incidence::Vector{Float64},
                                        pred_hospitalisations::Vector{Float64},
                                        pred_deaths::Vector{Float64},
										symptomaticrate::Float64,
										prop_critical::Float64,
										p_IS::Vector{Float64},
										Q_R::Vector{Float64},
                                        Q_ICUH::Vector{Float64},
                                        Q_HR::Vector{Float64},
                                        Q_ICUR::Vector{Float64},
										F_R::Vector{Float64})
	asymp_incidence = (1 - symptomaticrate) * simple_conv(pred_incidence, p_IS)
	symp_incidence = symptomaticrate * simple_conv(pred_incidence, p_IS)
	severe_incidence = (1 - prop_critical) * pred_hospitalisations
	crit_incidence = prop_critical * pred_hospitalisations
	death_incidence = pred_deaths
	asymp_prev = simple_conv(asymp_incidence, Q_R)
	critical_prev = simple_conv(crit_incidence, Q_ICUH)
	severe_prev = simple_conv(severe_incidence, Q_HR) .+ simple_conv(crit_incidence, Q_ICUR) .- critical_prev
	mild_prev = simple_conv(symp_incidence, Q_R) .- severe_prev .- critical_prev
	recovered_prev = simple_conv(asymp_incidence .+ symp_incidence, F_R)
	return (asymp_incidence = asymp_incidence,
            symp_incidence = symp_incidence,
            severe_incidence = severe_incidence,
            crit_incidence = crit_incidence,
            death_incidence = death_incidence,
            asymp_prev = asymp_prev,
            critical_prev = critical_prev,
            severe_prev = severe_prev,
            mild_prev = mild_prev,
            recovered_prev = recovered_prev)
end



"""
    predict_incidence_and_prev(modelfit;
                                    pred_incidence::Vector{Float64},
                                    pred_hospitalisations::Vector{Float64},
                                    pred_deaths::Vector{Float64},
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})

Convert a collected modelfit element into a hospitalisation prediction
"""
function predict_incidence_and_prev(modelfit;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})
	incidence = modelfit.mean_pred_incidence
	hospitalisations = modelfit.pred_hosps_ftc
	deaths = modelfit.pred_deaths_ftc
    name = modelfit.area
    prediction = predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
	return prediction, name
end



"""
    predict_incidence_and_prev(incidence_array,p_IH,p_ID,days;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})

Apply outcome predictions to each of a set of posterior draws for the incidence curves.
"""
function predict_incidence_and_prev(incidence_array,p_IH,p_ID,days;
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IS::Vector{Float64},
                                    Q_R::Vector{Float64},
                                    Q_ICUH::Vector{Float64},
                                    Q_HR::Vector{Float64},
                                    Q_ICUR::Vector{Float64},
                                    F_R::Vector{Float64})
        asymp_incidence = zeros(days,size(incidence_array,2))
        symp_incidence = similar(asymp_incidence)
        severe_incidence =similar(asymp_incidence)
        crit_incidence = similar(asymp_incidence)
        death_incidence = similar(asymp_incidence)
        asymp_prev = similar(asymp_incidence)
        critical_prev = similar(asymp_incidence)
        severe_prev = similar(asymp_incidence)
        mild_prev = similar(asymp_incidence)
        recovered_prev = similar(asymp_incidence)
        for k = 1:size(incidence_array,2)
                incidence = incidence_array[:,k]
                hospitalisations = KenyaSerology.simple_conv(incidence,0.001*p_IH)
                deaths = KenyaSerology.simple_conv(incidence,0.001*p_ID)
                prediction = KenyaSerology.predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
                asymp_incidence[:,k] = prediction.asymp_incidence
                symp_incidence[:,k]  = prediction.symp_incidence
                severe_incidence[:,k]  = prediction.severe_incidence
                crit_incidence[:,k]  = prediction.crit_incidence
                death_incidence[:,k]  = prediction.death_incidence
                asymp_prev[:,k]  = prediction.asymp_prev
                critical_prev[:,k]  = prediction.critical_prev
                severe_prev[:,k]  = prediction.severe_prev
                mild_prev[:,k]  = prediction.mild_prev
                recovered_prev[:,k]  = prediction.recovered_prev
        end
        return (asymp_incidence = asymp_incidence,
            symp_incidence = symp_incidence,
            severe_incidence = severe_incidence,
            crit_incidence = crit_incidence,
            death_incidence = death_incidence,
            asymp_prev = asymp_prev,
            critical_prev = critical_prev,
            severe_prev = severe_prev,
            mild_prev = mild_prev,
            recovered_prev = recovered_prev)
end

"""
	predict_incidence_and_prev(model::KenyaSerology.CoVAreaModel;
                                    proj_date::Date,
                                    cum_deaths::Integer,cum_deaths_date::Date,
                                    cum_hosps::Integer,cum_hosps_date::Date,
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IH,p_ID,
                                    d_incubation = LogNormal(1.644,0.363),#Lauer estimate
                                    d_duration_in_hosp = Weibull(1.572,17.819),#Haw et al
                                    d_ICUstay = Gamma(2.66,3.42),#Fit to IQR
                                    d_recovery = Exponential(2.4))

Output projections for relevant outcomes from a `KenyaSerology.CoVAreaModel` object.
"""
function predict_incidence_and_prev(model::KenyaSerology.CoVAreaModel;
                                    proj_date::Date,
                                    cum_deaths::Integer,cum_deaths_date::Date,
                                    cum_hosps::Integer,cum_hosps_date::Date,
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IH,p_ID,
                                    d_incubation = LogNormal(1.644,0.363),#Lauer estimate
                                    d_duration_in_hosp = Weibull(1.572,17.819),#Haw et al
                                    d_ICUstay = Gamma(2.66,3.42),#Fit to IQR
                                    d_recovery = Exponential(2.4))
    #Calculate lag Distributions
    #lag distributions
    p_IS = [cdf(d_incubation,t) - cdf(d_incubation,t-1) for t = 1:100]#infection to symptoms
    p_ICU = [cdf(d_ICUstay,t) - cdf(d_ICUstay,t-1) for t = 1:100]#ICU to leave ICU
    p_HR = [cdf(d_duration_in_hosp,t) - cdf(d_duration_in_hosp,t-1) for t = 1:100]#Hospital to discharge
    p_ICUR = KenyaSerology.simple_conv(p_ICU,p_HR)#ICU to discharge assuming its the sum of the two
    p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
    p_R = [cdf(d_recovery,t) - cdf(d_recovery,t-1) for t = 1:1000]
    #Upper tail functions
    Q_HR = vcat([1.],[1 - cumsum(p_HR)[t] for t = 1:100])
    Q_ICUH = vcat([1.],[1 - cumsum(p_ICU)[t] for t = 1:100])
    Q_ICUR = vcat([1.],[1 - cumsum(p_ICUR)[t] for t = 1:100])
    Q_R = vcat([1.],[1 - cumsum(p_R)[t] for t = 1:1000])
    F_R = 1 .- Q_R
    #Calculate fitted posterior mean daily incidence
    inc = KenyaSerology.incidence_across_samples(model,(proj_date-Date(2020,2,20)).value)
    incidence = KenyaSerology.create_credible_intervals(inc.true_incidence).mean_pred
    #Calculate unnormalised hosp and death prediction
    hospitalisations = KenyaSerology.simple_conv(incidence,p_IH)
    deaths = KenyaSerology.simple_conv(incidence,p_ID)
    #Calculate crude probability of hosp and death
    prob_hosp = cum_hosps/sum(hospitalisations[1:(cum_hosps_date - Date(2020,2,20)).value])
    prob_deaths = cum_deaths/sum(deaths[1:(cum_deaths_date - Date(2020,2,20)).value])

    hospitalisations .= prob_hosp.*hospitalisations
	deaths .= prob_deaths.*deaths
    name = model.areaname
    prediction = KenyaSerology.predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
	return prediction, name
end

function predict_incidence_and_prev(incidence::Vector;
                                    proj_date::Date,
                                    cum_deaths::Integer,cum_deaths_date::Date,
                                    cum_hosps::Integer,cum_hosps_date::Date,
                                    symptomaticrate::Float64,
                                    prop_critical::Float64,
                                    p_IH,p_ID,
                                    d_incubation = LogNormal(1.644,0.363),#Lauer estimate
                                    d_duration_in_hosp = Weibull(1.572,17.819),#Haw et al
                                    d_ICUstay = Gamma(2.66,3.42),#Fit to IQR
                                    d_recovery = Exponential(2.4))
    #Calculate lag Distributions
    #lag distributions
    p_IS = [cdf(d_incubation,t) - cdf(d_incubation,t-1) for t = 1:100]#infection to symptoms
    p_ICU = [cdf(d_ICUstay,t) - cdf(d_ICUstay,t-1) for t = 1:100]#ICU to leave ICU
    p_HR = [cdf(d_duration_in_hosp,t) - cdf(d_duration_in_hosp,t-1) for t = 1:100]#Hospital to discharge
    p_ICUR = KenyaSerology.simple_conv(p_ICU,p_HR)#ICU to discharge assuming its the sum of the two
    p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
    p_R = [cdf(d_recovery,t) - cdf(d_recovery,t-1) for t = 1:1000]
    #Upper tail functions
    Q_HR = vcat([1.],[1 - cumsum(p_HR)[t] for t = 1:100])
    Q_ICUH = vcat([1.],[1 - cumsum(p_ICU)[t] for t = 1:100])
    Q_ICUR = vcat([1.],[1 - cumsum(p_ICUR)[t] for t = 1:100])
    Q_R = vcat([1.],[1 - cumsum(p_R)[t] for t = 1:1000])
    F_R = 1 .- Q_R


    #Calculate unnormalised hosp and death prediction
    hospitalisations = KenyaSerology.simple_conv(incidence,p_IH)
    deaths = KenyaSerology.simple_conv(incidence,p_ID)
    #Calculate crude probability of hosp and death
    prob_hosp = cum_hosps/sum(hospitalisations[1:(cum_hosps_date - Date(2020,2,20)).value])
    prob_deaths = cum_deaths/sum(deaths[1:(cum_deaths_date - Date(2020,2,20)).value])

    hospitalisations .= prob_hosp.*hospitalisations
	deaths .= prob_deaths.*deaths
    prediction = KenyaSerology.predict_incidence_and_prev(pred_incidence=incidence,
                                        pred_hospitalisations=hospitalisations,
                                        pred_deaths=deaths,
                                        symptomaticrate=symptomaticrate,
                                        prop_critical=prop_critical,
                                        p_IS=p_IS,
                                        Q_R=Q_R,
                                        Q_ICUH=Q_ICUH,
                                        Q_HR=Q_HR,
                                        Q_ICUR=Q_ICUR,
                                        F_R=F_R)
	return prediction
end
