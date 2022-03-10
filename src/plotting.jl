
"""
    function plot_PCR_point_est_fit(θ_fitted::NamedTuple,model::CoVModel,k::Integer)

Plot the fit of a point estimate with and without vaccination. Incidence is not adjusted for number of test and probability of test positive among cases

"""
function plot_PCR_point_prediction(params,init,model,model_novacc,PCR_array,θ_fitted,all_IHRIFR_fits,p_ID,p_IH,scenario,startdate,name)

    #Test per day to scale prediction to match numbers of tests
    tests_per_day = sum(model.case_data[:,:,2],dims = 2)[:] #sum(model.case_data[:,:,:],dims = [2,3])[:]
    mean_tests_per_day_after_sept = mean(tests_per_day[300:end]) #+2000
    tests_per_day = tests_per_day[(startdate-Date(2020,2,20)).value:end]  # for when data is available after start date
    tests_per_day_new =  vcat(tests_per_day,ones(366-length(tests_per_day)).*mean_tests_per_day_after_sept)
    #tests_per_day_new = ones(366).* mean(tests_per_day)  # use average number of tests in the past for the year 2021

    immune_escape_variant_time = Float64((Date(2021,9,1) - Date(2021,1,31)).value) #  September 1st, time to implement variant immunity escape
    function immune_escape_variant!(integrator)
        # Reduce vaccine efficacies by 50%
        integrator.p[45:59] *= 0.5

        # Start with 'variant_E' exposed number of individuals
        variant_E = 50
        prop_vacc = sum(integrator.u[:,1:7,:],dims=2)[:,1,:] ./ sum(integrator.u[:,1:7,:],dims=[2,3])[:]
        integrator.u[:,2,:] .= variant_E .* prop_vacc # number exposed in each age group and vaccinatin status
        integrator.u[:,1,:] .= integrator.u[:,1,:] .- integrator.u[:,2,:] #substract E from S

        # move 50% from R -> S
        loss_susceptible = 0.5 .* integrator.u[:,7,:]
        integrator.u[:,1,:] .+= loss_susceptible
        integrator.u[:,7,:] .-= loss_susceptible
    end
    immune_escape_variant_cb = PresetTimeCallback([immune_escape_variant_time],immune_escape_variant!,save_positions=(false,false))

   # Update callbacks based on scenario
    if  scenario>2  # For Scenario 3 and 4
     model.cb = CallbackSet(model.cb, immune_escape_variant_cb)
     model_novacc.cb = CallbackSet(model_novacc.cb, immune_escape_variant_cb)
    end

    #Solve dynamical system
    p = copy(params)
    sol= solve(model.prob, Vern9(); callback = model.cb, u0=init, p=p, saveat = 1)
    p = copy(params)
    sol_novacc= solve(model_novacc.prob, Vern9(); callback = model_novacc.cb, u0=init, p=p, saveat = 1)

    incidence = KenyaCoVaccines.get_incidence_time_array(sol)
    PCR_pred = VectorOfArray([KenyaCoVaccines.simple_conv(incidence[:,j],PCR_array) for j = 1:size(incidence,2)])
    total_PCR_pred = sum(PCR_pred[:,:],dims = 2)[:]
    beta_bin_pred = θ_fitted.χ_PCR_array.*total_PCR_pred./((θ_fitted.χ_PCR_array-1).*total_PCR_pred .+ sum(model.N))
    pred = beta_bin_pred.*tests_per_day_new[1:size(PCR_pred,1)]

    incidence = KenyaCoVaccines.get_incidence_time_array(sol_novacc)
    PCR_pred = VectorOfArray([KenyaCoVaccines.simple_conv(incidence[:,j],PCR_array) for j = 1:size(incidence,2)])
    total_PCR_pred= sum(PCR_pred[:,:],dims = 2)[:]
    beta_bin_pred = θ_fitted.χ_PCR_array.*total_PCR_pred./((θ_fitted.χ_PCR_array-1).*total_PCR_pred .+ sum(model.N))
    pred_novacc = beta_bin_pred.*tests_per_day_new[1:size(PCR_pred,1)]

    cum_inc_vac_children = VectorOfArray([sum(u[1:4,8,:]) for u in sol.u])
    cum_inc_novac_children = VectorOfArray([sum(u[1:4,8,:]) for u in sol_novacc.u])
    cum_inc_vac_adults = VectorOfArray([sum(u[5:14,8,:]) for u in sol.u])
    cum_inc_novac_adults = VectorOfArray([sum(u[5:14,8,:]) for u in sol_novacc.u])
    cum_inc_vac_elderly = VectorOfArray([sum(u[15:end,8,:]) for u in sol.u])
    cum_inc_novac_elderly = VectorOfArray([sum(u[15:end,8,:]) for u in sol_novacc.u])

    #xticks_monthly = [(Date(y,m,1) - Date(2020,2,20)).value for m = 1:12,y = 2020:2021][1:end]
    #xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2020:2021][1:end]
    xticks_monthly = [(Date(y,m,1) - Date(2021,1,1)).value for m = 1:12,y = 2021][1:end]
    xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2021][1:end]

    pcr_pos = sum(model.case_data[:,:,1],dims = 2)[(Date(2021,1,31)-Date(2020,2,20)).value : end]
    pcr_pos_ma = KenyaCoVaccines.sma(pcr_pos,7)  # add 7-day movinf average to the plot
    plt = plot(pred,
            lab = "With vaccination",
            lw=1,
            legend = :topright,
            title = "PCR prevalence $(name)",
            xticks = (xticks_monthly,xticklabs_monthly),
            ylabel = "Daily cases",
            color=:blue,
            size = (800,600),dpi = 300);
    plot!(plt,pred_novacc,color=:red,lab = "Without vaccination",lw=1);
    plot!(plt,pcr_pos_ma,lw=1,color=:green,lab="7-day MA PCR prevalence");
    scatter!(plt,pcr_pos,color=:green,lab="Observed PCR prevalence")



    plt2=plot(cum_inc_vac_children,lab = "vaccines: 0-19 yrs old",legend = :bottomright,title = "Cumulative incidence $(name)",xlabel = "days",lw=3);
    plot!(plt2,cum_inc_novac_children,lab = "no vaccines: 0-19 yrs old",lw=3)

    plt3=plot(cum_inc_vac_adults,lab = "vaccines: 20-64 yrs old",legend = :bottomright,title = "Cumulative incidence $(name)",xlabel = "days",lw=3);
    plot!(plt3,cum_inc_novac_adults,lab = "no vaccines: 20-64 yrs old",lw=3)

    plt4=plot(cum_inc_vac_elderly,lab = "vaccines: 65+ yrs old",legend = :bottomright,title = "Cumulative incidence $(name)",xlabel = "days",lw=3);
    plot!(plt4,cum_inc_novac_elderly,lab = "no vaccines: 65+ yrs old",lw=3)

    p = copy(params)
    vacc_efficacy_severe_disease_risk = p[55:59]

    IHR = all_IHRIFR_fits[[nm[1] .== name for nm in all_IHRIFR_fits]][1].IHR

    υ = (1 .- vacc_efficacy_severe_disease_risk') .* IHR # efficacy multiplied by risk of severe disease by age
    υ_nv = ones(5)' .* IHR # No vaccine efficacy

    incidence = diff([u[:,8,:] for u in sol.u])
    dis_inc =  [inc .* υ for inc in incidence]
    disease_inc =  [sum(inc) for inc in dis_inc]  #sum over age and dose
    disease_inc =  KenyaCoVaccines.simple_conv(disease_inc,p_IH)  # convolve with time to hsop function to get right timing of Hospitalisation

    incidence_nv = diff([u[:,8,:] for u in sol_novacc.u])
    dis_inc_nv =  [inc .* υ_nv for inc in incidence_nv]
    disease_inc_nv = [sum(inc) for inc in dis_inc_nv]  #sum over age and dose
    disease_inc_nv =  KenyaCoVaccines.simple_conv(disease_inc_nv,p_IH)  # convolve with time to hsop function to get right timing of Hospitalisation

    l = min(length(disease_inc_nv),length(disease_inc))
    severe_disease_averted = cumsum(disease_inc_nv[1:l] .- disease_inc[1:l])  # cummulative

    cum_incidence =  VectorOfArray([sum(u[:,8,:]) for u in sol.u])
    cum_incidence_nv =  VectorOfArray([sum(u[:,8,:]) for u in sol_novacc.u])
    cases_averted  = cum_incidence_nv[1:(l-1)] .- cum_incidence[1:(l-1)] # cummulative

    # deaths averted
    IFR = all_IHRIFR_fits[[nm[1] .== name for nm in all_IHRIFR_fits]][1].IFR
    υ = (1 .- vacc_efficacy_severe_disease_risk') .* IFR # efficacy multiplied by risk of death by age, assumes only those with severe disease die
    υ_nv = ones(5)' .* IFR # No vaccine efficacy

    incidence = diff([u[:,8,:] for u in sol.u])
    d_inc =  [inc .* υ for inc in incidence]
    death_inc =  [sum(inc) for inc in d_inc]  #sum over age and dose
    deaths =  KenyaCoVaccines.simple_conv(death_inc,p_ID)  # convolve with time to death function to get right timing of deaths
    cumm_deaths = cumsum(deaths)

    incidence = diff([u[:,8,:] for u in sol_novacc.u])
    d_inc =  [inc .* υ_nv for inc in incidence]
    death_inc =  [sum(inc) for inc in d_inc]  #sum over age and dose
    deaths =  KenyaCoVaccines.simple_conv(death_inc,p_ID)  # convolve with time to death function to get right timing of deaths
    cumm_deaths_novacc = cumsum(deaths)

    deaths_averted = cumm_deaths_novacc .- cumm_deaths

    plot(severe_disease_averted,dpi = 300,lab="Severe disease", title = "Cummulative Cases/Disease/Deaths averted $(name)",
     lw=3,color=:red,legend =:topleft,ylabel= "Disease/Deaths",xticks=nothing)
    plot!(cases_averted,lab="Cases", ylabel= "Cases",xticks=nothing,legend =:topleft, lw=3,color=:blue)
    plt5=plot!(twinx(),deaths_averted,lab="Deaths",size=(800,600),xticks=(xticks_monthly,xticklabs_monthly), lw=3, color=:black,legend =:topright)

    @save("scenario_$(scenario)/predicted_states/solutions_$(name).jld2",sol)
    @save("scenario_$(scenario)/deaths_averted/deaths_averted_$(name).jld2",deaths_averted)
    @save("scenario_$(scenario)/cases_averted/cases_averted_$(name).jld2", cases_averted)
    @save("scenario_$(scenario)/severe_disease_averted/severe_disease_averted_$(name).jld2",severe_disease_averted)

    @save("scenario_$(scenario)/deaths/deaths_$(name).jld2",cumm_deaths)
    @save("scenario_$(scenario)/cases/cases_$(name).jld2", cum_incidence)
    @save("scenario_$(scenario)/severe_disease/severe_disease_$(name).jld2",disease_inc)

    @save("scenario_$(scenario)/deaths/deaths_$(name).jld2",cumm_deaths_novacc)
    @save("scenario_$(scenario)/cases/cases_$(name).jld2", cum_incidence_nv)
    @save("scenario_$(scenario)/severe_disease/severe_disease_$(name).jld2",disease_inc_nv)


    return plt, plt2, plt3, plt4, plt5
end

"""
Function to plot vaccine uptake for a given county
    Plots vaccine doses over time (daily) and the cumulative number of doses at a given time by age group
"""

function plot_vaccination_uptake(time,dose1_data,dose2_data)
    xticks_monthly = [(Date(y,m,1) - Date(2021,1,31)).value for m = 1:12,y = 2021][1:end]
    xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2021][1:end]

    vacc_data = dose1_data .+ dose2_data

    plt1= plot(sum(vacc_data, dims=2), lw=3, lab="Dose 1 + Dose 2",xticks = (xticks_monthly,xticklabs_monthly),legend = :topleft, ylims=[0,maximum(sum(vacc_data, dims=2)[:])+500])
    plot!(plt1,sum(dose1_data,dims=2),ylabel="Number of doses", xlabel = "Time",lab="Dose 1", lw=3)

    cum_doses = zeros(size(dose1_data))
    for  age in 1:17
     cum_doses[:,age] = cumsum(vacc_data[:,age])
    end

    xticks_age = 1:17
    xticklabs_age= vcat([string(a*5)*"-"*string(a*5 + 4) for a = 0:15],"80+")

    plt2=bar(1:17,cum_doses[time,:], xticks = (xticks_age,xticklabs_age),
        legend= :none,
        xlabel = "Age group",
        title = "Number of doses administered as at $(Date(2021,1,31) + Dates.Day(time))",
        ylabel = "Cummulative number of doses",
        size = (800,600))

    return plt1, plt2
end


"""
    function gather_monthly_sero_data(sero_data,years,date_first_day::Date)
Calculate the all age group sero-positivity with Jeffery 95% CIs for the point estimate. X tick positions are given
relative to the start date of the data set: `date_first_day`.
"""


function gather_monthly_sero_data(sero_data,years,date_first_day::Date)
        data_dates = [date_first_day + Day(k-1) for k=1:size(sero_data,1)]
        monthly_x_pos = Float64[]
        monthly_seropos = Float64[]
        lb_seropos = Float64[]
        ub_seropos = Float64[]


        for yr in years
                for mnth in 1:12
                        daysinmonthyear = findall((month.(data_dates).== mnth).&(year.(data_dates).== yr) )
                        pos_sero = sum(sero_data[daysinmonthyear,:,1])
                        neg_sero = sum(sero_data[daysinmonthyear,:,2])
                        if pos_sero + neg_sero > 0
                                p_sero = pos_sero/(pos_sero+neg_sero)
                                push!(monthly_x_pos,mean(daysinmonthyear))
                                push!(monthly_seropos,p_sero)
                                d = Beta(pos_sero +0.5,neg_sero +0.5)
                                push!(lb_seropos,p_sero-invlogcdf(d,log(0.025)))
                                push!(ub_seropos,invlogcdf(d,log(0.975)) - p_sero )
                        end
                end
        end
        return (monthly_x_pos=monthly_x_pos,monthly_seropos=monthly_seropos,lb_seropos=lb_seropos,ub_seropos=ub_seropos)
end

"""
    function plot_PCR_point_est_fit(θ_fitted::NamedTuple,model::CoVModel,k::Integer)
Plot the fit of a point estimate `θ_fitted` of the underlying parameters to the PCR data for the `k`th
area of the `CoVModel` type `model`.
"""
function plot_PCR_point_est_fit(θ_fitted::NamedTuple,model,k;startdate = Date(2020,2,20),enddate = Date(2021,3,1))
        #  Get parameters
        u0 = KenyaCoVaccines.create_initial_conditions(θ_fitted.E₀_array[k],model.subareas[k].N)
        p = deepcopy(model.subareas[k].p)
        p[1] = θ_fitted.β₀
        p[41] = θ_fitted.β_home
        p[42] = θ_fitted.β_school
        p[43] = θ_fitted.β_other
        p[44] = θ_fitted.β_work

        #Solve dynamical system
        l_solve = (enddate - startdate).value
        sol = solve(model.subareas[k].prob, BS3();callback = model.subareas[k].cb, u0=u0, p=p, saveat = 1,tspan = (0,l_solve))
        incidence = KenyaCoVaccines.get_incidence_time_array(sol)
        PCR_pred = VectorOfArray([KenyaCoVaccines.simple_conv(incidence[:,j],model.PCR_array) for j = 1:size(incidence,2)])
        total_PCR_pred = sum(PCR_pred[:,:],dims = 2)[:]

        #Find days with and without neg. PCR test result data and scale prediction to match numbers of tests
        total_days_with_neg_test = trues(length(total_PCR_pred))
        days_without_neg_test = findall(model.subareas[k].case_data[:,1,2] .== -1)
        days_with_neg_test = findall(model.subareas[k].case_data[:,1,2] .>= 0)
        total_days_with_neg_test[days_without_neg_test] .= false
        tests_per_day_from_data = sum(model.subareas[k].case_data[:,:,2],dims = 2)[:]
        mean_tests_per_day_after_sept = mean(tests_per_day_from_data[194:end])

        tests_per_day_for_pred = fill(mean_tests_per_day_after_sept,length(total_PCR_pred))
        tests_per_day_for_pred[days_with_neg_test] .= tests_per_day_from_data[days_with_neg_test]


        neg_bin_pred = total_PCR_pred.*θ_fitted.p_test_array[k].*model.relative_testing_rate[1:size(PCR_pred,1)]
        beta_bin_pred = θ_fitted.χ_PCR_array[k].*total_PCR_pred./((θ_fitted.χ_PCR_array[k]-1).*total_PCR_pred .+ sum(model.subareas[k].N))
        beta_bin_pred = beta_bin_pred.*tests_per_day_for_pred
        pred = neg_bin_pred.*(.~total_days_with_neg_test) .+ beta_bin_pred.*total_days_with_neg_test

        #Plot

        xticks_monthly = [(Date(y,m,1) - Date(2020,2,20)).value for m = 1:12,y = 2020:2021][3:end]
        xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2020:2021][3:end]
        plt = scatter(sum(model.subareas[k].case_data[:,:,1],dims = 2),
                lab = "Data",
                legend = :topleft,
                title = "KenyaCoV: fit to $(model.subareas[k].name) PCR data",
                xticks = (xticks_monthly,xticklabs_monthly),
                ylabel = "Daily cases",
                size = (800,600),dpi = 250);
        plot!(plt,pred,lab = "KenyaCoV fit",lw=2);
        return plt
end

# function plot_PCR_point_est_fit(E₀_var,θ_fitted::NamedTuple,model,k;startdate = Date(2020,2,20),enddate = Date(2021,3,1))
#         #  Get parameters
#         u0 = KenyaCoVaccines.create_initial_conditions(θ_fitted.E₀_array[k],model.subareas[k].N)
#         p = deepcopy(model.subareas[k].p)
#         p[1] = θ_fitted.β₀
#         p[41] = θ_fitted.β_home
#         p[42] = θ_fitted.β_school
#         p[43] = θ_fitted.β_other
#         p[44] = θ_fitted.β_work
#         #New variant callback
#         newvariant_time = Float64((Date(2021,1,31) - Date(2020,2,20)).value)
#     function affect_newvariant!(integrator)
#         loss_susceptible = 0.75.*integrator.u[:,7,1]
#         integrator.u[:,1,1] .+= loss_susceptible
#         integrator.u[:,7,1] .-= loss_susceptible
#         u0_nv = KenyaCoVaccines.create_initial_conditions(E₀_var,model.subareas[k].N)
#         integrator.u[:,2,1] .= u0_nv[:,2,1]
#         integrator.u[:,3:6,1] .= 0.
#     end
#     newvariant_cb = PresetTimeCallback([newvariant_time],affect_newvariant!,save_positions=(false,false))

#     cb = CallbackSet(model.subareas[k].cb,newvariant_cb)

#         #Solve dynamical system
#         l_solve = (enddate - startdate).value
#         sol = solve(model.subareas[k].prob, BS3();callback = cb, u0=u0, p=p, saveat = 1,tspan = (0,l_solve))
#         incidence = KenyaCoVaccines.get_incidence_time_array(sol)
#         PCR_pred = VectorOfArray([KenyaCoVaccines.simple_conv(incidence[:,j],model.PCR_array) for j = 1:size(incidence,2)])
#         total_PCR_pred = sum(PCR_pred[:,:],dims = 2)[:]

#         #Find days with and without neg. PCR test result data and scale prediction to match numbers of tests
#         total_days_with_neg_test = trues(length(total_PCR_pred))
#         days_without_neg_test = findall(model.subareas[k].case_data[:,1,2] .== -1)
#         days_with_neg_test = findall(model.subareas[k].case_data[:,1,2] .>= 0)
#         total_days_with_neg_test[days_without_neg_test] .= false
#         tests_per_day_from_data = sum(model.subareas[k].case_data[:,:,2],dims = 2)[:]
#         mean_tests_per_day_after_sept = mean(tests_per_day_from_data[194:end])

#         tests_per_day_for_pred = fill(mean_tests_per_day_after_sept,length(total_PCR_pred))
#         tests_per_day_for_pred[days_with_neg_test] .= tests_per_day_from_data[days_with_neg_test]


#         neg_bin_pred = total_PCR_pred.*θ_fitted.p_test_array[k].*model.relative_testing_rate[1:size(PCR_pred,1)]
#         beta_bin_pred = θ_fitted.χ_PCR_array[k].*total_PCR_pred./((θ_fitted.χ_PCR_array[k]-1).*total_PCR_pred .+ sum(model.subareas[k].N))
#         beta_bin_pred = beta_bin_pred.*tests_per_day_for_pred
#         pred = neg_bin_pred.*(.~total_days_with_neg_test) .+ beta_bin_pred.*total_days_with_neg_test

#         #Plot

#         xticks_monthly = [(Date(y,m,1) - Date(2020,2,20)).value for m = 1:12,y = 2020:2021][3:end]
#         xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2020:2021][3:end]
#         plt = scatter(sum(model.subareas[k].case_data[:,:,1],dims = 2),
#                 lab = "Data",
#                 legend = :topleft,
#                 title = "KenyaCoV: fit to $(model.subareas[k].name) PCR data",
#                 xticks = (xticks_monthly,xticklabs_monthly),
#                 ylabel = "Daily cases");
#         plot!(plt,pred,lab = "KenyaCoV fit",lw=2);
#         return plt
# end


"""
    function plot_sero_point_est_fit(θ_fitted::NamedTuple,model::CoVModel,k::Integer)
Plot the fit of a point estimate `θ_fitted` of the underlying parameters to the serological data for the `k`th
area of the `CoVModel` type `model`.
"""

function plot_sero_point_est_fit(θ_fitted::NamedTuple,model,k;startdate = Date(2020,2,20),enddate = Date(2021,3,1))
        #  Get parameters
        u0 = KenyaCoVaccines.create_initial_conditions(θ_fitted.E₀_array[k],model.subareas[k].N)
        p = deepcopy(model.subareas[k].p)
        p[1] = θ_fitted.β₀
        p[41] = θ_fitted.β_home
        p[42] = θ_fitted.β_school
        p[43] = θ_fitted.β_other
        p[44] = θ_fitted.β_work

        #Find the best fit sero-detection curve
        sero_array = vcat(model.baseline_sero_array[1:30],[(1-θ_fitted.seroreversionrate)^day for day in 1:length(model.baseline_sero_array[31:end])])

        #Solve dynamical system
        l_solve = (enddate - startdate).value
        sol = solve(model.subareas[k].prob, BS3();callback = model.subareas[k].cb, u0=u0, p=p, saveat = 1,tspan = (0,l_solve))
        incidence = KenyaCoVaccines.get_incidence_time_array(sol)
        sero_pred = VectorOfArray(model.sero_sensitivity*[KenyaCoVaccines.simple_conv(incidence[:,j],sero_array) for j = 1:size(incidence,2)])
        total_prop_inf = [sum(u[:,8,:]) for u in sol.u]./sum(model.subareas[k].N)

        #Rescale sero_pred numbers by pop. size in each age group
        pred_sero_prop_by_age = sero_pred./repeat(model.subareas[k].N',size(sero_pred,1),1)

        if sum(model.subareas[k].sero_data) > 0
                #Rescale the sero proportion prediction by the observed testing rate per age group
                num_tests_by_age = sum(model.subareas[k].sero_data[:,:,:],dims = [1,3])[:]
                num_tests = sum(model.subareas[k].sero_data[:,:,:])
                pred_sero_prop_overall = pred_sero_prop_by_age*(num_tests_by_age./num_tests)

                #Gather monthly sero-pos data
                monthly_sero_data = gather_monthly_sero_data(model.subareas[k].sero_data,[2020,2021],Date(2020,2,21))

                xticks_monthly = [(Date(y,m,1) - Date(2020,2,20)).value for m = 1:12,y = 2020:2021][3:end]
                xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2020:2021][3:end]

                plt = scatter(monthly_sero_data.monthly_x_pos,monthly_sero_data.monthly_seropos,
                                yerror = (monthly_sero_data.lb_seropos,monthly_sero_data.ub_seropos),
                                ms = 8,
                                color = :green,
                                xticks = (xticks_monthly,xticklabs_monthly),
                                lab = "Data",
                                legend = :topleft,
                                ylims = (-0.02,1),
                                ylabel = "Prop. of population",
                                title = "KenyaCoV: fit to $(model.subareas[k].name) sero data",
                                size = (800,600),dpi = 250);
                plot!(plt,pred_sero_prop_overall,
                        lab = "KenyaCoV sero. pred.",
                        lw = 2,
                        color = :green);

                plot!(plt,total_prop_inf,
                lab = "KenyaCoV cum. inf. pred.",
                lw = 2,
                color = :red);


                return plt
        else
                xticks_monthly = [(Date(y,m,1) - Date(2020,2,20)).value for m = 1:12,y = 2020:2021][3:end]
                xticklabs_monthly = [monthname(m)[1:3]*"-"*string(y-2000) for m = 1:12,y = 2020:2021][3:end]

                plt = plot(total_prop_inf,
                        xticks = (xticks_monthly,xticklabs_monthly),
                        lab = "KenyaCoV cum. inf. pred.",
                        legend = :topleft,
                        ylabel = "Prop. of population",
                        lw = 2,
                        title = "KenyaCoV: fit to $(model.subareas[k].name), no sero data",
                        color = :red,
                        size = (800,600),dpi = 250);


                return plt
        end
end

"""
    function plot_sero_age_for_monthyear_point_est(θ_fitted::NamedTuple,model,k,mnth,yr,date_first_day::Date)
Plot the fit of a point estimate `θ_fitted` of the underlying parameters to the age profile of serological data for the `k`th
area of the `CoVModel` type `model` in month `mnth`, year `yr`.
"""

function plot_sero_age_for_monthyear_point_est(θ_fitted::NamedTuple,model,k,mnth,yr,date_first_day::Date;startdate = Date(2020,2,20),enddate = Date(2021,3,1))
        #  Get parameters
        u0 = KenyaCoVaccines.create_initial_conditions(θ_fitted.E₀_array[k],model.subareas[k].N)
        p = deepcopy(model.subareas[k].p)
        p[1] = θ_fitted.β₀
        p[41] = θ_fitted.β_home
        p[42] = θ_fitted.β_school
        p[43] = θ_fitted.β_other
        p[44] = θ_fitted.β_work

        #Find the best fit sero-detection curve
        sero_array = vcat(model.baseline_sero_array[1:30],[(1-θ_fitted.seroreversionrate)^k for k in 1:length(model.baseline_sero_array[31:end])])

        #Solve dynamical system
        l_solve = (enddate - startdate).value
        sol = solve(model.subareas[k].prob, BS3();callback = model.subareas[k].cb, u0=u0, p=p, saveat = 1,tspan = (0,l_solve))
        incidence = KenyaCoVaccines.get_incidence_time_array(sol)
        sero_pred = VectorOfArray(model.sero_sensitivity*[KenyaCoVaccines.simple_conv(incidence[:,j],sero_array) for j = 1:size(incidence,2)])
        total_prop_inf = [sum(u[:,8,1]) for u in sol.u]./sum(model.subareas[k].N)

        #Rescale sero_pred numbers by pop. size in each age group and select the month-year
        data_dates = [date_first_day + Day(k-1) for k=1:size(model.subareas[k].sero_data,1)]
        daysinmonthyear = findall((month.(data_dates).== mnth).&(year.(data_dates).== yr) )

        pred_sero_prop_by_age = sero_pred./repeat(model.subareas[k].N',size(sero_pred,1),1)
        #Find average predicted sero-pos by age group in month
        pred_sero_prop_by_age = mean(pred_sero_prop_by_age[daysinmonthyear,:],dims = 1)[:]
        #Gather data for selected month-year

        pos_sero_by_age = sum(model.subareas[k].sero_data[daysinmonthyear,:,1],dims = 1)[:]
        neg_sero_by_age = sum(model.subareas[k].sero_data[daysinmonthyear,:,2],dims = 1)[:]
        prop_sero_by_age = pos_sero_by_age./(pos_sero_by_age .+ neg_sero_by_age)

        data_err = zeros(17)
        for a in 1:length(data_err)
                if pos_sero_by_age[a] + neg_sero_by_age[a] > 0
                        d = Beta(pos_sero_by_age[a]+0.5,neg_sero_by_age[a]+0.5)
                        P_a = pos_sero_by_age[a]/(pos_sero_by_age[a]+neg_sero_by_age[a])
                        data_err[a] = P_a - invlogcdf(d,log(0.025))
                end
        end


        plt = groupedbar(hcat(prop_sero_by_age,pred_sero_prop_by_age),
                        yerr = hcat(data_err,zeros(17)),
                        lab = ["KNBS data" "KenyaCoV Pred."],
                        xticks = (1:17,vcat([string((a-1)*5)*"-"*string(a*5-1) for a = 1:16 ],"80+") ),
                        legend = :topleft,
                        size = (700,400),dpi =250,
                        ylabel = "Prop. of age group",
                        title = "KenyaCoV: fit to $(model.subareas[k].name) sero data in $(monthname(mnth)), $(yr)");
        return plt
end

"""
function plot_sero_age_for_monthyear_point_est_fit_with_comparison(θ_fitted::NamedTuple,model,k,mnth,yr,date_first_day::Date;
                                                    comparison_sero_array,comparison_sero_label::String)
Plot the fit of a point estimate `θ_fitted` of the underlying parameters to the age profile of serological data for the `k`th
area of the `CoVModel` type `model` in month `mnth`, year `yr`. Also, point an additional comparison array `comparison_sero_array`.
"""
function plot_sero_age_for_monthyear_point_est_fit_with_comparison(θ_fitted::NamedTuple,model,k,mnth,yr,date_first_day::Date;
                                                    comparison_sero_array,comparison_sero_label::String,startdate = Date(2020,2,20),enddate = Date(2021,3,1))
        #  Get parameters
        u0 = KenyaCoVaccines.create_initial_conditions(θ_fitted.E₀_array[k],model.subareas[k].N)
        p = deepcopy(model.subareas[k].p)
        p[1] = θ_fitted.β₀
        p[41] = θ_fitted.β_home
        p[42] = θ_fitted.β_school
        p[43] = θ_fitted.β_other
        p[44] = θ_fitted.β_work

        #Find the best fit sero-detection curve
        sero_array = vcat(model.baseline_sero_array[1:30],[(1-θ_fitted.seroreversionrate)^k for k in 1:length(model.baseline_sero_array[31:end])])

        #Solve dynamical system
        l_solve = (enddate - startdate).value
        sol = solve(model.subareas[k].prob, BS3();callback = model.subareas[k].cb, u0=u0, p=p, saveat = 1, tspan = (0,l_solve))
        incidence = KenyaCoVaccines.get_incidence_time_array(sol)
        sero_pred = VectorOfArray(model.sero_sensitivity*[KenyaCoVaccines.simple_conv(incidence[:,j],sero_array) for j = 1:size(incidence,2)])
        total_prop_inf = [sum(u[:,8,:]) for u in sol.u]./sum(model.subareas[k].N)

        #Rescale sero_pred numbers by pop. size in each age group and select the month-year
        data_dates = [date_first_day + Day(k-1) for k=1:size(model.subareas[k].sero_data,1)]
        daysinmonthyear = findall((month.(data_dates).== mnth).&(year.(data_dates).== yr) )

        pred_sero_prop_by_age = sero_pred./repeat(model.subareas[k].N',size(sero_pred,1),1)
        #Find average predicted sero-pos by age group in month
        pred_sero_prop_by_age = mean(pred_sero_prop_by_age[daysinmonthyear,:],dims = 1)[:]


        plt = groupedbar(hcat(pred_sero_prop_by_age,comparison_sero_array),
                        lab = ["KenyaCoV Pred." comparison_sero_label],
                        xticks = (1:17,vcat([string((a-1)*5)*"-"*string(a*5-1) for a = 1:16 ],"80+") ),
                        legend = :topleft,
                        size = (700,400),dpi =250,                        ylabel = "Prop. of age group",
                        title = "KenyaCoV: fit to $(model.subareas[k].name) sero data in $(monthname(mnth)), $(yr)");
        return plt
end


function movingaverage(X::Vector,numofele::Int)
    BackDelta = div(numofele,2)
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(X)
    Y = similar(X)
    for n = 1:len
        lo = max(1,n - BackDelta)
        hi = min(len,n + ForwardDelta)
        Y[n] = mean(X[lo:hi])
    end
    return Y
end
