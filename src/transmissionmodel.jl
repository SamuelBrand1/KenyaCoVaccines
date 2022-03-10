

"""
        function build_vt(A::AbstractArray,extrapolation_bc::Union{Interpolations.Line,Interpolations.Flat})
Get a tuple of functions: 1) the cublic spline interpolation over the columns of the rate array `A`,  `extrapolation_bc` gives whether extrapolation is flat or linear.
"""
function build_vt(A::AbstractArray, extrapolation_bc::Union{Interpolations.Line,Interpolations.Flat})
        xs = 1:size(A, 1)
        extrap_array = [CubicSplineInterpolation(xs, A[:, a], extrapolation_bc = extrapolation_bc) for a = 1:size(A, 2)]
        return t -> [extrap(t) for extrap in extrap_array]
end

"""
        function build_vt_onedim(A::AbstractArray,extrapolation_bc::Union{Interpolations.Line,Interpolations.Flat})
Get a tuple of functions: 1) the cublic spline interpolation over the columns of the rate array `A`,  `extrapolation_bc` gives whether extrapolation is flat or linear.
"""
function build_vt_onedim(A::AbstractVector, extrapolation_bc::Union{Interpolations.Line,Interpolations.Flat})
        xs = 1:length(A)
        extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = extrapolation_bc)
        return t -> max(extrap(t),0.0)
end

"""
function build_prob_vaccine_modelling(N,contact_data,extrapolation_bc,T_ho,T_other,T_school;
                                                startdate::Date = Date(2020,2,20),
                                                enddate::Date = Date(2020,9,30))

Create an `ODEProblem`, which represents the age structured model dynamics for a single county
With the following vaccination states:
0d = no vaccine
1df = 1st dose at maximum protection
2df = maximum protection from 2 doses
"""
function build_prob_vaccine_modelling(N, vacc_rate_1, vacc_rate_2, T_ho, T_other, T_school, T_work;
        time_to_full_efficacy::Float64 = 14.0, # Time to full vaccine vacc_efficacy
        hosp_rate_by_age::Vector{Float64} = zeros(6), # Hospitalisation rate by age
        p_symp::Vector{Float64} = zeros(6), # probability of developing symptoms
        p_sus::Vector{Float64} = zeros(6),  # relative susceptibility by age
        p_sus_w::Vector{Float64} = zeros(6),  # relative susceptibility by age following first epissode
        startdate::Date = Date(2020, 2, 20),
        enddate::Date = Date(2020, 9, 30))

        vacc_rate_1 = vacc_rate_1[startdate.value:end, :]  # Trim the vacc data to cover only the days being fitted/modeled relative to time zero being 20th Feb 2020
        vacc_rate_2 = vacc_rate_2[startdate.value:end, :]  # Trim the vacc data to cover only the days being fitted/modeled relative to time zero being 20th Feb 2020

        extrapolation_bc = Interpolations.Flat()
        vt = KenyaCoVaccines.build_vt(vacc_rate_1, extrapolation_bc)  #first dose
        vt2 = KenyaCoVaccines.build_vt(vacc_rate_2, extrapolation_bc) #second dose

        agegroup_strs = vcat("0-19", "20-49", "50-59", "60-69", "70-79", "80_")
        disease_status_strs = ["S", "E", "A", "P", "M", "V", "R", "W", "F", "C"]
        vaccine_status = ["0d", "1d", "2d"]
        symbols = [Symbol(dis_str * "_" * age_str * "_" * vacc_str) for age_str in agegroup_strs, dis_str in disease_status_strs, vacc_str in vaccine_status]

        tspan = (0.0, (enddate - startdate).value)

        function kenyacov_ode(du, u, p, t)
                #Gather parameters
        
                #Basic parameters
                β₀, β_home, β_school, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3 = p
        
                #  Update efficacies in case of am integrator call back affecting them
                vacc_efficacy_acquisition = vcat(ve_ac1, ve_ac2, ve_ac3)  #Vaccine efficacy against acquisition by dose
                vacc_efficacy_infectiousness = vcat(ve_inf1, ve_inf2, ve_inf3)       #Vaccine efficacy against transmission by dose
                vacc_efficacy_severe_disease_risk = vcat(ve_dis1, ve_dis2, ve_dis3)
        
                # edit transmission to include changes due to variants
                β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_αvar, κ = time_scale_αvar, x₀ = mid_point_αvar))
                β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_δvar, κ = time_scale_δvar, x₀ = mid_point_δvar))
        
                # Contact rates
                ct = 1.0  # modelling begins after contact rates are back to normal.
        
                #Age dependent rates
                σω = ones(3)' .* p_sus_w  # relative susceptibility following first episode
                δ = ones(3)' .* p_symp   # probability of developing symptoms
                σ = ones(3)' .* p_sus    # relative susceptibility by age
        
                #Vaccine effectiveness by dose stage
                υ = (1 .- vacc_efficacy_severe_disease_risk') .* hosp_rate_by_age # efficacy multiplied by risk of severe disease by age
                # rate_vacc_progress_1 = 1 / time_to_full_efficacy
                # rate_vacc_progress_2 = 1 / time_to_full_efficacy
                rate_vacc_progress_1 = 1.0
                rate_vacc_progress_2 = 1.0
                γV = γM
        
                #Gather state
                S = @view u[:, 1, :]
                E = @view u[:, 2, :]
                A = @view u[:, 3, :]
                P = @view u[:, 4, :]
                M = @view u[:, 5, :]
                V = @view u[:, 6, :]
                R = @view u[:, 7, :]
                W = @view u[:, 8, :]
                F = @view u[:, 9, :]
                C = @view u[:, 10, :]
        
                # I = zeros(6)
        
                I = (1 - vacc_efficacy_infectiousness[1]) .* (ϵ .* A[:, 1] .+ P[:, 1] .+ M[:, 1] .+ V[:, 1])
                for dose = 2:3
                        I .+= (1 - vacc_efficacy_infectiousness[dose]) .* (ϵ .* A[:, dose] .+ P[:, dose] .+ M[:, dose] .+ V[:, dose])
                end
        
                λ = β₀ .* (β_home .* T_ho * I .+ β_school .* T_school * I .+ ct .* (β_other .* T_other * I .+ β_work .* T_work * I)) ./ N
        
        
                # Transmission events
                du[:, 1, :] .= -σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)  # Reduce force of infection by vaccine efficacy against acquisition
                du[:, 2, :] .= σ .* (S .+ σω .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ) .- α .* E
                du[:, 3, :] .= α .* E .* (1 .- δ) .- γA .* A
                du[:, 4, :] .= α .* E .* δ .- αP .* P
                du[:, 5, :] .= αP .* P .* (1 .- υ) .- γM .* M
                du[:, 6, :] .= αP .* P .* υ .- γV .* V
                du[:, 7, :] .= γA .* A .+ γM .* M .+ γV .* V .- ω .* R
                du[:, 8, :] .= ω .* R .- σω .* σ .* W .* ((1 .- vacc_efficacy_acquisition') .* λ)
                du[:, 9, :] .= σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)
                du[:, 10, :] .= σ .* (S .+ σω .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ)
        
                # Vaccination events
                du[:, 1:8, 1] .+= -u[:, 1:8, 1] .* (vt(t) ./ N) * rate_vacc_progress_1 # Vaccination independent of exposure status
                du[:, 1:8, 2] .+= u[:, 1:8, 1] .* (vt(t) ./ N) * rate_vacc_progress_1 .- u[:, 1:8, 2] .* (vt2(t) ./ N) * rate_vacc_progress_2
                du[:, 1:8, 3] .+= u[:, 1:8, 2] .* (vt2(t) ./ N) * rate_vacc_progress_2
        
                return nothing
        end

        ff = ODEFunction(kenyacov_ode; syms = symbols[:])
        u0 = zeros(6, 10, 3) # place holder
        return ODEProblem(ff, u0, tspan)
end
"""
function build_prob_vaccine_modelling_simplified(N, fully_vacc_rate::AbstractVector, T_ho, T_other, T_school, T_work;
        maximum_uptake = 0.8;#No age group has higher percentage uptake
        hosp_rate_by_age::Vector{Float64} = zeros(6), # Hospitalisation rate by age
        prop_symp_17::Vector{Float64} = [0.0638201, 0.0139376, 0.0195662, 0.0242427, 0.0631479, 0.093603, 0.0992078, 0.0895518, 0.0963006, 0.120111, 0.18606, 0.253572, 0.286403, 0.532203, 0.628416, 0.532232, 0.8299886], # probability of developing symptoms
        prob_sus_w_17::Vector{Float64} = [0.636022, 0.37745, 0.424016, 0.456355, 0.633717, 0.725298, 0.739908, 0.714375, 0.7324, 0.790052, 0.917996, 1.02082, 1.06435, 1.31636, 1.39356, 1.31638, 1.531594],  # relative susceptibility by age
        p_sus_w::Vector{Float64} = fill(0.16, 6),  # relative susceptibility by age following first epissode
        startdate::Date = Date(2020, 2, 20),
        enddate::Date = Date(2020, 9, 30))

Simplified version of building an ODEProblem with vaccination rates. The baseline assumption here is that strategy is always to offer to over 50s 
until some threshold `maximum_uptake` is reached (and people stop coming forwards from that group). Then vaccines are offered to over 18s.
"""
function build_prob_vaccine_modelling_simplified(N, fully_vacc_rate::AbstractVector, T_ho, T_other, T_school, T_work;
        maximum_uptake = 0.8,#No age group has higher percentage uptake
        hosp_rate_by_age::Vector{Float64} = zeros(6), # Hospitalisation rate by age
        prop_symp_17::Vector{Float64} = [0.0638201, 0.0139376, 0.0195662, 0.0242427, 0.0631479, 0.093603, 0.0992078, 0.0895518, 0.0963006, 0.120111, 0.18606, 0.253572, 0.286403, 0.532203, 0.628416, 0.532232, 0.8299886], # probability of developing symptoms
        prob_sus_w_17::Vector{Float64} = [0.636022, 0.37745, 0.424016, 0.456355, 0.633717, 0.725298, 0.739908, 0.714375, 0.7324, 0.790052, 0.917996, 1.02082, 1.06435, 1.31636, 1.39356, 1.31638, 1.531594],  # relative susceptibility by age
        p_sus_w::Vector{Float64} = fill(0.16, 6),  # relative susceptibility by age following first epissode
        startdate::Date = Date(2020, 2, 20),
        enddate::Date = Date(2020, 9, 30))

        #reduce susceptibility and symptom risk vectors to six age groups
        p_sus = vcat(mean(prob_sus_w_17[1:4]), mean(prob_sus_w_17[5:10]), mean(prob_sus_w_17[11:12]), mean(prob_sus_w_17[13:14]), mean(prob_sus_w_17[15:16]), prob_sus_w_17[17])
        p_symp = vcat(mean(prop_symp_17[1:4]), mean(prop_symp_17[5:10]), mean(prop_symp_17[11:12]), mean(prop_symp_17[13:14]), mean(prop_symp_17[15:16]), prop_symp_17[17])

        #Generate vaccination rate functions
        extrapolation_bc = Interpolations.Flat()
        t_0 = (startdate - Date(2020, 2, 20)).value
        vt2 = build_vt_onedim(fully_vacc_rate[t_0:end], extrapolation_bc) #second dose
        vt = t -> vt2(t + 56.0) #First dose became effective 8 weeks before second dose becomes effective

        #Names of states
        agegroup_strs = vcat("0-19", "20-49", "50-59", "60-69", "70-79", "80_")
        disease_status_strs = ["S", "E", "A", "P", "M", "V", "R", "W", "F", "C"]
        vaccine_status = ["0d", "1d", "2d"]
        symbols = [Symbol(dis_str * "_" * age_str * "_" * vacc_str) for age_str in agegroup_strs, dis_str in disease_status_strs, vacc_str in vaccine_status]

        #Simulation period
        tspan = (0.0, (enddate - startdate).value)

        #Age dependent rates
        σω = ones(3)' .* p_sus_w  # relative susceptibility following first episode
        δ = ones(3)' .* p_symp   # probability of developing symptoms
        σ = ones(3)' .* p_sus    # relative susceptibility by age


function kenyacov_ode(du, u, p, t)

        #Basic parameters
        β₀, β_home, β_school, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3 = p
        γV = γM
        ct = 1.0 #Assumed that contact rate back to baseline
        #  Update efficacies :CAREFUL in case of an integrator call back affecting them
        vacc_efficacy_acquisition = vcat(ve_ac1, ve_ac2, ve_ac3)  #Vaccine efficacy against acquisition by dose
        vacc_efficacy_infectiousness = vcat(ve_inf1, ve_inf2, ve_inf3)       #Vaccine efficacy against transmission by dose
        vacc_efficacy_severe_disease_risk = vcat(ve_dis1, ve_dis2, ve_dis3)

        #Vaccine effectiveness by dose stage
        υ = (1 .- vacc_efficacy_severe_disease_risk') .* hosp_rate_by_age # efficacy multiplied by risk of severe disease by age


        # edit transmission to include changes due to variants
        β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_αvar, κ = time_scale_αvar, x₀ = mid_point_αvar))
        β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_δvar, κ = time_scale_δvar, x₀ = mid_point_δvar))

        #Gather state
        S = @view u[:, 1, :]
        E = @view u[:, 2, :]
        A = @view u[:, 3, :]
        P = @view u[:, 4, :]
        M = @view u[:, 5, :]
        V = @view u[:, 6, :]
        R = @view u[:, 7, :]
        W = @view u[:, 8, :]
        F = @view u[:, 9, :]
        C = @view u[:, 10, :]

        #Proportion of each age-vaccination group in each disease state
        N_overfifty_unvac = sum(u[3:end, 1:8, 1])
        N_underfifty_unvac = sum(u[2, 1:8, 1])
        N_overfifty_1vac = sum(u[3:end, 1:8, 2])
        N_underfifty_1vac = sum(u[2, 1:8, 2])

        overfifty_unvac_prop = u[3:end, 1:8, 1] ./ N_overfifty_unvac
        overfifty_1vac_prop = N_overfifty_1vac ≈ 0.0 ? zeros(size(u[3:end, 1:8, 2])) : u[3:end, 1:8, 2] ./ N_overfifty_1vac
        underfifty_unvac_prop = u[2, 1:8, 1] ./ N_underfifty_unvac
        underfifty_1vac_prop = N_underfifty_1vac ≈ 0.0 ? zeros(size(u[2, 1:8, 2])) : u[2, 1:8, 2] ./ N_underfifty_1vac

        I = (1 - vacc_efficacy_infectiousness[1]) .* (ϵ .* A[:, 1] .+ P[:, 1] .+ M[:, 1] .+ V[:, 1])
        for dose = 2:3
                I .+= (1 - vacc_efficacy_infectiousness[dose]) .* (ϵ .* A[:, dose] .+ P[:, dose] .+ M[:, dose] .+ V[:, dose])
        end

        λ = β₀ .* (β_home .* T_ho * I .+ β_school .* T_school * I .+ ct .* (β_other .* T_other * I .+ β_work .* T_work * I)) ./ N

        # Transmission events
        du[:, 1, :] .= -σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)  # Reduce force of infection by vaccine efficacy against acquisition
        du[:, 2, :] .= σ .* (S .+ σω .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ) .- α .* E
        du[:, 3, :] .= α .* E .* (1 .- δ) .- γA .* A
        du[:, 4, :] .= α .* E .* δ .- αP .* P
        du[:, 5, :] .= αP .* P .* (1 .- υ) .- γM .* M
        du[:, 6, :] .= αP .* P .* υ .- γV .* V
        du[:, 7, :] .= γA .* A .+ γM .* M .+ γV .* V .- ω .* R
        du[:, 8, :] .= ω .* R .- σω .* σ .* W .* ((1 .- vacc_efficacy_acquisition') .* λ)
        du[:, 9, :] .= σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)
        du[:, 10, :] .= σ .* (S .+ σω .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ)

        #Determine if over 50s have reached their threshold uptake
        over_fifty_prop_vac_below_threshold = sum(u[3:end, 1:8, 2:3]) / sum(N[3:end]) < maximum_uptake
        over_fifty_prop_vac_below_2vac_threshold = sum(u[3:end, 1:8, 3]) / sum(N[3:end]) < maximum_uptake
        over_fifty_prop_vac_above_threshold_under_fifty_below_threshold = (~over_fifty_prop_vac_below_threshold) & (sum(u[2, 1:8, 2:3]) / N[2] < maximum_uptake)
        over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold = (~over_fifty_prop_vac_below_2vac_threshold) & (sum(u[2, 1:8, 3]) / N[2] < maximum_uptake)

        # Vaccination rates for over 50s 
        du[3:end, 1:8, 1] .+= -overfifty_unvac_prop .* vt(t) .* over_fifty_prop_vac_below_threshold# Vaccination independent of exposure status
        du[3:end, 1:8, 2] .+= overfifty_unvac_prop .* vt(t) .* over_fifty_prop_vac_below_threshold .- overfifty_1vac_prop .* vt2(t) .* over_fifty_prop_vac_below_2vac_threshold
        du[3:end, 1:8, 3] .+= overfifty_1vac_prop .* vt2(t) .* over_fifty_prop_vac_below_2vac_threshold
        # # Vaccination rates for under 50s 
        du[2, 1:8, 1] .+= -underfifty_unvac_prop .* vt(t) .* over_fifty_prop_vac_above_threshold_under_fifty_below_threshold  #Vaccination independent of exposure status
        du[2, 1:8, 2] .+= underfifty_unvac_prop .* vt(t) .* over_fifty_prop_vac_above_threshold_under_fifty_below_threshold .- underfifty_1vac_prop .* vt2(t) .* over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold
        du[2, 1:8, 3] .+= underfifty_1vac_prop .* vt2(t) .* over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold

        return nothing
end

        ff = ODEFunction(kenyacov_ode; syms = symbols[:])
        u0 = zeros(6, 10, 3) # place holder
        return ODEProblem(ff, u0, tspan)
end

function build_prob_vaccine_modelling_vaccine_escape(N,prop_vaccines, T_ho, T_other, T_school, T_work;
        maximum_uptake = 0.8,#No age group has higher percentage uptake
        hosp_rate_by_age::Vector{Float64} = zeros(6), # Hospitalisation rate by age
        prop_symp_17::Vector{Float64} = [0.0638201, 0.0139376, 0.0195662, 0.0242427, 0.0631479, 0.093603, 0.0992078, 0.0895518, 0.0963006, 0.120111, 0.18606, 0.253572, 0.286403, 0.532203, 0.628416, 0.532232, 0.8299886], # probability of developing symptoms
        prob_sus_w_17::Vector{Float64} = [0.636022, 0.37745, 0.424016, 0.456355, 0.633717, 0.725298, 0.739908, 0.714375, 0.7324, 0.790052, 0.917996, 1.02082, 1.06435, 1.31636, 1.39356, 1.31638, 1.531594],  # relative susceptibility by age
        p_sus_w::Vector{Float64} = fill(0.16, 6),  # relative susceptibility by age following first epissode
        startdate::Date = Date(2020, 2, 20),
        enddate::Date = Date(2020, 9, 30))

        #reduce susceptibility and symptom risk vectors to six age groups
        p_sus = vcat(mean(prob_sus_w_17[1:4]), mean(prob_sus_w_17[5:10]), mean(prob_sus_w_17[11:12]), mean(prob_sus_w_17[13:14]), mean(prob_sus_w_17[15:16]), prob_sus_w_17[17])
        p_symp = vcat(mean(prop_symp_17[1:4]), mean(prop_symp_17[5:10]), mean(prop_symp_17[11:12]), mean(prop_symp_17[13:14]), mean(prop_symp_17[15:16]), prop_symp_17[17])

        
        #Names of states
        agegroup_strs = vcat("0-19", "20-49", "50-59", "60-69", "70-79", "80_")
        disease_status_strs = ["S", "E", "A", "P", "M", "V", "R", "W", "F", "C"]
        vaccine_status = ["0d", "1d", "2d"]
        symbols = [Symbol(dis_str * "_" * age_str * "_" * vacc_str) for age_str in agegroup_strs, dis_str in disease_status_strs, vacc_str in vaccine_status]

        #Simulation period
        tspan = (0.0, (enddate - startdate).value)

        #Age dependent rates
        σω = ones(3)' .* p_sus_w  # relative susceptibility following first episode
        δ = ones(3)' .* p_symp   # probability of developing symptoms
        σ = ones(3)' .* p_sus    # relative susceptibility by age

        υ = zeros(6,3)
        
        function kenyacov_ode(du, u, p, t)
        
                #Basic parameters
                β₀, β_home, β_school, β_other, β_work, α, ϵ, αP, γA, γM, ω, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale, ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3, ve_dis1, ve_dis2, ve_dis3, red_prot, vt1, vt2 = p
                γV = γM
                ct = 1.0 #Assumed that contact rate back to baseline
                #  Update efficacies :CAREFUL in case of an integrator call back affecting them
                vacc_efficacy_acquisition = vcat(ve_ac1, ve_ac2, ve_ac3)  #Vaccine efficacy against acquisition by dose
                vacc_efficacy_infectiousness = vcat(ve_inf1, ve_inf2, ve_inf3)       #Vaccine efficacy against transmission by dose
                vacc_efficacy_severe_disease_risk = vcat(ve_dis1, ve_dis2, ve_dis3)
        
                #Vaccine effectiveness by dose stage
                # υ = (1 .- vacc_efficacy_severe_disease_risk') .* hosp_rate_by_age # efficacy multiplied by risk of severe disease by age
        
        
                # edit transmission to include changes due to variants
                β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_αvar, κ = time_scale_αvar, x₀ = mid_point_αvar))
                β₀ = β₀ * (1 + KenyaCoVaccines.gen_logistic(t; L = inc_R_δvar, κ = time_scale_δvar, x₀ = mid_point_δvar))
        
                #Gather state
                S = @view u[:, 1, :]
                E = @view u[:, 2, :]
                A = @view u[:, 3, :]
                P = @view u[:, 4, :]
                M = @view u[:, 5, :]
                V = @view u[:, 6, :]
                R = @view u[:, 7, :]
                W = @view u[:, 8, :]
                F = @view u[:, 9, :]
                C = @view u[:, 10, :]
        
                #Proportion of each age-vaccination group in each disease state
                N_overfifty_unvac = sum(u[3:end, 1:8, 1])
                N_underfifty_unvac = sum(u[2, 1:8, 1])
                N_overfifty_1vac = sum(u[3:end, 1:8, 2])
                N_underfifty_1vac = sum(u[2, 1:8, 2])
        
                overfifty_unvac_prop = u[3:end, 1:8, 1] ./ N_overfifty_unvac
                overfifty_1vac_prop = N_overfifty_1vac ≈ 0.0 ? zeros(size(u[3:end, 1:8, 2])) : u[3:end, 1:8, 2] ./ N_overfifty_1vac
                underfifty_unvac_prop = u[2, 1:8, 1] ./ N_underfifty_unvac
                underfifty_1vac_prop = N_underfifty_1vac ≈ 0.0 ? zeros(size(u[2, 1:8, 2])) : u[2, 1:8, 2] ./ N_underfifty_1vac
        
                I = (1 - vacc_efficacy_infectiousness[1]) .* (ϵ .* A[:, 1] .+ P[:, 1] .+ M[:, 1] .+ V[:, 1])
                for dose = 2:3
                        I .+= (1 - vacc_efficacy_infectiousness[dose]) .* (ϵ .* A[:, dose] .+ P[:, dose] .+ M[:, dose] .+ V[:, dose])
                end
        
                λ = β₀ .* (β_home .* T_ho * I .+ β_school .* T_school * I .+ ct .* (β_other .* T_other * I .+ β_work .* T_work * I)) ./ N
        
                # Transmission events
                du[:, 1, :] .= -σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)  # Reduce force of infection by vaccine efficacy against acquisition
                du[:, 2, :] .= σ .* (S .+ (red_prot .* σω) .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ) .- α .* E
                du[:, 3, :] .= α .* E .* (1 .- δ) .- γA .* A
                du[:, 4, :] .= α .* E .* δ .- αP .* P
                du[:, 5, :] .= αP .* P .* (1 .- υ) .- γM .* M
                du[:, 6, :] .= αP .* P .* υ .- γV .* V
                du[:, 7, :] .= γA .* A .+ γM .* M .+ γV .* V .- ω .* R
                du[:, 8, :] .= ω .* R .- (red_prot .* σω) .* σ .* W .* ((1 .- vacc_efficacy_acquisition') .* λ)
                du[:, 9, :] .= σ .* S .* ((1 .- vacc_efficacy_acquisition') .* λ)
                du[:, 10, :] .= σ .* (S .+ (red_prot .* σω) .* W) .* ((1 .- vacc_efficacy_acquisition') .* λ)
        
                #Determine if over 50s have reached their threshold uptake
                over_fifty_prop_vac_below_threshold = sum(u[3:end, 1:8, 2:3]) / sum(N[3:end]) < maximum_uptake
                over_fifty_prop_vac_below_2vac_threshold = sum(u[3:end, 1:8, 3]) / sum(N[3:end]) < maximum_uptake
                over_fifty_prop_vac_above_threshold_under_fifty_below_threshold = (~over_fifty_prop_vac_below_threshold) & (sum(u[2, 1:8, 2:3]) / N[2] < maximum_uptake)
                over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold = (~over_fifty_prop_vac_below_2vac_threshold) & (sum(u[2, 1:8, 3]) / N[2] < maximum_uptake)
        
                # Vaccination rates for over 50s 
                du[3:end, 1:8, 1] .+= -overfifty_unvac_prop .* prop_vaccines .* vt1 .* over_fifty_prop_vac_below_threshold# Vaccination independent of exposure status
                du[3:end, 1:8, 2] .+= overfifty_unvac_prop .* prop_vaccines .* vt1 .* over_fifty_prop_vac_below_threshold .- overfifty_1vac_prop .* prop_vaccines .* vt2 .* over_fifty_prop_vac_below_2vac_threshold
                du[3:end, 1:8, 3] .+= overfifty_1vac_prop .* prop_vaccines .* vt2 .* over_fifty_prop_vac_below_2vac_threshold
                # # Vaccination rates for under 50s 
                du[2, 1:8, 1] .+= -underfifty_unvac_prop .* prop_vaccines .* vt1 .* over_fifty_prop_vac_above_threshold_under_fifty_below_threshold  #Vaccination independent of exposure status
                du[2, 1:8, 2] .+= underfifty_unvac_prop .* prop_vaccines .* vt1 .* over_fifty_prop_vac_above_threshold_under_fifty_below_threshold .- underfifty_1vac_prop .* prop_vaccines .* vt2 .* over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold
                du[2, 1:8, 3] .+= underfifty_1vac_prop .* prop_vaccines .* vt2 .* over_fifty_prop_vac_above_threshold_under_fifty_below_2vac_threshold
        
                return nothing
        end

        ff = ODEFunction(kenyacov_ode; syms = symbols[:])
        u0 = zeros(6, 10, 3) # place holder
        return ODEProblem(ff, u0, tspan)
end



"""
  function get_incidence_time_array(sol)

Convert an ODESolution `sol` into a prediction of daily incidence by age group.
"""
function get_incidence_time_array_by_age(sol, state)
        Matrix(VectorOfArray(diff([sum(u[:, state, :], dims = [2, 3])[:] for u in sol.u]))')
end

function get_incidence_time_array_by_age_dose(sol, state)
        permutedims(Array(VectorOfArray(diff([u[:, state, :][:, :] for u in sol.u]))), [3, 1, 2])
end

function get_incidence_time_array(sol, state) # sum over age and doses
        inc = Matrix(VectorOfArray(diff([sum(u[:, state, :], dims = [2, 3])[:] for u in sol.u]))')
        tinc = sum(inc, dims = [2])[:]
        return tinc
end

function get_cummulative_incidence_time_array(sol, state) # sum over age and doses
        inc = Matrix(VectorOfArray([sum(u[:, state, :], dims = [2, 3])[:] for u in sol.u])')
        tinc = sum(inc, dims = [2])[:]
        return tinc
end
