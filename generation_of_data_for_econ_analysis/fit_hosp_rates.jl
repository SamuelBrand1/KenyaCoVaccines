using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks, LogExpFunctions
using RecursiveArrayTools, Interpolations, MCMCChains, Suppressor, Interpolations
import KenyaCoVaccines

##Load methods and group data into correct age groups

include("generation_methods.jl");
include("inference_methods.jl");
@load("data/p_ID.jld2")
@load("data/p_IH.jld2")
@load("data/cleaned_linelist20211109_deaths_c_age__date_of_lab_confirmation.jld2")

##Gather infection predictions by county

fitfiles = readdir("modelfits", join = true)
inf_predictions = Vector{NamedTuple}(undef, 47)
prev_inf_predictions = Vector{Vector}(undef, 47)

for (i, filename) in enumerate(fitfiles)
    model = @suppress_err load(filename)["model"]
    KenyaCoVaccines.change_prob!(model, zeros(500); startdate = Date(2020, 12, 1), enddate = Date(2021, 9, 27))
    println("Making infection predictions for county $(model.areaname)")
    inf_predictions[i] = gather_post_mean_infection_rates_by_type_novac(model, Date(2021, 9, 27))
    prev_inf_predictions[i] = model.init_sero .* mean(model.MCMC_results.chain[:init_scale]) .* model.N
end

##

@load("data/fitted_mort_by_county_age_variant.jld2", mort_scales, mort_age, mort_var)

##unscaled hosp pred infections

unscaled_hosp_pred_inf = zeros(300)
unscaled_hosp_pred_inf_wt = zeros(300)
unscaled_hosp_pred_inf_ab = zeros(300)
unscaled_hosp_pred_inf_delta = zeros(300)

temp_hosp_pred = zeros(300)
for i = 1:47
    for a = 1:6
        temp_hosp_pred .= inf_predictions[i].infections_by_age_first_pred[a] .* mort_scales[i] .* mort_age[a]
        temp_hosp_pred .+= 0.15 .* inf_predictions[i].infections_by_age_reinf_pred[a] .* mort_scales[i] .* mort_age[a]
        unscaled_hosp_pred_inf_wt .+= temp_hosp_pred .* inf_predictions[i].wt_frq_pred
        unscaled_hosp_pred_inf_ab .+= (mort_var[1]) .* temp_hosp_pred .* inf_predictions[i].alpha_beta_frq_pred
        unscaled_hosp_pred_inf_delta .+= (mort_var[2]) .* temp_hosp_pred .* inf_predictions[i].delta_frq_pred

        temp_hosp_pred .*= (inf_predictions[i].wt_frq_pred .+ (mort_var[1]) .* inf_predictions[i].alpha_beta_frq_pred .+ (mort_var[2]) .* inf_predictions[i].delta_frq_pred)
        unscaled_hosp_pred_inf .+= temp_hosp_pred

    end
end

##Plots

scatter(sum(deaths_data.deaths[deaths_data.dates.>=Date(2020, 12, 25), :, :], dims = [2, 3])[:], lab = "")
plot!(unscaled_hosp_pred_inf, lab = "", lw = 3)

##Gather hospitalisation data

kenya_MoH = CSV.File("data/2022-01-10_MOH_Kenya_user_tweets__EXTRACTED_cleaned_dates.csv") |> DataFrame
kenya_MoH.dates

kenya_MoH.rel_time = [(d - Date(2020, 12, 1)).value for d in kenya_MoH.dates]
idxs_hosp_valid = kenya_MoH.in_hospital .> 0
idxs_ICU_valid = (kenya_MoH.intensive_care_unit .> 0) .& (kenya_MoH.intensive_care_unit .< 4000)
idxs_suppoxy_valid = kenya_MoH.supplementary_oxygen .> 0

##

###
### CONSTRUCT DELAY DISTRIBUTIONS
###

d_incubation = LogNormal(1.644, 0.363)#Lauer estimate
# d_duration_in_hosp = Weibull(1.572, 17.819)#Haw et al
d_duration_in_hosp = Gamma(8.0472378052287681, 0.0017022772867586302 * 365)
mean(d_duration_in_hosp)
std(d_duration_in_hosp)
# shape = 53.273947475158401, scale = 0.0003599898843749761

# d_ICUstay = Gamma(2.66, 3.42)#Fit to IQR
d_ICUstay = Gamma(53.273947475158401, 0.0003599898843749761 * 365) #From Cameline
mean(d_ICUstay)
std(d_ICUstay)
d_recovery = Exponential(2.4)
#lag distributions
p_IS = [cdf(d_incubation, t) - cdf(d_incubation, t - 1) for t = 1:100]#infection to symptoms
p_ICU = [cdf(d_ICUstay, t) - cdf(d_ICUstay, t - 1) for t = 1:100]#ICU to leave ICU
p_HR = [cdf(d_duration_in_hosp, t) - cdf(d_duration_in_hosp, t - 1) for t = 1:100]#Hospital to discharge
p_ICUR = KenyaCoVaccines.simple_conv(p_ICU, p_HR)#ICU to discharge assuming its the sum of the two
p_SH = [0.2 for t = 1:5] #Moghadas et al uniform 1-5 day estimate
p_R = [cdf(d_recovery, t) - cdf(d_recovery, t - 1) for t = 1:1000]
#Upper tail functions
Q_HR = vcat([1.0], [1 - cumsum(p_HR)[t] for t = 1:100])
Q_ICUH = vcat([1.0], [1 - cumsum(p_ICU)[t] for t = 1:100])
Q_ICUR = vcat([1.0], [1 - cumsum(p_ICUR)[t] for t = 1:100])
Q_R = vcat([1.0], [1 - cumsum(p_R)[t] for t = 1:1000])
F_R = 1 .- Q_R

##Functions for creating fit
"""
    function pred_ICU_incidence(H_ICU, η_hosp, unscaled_pred_inf, p_IH)

Predict daily incidence at ICU
"""
function pred_ICU_incidence(H_ICU, η_hosp, unscaled_pred_inf, p_IH)
    T = eltype(H_ICU)
    _p_IH = normalize([exp(η_hosp * t) * p_IH[t] for t = 1:50], 1)
    pred_ICU_inc = zeros(T, size(unscaled_pred_inf))
    for t = 1:length(pred_ICU_inc)
        for s = 1:min(t - 1, 50)
            pred_ICU_inc[t] += H_ICU * unscaled_pred_inf[t-s] * _p_IH[s]
        end
    end
    return pred_ICU_inc
end

pred_ICU_inc = pred_ICU_incidence(3.0, 0.0, unscaled_hosp_pred_inf, p_IH)

plot(pred_ICU_inc)
scatter!(kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid])

"""
    function pred_ICU_occupancy(pred_ICU_inc, Q_ICUH)

Predict the number of people in ICU on any day.        
"""
function pred_ICU_occupancy(pred_ICU_inc, Q_ICUH)
    pred_ICU_occupancy = zeros(eltype(pred_ICU_inc), length(pred_ICU_inc))
    for t = 1:length(pred_ICU_occupancy)
        for s = 1:min(t - 1, 100)
            pred_ICU_occupancy[t] += pred_ICU_inc[t-s] * Q_ICUH[s]
        end
    end
    return pred_ICU_occupancy
end
pred_ICU_inc = pred_ICU_incidence(1.5, 0.5, unscaled_hosp_pred_inf, p_IH)
pred_ICU_occ = pred_ICU_occupancy(pred_ICU_incidence(1.5, 0.5, unscaled_hosp_pred_inf, p_IH), Q_ICUH)


plot(pred_ICU_occ)
scatter!(kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid])

function regularization_ICU(H_ICU, η_hosp, α)
    reg_loss = -logpdf(Normal(log(1.5), 0.2), H_ICU)
    reg_loss += -logpdf(Normal(0.0, 0.2), η_hosp)
    reg_loss += -logpdf(Normal(log(0.25), 0.5), α)

    return reg_loss
end

function logpdf_mean_alpha_param_NB(data, mean, alpha)
    σ² = clamp(mean + (alpha * mean^2), 0.00001, 1e6)
    p_negbin = clamp(1 - (alpha * mean^2 / σ²), 1e-10, 0.9999999999)
    r_negbin = clamp(1 / alpha, 0, 1000000)
    return logpdf(NegativeBinomial(r_negbin, p_negbin), data)
end

function loss_ICU_occupancy(H_ICU, η_hosp, α, unscaled_pred_inf, p_IH, Q_ICUH, ts, ICU_occ_data)
    T = eltype(H_ICU)
    loss = T(0.0)
    loss += regularization_ICU(H_ICU, η_hosp, α)
    pred_ICU_occ = pred_ICU_occupancy(pred_ICU_incidence(H_ICU, η_hosp, unscaled_pred_inf, p_IH), Q_ICUH)
    for (k, t) in enumerate(ts)
        if t >= 90 && t <= 300
            loss -= logpdf_mean_alpha_param_NB(ICU_occ_data[k], pred_ICU_occ[t], α)
        end
    end
    return loss, pred_ICU_occ
end

@time loss, pred = loss_ICU_occupancy(1.5, 0.3, 0.25, unscaled_hosp_pred_inf, p_IH, Q_ICUH, kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid])

##Use GalacticOptim

ICU_fixed_data = (unscaled_pred_inf = unscaled_hosp_pred_inf, p_IH = p_IH, Q_ICUH = Q_ICUH, ts = kenya_MoH.rel_time[idxs_ICU_valid], ICU_occ_data = kenya_MoH.intensive_care_unit[idxs_ICU_valid])

function ICU_loss(x, p)
    H_ICU = exp(x[1])
    η_hosp = x[2]
    α = exp(x[3])
    @unpack unscaled_pred_inf, p_IH, Q_ICUH, ts, ICU_occ_data = p
    loss, pred = loss_ICU_occupancy(H_ICU, η_hosp, α, unscaled_pred_inf, p_IH, Q_ICUH, ts, ICU_occ_data)
    return loss, pred
end

_loss, _pred = ICU_loss([log(1.5), 0.3, log(0.25)], ICU_fixed_data)

recorded_losses = []
function plot_ICU_occ(x, loss, pred)
    println("loss = $(loss), alpha = $(exp(x[3]))")
    push!(recorded_losses, loss)
    plt = scatter(kenya_MoH.rel_time[idxs_ICU_valid], kenya_MoH.intensive_care_unit[idxs_ICU_valid], lab = "ICU occ data")
    plot!(plt, pred, lw = 3, lab = "pred")
    display(plt)
    return false
end

##Fit
x0 = [log(1.5), 0.3, log(0.25)]

f = OptimizationFunction(ICU_loss, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(f, x0, ICU_fixed_data)
sol = solve(prob, Flux.ADAM(0.001); cb = plot_ICU_occ, maxiters = 3000)

##
x = sol.u
H_ICU_to_death_fit = exp(x[1])
η_ICU = x[2]
alpha_ICU = exp(x[3])
# @save("data/fitted_rel_risk_ICU.jld2", H_ICU_to_death_fit, η_ICU, alpha_ICU)

## HOSPITALISATION
best_loss_ICU, best_pred_ICU = ICU_loss(sol.u, ICU_fixed_data)
best_pred_ICU_inc = pred_ICU_incidence(H_ICU_to_death_fit, η_ICU, ICU_fixed_data.unscaled_pred_inf, ICU_fixed_data.p_IH)
best_pred_H_inc_from_ICU = KenyaCoVaccines.simple_conv(best_pred_ICU_inc, p_ICU)

function regularization_hosp(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, α)
    reg_loss = -logpdf(Normal(log(15.0), 5.0), H_hosp)
    reg_loss = -logpdf(Normal(log(15.0), 5.0), H_var_ab)
    reg_loss = -logpdf(Normal(log(17.0), 5.0), H_var_delta)

    reg_loss += -logpdf(Normal(0.0, 0.2), η_hosp)
    reg_loss += -logpdf(Normal(0.0, 0.2), η_hosp_delta)

    reg_loss += -logpdf(Normal(log(0.25), 0.5), α)

    return reg_loss
end

function pred_hosp_incidence(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, unscaled_pred_inf_wt, unscaled_pred_inf_ab, unscaled_pred_inf_delta, p_IH)
    T = eltype(H_hosp)
    _p_IH = normalize([exp(η_hosp * t) * p_IH[t] for t = 1:50], 1)
    _p_IH_delta = normalize([exp(η_hosp_delta * t) * p_IH[t] for t = 1:50], 1)

    pred_hosp_inc = zeros(T, size(unscaled_pred_inf_wt))
    for t = 1:length(pred_ICU_inc)
        for s = 1:min(t - 1, 50)
            pred_hosp_inc[t] += (H_hosp * unscaled_pred_inf_wt[t-s] + H_var_ab * unscaled_pred_inf_ab[t-s]) * _p_IH[s]
            pred_hosp_inc[t] += (H_var_delta * unscaled_pred_inf_delta[t-s]) * _p_IH_delta[s]

        end
    end
    return pred_hosp_inc
end

pred_hosp_inc = pred_hosp_incidence(10.0, 10.0, 15.0, 0.3, 0.6, unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, ICU_fixed_data.p_IH)
plot(pred_hosp_inc)

function pred_hosp_occupancy(pred_hosp_inc, best_pred_H_inc_from_ICU, Q_HR)
    pred_hosp_occupancy = zeros(eltype(pred_hosp_inc), length(pred_hosp_inc))
    for t = 1:length(pred_hosp_occupancy)
        for s = 1:min(t - 1, 100)
            # pred_hosp_occupancy[t] += (pred_hosp_inc[t-s] + best_pred_H_inc_from_ICU[t-s]) * Q_ICUH[s]
            pred_hosp_occupancy[t] += (pred_hosp_inc[t-s] + best_pred_H_inc_from_ICU[t-s]) * Q_HR[s]
        end
    end
    return pred_hosp_occupancy
end

pred_hosp_occ = pred_hosp_occupancy(pred_hosp_inc, best_pred_H_inc_from_ICU, Q_HR)
plot(pred_hosp_occ)
scatter!(kenya_MoH.rel_time[idxs_hosp_valid], kenya_MoH.in_hospital[idxs_hosp_valid])


function loss_hosp_occupancy(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, α, unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, p_IH, Q_HR, ts, hosp_occ_data)
    T = eltype(H_hosp)
    loss = T(0.0)
    loss += regularization_hosp(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, α)
    pred_hosp_occ = pred_hosp_occupancy(pred_hosp_incidence(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, p_IH), best_pred_H_inc_from_ICU, Q_HR)
    for (k, t) in enumerate(ts)
        if t >= 90 && t <= 300
            loss -= logpdf_mean_alpha_param_NB(hosp_occ_data[k], pred_hosp_occ[t], α)
        end
    end
    return loss, pred_hosp_occ
end

loss, pred = loss_hosp_occupancy(10.0, 10.0, 15.0, 0.4, 0.25, 0.4, unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, ICU_fixed_data.p_IH, Q_HR, kenya_MoH.rel_time[idxs_hosp_valid], kenya_MoH.in_hospital[idxs_hosp_valid])
plot!(pred)

hosp_fixed_data = (unscaled_hosp_pred_inf_wt = unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab = unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta = unscaled_hosp_pred_inf_delta,
    p_IH = p_IH, Q_HR = Q_HR, ts = kenya_MoH.rel_time[idxs_hosp_valid], hosp_occ_data = kenya_MoH.in_hospital[idxs_hosp_valid])

function hosp_loss(x, p)
    H_hosp = exp(x[1])
    H_var_ab = exp(x[2])
    H_var_delta = exp(x[3])
    η_hosp = x[4]
    η_hosp_delta = x[5]
    α = exp(x[6])
    @unpack unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, p_IH, Q_HR, ts, hosp_occ_data = p
    loss, pred = loss_hosp_occupancy(H_hosp, H_var_ab, H_var_delta, η_hosp, η_hosp_delta, α, unscaled_hosp_pred_inf_wt, unscaled_hosp_pred_inf_ab, unscaled_hosp_pred_inf_delta, p_IH, Q_HR, ts, hosp_occ_data)
    return loss, pred
end
_loss, _pred = hosp_loss([log(10.0), log(10.0), log(15.0), 0.4, 0.4, log(0.25)], hosp_fixed_data)
plot(_pred)
scatter!(kenya_MoH.rel_time[idxs_hosp_valid], kenya_MoH.in_hospital[idxs_hosp_valid])

recorded_losses_hosp = []
function plot_hosp_occ(x, loss, pred)
    println("loss = $(loss), alpha = $(exp(x[6]))")
    push!(recorded_losses_hosp, loss)
    plt = scatter(kenya_MoH.rel_time[idxs_hosp_valid], kenya_MoH.in_hospital[idxs_hosp_valid], lab = "Hosp occ data")
    plot!(plt, pred, lw = 3, lab = "pred")
    display(plt)
    return false
end


##
x0_h = [log(10.0), log(10.0), log(15.0), 0.4, 0.4, log(0.25)]

f_h = OptimizationFunction(hosp_loss, GalacticOptim.AutoForwardDiff())
prob_h = OptimizationProblem(f_h, x0_h, hosp_fixed_data)
sol_h = solve(prob_h, Flux.ADAM(0.001); cb = plot_hosp_occ, maxiters = 3000)
scatter(recorded_losses)

x = sol_h.u

H_hosp = exp(x[1])
H_var_ab = exp(x[2])
H_var_delta = exp(x[3])
η_hosp = x[4]
η_hosp_delta = x[5]
α = exp(x[6])
alpha_hosp = α
@save("data/fitted_rel_risk_hosp.jld2", H_hosp, H_var_ab, H_var_delta, η_hosp, alpha_hosp)


