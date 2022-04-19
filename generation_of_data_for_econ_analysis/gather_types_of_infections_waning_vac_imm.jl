using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra, LogExpFunctions
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks, Sundials, Tullio
using RecursiveArrayTools, Interpolations, MCMCChains, Suppressor, Interpolations
using RCall
import KenyaCoVaccines

##

@load("data/p_ID.jld2")
@load("data/N_kenya.jld2")
model = @suppress_err load("modelfits/Nairobi_model.jld2")["model"]
age_symptoms = model.p_symp
veff1 = 1 - model.VE_severe_disease_risk[2]
veff2 = 1 - model.VE_severe_disease_risk[3]
veff_death1 = 0.2
veff_death2 = 0.1

include("generation_methods.jl");
include("inference_methods.jl");

@load("data/fitted_mort_by_county_age_variant.jld2")
@load("data/fitted_rel_risk_ICU.jld2")
@load("data/fitted_rel_risk_hosp.jld2")

## Get vaccine scenarios

t_0 = (Date(2021, 9, 27) - Date(2020, 2, 20)).value #Vaccines start Sept 2021
t_end = (Date(2021, 9, 27) + Month(18) - Date(2020, 2, 20)).value #Option where vaccines finish in 18 months
t_end_rapid = (Date(2021, 9, 27) + Month(6) - Date(2020, 2, 20)).value #Option where vaccines finish in 6 months
ts = [0; t_0; t_end; t_end + 100] .+ 14 # Time points for 18 month deployment lagged 14 days to account for delay between jab and maximum effectiveness
ts_rapid = [0; t_0; t_end_rapid; t_end_rapid + 100] .+ 14 # Time points for 6 month deployment

# Number of people DOUBLE jabbed by time points under three coverages
fully_vaccinated_30perc = [0; 0; 15.9e6 / 2; 15.9e6 / 2]
fully_vaccinated_50perc = [0; 0; 25e6 / 2; 25e6 / 2]
fully_vaccinated_70perc = [0; 0; 35e6 / 2; 35e6 / 2]

## Create vaccine callbacks

include("create_vac_callbacks.jl");

## Get basic models for each county

fitfiles = readdir("modelfits", join=true)
model_vect = [@suppress_err load(filename)["model"] for filename in fitfiles]

##Test model

model = model_vect[30]
perc_over_18 = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])
KenyaCoVaccines.change_prob_immune_escape_with_waning!(model, perc_over_18; startdate=Date(2020, 12, 1), enddate=Date(2021, 12, 1))
VE_acq_base = model.VE_acquisition
VE_inf_base = model.VE_infectiousness
VE_sev_base = model.VE_severe_disease_risk

##

θs = model |> flatten_chain
# function immune_escape!(integrator)
#     println(size(integrator.u))
#     R = @view integrator.u[:, 7, :]
#     #Reduce generation time by 30%
#     gen_time_scale = 0.7
#     integrator.p[1] = 1 * (1 / gen_time_scale) * integrator.p[1]# β₀
#     integrator.p[6] = (1 / gen_time_scale) * integrator.p[6] # α
#     integrator.p[8:10] = (1 / gen_time_scale) * integrator.p[8:10]# αP, γA, γM

#     #Reduce vaccine protection against acquiring and transmitting infection by 50%
#     integrator.p[19:24] = 0.5 * integrator.p[19:24] #ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3
#     #Reduce protection from reinfection by 50%
#     reinf_red = 0.5
#     integrator.p[28] = (1.0 * reinf_red + 0.16 * (1 - reinf_red)) / 0.16
#     #50% of R -> W at emergence of immune escape variant but set loss of complete immunity to 0
#     integrator.u[:, 7, :] .-= 0.5 .* R
#     integrator.u[:, 8, :] .+= 0.5 .* R
#     integrator.p[11] = 0.0 # ω

# end

fs = falses(6, 11, 5)
fs[:, 6, :] .= true
idxs_R_cvode = fs[:]
fs = falses(6, 11, 5)
fs[:, 7, :] .= true
idxs_W1_cvode = fs[:]
fs = falses(6, 11, 5)
fs[:, 8, :] .= true
idxs_W2_cvode = fs[:]

function immune_escape_CVODE!(integrator)

    #Reduce generation time by 30%
    gen_time_scale = 0.7
    integrator.p[1] = 1 * (1 / gen_time_scale) * integrator.p[1]# β₀
    integrator.p[6] = (1 / gen_time_scale) * integrator.p[6] # α
    integrator.p[8:10] = (1 / gen_time_scale) * integrator.p[8:10]# αP, γA, γM

    #Reduce vaccine protection against acquiring and transmitting infection by 50%
    integrator.p[19:24] = 0.5 * integrator.p[19:24] #ve_ac1, ve_ac2, ve_ac3, ve_inf1, ve_inf2, ve_inf3
    #Reduce protection from reinfection in W1 by 0%
    reinf_red = 0.0
    integrator.p[28] = (1.0 * reinf_red + 0.16 * (1 - reinf_red)) / 0.16
    #50% of R -> W2 and W1 -> W2 at emergence of immune escape variant but set loss of complete immunity to 0
    integrator.u[idxs_W2_cvode] .+= 0.5 .* integrator.u[idxs_R_cvode]
    integrator.u[idxs_W2_cvode] .+= 0.5 .* integrator.u[idxs_W1_cvode]
    integrator.u[idxs_R_cvode] .-= 0.5 .* integrator.u[idxs_R_cvode]
    integrator.u[idxs_W1_cvode] .-= 0.5 .* integrator.u[idxs_W1_cvode]

    # integrator.p[11] = 0.0 # ω

end

cb_imm_esc_cvode = PresetTimeCallback([349.0], immune_escape_CVODE!, save_positions=(false, false))

# cb_imm_esc = CallbackSet()
solver = CVODE_BDF(linear_solver=:GMRES)
# solver = Rosenbrock23(autodiff = false)
##Generate random vaccine effectivenesses
VE_acqs = [[0.0, rand(Uniform(0.5, 0.6)), rand(Uniform(0.6, 0.85))] for i = 1:2000]
VE_infs = [[0.0, 0.0, rand(Uniform(0.0, 0.69))] for i = 1:2000]
VE_sevs = [[0.0, rand(Uniform(0.8, 0.9)), rand(Uniform(0.95, 0.99))] for i = 1:2000]
VE_deaths = [[0.0, rand(Uniform(0.9, 0.95)), rand(Uniform(0.95, 0.99))] for i = 1:2000]

@save("vac_effs.jld2", VE_acqs, VE_infs, VE_sevs, VE_deaths)
#Matrix versions
VE_waning = 0.0
VE_sevs_mat = 1 .- [vecvec_to_mat(VE_sevs) vecvec_to_mat(VE_sevs)[:, end] VE_waning .* vecvec_to_mat(VE_sevs)[:, end]]
VE_deaths_mat = 1 .- [vecvec_to_mat(VE_deaths) vecvec_to_mat(VE_deaths)[:, end] VE_waning .* vecvec_to_mat(VE_deaths)[:, end]]

##
rep = 1062
@time sol_vw = solve_model_waning_vac(θs[rep], model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[7], solver, VE_acqs[rep], VE_infs[rep], VE_sevs[rep])
# ι_all = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol_vw, 10), 0.0, Inf) #All infections
# ι_dis_weighted = clamp!(KenyaCoVaccines.get_incidence_time_array_by_age_dose(sol_vw, 11), 0.0, Inf) #Weighted by infection risk
E = [sum(u[:, 2, :]) for u in sol_vw.u]
A = [sum(u[:, 3, :]) for u in sol_vw.u]
R = [sum(u[:, 6, :]) for u in sol_vw.u]
W = [sum(u[:, 7, :]) for u in sol_vw.u]
W2 = [sum(u[:, 8, :]) for u in sol_vw.u]
F = [sum(u[:, 9, :]) for u in sol_vw.u]
C = [sum(u[:, 10, :]) for u in sol_vw.u]
C2 = [sum(u[:, 11, :]) for u in sol_vw.u]

V1 = [sum(u[:, 1:8, 2]) for u in sol_vw.u]
V2 = [sum(u[:, 1:8, 3]) for u in sol_vw.u]
WV1 = [sum(u[:, 1:8, 4]) for u in sol_vw.u]
WV2 = [sum(u[:, 1:8, 5]) for u in sol_vw.u]
##
plot(R ./ sum(model.N))
plot(W)
plot((R .+ W .+ W2) / sum(model.N), lab="")
plot(E)
plot(A, lab="")
plot(W2 ./ sum(model.N))

plot(diff(C))
plot(V1, lab="")
plot(V2, lab="")
plot(WV1, lab="")
plot(WV2, lab="")
plot((V1 .+ V2 .+ WV1 .+ WV2) ./ sum(model.N[2:6]))
plot((V2 .+ V1) ./ sum(model.N[2:6]))

##
# findall(sol_ie.u[end] .< 0)
##
# plot(diff(C), lab="All infections - no vaccine")
plot(C, lab="All infections - 50% rapid vaccine")
# plot!(diff(C2), lab="Dis risk weighted infections")
# plot!(diff(F), lab="First time infections")
plot!(C, lab="All infections - 70% rapid vaccine")
plot!(xlims=(30, 800),
    ylabel="Daily infections",
    xlabel="Days after 1/12/20",
    title="Nairobi infections with waning immunity",
    legend=:bottomright)
plot!(sum(ι_all, dims=[2, 3])[:])

plot!(sum(ι_dis_weighted, dims=[2, 3])[:])
plot!(sum(ι_all2, dims=[2, 3])[:])

vline!([349], lab="Variant introduction")
vline!([427], lab="")
vline!([230], lab="First jabs start")
##
@time infs = project_daily_infections_by_age_and_vaccine_waning_vac(θs[2], model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[7], solver)
ι_all2 = sum(infs[1], dims=[2, 3])[:]
ι_dis_weighted2 = sum(infs[2], dims=[2, 3])[:]
plot!(ι_all2)
plot!(ι_dis_weighted2)

##
θs = model |> flatten_chain
infections_by_age_dose = map((θ, VE_acq, VE_inf, VE_sev) -> project_daily_infections_by_age_and_vaccine_waning_vac(θ, model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[7], solver, VE_acq, VE_inf, VE_sev),
    θs, VE_acqs, VE_infs, VE_sevs)
# infections_by_age_dose = θs .|> θ -> project_daily_infections_by_age_and_vaccine_waning_vac(θ, model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[7], solver)

##

infs_gathered = gather_infection_type_for_econ_analysis_waning_vac(model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[1], solver, VE_acqs, VE_infs, VE_sevs, VE_deaths;
    mort_scale=mort_scales[30], mort_age=mort_age, mort_var=mort_var,
    H_ICU_to_death_fit=H_ICU_to_death_fit, H_hosp=H_hosp, H_var_ab=H_var_ab, H_var_delta=H_var_delta)

##
A_pred = get_credible_intervals(sum(infs_gathered.asymptomatic_infs, dims=2)[:, 1, :])
M_pred = get_credible_intervals(sum(infs_gathered.mild_infections, dims=2)[:, 1, :])
sev_pred = get_credible_intervals(sum(infs_gathered.severe_infs, dims=2)[:, 1, :])
crit_pred = get_credible_intervals(sum(infs_gathered.critical_infs, dims=2)[:, 1, :])
deaths_pred = get_credible_intervals(sum(infs_gathered.deadly_infs, dims=2)[:, 1, :])

plot!(A_pred.pred, ribbon=(A_pred.lb, A_pred.ub))
plot(deaths_pred.pred, ribbon=(deaths_pred.lb, deaths_pred.ub))

##
scenarios_names = ["no_vac", "30_perc", "50_perc", "70_perc", "30_perc_rapid", "50_perc_rapid", "70_perc_rapid"]


for scen_num = 1:7
    println("Beginning waning vac. eff. scenario $(scen_num)")
    for k = 1:length(model_vect)
        perc_over_18 = sum(N_kenya[5:end, model.areaname]) / sum(N_kenya[5:end, :])
        KenyaCoVaccines.change_prob_immune_escape_with_waning!(model_vect[k], perc_over_18; startdate=Date(2020, 12, 1), enddate=Date(2023, 12, 1))
    end

    infs_types = mapreduce((model, mort_scale) -> gather_infection_type_for_econ_analysis_waning_vac(model, Date(2023, 6, 1), cb_imm_esc_cvode, vac_cbs[scen_num], solver, VE_acqs, VE_infs, VE_sevs, VE_deaths;
            mort_scale=mort_scale, mort_age=mort_age, mort_var=mort_var,
            H_ICU_to_death_fit=H_ICU_to_death_fit, H_hosp=H_hosp, H_var_ab=H_var_ab, H_var_delta=H_var_delta),
        add_infs, model_vect, mort_scales)

    kenya_inc_A = infs_types.asymptomatic_infs
    kenya_inc_M = infs_types.mild_infections
    kenya_inc_sev = infs_types.severe_infs
    kenya_inc_crit = infs_types.critical_infs
    kenya_inc_deaths = infs_types.deadly_infs

    savefilename_inc_A = "generation_of_data_for_econ_analysis/waning_imm_scenario_$(scen_num)/kenya_inc_A_" * scenarios_names[scen_num] * "_w_vac.rda"
    savefilename_inc_M = "generation_of_data_for_econ_analysis/waning_imm_scenario_$(scen_num)/kenya_inc_M_" * scenarios_names[scen_num] * "_w_vac.rda"
    savefilename_inc_sev = "generation_of_data_for_econ_analysis/waning_imm_scenario_$(scen_num)/kenya_inc_sev_" * scenarios_names[scen_num] * "_w_vac.rda"
    savefilename_inc_crit = "generation_of_data_for_econ_analysis/waning_imm_scenario_$(scen_num)/kenya_inc_crit_" * scenarios_names[scen_num] * "_w_vac.rda"
    savefilename_inc_deaths = "generation_of_data_for_econ_analysis/waning_imm_scenario_$(scen_num)/kenya_inc_deaths_" * scenarios_names[scen_num] * "_w_vac.rda"


    @rput kenya_inc_A kenya_inc_M kenya_inc_sev kenya_inc_crit kenya_inc_deaths
    @rput savefilename_inc_A savefilename_inc_M savefilename_inc_sev savefilename_inc_crit savefilename_inc_deaths

    R"""
        save(kenya_inc_A,file = savefilename_inc_A)
        save(kenya_inc_M,file = savefilename_inc_M)
        save(kenya_inc_sev,file = savefilename_inc_sev)
        save(kenya_inc_crit,file = savefilename_inc_crit)
        save(kenya_inc_deaths,file = savefilename_inc_deaths)
    """
end

##Test randomized draws


R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_6/kenya_inc_A_50_perc_rapid_w_vac.rda")
"""
kenya_inc_A_sc_6 = @rget kenya_inc_A
kenya_inc_A_sc_6 = sum(kenya_inc_A_sc_6, dims=2)[:, 1, :]
kenya_inc_A = 0
R"""
kenya_inc_A = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/waning_imm_scenario_7/kenya_inc_A_70_perc_rapid_w_vac.rda")
"""
kenya_inc_A_sc_7 = @rget kenya_inc_A
kenya_inc_A_sc_7 = sum(kenya_inc_A_sc_7, dims=2)[:, 1, :]
kenya_inc_A = 0
R"""
kenya_inc_A = 0 #Free memory
"""

diff_A = (sum(kenya_inc_A_sc_6, dims=1).-sum(kenya_inc_A_sc_7, dims=1))[:]
boxplot(diff_A)

pred_A_sc_6 = get_credible_intervals(kenya_inc_A_sc_6)
pred_A_sc_7 = get_credible_intervals(kenya_inc_A_sc_7)
pred_A_diff_cum = get_credible_intervals(cumsum(kenya_inc_A_sc_6 .- kenya_inc_A_sc_7, dims=1))
pred_A_diff = get_credible_intervals(kenya_inc_A_sc_6 .- kenya_inc_A_sc_7)

plot(pred_A_sc_6.pred, ribbon=(pred_A_sc_6.lb, pred_A_sc_6.ub),
    lab="50% rapid", title="Asymptomatic infections", fillalpha=0.1)
plot!(pred_A_sc_7.pred, ribbon=(pred_A_sc_7.lb, pred_A_sc_7.ub),
    lab="70% rapid", fillalpha=0.1)

plot(pred_A_diff.pred, ribbon=(pred_A_diff.lb, pred_A_diff.ub),
    lab="", title="Daily Difference between 50% and 70%", fillalpha=0.2, lw=3)

## No waning

R"""
load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_6/kenya_inc_A_50_perc_rapid_imm_esc.rda")
"""
kenya_inc_A_sc_6_nw = @rget kenya_inc_A
kenya_inc_A_sc_6_nw = sum(kenya_inc_A_sc_6_nw, dims=2)[:, 1, :]
kenya_inc_A = 0
R"""
kenya_inc_A = 0 #Free memory
"""
R"""
load(file = "generation_of_data_for_econ_analysis/immune_escape_scenario_7/kenya_inc_A_70_perc_rapid_imm_esc.rda")
"""
kenya_inc_A_sc_7_nw = @rget kenya_inc_A
kenya_inc_A_sc_7_nw = sum(kenya_inc_A_sc_7_nw, dims=2)[:, 1, :]
kenya_inc_A = 0
R"""
kenya_inc_A = 0 #Free memory
"""

pred_A_sc_6_nw = get_credible_intervals(kenya_inc_A_sc_6_nw)
pred_A_sc_7_nw = get_credible_intervals(kenya_inc_A_sc_7_nw)
pred_A_diff_cum_nw = get_credible_intervals(cumsum(kenya_inc_A_sc_6_nw .- kenya_inc_A_sc_7_nw, dims=1))
pred_A_diff_nw = get_credible_intervals(kenya_inc_A_sc_6_nw .- kenya_inc_A_sc_7_nw)

plot(pred_A_sc_6_nw.pred, ribbon=(pred_A_sc_6_nw.lb, pred_A_sc_6_nw.ub),
    lab="50% rapid", title="Asymptomatic infections (no vac waning)", fillalpha=0.1)
plot!(pred_A_sc_7_nw.pred, ribbon=(pred_A_sc_7_nw.lb, pred_A_sc_7_nw.ub),
    lab="70% rapid", fillalpha=0.1)

plot(pred_A_diff_nw.pred, ribbon=(pred_A_diff_nw.lb, pred_A_diff_nw.ub),
    lab="", title="Daily Difference between 50% and 70% (no vac waning)", fillalpha=0.2, lw=3)

plot(pred_A_diff_cum_nw.pred, ribbon=(pred_A_diff_cum_nw.lb, pred_A_diff_cum_nw.ub),
    lab="", title="Cum. Difference between 50% and 70% (no vac waning)", fillalpha=0.2, lw=3)


