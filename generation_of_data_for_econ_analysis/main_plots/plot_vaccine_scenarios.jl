
using StatsPlots, Dates, CSV, DataFrames, JLD2, NamedArrays

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

#Daily rate of vaccination
scen_str = ["no vac", "30 perc", "50 perc", "70 perc", "30 perc rapid", "50 perc rapid", "70 perc rapid"]

scen_vac_rates = [0.0, (fully_vaccinated_30perc[3] - fully_vaccinated_30perc[2]) / (ts[3] - ts[2]),
    (fully_vaccinated_50perc[3] - fully_vaccinated_50perc[2]) / (ts[3] - ts[2]),
    (fully_vaccinated_70perc[3] - fully_vaccinated_70perc[2]) / (ts[3] - ts[2]),
    (fully_vaccinated_30perc[3] - fully_vaccinated_30perc[2]) / (ts_rapid[3] - ts_rapid[2]),
    (fully_vaccinated_50perc[3] - fully_vaccinated_50perc[2]) / (ts_rapid[3] - ts_rapid[2]),
    (fully_vaccinated_70perc[3] - fully_vaccinated_70perc[2]) / (ts_rapid[3] - ts_rapid[2])]

N_ovr_50 = sum(N_kenya[11:end, :])
N_18_49 = sum(N_kenya[5:10, :]) + 0.4 * sum(N_kenya[4, :])
max_level_ovr_50s = N_ovr_50 .* 0.8
max_level_18_49 = [0, fully_vaccinated_30perc[end], fully_vaccinated_50perc[end], fully_vaccinated_70perc[end],
    fully_vaccinated_30perc[end], fully_vaccinated_50perc[end], fully_vaccinated_70perc[end]] .- max_level_ovr_50s

max_level_18_49_perc = 100 .* max_level_18_49 ./ N_18_49

t_hit_max_lvl_ovr_50 = map(rate -> max_level_ovr_50s / rate, scen_vac_rates)
t_hit_max_lvl_18_49 = map((rate, max_lvl) -> max_lvl / rate, scen_vac_rates, max_level_18_49) .+ t_hit_max_lvl_ovr_50


## Get population sizes

@load("data/N_kenya.jld2")

## X ticks and dates

xtickdates = [Date(2020, 2, 1) + Month(k) for k = 19:29]
xtickpos = [(d - Date(2020, 2, 19)).value for d in xtickdates]
xticklab = [monthname(d)[1:3] * "-" * string(year(d) - 2000) for d in xtickdates]
t_sept_model = (Date(2021, 9, 1) - Date(2021, 1, 1)).value
t_sept_data = (Date(2021, 9, 1) - Date(2020, 2, 19)).value
t_end_plot = (Date(2022, 7, 1) - Date(2020, 2, 19)).value

t_start_model = (Date(2021, 1, 1) - Date(2020, 2, 19)).value
t_end_model = (Date(2022, 7, 1) - Date(2020, 2, 19)).value
t_end_model_fit = (Date(2021, 9, 1) - Date(2021, 1, 1)).value
t_end_fit = (Date(2021, 9, 1) - Date(2020, 2, 19)).value

xs_model_proj = (t_end_fit+1):t_end_model
xs_model_fitted_to_data = t_start_model:t_end_fit

n_model_proj = length(xs_model_proj)
n_model_fit = length(xs_model_fitted_to_data)

## Plot cumulative vaccines over 50s

plt_vac_ovr_50 = plot(xticks = (xtickpos, xticklab),
    size = (800, 400),
    left_margin = 5mm, right_margin = 2mm,
    xlims = (t_end_fit, t_end_plot + 5),
    ylims = (0, 100),
    legend = :bottomright,
    ylabel = "Fully vaccinated (%)",
    title = "Scenarios: proportion fully vaccinated (over 50 yo)",
    tickfont = 11, guidefont = 13, titlefont = 17)

plot!(plt_vac_ovr_50, [ts[1:2]; [10000, 100000]], [0, 0, 0, 0],
    lw = 2, color = :grey, lab = "No vaccination")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[2] + ts[2], t_hit_max_lvl_ovr_50[2] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, color = 1, lab = "30% coverage")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[3] + ts[2], t_hit_max_lvl_ovr_50[3] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, color = 2, lab = "50% coverage")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[4] + ts[2], t_hit_max_lvl_ovr_50[4] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, color = 3, lab = "70% coverage")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[5] + ts[2], t_hit_max_lvl_ovr_50[5] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, ls = :dash, color = 1, lab = "30% coverage, rapid rollout")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[6] + ts[2], t_hit_max_lvl_ovr_50[6] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, ls = :dash, color = 2, lab = "50% coverage, rapid rollout")
plot!(plt_vac_ovr_50, [ts[1:2]; [t_hit_max_lvl_ovr_50[7] + ts[2], t_hit_max_lvl_ovr_50[7] + ts[2] + 1000]], [0, 0, 80, 80],
    lw = 2, ls = :dash, color = 3, lab = "70% coverage, rapid rollout")

# savefig(plt_vac_ovr_50,"generation_of_data_for_econ_analysis/main_plots/vaccine_rate_ovr_50.png")

## Plot cumulative vaccines 18-49

plt_vac_18_49 = plot(xticks = (xtickpos, xticklab),
    size = (800, 400),
    left_margin = 5mm, right_margin = 2mm,
    xlims = (t_end_fit, t_end_plot + 5),
    ylims = (0, 100),
    legend = :topleft,
    ylabel = "Fully vaccinated (%)",
    title = "Scenarios: proportion fully vaccinated (18-49 yo)",
    tickfont = 11, guidefont = 13, titlefont = 17)

plot!(plt_vac_18_49, [ts[1:2]; [10000, 100000]], [0, 0, 0, 0],
    lw = 2, color = :grey, lab = "No vaccination")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[2] + ts[2], t_hit_max_lvl_18_49[2] + ts[2], t_hit_max_lvl_18_49[2] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[2], max_level_18_49_perc[2]],
    lw = 2, color = 1, lab = "30% coverage")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[3] + ts[2], t_hit_max_lvl_18_49[3] + ts[2], t_hit_max_lvl_18_49[3] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[3], max_level_18_49_perc[3]],
    lw = 2, color = 2, lab = "50% coverage")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[4] + ts[2], t_hit_max_lvl_18_49[4] + ts[2], t_hit_max_lvl_18_49[4] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[4], max_level_18_49_perc[4]],
    lw = 2, color = 3, lab = "70% coverage")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[5] + ts[2], t_hit_max_lvl_18_49[5] + ts[2], t_hit_max_lvl_18_49[5] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[5], max_level_18_49_perc[5]],
    lw = 2, ls = :dash, color = 1, lab = "30% coverage, rapid rollout")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[6] + ts[2], t_hit_max_lvl_18_49[6] + ts[2], t_hit_max_lvl_18_49[6] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[6], max_level_18_49_perc[6]],
    lw = 2, ls = :dash, color = 2, lab = "50% coverage, rapid rollout")
plot!(plt_vac_18_49, [ts[1:2]; [t_hit_max_lvl_ovr_50[7] + ts[2], t_hit_max_lvl_18_49[7] + ts[2], t_hit_max_lvl_18_49[7] + ts[2] + 1000]], [0, 0, 0, max_level_18_49_perc[7], max_level_18_49_perc[7]],
    lw = 2, ls = :dash, color = 3, lab = "70% coverage, rapid rollout")

# savefig(plt_vac_18_49, "generation_of_data_for_econ_analysis/main_plots/vaccine_rate_18_49.png")
