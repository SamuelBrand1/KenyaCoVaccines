## Script for plots of fixed parameters

using Dates, Statistics
using StatsPlots, JLD2, Plots.PlotMeasures, CSV, DataFrames
using Parameters, Distributions, LinearAlgebra
using RecursiveArrayTools, MCMCChains, Interpolations, Suppressor
import KenyaCoVaccines


## Age dependent rates fixed in model
# 0-19, 20-49, 50-59, 60-69, 70-79 and 80+ year olds

model = load("modelfits/Nairobi_model.jld2")["model"]
age_strs = ["0-19", "20-49", "50-59", "60-69", "70-79", "80+"]
rel_sus_plt = bar(model.p_sus,
    lab = "",
    ylabel = "Rel. susceptibility",
    xlabel = "Age group",
    xticks = (1:6, age_strs),
    title = "Relative susceptibility in transmission model",
    tickfont = 11, guidefont = 13, titlefont = 16)

symptom_plt = bar(model.p_symp,
    lab = "",
    ylabel = "Prob. symptoms",
    xlabel = "Age group",
    xticks = (1:6, age_strs),
    title = "Probability of symptomatic disease in transmission model",
    tickfont = 11, guidefont = 13, titlefont = 12)

# savefig(rel_sus_plt, "generation_of_data_for_econ_analysis/SI_plots/rel_sus_plt.png")
# savefig(symptom_plt, "generation_of_data_for_econ_analysis/SI_plots/symp_plt.png")

## Age dependent contact rates fixed in model
max_mixing_rate = maximum([maximum(model.M_county_ho), maximum(model.M_county_other), maximum(model.M_county_school), maximum(model.M_county_work)])

plt_M_home = heatmap(model.M_county_ho,
    clims = (0, max_mixing_rate),
    yticks = (1:6, age_strs),
    xticks = (1:6, age_strs),
    ylabel = "Age group (from)",
    xlabel = "Age group (to)",
    title = "Scaled contact rates from person to group (Home)")

plt_M_other = heatmap(model.M_county_other,
    clims = (0, max_mixing_rate),
    yticks = (1:6, age_strs),
    xticks = (1:6, age_strs),
    ylabel = "Age group (from)",
    xlabel = "Age group (to)",
    title = "Scaled contact rates from person to group (Other)")

plt_M_school = heatmap(model.M_county_school,
    clims = (0, max_mixing_rate),
    yticks = (1:6, age_strs),
    xticks = (1:6, age_strs),
    ylabel = "Age group (from)",
    xlabel = "Age group (to)",
    title = "Scaled contact rates from person to group (School)")

plt_M_work = heatmap(model.M_county_work,
    clims = (0, max_mixing_rate),
    yticks = (1:6, age_strs),
    xticks = (1:6, age_strs),
    ylabel = "Age group (from)",
    xlabel = "Age group (to)",
    title = "Scaled contact rates from person to group (Work)")

# savefig(plt_M_home, "generation_of_data_for_econ_analysis/SI_plots/home_mixing.png")
# savefig(plt_M_other, "generation_of_data_for_econ_analysis/SI_plots/other_mixing.png")
# savefig(plt_M_school, "generation_of_data_for_econ_analysis/SI_plots/school_mixing.png")
# savefig(plt_M_work, "generation_of_data_for_econ_analysis/SI_plots/work_mixing.png")


## Plots for the infection model

@load("data/fitted_mort_by_county_age_variant.jld2")
@load("data/fitted_rel_risk_ICU.jld2")
@load("data/fitted_rel_risk_hosp.jld2")



# countynames
counties_with_serodata = ["Nairobi", "Mombasa",
    "Nakuru", "Uasin Gishu", "Embu", "Kisumu", "Siaya", "Kisii", "Nyeri",
    "Kilifi", "Kwale"]

semiurbancounties_lowserodata = ["Kajiado", "Kiambu", "Machakos", "Isiolo", "Kirinyaga", "Murang'a",
    "Trans Nzoia", "Nyandarua", "Meru", "Bungoma", "Nyamira", "Kakamega",
    "Laikipia", "Taita Taveta", "Makueni", "Elgeyo-Marakwet",
    "Tharaka-Nithi", "Kericho", "Bomet", "Homa Bay", "Busia", "Migori",
    "Nandi", "Vihiga"]

ruralcounties_lowserodata = ["Garissa", "Kitui",
    "Lamu", "Marsabit", "Narok", "Tana River", "Turkana", "Wajir",
    "Baringo", "Mandera", "Samburu", "West Pokot"]

countynames = sort([counties_with_serodata; semiurbancounties_lowserodata; ruralcounties_lowserodata])

plt_rel_report = bar(mort_scales,
    orientation = :horizontal,
    lab = "",
    title = "Relative reporting rate for clinical outcomes",
    yticks = (1:47, countynames),
    xlabel = "Rel. reporting rate",
    xticks = 0:0.1:1, xlims = (0, 1.0),
    tickfont = 9, guidefont = 13, titlefont = 12, dpi = 250,
    size = (500, 900), left_margin = 20mm, right_margin = 20mm)

# savefig(plt_rel_report, "generation_of_data_for_econ_analysis/SI_plots/rel_reporting.png")

plt_mort_age = bar(mort_age,
    lab = "",
    xticks = (1:6, age_strs),
    ylabel = "Risk per infection",
    title = "Baseline fatality risk per infection",
    xlabel = "Age group",
    yticks = 0.0:0.005:0.04,
    tickfont = 9, guidefont = 13, titlefont = 14, dpi = 250)

# savefig(plt_mort_age, "generation_of_data_for_econ_analysis/SI_plots/mort_age.png")
var_str = ["wild-type", "Alpha/Beta variant", "Delta variant"]
plt_mort_var = bar([1.0; mort_var],
    lab = "",
    xticks = (1:3, var_str),
    ylabel = "Rel. risk due to variant",
    title = "Rel. fatality risk by variant",
    xlabel = "Variant",
    tickfont = 9, guidefont = 13, titlefont = 14, dpi = 250)

# savefig(plt_mort_var, "generation_of_data_for_econ_analysis/SI_plots/mort_var.png")

##
plt_mort_age_outcome = groupedbar(mort_age .* [1.0; H_ICU_to_death_fit; H_hosp]',
    lab = ["Deadly" "Critical" "Severe"], legend = :topleft,
    xticks = (1:6, age_strs),
    ylabel = "Risk per infection",
    title = "Baseline risk per infection",
    xlabel = "Age group",
    # yticks = 0.0:0.005:0.04,
    tickfont = 9, guidefont = 13, titlefont = 14, dpi = 250)
savefig(plt_mort_age_outcome, "generation_of_data_for_econ_analysis/SI_plots/mort_age_outcome.png")

plt_mort_var_outcome = groupedbar(hcat([1.0; mort_var], [1.0; mort_var]),
    lab = ["Deadly" "Critical" "Severe"], legend = :topleft,
    # xticks = (1:6, age_strs),
    ylabel = "Risk per infection",
    title = "Baseline risk per infection",
    xlabel = "Age group",
    # yticks = 0.0:0.005:0.04,
    tickfont = 9, guidefont = 13, titlefont = 14, dpi = 250)







