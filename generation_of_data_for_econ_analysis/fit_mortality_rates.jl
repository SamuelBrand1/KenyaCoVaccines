using DataInterpolations, FileIO, CSV, DataFrames, Dates, Statistics
using StatsPlots, JLD2, Flux, ForwardDiff, GalacticOptim, LinearAlgebra
using Parameters, Distributions, OrdinaryDiffEq, DiffEqCallbacks, LogExpFunctions
using RecursiveArrayTools, Interpolations, MCMCChains, Suppressor, Interpolations
import KenyaCoVaccines

##Load methods and group data into correct age groups

include("generation_methods.jl");
include("inference_methods.jl");
@load("data/p_ID.jld2")
@load("data/cleaned_linelist20211109_deaths_c_age__date_of_lab_confirmation.jld2")

t0 = findfirst(deaths_data.dates .== Date(2020, 12, 2))
t1 = findfirst(deaths_data.dates .== Date(2021, 9, 27))

deathsdata = cat(sum(deaths_data.deaths[t0:t1, :, 1:4], dims = 3),
    sum(deaths_data.deaths[t0:t1, :, 5:10], dims = 3),
    sum(deaths_data.deaths[t0:t1, :, 11:12], dims = 3),
    sum(deaths_data.deaths[t0:t1, :, 13:14], dims = 3),
    sum(deaths_data.deaths[t0:t1, :, 15:16], dims = 3),
    sum(deaths_data.deaths[t0:t1, :, 17], dims = 3), dims = 3)

deathsdata_vect = [Matrix(deathsdata[:, k, :][:, :]') for k = 1:47]

prev_deaths_vect = [[sum(deaths_data.deaths[1:(t0+30), k, 1:4])
    sum(deaths_data.deaths[1:(t0+30), k, 5:10])
    sum(deaths_data.deaths[1:(t0+30), k, 11:12])
    sum(deaths_data.deaths[1:(t0+30), k, 13:14])
    sum(deaths_data.deaths[1:(t0+30), k, 15:16])
    sum(deaths_data.deaths[1:(t0+30), k, 17])] for k = 1:47]


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
function regularization(mort_scales, mort_age, mort_var, η)
    T = eltype(mort_scales)
    total_prior = T(0.0)
    for k = 1:29
        total_prior -= logpdf(Beta(10, 90), mort_scales[k])
    end
    for k = 31:47
        total_prior -= logpdf(Beta(10, 90), mort_scales[k])
    end
    total_prior -= logpdf(Normal(0, 0.01), log(mort_scales[30]))
    total_prior -= logpdf(Beta(0.0001 * 50, (1 - 0.0001) * 50), mort_age[1])
    total_prior -= logpdf(Beta(0.00002 * 50, (1 - 0.00002) * 50), mort_age[2])
    total_prior -= logpdf(Beta(0.00002 * 50, (1 - 0.00002) * 50), mort_age[3])
    total_prior -= logpdf(Beta(0.001 * 50, (1 - 0.001) * 50), mort_age[4])
    total_prior -= logpdf(Beta(0.01 * 20, (1 - 0.01) * 20), mort_age[5])
    total_prior -= logpdf(Beta(0.01 * 20, (1 - 0.01) * 20), mort_age[6])
    total_prior -= logpdf(Normal(log(1.2), 0.2), mort_var[1])
    total_prior -= logpdf(Normal(log(1.2), 0.2), mort_var[2])
    total_prior -= logpdf(Normal(0.0, 0.05), η)
    return total_prior
end



function poisson_loss_death_pred(mort_scales, mort_age, mort_var, mort_reinf, deathsdata_vect, prev_deaths_vect, prev_inf_predictions, predictions_vect, η, p_ID)
    T = eltype(mort_scales)
    total_loss = T(0.0)
    n = length(inf_predictions[1].infections_by_age_first_pred[1])
    overall_pred = zeros(T,n)

    for k = 1:47
        for a = 1:6
            total_loss -= logpdf(Poisson(mort_age[a] * mort_scales[k] * prev_inf_predictions[k][a]), prev_deaths_vect[k][a]) #Loss due to previous infections
        end
        _loss, _pred = poisson_loss_death_pred_county(deathsdata_vect[k], mort_scales[k], mort_age, mort_var, mort_reinf, predictions_vect[k], η, p_ID) #loss due to 
        total_loss += _loss
        overall_pred .+= _pred
    end
    total_loss += regularization(mort_scales, mort_age, mort_var, η)
    return total_loss, overall_pred
end

x = ones(47)
x[30] = 0.0
y = zeros(47)
y[30] = 1.0
scales_guess = 0.1 * x + y

deaths_loss_data = (; deathsdata_vect, prev_deaths_vect, prev_inf_predictions, inf_predictions, p_ID)

function deaths_loss(x, p)
    @unpack deathsdata_vect, prev_deaths_vect, prev_inf_predictions, inf_predictions, p_ID = deaths_loss_data
    mort_scales = [logistic.(x[1:29]); exp(x[30]); logistic.(x[31:47])]
    mort_age = logistic.(x[(47+1):47+6])
    mort_var = exp.(x[54:55])
    η = x[56]
    return poisson_loss_death_pred(mort_scales, mort_age, mort_var, 0.15, deathsdata_vect, prev_deaths_vect, prev_inf_predictions, inf_predictions, η, p_ID)
end
##

recorded_losses = Any[]


function plot_deaths(x, loss, pred)
    println("loss = $(loss)")
    push!(recorded_losses,loss)
    plt = scatter(sum(deathsdata, dims = [2, 3])[:], lab = "data")
    plot!(plt, pred, lw = 3, lab = "pred")
    display(plt)
    return false
end

plot_deaths(loss, pred)


x = ones(47)
x[30] = 0.0
y = zeros(47)
y[30] = 1.0
scales = 0.1 * x + y
x₀ = [logit.(scales[1:29])
    log.(scales[30])
    logit.(scales[31:47])
    logit.(fill(0.001, 6))
    log.([1.2, 1.2])
    0.0]

loss, pred = deaths_loss(x₀, deaths_loss_data)

##


f = OptimizationFunction(deaths_loss, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(f, x₀, deaths_loss_data)
sol = solve(prob, Flux.ADAM(0.01); cb = plot_deaths,maxiters = 200)



##

loss, pred = poisson_loss_death_pred(scales, fill(0.001, 6), [1.2, 1.2], 0.15, deathsdata_vect, prev_deaths_vect, prev_inf_predictions, inf_predictions, 0.0, p_ID)

plot!(pred)
scatter!(sum(deathsdata, dims = [2, 3])[:])
vline!([30])