using JLD2, MCMCChains
import KenyaCoVaccines
using StatsPlots, Distributions, Optim

## Group fitted models by urban (Nairobi and Mombasa) and non-urban

non_urban_filenames = ["modelfits/Nakuru_model.jld2",
    "modelfits/Uasin Gishu_model.jld2",
    "modelfits/Embu_model.jld2",
    "modelfits/Kisumu_model.jld2",
    "modelfits/Nyeri_model.jld2",
    "modelfits/Siaya_model.jld2",
    "modelfits/Kisii_model.jld2",
    "modelfits/Kilifi_model.jld2",
    "modelfits/Kwale_model.jld2"]

## Combine posterior distribution draws of detection parameters across groups

non_urban_parameter_chains = []

for filename in non_urban_filenames
    fit_dict = load(filename)
    model = fit_dict[first(keys(fit_dict))]
    push!(non_urban_parameter_chains, model.MCMC_results.chain)
end

## Define a posterior density for the hyperparameters of an underlying Gamma distribution
#The assumption is that the non-cities have specific parameters drawn from a hierarchy (See Gelman et al)

# log-posterior density
function log_post_density_for_detection_parameters(k, θ, post_array::AbstractArray, param_sym::Symbol, prior_k, prior_θ)
    # LPD=0.0
    LPD = logpdf(prior_k, k) + logpdf(prior_θ, θ)
    for chn in post_array
        param_draws = get(chn, param_sym)
        LPD += log(mean([pdf(Gamma(k, θ), p) for p in param_draws[param_sym]]))
    end
    LPD
end


###
### Fit the population distribution of detection parameters 
### 
obj_p_test(x) = -log_post_density_for_detection_parameters(exp(x[1]), exp(x[2]), non_urban_parameter_chains, :p_test, Exponential(10), Exponential(0.5 / 10))
obj_χ(x) = -log_post_density_for_detection_parameters(exp(x[1]), exp(x[2]), non_urban_parameter_chains, :χ, Exponential(10), Exponential(1.5 / 10))
## 


fit1 = optimize(obj_p_test, log.([10.0, 0.1]), BFGS())
fit2 = optimize(obj_χ, log.([10.0, 0.1]), BFGS())

fitted_pop_paramsptest = exp.(fit1.minimizer)
fitted_pop_paramschi = exp.(fit2.minimizer)


## Save fitted priors --
## NB: Prior variance is doubled compared to fit
fitted_priors = (d_ptest = Gamma(fitted_pop_paramsptest[1] / 2, fitted_pop_paramsptest[2] * 2),
    d_χ = Gamma(fitted_pop_paramschi[1] / 2, fitted_pop_paramschi[2] * 2))

@save("data/non_urban_detection_fits.jld2", fitted_priors)


##Plots to examine the inferred priors

plot(fitted_priors.d_ptest, title = "Fitted prior for p_test (doubled variance)", lab = "")
plot(fitted_priors.d_χ, title = "Fitted prior for chi (doubled variance)", lab = "")

