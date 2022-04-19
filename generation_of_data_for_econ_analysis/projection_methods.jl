###Projection methods aim at post-processing the 
### saved scenario projections rather than
### directly operating on the saved model fits.

"""
    function create_NB_dist(μ,α)

Create a `NegativeBinomial` distribution with the correct `r,p` parameterisation from a `μ,α` parameterisation.        
"""
function create_NB_dist(μ, α)
    μ += 0.001*(μ == 0.0)

    σ² = μ + α * μ^2
    p_negbin = 1 - (α * μ^2 / σ²)
    r_negbin = 1 / α
    d = NegativeBinomial(r_negbin, p_negbin)
    return d
end


