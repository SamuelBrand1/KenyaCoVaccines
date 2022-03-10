

"""
        function find_initial_MAP(model,trans;serowaningrate = 0.0)

This finds a MAP estimator for the `model` object with variable transformation `trans`. Step 1 is BlackBox optimisation
        using a differential evolution optimizer, which is then refined in step 2 by a local Nelder-Mead optimisation.
"""
function find_initial_MAP(model::CoVAreaModel,trans;serowaningrate = 0.0,VE_acquisition=zeros(3),VE_infectiousness=zeros(3),VE_severe_disease_risk=zeros(3))
	println("Searching for a good initial condition for $(model.areaname)")
	l_area(x) = model.log_likelihood(x,model,serowaningrate,VE_acquisition,VE_infectiousness,VE_severe_disease_risk) + model.log_priors(x)
	f(x) = -transform_logdensity(trans, l_area,x)
	searchrange = fill((-3.,3.),TransformVariables.dimension(trans))
	res = bboptimize(f; SearchRange = searchrange,PopulationSize=100,MaxSteps=10000)
	q₀ = best_candidate(res)
	res_nm = optimize(f,q₀,NelderMead(),Optim.Options(iterations = 1000,show_trace = true))
	return res_nm.minimizer
end


"""
        function PSGLD!(model,samples,trans,stepsize_params,q₀;serowaningrate = 0.0,beta = 0.99,λ =1e-5)

This draws a chain of parameter samples of length `samples` from the log-posterior defined by the `model` object using pre-conditioned Stochastic
Gradient Langevin Dynamics (PSGLD). The step size at time step `t` is `a*(b + t)^(-γ)`. The `a`, `b` and `γ` values must be given
in the `NamedTuple` object `stepsize_params`. The starting point is at `q₀`. The method is based on _C. Li, et al, Preconditioned stochastic gradient Langevin dynamics for deep neural networks. Proceedings of the AAAI Conference on Artificial Intelligence. 30, 1788–1794 (2016)_
As in the paper, the variable `beta` is the learning rate for the pre-conditioner and `λ` sets the maximum implied curvature of the parameter space. These are set to defaults from the paper.

The MCMC results are stored in the model object as per usual.
"""
function PSGLD!(model::CoVAreaModel,samples,trans,stepsize_params,q₀;serowaningrate = 0.0,VE_acquisition=zeros(3),VE_infectiousness=zeros(3),VE_severe_disease_risk=zeros(3),beta = 0.99,λ =1e-5)
	θ = copy(q₀)
	l_area(x) = model.log_likelihood(x,model,serowaningrate,VE_acquisition,VE_infectiousness,VE_severe_disease_risk) + model.log_priors(x)
	f(x) = transform_logdensity(trans, l_area,x)
	#ℓ = TransformedLogDensity(trans, l_area)#transformed log-likelihood
	#∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters

	log_post_den_vals = zeros(samples)
	param_draws = zeros(samples/10,length(θ)) # save every 10th draw
	precond = zeros(length(θ))
	G_diag = zeros(length(θ))

	for t in 1:samples
		################GRADIENT CALCULATIONS
		#lpd,∇L = LogDensityProblems.logdensity_and_gradient(∇ℓ,θ)
		lpd = f(θ)
		∇L =  ForwardDiff.gradient(f,θ)
		#Step size
		ϵ = stepsize_params.a*(stepsize_params.b + t)^(-stepsize_params.γ)

		###############PRECONDITIONING#####################
		#Learns the curvature of the space
		precond *= beta
		precond += (1-beta)*(∇L .*∇L)
		G_diag .= 1.0./(λ .+ sqrt.(precond))

		###############LOG-POSTERIOR ASCENT WITH RANDOM WALK COMPONENT#############
		θ .+= 0.5.*ϵ.*G_diag.*∇L  .+ sqrt.(ϵ.*G_diag).*randn(length(∇L))

		#################BOOKKEEPING############################
		if t % 10 == 0  # edited to save every 10the draw
			param_draws[t,:] = θ
			log_post_den_vals[t] = lpd
			println(θ[1])
			println("On step $(t) log-posterior density is $(lpd). Step size is $(ϵ*1e4)e-4. Max component of G is $(maximum(G_diag)).")
		end
	end
	transformed_samples = [TransformVariables.transform(trans,param_draws[t,:]) for t in 1:size(param_draws,1)]
	val = zeros(length(transformed_samples),length(transformed_samples[1]),1)
	for i = 1:size(val,1),j = 1:size(val,2)
		val[i,j,1] = transformed_samples[i][j]
	end
	chn = Chains(val,[String(k) for k in keys(transformed_samples[1])])
	model.MCMC_results = KenyaCoVaccines.MCMCResults(chn,
                                        [l_area(transformed_samples[i]) - model.log_priors(transformed_samples[i]) for i = 1:length(transformed_samples)],
                                        DynamicHMC.TreeStatisticsNUTS[])

	return param_draws,log_post_den_vals
end

"""
        function inferparameters!(model,samples,trans,psgld_stepsize::NamedTuple;serowaningrate = 0.0,beta = 0.99,λ = 1e-5)

                Infer the unknown transmission and observation parameters for the `model` county. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search has a pre-optimisation using the default adaptive differential evolution algorithm implemented by the `BlackBoxOptim` package,
MCMC draws are done via PSGLD. Inference output is stored in-place.
"""
function inferparameters!(model::CoVAreaModel,samples,trans,psgld_stepsize::NamedTuple;serowaningrate = 0.0,VE_acquisition=zeros(3),VE_infectiousness=zeros(3),VE_severe_disease_risk=zeros(3),beta = 0.99,λ = 1e-5)
	q₀ = find_initial_MAP(model,trans;serowaningrate =serowaningrate,VE_acquisition=VE_acquisition,VE_infectiousness=VE_infectiousness,VE_severe_disease_risk=VE_severe_disease_risk)
	PSGLD!(model,samples,trans,psgld_stepsize,q₀;serowaningrate =serowaningrate,VE_acquisition=VE_acquisition,VE_infectiousness=VE_infectiousness,VE_severe_disease_risk=VE_severe_disease_risk,beta = beta,λ = λ)
	return nothing
end


"""
        function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/180,num_chains = 1)
Infer the unknown transmission and observation parameters for the `k`th county in by `country_model`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search starts at `q₀` in the transformed domain defined by `trans`.
    whilst the initial kinetic energy search begins at `D` and the HMC step size is fixed to be `stepsize`. Inference output is stored in-place.
    Serowaning rate can be added, as well as additional chains for cross-comparison and validation.
"""
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,q₀;serowaningrate = 1/365,num_chains = 1)
        println("Starting MCMC parameter inference for $(model.areaname)")
        n = length(q₀)
        l_area(x) = model.log_likelihood(x,model,serowaningrate) + model.log_priors(x)
		ℓ = TransformedLogDensity(trans, l_area)#transformed log-likelihood
        ∇ℓ = LogDensityProblems.ADgradient(:ForwardDiff, ℓ)#transformed log-likelihood gradient wrt the parameters
		#θ = copy(q₀)
		#f(x) = - transform_logdensity(trans, l_area,x)
		#∇ℓ =  ForwardDiff.gradient(f,θ)

        results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                initialization = (q = q₀,κ=GaussianKineticEnergy(D),ϵ = stepsize),
                                warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric),
                                reporter = NoProgressReport())

        transformed_results = TransformVariables.transform.(trans,results.chain)
        val = zeros(length(transformed_results),length(transformed_results[1]),1)

        for i = 1:size(val,1),j = 1:size(val,2)
                val[i,j,1] = transformed_results[i][j]
        end

        for chain_reps = 2:num_chains
                println("Running additional chain $(chain_reps) for $(model.areaname)")
                results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇ℓ, samples,
                                        initialization = (q = q₀,κ=GaussianKineticEnergy(D),ϵ = stepsize),
                                        warmup_stages = fixed_stepsize_warmup_stages(M=Symmetric),
                                        reporter = NoProgressReport())
                transformed_results = TransformVariables.transform.(trans,results.chain)
                for i = 1:size(val,1),j = 1:size(val,2)
                        val[i,j,chain_reps] = transformed_results[i][j]
                end
        end

        chn = Chains(val,[String(k) for k in keys(transformed_results[1])])

        model.MCMC_results = MCMCResults(chn,
                                        [l_area(transformed_results[i]) - model.log_priors(transformed_results[i]) for i = 1:length(transformed_results)],
                                        results.tree_statistics)

        return nothing
end


"""
        function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal,searchrange;serowaningrate = 1/180,num_chains = 1)
Infer the unknown transmission and observation parameters for the `k`th county in by `country_model`. Inference is performed by drawing `samples` number of replicates from the posterior distribution with log-density
    given in the model with sero-waning rate fixed. The parameters and their transformation are encoded in `trans`. Initial parameter search has a pre-optimisation using the default adaptive differential evolution algorithm implemented by the `BlackBoxOptim` package,
    whilst the initial kinetic energy search begins at `D` and the HMC step size is fixed to be `stepsize`. Inference output is stored in-place.
    serowaning rate can be added, as well as additional chains for cross-comparison and validation.
"""
function inferparameters!(model,samples,trans,stepsize::Float64,D::Diagonal;serowaningrate = 1/365,num_chains = 1)
        println("Searching for a good initial condition for $(model.areaname)")
        l_area(x) = model.log_likelihood(x,model,serowaningrate) + model.log_priors(x)
        f(x) = -transform_logdensity(trans, l_area,x)
        searchrange = fill((-3.,3.),TransformVariables.dimension(trans))
        res = bboptimize(f; SearchRange = searchrange,PopulationSize=100,MaxSteps=30000,TraceMode = :silent)
        q₀ = best_candidate(res)
        inferparameters!(model,samples,trans,stepsize,D,q₀;serowaningrate=serowaningrate,num_chains=num_chains)
        return nothing
end
