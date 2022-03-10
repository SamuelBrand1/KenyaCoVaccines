"""
	function fit_psgld_step!(θ,∇L,precond,t, a, b, γ)

Perform one inplace draw of parameters using preconditioned stochastic gradient
Langevin dynamics.
"""
function fit_psgld_step!(θ,∇L,precond,step, a, b, γ)
    θ = θ₀
	##########DEFINE########################
	beta = 0.9;
	λ =1e-8;

	################GRADIENT CALCULATIONS
	 ϵ = a*(b + step)^-γ

	###############PRECONDITIONING#####################
	 if step == 1
	   precond[:] .= ∇L.*∇L
     else
	   precond .*= beta
	   precond .+= (1-beta).*(∇L .*∇L)
     end
	 #m = λ .+ sqrt.(precond/((1-(beta)^t)))
	 G_inv_diag = λ .+ sqrt.(precond)

	 ###############ASCENT###############################
	 for i in 1:length(∇L)
		 noise = sqrt(ϵ/G_inv_diag[i])*randn()
		 θ[i] += (0.5*ϵ*∇L[i]/G_inv_diag[i] +  noise)
	 end
	 return nothing
end
