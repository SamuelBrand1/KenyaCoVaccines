
function priors_onegroup_alpha_delta_variant_cities(θ)
	@unpack β₀,β_school,β_home,β_other,β_work,ϵ,χ,p_test,E₀,inc_R_αvar,time_scale_αvar,mid_point_αvar,inc_R_δvar,time_scale_δvar,mid_point_δvar,init_scale= θ #,S₁,S₂,S₃,S₄,S₅,S₆

    LP = 0.
    LP += logpdf(Beta(50,50),ϵ)
    LP += logpdf(Gamma(3,4.5/3),χ)
    LP += logpdf(Gamma(3,5/3),p_test) ##units relative to 1% detection chance

	LP += logpdf(Gamma(3,1000/3),E₀)

    LP += logpdf(Gamma(10,1.5/10),β₀)
    LP += logpdf(Gamma(10,1.5/10),β_other)
    LP += logpdf(Gamma(10,1.5/10),β_work)
    LP += logpdf(Gamma(10,1.5/10),β_school)
    LP += logpdf(Gamma(10,1.5/10),β_home)

	LP += logpdf(Gamma(15,0.4/15),inc_R_αvar)
	LP += logpdf(Gamma(15,0.15/15),time_scale_αvar)
	LP += logpdf(Gamma(5,0.5/5),mid_point_αvar) #In units of 30 days after Feb 1st
	LP += logpdf(Gamma(15,0.6/15),inc_R_δvar)
	LP += logpdf(Gamma(15,0.2/15),time_scale_δvar)
	LP += logpdf(Gamma(5,0.5/5),mid_point_δvar)
	LP += logpdf(Gamma(5,1.5/5),init_scale)

    return LP
end

function priors_onegroup_alpha_delta_variant_outside_cities(θ)
	@unpack β₀,β_school,β_home,β_other,β_work,ϵ,χ,p_test,E₀,inc_R_αvar,time_scale_αvar,mid_point_αvar,inc_R_δvar,time_scale_δvar,mid_point_δvar,init_scale = θ
	LP = 0.
    LP += logpdf(Beta(50,50),ϵ)
    LP += logpdf(Gamma(3,4.5/3),χ)
    LP += logpdf(Gamma(3,2/3),p_test) #units relative to 1% detection chance

	LP += logpdf(Gamma(3,1000/3),E₀)

    LP += logpdf(Gamma(10,1.5/10),β₀)
    LP += logpdf(Gamma(10,1.5/10),β_other)
    LP += logpdf(Gamma(10,1.5/10),β_work)
    LP += logpdf(Gamma(10,1.5/10),β_school)
    LP += logpdf(Gamma(10,1.5/10),β_home)

	LP += logpdf(Gamma(15,0.4/15),inc_R_αvar)
	LP += logpdf(Gamma(15,0.15/15),time_scale_αvar)
	LP += logpdf(Gamma(5,0.5/5),mid_point_αvar) #In units of 30 days after Feb 1st
	LP += logpdf(Gamma(15,0.6/15),inc_R_δvar)
	LP += logpdf(Gamma(15,0.2/15),time_scale_δvar)
	LP += logpdf(Gamma(5,0.5/5),mid_point_δvar)
	LP += logpdf(Gamma(5,1.5/5),init_scale)

    return LP
end


function priors_onegroup_alpha_delta_variant_noncities_fitted(θ)
    @unpack β₀, β_school, β_home, β_other, β_work, ϵ, χ, p_test, E₀, inc_R_αvar, time_scale_αvar, mid_point_αvar, inc_R_δvar, time_scale_δvar, mid_point_δvar, init_scale = θ


    LP = 0.0
    LP += logpdf(Beta(10, 10), ϵ)

    LP += logpdf(fitted_priors.d_χ, χ)  #
    LP += logpdf(fitted_priors.d_ptest, p_test) # units relative to 1% detection chance

    LP += logpdf(Gamma(3, 1000 / 3), E₀)

    LP += logpdf(Gamma(5, 1.5 / 5), β₀)
    LP += logpdf(Gamma(5, 1.5 / 5), β_other)
    LP += logpdf(Gamma(5, 1.5 / 5), β_work)
    LP += logpdf(Gamma(5, 1.5 / 5), β_school)
    LP += logpdf(Gamma(5, 1.5 / 5), β_home)

    LP += logpdf(Gamma(15, 0.4 / 15), inc_R_αvar)
    LP += logpdf(Gamma(15, 0.15 / 15), time_scale_αvar)
    LP += logpdf(Gamma(5, 0.5 / 5), mid_point_αvar) #In units of 30 days after Feb 1st
    LP += logpdf(Gamma(15, 0.6 / 15), inc_R_δvar)
    LP += logpdf(Gamma(15, 0.2 / 15), time_scale_δvar)
    LP += logpdf(Gamma(5, 0.5 / 5), mid_point_δvar)

	LP += logpdf(Gamma(3, 1.2 / 3), init_scale)

    return LP
end
