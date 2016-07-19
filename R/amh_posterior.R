amh_posterior <- function( pars , model , prior , data){
	#*** evaluate likelihood
	ll <- amh_loglike( model=model, data=data , pars=pars)
	#*** evaluate priors
	priorval <- amh_eval_priors( pars= pars , prior=prior )
	#*** compute posterior
	posteriorval <- ll + priorval
	#--- output
	res <- base::list( ll = ll , priorval = priorval ,
				posteriorval = posteriorval	, pars = pars )
	base::return(res)
}