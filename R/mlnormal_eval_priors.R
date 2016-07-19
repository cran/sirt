
mlnormal_eval_priors <- function( pars , prior, sum_all = FALSE ){
	NP <- base::length(pars)
	priorval0 <- base::rep(NA,NP)
	eps <- 1E-100	
	for (pp in 1:NP){
		#*** evaluate
		pars_pp <- base::names(pars)[pp]
		prior_pp <- prior[[ pars_pp ]]
		prior_arg_pp <- prior_pp[[2]]
		prior_arg_pp[[1]] <- pars[pp]
		priorval_pp <- base::log( base::do.call( prior_pp[[1]], prior_arg_pp ) + eps )
		priorval0[pp] <- priorval_pp
	}
	if ( sum_all ){
		priorval0 <- base::sum(priorval0)
	}
	base::return(priorval0)
}