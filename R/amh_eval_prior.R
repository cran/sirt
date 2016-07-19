

amh_eval_priors <- function( pars , prior ){
	NP <- base::length(pars)
	priorval0 <- 0
	eps <- 1E-100
	for (pp in 1:NP){
		#*** evaluate
		pars_pp <- base::names(pars)[pp]
		prior_arg_pp <- prior[[pp]][[2]]
		prior_arg_pp[[1]] <- pars[pp]
		priorval_pp <- base::log( base::do.call( prior[[pp]][[1]], prior_arg_pp ) + eps )
		priorval0 <- priorval0 + priorval_pp
	}
	base::return(priorval0)
}