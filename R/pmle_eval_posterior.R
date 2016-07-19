

##################################################
# evaluate posterior
pmle_eval_posterior <- function( data , model , prior ,
				pars ,  eps = 1E-100 ){
	NP <- base::length(pars)
	#--- evaluate log-likelihood function
	ll <- base::do.call( model , base::list( pars = pars , data = data )  )	
	#--- evaluate prior distributions
	prior1 <- 0
	for (pp in 1:NP){
		prior_arg_pp <- prior[[pp]][[2]]
		prior_arg_pp[[1]] <- pars[pp]
		priorval <- base::log( base::do.call( prior[[pp]][[1]] , prior_arg_pp )  + eps )
		prior1 <- prior1 + priorval			
	}
	base::names(prior1) <- NULL
	#--- compute objective functions (posterior)
	post <- ll + prior1 
	#--- output
	res <- base::list( ll = ll , prior = prior1 , post = post )
	base::return(res)
}
##############################################################################			
