

###################################################################
amh_ic <- function( dev , N , pars , amh_summary , model , data , priorval ){	   	   
	# Information criteria
	ic <- base::list( "deviance" = dev , "n" = N )		
	ic$np <- base::length(pars)
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( base::log(ic$n) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( base::log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	#---- DIC
	p1 <- amh_summary[ 1 , "Mean"]	
	pars_MAP <- amh_summary$Mean[-1]
	base::names(pars_MAP) <- base::names(pars)
	ll <- base::do.call( model , base::list( pars = pars_MAP , data = data ) )
    ic$pD <-  p1 - ( -2*ll)	
	ic$DIC <- dev + 2*ic$pD	
	ic$loglik <- -dev/2
	ic$logprior <- priorval
	ic$logpost <- ic$loglik + ic$logprior
	base::return(ic)	
}	
###################################################################		