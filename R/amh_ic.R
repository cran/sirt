

###################################################################
amh_ic <- function( dev , N , pars , amh_summary , model , data ){	   	   
	# Information criteria
	ic <- list( "deviance" = dev , "n" = N )		
	ic$np <- length(pars)
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (consistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	
	#---- DIC
	p1 <- amh_summary[ 1 , "Mean"]	
	pars_MAP <- amh_summary$Mean[-1]
	names(pars_MAP) <- names(pars)
	ll <- do.call( model , list( pars = pars_MAP , data = data ) )
    ic$pD <-  p1 - ( -2*ll)	
	ic$DIC <- dev + 2*ic$pD	
	return(ic)	
		}	
###################################################################		