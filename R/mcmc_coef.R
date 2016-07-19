
###########################################
# coefficients from one MCMC chain
mcmc_coef <- function( mcmcobj , exclude="deviance" ){	
	mcmcobj <- mcmcobj[ , ! ( colnames(mcmcobj) %in% exclude ) ]
	res <- base::colMeans(mcmcobj)	
	colnames(mcmcobj) -> names(res)
	base::return(res)
		}