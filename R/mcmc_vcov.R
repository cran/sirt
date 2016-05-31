

###########################################
# variance covariance matrix
mcmc_vcov <- function( mcmcobj , exclude = "deviance" ){	
	mcmcobj <- mcmcobj[ , ! ( colnames(mcmcobj) %in% exclude ) ]
	res <- var(mcmcobj)
	colnames(mcmcobj) -> colnames(res) -> rownames(res)
	return(res)
		}