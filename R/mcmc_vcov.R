

###########################################
# variance covariance matrix
mcmc_vcov <- function( mcmcobj , exclude = "deviance" ){	
	mcmcobj <- mcmcobj[ , ! ( base::colnames(mcmcobj) %in% exclude ) ]
	res <- stats::var(mcmcobj)
	base::colnames(mcmcobj) -> base::colnames(res) -> base::rownames(res)
	base::return(res)
		}