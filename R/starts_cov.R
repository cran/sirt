
starts_cov <- function(W , var_trait , var_ar , var_state , a )
{
	covmat <- matrix( var_trait , nrow=W , ncol=W )  
	for (ii in 1:W){
		covmat[ii,ii] <- covmat[ii,ii] + var_state
	}
	matr <- 0*covmat
	for (ii in 1:W){
		for (jj in 1:W){
			matr[ii,jj] <- a^( abs(ii-jj) ) * var_ar 
		}
	}
	covmat <- covmat + matr
	return(covmat)
}
