
mcmc_Rhat <- function( mcmc_object , n_splits = 3 )
{
	n_samples <- base::nrow(mcmc_object)
	n_pars <- base::ncol(mcmc_object)
	n_within <- base::floor( n_samples / n_splits )
	rhat_vec <- base::rep(NA , n_pars)
	base::names(rhat_vec) <- base::colnames(mcmc_object)
	for (pp in 1:n_pars){
		# pp <- 1
		matr <- base::matrix( NA , nrow= n_within , ncol=n_splits)
		for (ss in 1:n_splits){
			matr[,ss] <- mcmc_object[ (ss-1)* n_within + 1:n_within , pp ]
		}
		rhat_vec[pp] <- Rhat1(matr)
	}
	base::return(rhat_vec)
}