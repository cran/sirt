
mlnormal_covmat_add_ridge <- function( covmat , eps = 1E-3)
{
	base::diag(covmat) <- base::diag(covmat)*( 1 + eps )
	base::return(covmat)
}