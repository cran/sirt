
mlnormal_covmat_add_ridge <- function( covmat , eps = 1E-3)
{
	diag(covmat) <- diag(covmat)*( 1 + eps )
	return(covmat)
}
