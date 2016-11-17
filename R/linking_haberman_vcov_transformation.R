
linking_haberman_vcov_transformation <- function( H1 , aj_vcov )
{
	aj_vcov <- H1 %*% aj_vcov %*% base::t(H1)
	aj_se <- base::c( base::sqrt( base::diag( aj_vcov ) ) )
	res <- list( vcov = aj_vcov , se = aj_se )
	return(res)
}