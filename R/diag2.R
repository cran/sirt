
diag2 <- function( vec)
{
	if ( base::length(vec) > 1){
		res <- base::diag(vec)
	} else {
		res <- matrix(vec, nrow=1,ncol=1)
	}
	base::return(res)
}