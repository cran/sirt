
vcov.xxirt <- function( object , ...){
		res <- xxirt_hessian( object )
		base::return( base::solve(-res) )
}