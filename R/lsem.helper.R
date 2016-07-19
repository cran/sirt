

sirt.wtdSD <- function( x , w ){
	res1 <- base::sum( x*w )
	res2 <- base::sum( x^2*w)
	res <- base::sqrt( res2 - res1^2 )
	base::return(res)
		}