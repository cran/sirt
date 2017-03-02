

###################################################################
# inverse gamma distribution for variance
rinvgamma2 <- function( n , n0 , var0 ){
    # INPUT:
    # N ... number of random draws
    # n0 ... sample size prior
    # var0 ... prior variance	
#	res <- 1/ stats::rgamma( N , n0 / 2 ,  n0 * var0 / 2 )
	res <- MCMCpack::rinvgamma( n , shape=n0 / 2 ,  scale=n0 * var0 / 2 )
    return(res) 
        }
#####################################################################
dinvgamma2 <- function( x , n0 , var0 ){
	res <- MCMCpack::dinvgamma( x , shape=n0 / 2 ,  scale=n0 * var0 / 2 )
    return(res) 
        }
