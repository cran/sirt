
###########################################
# simulate from a Dirichlet distribution
dirichlet.simul <- function( alpha ){
    # alpha is a matrix input
    #-----
    N <- nrow(alpha)
    K <- ncol(alpha)
    ygamma <- 0*alpha
    for (ii in 1:K){   # ii <- 1
        ygamma[,ii] <- rgamma( n=N , shape=alpha[,ii] )
                }
    x <- ygamma / rowSums(ygamma)
    return(x)
            }
#################################################
# derivative of digamma function
digamma1 <- function(x,h=.001){
    ( digamma(x+h) - digamma(x-h) ) / (2*h)
                }
##################################################
# Maximum likelihood estimation of distribution parameters
dirichlet.mle <- function( x , eps=10^(-16),convcrit=.00001 , maxit=1000,
		progress=FALSE){
    N <- nrow(x)
    K <- ncol(x)
    # compute log pbar
	x <- x+eps
	x <- x / rowSums(x)
    log.pbar <- colMeans( log( x+eps ) )
    # compute inits
    alphaprob <- colMeans( x )
    p2 <- mean( x[,1]^2 )
    xsi <- ( alphaprob[1] - p2 ) / ( p2 - ( alphaprob[1] )^2 )
    alpha <- xsi * alphaprob 
    K1 <- matrix(1,K,K)
    conv <- 1
	iter <- 1
	#******************************
    # BEGIN iterations
    while( ( conv > convcrit ) & (iter < maxit) ){
        alpha0 <- alpha
        g <- N * digamma( sum(alpha ) ) - N * digamma(alpha) + N * log.pbar
        z <- N * digamma1( sum(alpha ))
        H <- diag( -N*digamma1( alpha ) ) + z
        alpha <- alpha0 - solve(H , g )
		alpha[ alpha < 0 ] <- 10^(-10)
        conv <- max( abs( alpha0 - alpha ) )
		if (progress){       print( paste( iter , sum(alpha) , conv) ) }
		iter <- iter+1
		flush.console()
                }
    alpha0 <- sum(alpha)
    xsi <- alpha / alpha0
    res <- list( "alpha"=alpha , "alpha0" = alpha0 , "xsi" = xsi )
    return(res)
        }
##############################################################
    				