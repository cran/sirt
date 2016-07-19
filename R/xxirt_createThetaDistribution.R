
######################################################################
xxirt_createThetaDistribution <- function( par , est , P , prior=NULL,
		prior_par1 = NULL , prior_par2 = NULL	){
    res <- base::list()
    res$par <- par
    res$est <- est
    res$P <- P
	NP <- base::length(par)
#	if ( is.null(lower) ){
#		res$lower <- rep( - Inf , NP )
#			} else { 
#		res$lower <- lower
#				}
#	if ( is.null(upper) ){
#		res$upper <- rep( - Inf , NP )
#			} else { 
#		res$upper <- upper
#				}		
    res$prior <- prior
	res$prior_par1 <- prior_par1
	res$prior_par2 <- prior_par2
    base::class(res) <- "ThetaDistribution"
    base::return(res)
}
######################################################################				