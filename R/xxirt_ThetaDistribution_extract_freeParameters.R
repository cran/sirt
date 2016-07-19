

xxirt_ThetaDistribution_extract_freeParameters <- function( customTheta ){
		est <- customTheta$est
		if ( base::sum(est) == 0 ){
			par1 <- NULL 
		} else {
			par1 <- customTheta$par[ est ]
		}
		base::return(par1)
}