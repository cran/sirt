

xxirt_parTheta_extract_freeParameters <- function( customTheta ){		
		ind <- customTheta$est
		p1 <- customTheta$par[ ind ]
		base::names(p1) <- base::names(customTheta$par)[ind]
		base::return(p1)
}