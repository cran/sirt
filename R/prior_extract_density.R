

###########################################################
# extract density of all variables
prior_extract_density <- function( prior ){
	NP <- base::length(prior)
	dens <- base::rep("" ,NP)
	base::names(dens) <- base::names(prior)
	for (pp in 1:NP){
		# pp <- 1
		v1 <- prior[[pp]][[1]]
		v2 <- prior[[pp]][[2]]
		NV <- base::length(v2)
		v_pp <- paste0( v1 , "(" )
		if (NV > 1){
			for (vv in 2:NV){
				# vv <- 2
				h1 <- base::names(v2)[vv]
				if ( ! base::is.null(h1) ){	
					if ( ! base::is.na(h1) ){
						if ( !(  h1 %in% c("NA","") ) ){
							v_pp <- base::paste0( v_pp , h1 , "=" )
						}
					}
				}
				v_pp <- base::paste0( v_pp , v2[vv])
				if (vv < NV){ 
					v_pp <- base::paste0( v_pp , "," ) 
				}
			}
		}
		dens[pp] <- base::paste0( v_pp , ")" )
	}
	base::return(dens)
}