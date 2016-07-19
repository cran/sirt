
###############################################################
# pmle data processing
pmle_data_proc <- function( pars , pars_lower , pars_upper ){
	NP <- base::length(pars)
	if ( base::is.null( pars_lower) ){
		pars_lower <- base::rep( - Inf , NP )
	}	
	if ( base::is.null( pars_upper) ){
		pars_upper <- base::rep(  Inf , NP )
	}	
	res <- base::list( pars_lower = pars_lower ,
				pars_upper = pars_upper )
	base::return(res)
}