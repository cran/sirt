
###########################################################
pmle_process_prior <- function( prior , pars ){
	#*** convert string to prior list
	if ( base::is.character(prior) ){	
		prior <- prior_model_parse( prior_model = prior )		
	}
	#*** correct list specification
	parnames <- base::names(pars)
	NP <- base::length(parnames)
	prior1 <- base::as.list(1:NP)
	base::names(prior1) <- parnames

	#**** correct prior specification
	prior0 <- base::list( "dimproper" , base::list(NA) )
	parnames0 <- base::setdiff( parnames , base::names(prior) )
	NP0 <- base::length(parnames0)
	if (NP0 > 0 ){
		for (pp in 1:NP0){
			prior[[ parnames0[pp] ]] <- prior0
		}
	}	
	#**** write prior
	for (pp in 1:NP){
		prior1[[pp]] <- prior[[ parnames[pp] ]]
	}	
	dens1 <- prior_extract_density( prior=prior1 )
	dens1 <- base::data.frame("parameter" = base::names(pars) , "prior" = dens1 )
	base::rownames(dens1) <- NULL	
	res <- base::list("prior"=prior1 , "dens" = dens1 )					
	base::return(res)
}					
###########################################################					


