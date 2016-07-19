
		
#################################################
# clean string		
prior_model_pars_CleanString <- function( ps ){		
	ps <- base::gsub( " " , "" , ps )
	ps <- ps[ ps != "" ]
	NP <- base::length(ps)
	for (pp in 1:NP){		
		# pp <- 1
		ps_pp <- ps[pp]	
		# locate comment symbol
		h1 <- base::gregexpr(pattern ='#', text= ps_pp)
		if ( h1 > 0){
			ps[pp] <- base::substring( ps_pp , 1 , h1[[1]] -1)	
		}
	}
	ps <- ps[ ps != "" ]
	base::return(ps)
}
##################################################	