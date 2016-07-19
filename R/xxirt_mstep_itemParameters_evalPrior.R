
############################################################
# evaluate prior in M-step
xxirt_mstep_itemParameters_evalPrior <- function(partable, h = 0){
		
		eps <- 1E-300
        partable1 <- partable[ partable$parfree == 1 , ]
		NP <- base::nrow(partable1)		
		NP1 <- base::max( NP , 1)
		#*** evaluate prior distributions in partable
		pen <- base::rep(0,NP1)
        if (NP>0){			
			for (pp in 1:NP){
				# pp <- 1
				if ( ! base::is.na( partable1[pp,"prior"] ) ) {				
					prior_pp <- partable1[pp,"prior"]
					val <- partable1[pp,"value"] + h
					prior_args_pp <- base::list( val , partable1[pp,"prior_par1"] ,  
												 partable1[pp,"prior_par2"] )
					prior_val <- base::do.call( prior_pp , prior_args_pp )
					pen[pp] <- - base::log( prior_val + eps )
				}  # end if prior
			  }	# end pp	
		}  # end if NP > 0
		base::return(pen)
}
##########################################################			