
############################################################
# evaluate prior in M-step
xxirt_mstep_itemParameters_evalPrior <- function(partable, h = 0){

        partable1 <- partable[ partable$parfree == 1 , ]
		NP <- nrow(partable1)		
		NP1 <- max( NP , 1)
		#*** evaluate prior distributions in partable
		pen <- rep(0,NP1)
        if (NP>0){			
			for (pp in 1:NP){
				# pp <- 1
				if ( ! is.na( partable1[pp,"prior"] ) ) {				
						prior_pp <- partable1[pp,"prior"]
						val <- partable1[pp,"value"] + h
						prior_val <- do.call( prior_pp , list( val , partable1[pp,"prior_par1"] ,  
									partable1[pp,"prior_par2"] ) )											
						pen[pp] <- - log( prior_val + 1E-300 )
								}  # end if prior
						   }	# end pp	
# Revalpr("pen")						   
						}  # end if NP > 0
		return(pen)
			}
##########################################################			