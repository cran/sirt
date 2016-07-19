
#################################################
# process parameter table
xxirt_proc_ParTable <- function( itemtype , partable , items ){								
		#*** extract item types from partable
		itemtype <- base::unlist( base::sapply( items , FUN = function(ii){ 
						partable$type[ base::paste(partable$item) == ii ][1] 
						} )  )					
		I <- base::length(items)
		partable$rowindex <- base::seq( 1 , base::nrow(partable) )
		#*** parameter index
		partable[ ! partable$est , "parindex" ] <- 0
		# indices <- sort( unique( partable$parindex ) )
		indices <- base::unique( partable$parindex )
		indices <- base::c( 0 , base::setdiff( indices , base::c(0) ) )
		IN <- base::length(indices)
		partable$parindex <- base::match( partable$parindex , indices ) - 1									
		partable$parindex[ partable$parindex == 0 ] <- NA
		#*** set prior distributions of fixed parameters to NA
		partable[ ! partable$est , "prior" ] <- NA
		#*** list with parameter table indices		
		partable$parfree <- 1*partable$est
		partable_index <- base::as.list( 1:I )
		
		for (ii in 1:I){
			partable_index[[ii]] <- base::which( partable$itemnr == ii )
		}
		ind <- base::which( base::duplicated( partable$parindex ) & ( ! is.na( partable$prior ) ) )	
		if ( base::length(ind) > 0 ){
			partable[ind , c("prior","prior_par1","prior_par2", "parlabel") ] <- NA		
		}	
		ind <- base::which( base::duplicated( partable$parindex ) )	
		if ( base::length(ind) > 0 ){
			partable[ind,"parfree"] <- 0
		}									
														
		#*** extract ncat and maxK
		p1 <- partable[ ! base::duplicated( partable$item) , ]
		ncat <- p1$ncat	
		base::names(ncat) <- base::paste(p1$item)
		maxK <- base::max(ncat)
		#*** extract M-step method
		p1 <- partable[ partable$parfree == 1 , ]
        m1 <- TRUE
		if ( base::nrow(p1) > 0 ){
			m1 <- ( base::mean( p1$lower == - Inf ) < 1 ) | ( base::mean( p1$upper == Inf ) < 1 )
		}
		mstep_method <- if (m1){ "L-BFGS-B" } else { "BFGS" }
		
		#**** item indices per parameter
		NP <- -Inf
        if ( base::sum(partable$est) > 0 ){	
			NP <- base::max( partable$parindex , na.rm=TRUE )
		}
		
		if ( NP > -Inf){
			item_index <- base::as.list( 1:NP )
			for (pp in 1:NP){
				p1 <- partable[ paste0( partable$parindex) == pp , ]
				base::names(item_index)[pp] <- p1[1,"parlabel"]
				item_index[[pp]] <- p1$itemnr			
			}
		} else {
			item_index <- base::list()								
		}
				
		#----- output
		res <- base::list( itemtype = itemtype , partable = partable ,
					partable_index = partable_index , ncat = ncat ,
					maxK = maxK , mstep_method = mstep_method , 
					item_index=item_index )
		base::return(res)									
}
################################################						
