 

#-----------------------------------------------------------------------------------------------------
detect.index <- function( ccovtable , itemcluster ){
    # INPUT:
    # result from ccov.np 
    # itemcluster ... identifies an item cluster for each item
    #.............................
    # calculate delta
    ccovtable <- ccovtable$ccov.table
    ccovtable$delta <- ifelse( itemcluster[ ccovtable$item1ID ] == 
				itemcluster[ ccovtable$item2ID ] , 1 , -1 )
	
	#******************************************
    #--- calculate unweighted and weighted indizes
	ccov <- ccovtable$ccov			# conditional covariance
	delta <- ccovtable$delta		# indicator for partition
	N <- ccovtable$N
	sqrt_N <- base::sqrt(N)
	sign_ccov <- base::sign(ccov)
	abs_ccov <- base::abs( ccov )
	# number of parameters
	np <- 5
	parnames <- weighted.indizes <- indizes <- rep(NA,np)
	
	#--- DETECT
	ii <- 1
	indizes[ii] <- 100*base::mean(ccov*delta)
	weighted.indizes[ii] <- 100*stats::weighted.mean( ccov * delta , sqrt_N ) 
	parnames[ii] <- "DETECT"
	#--- ASSI
	ii <- 2
	indizes[ii] <- base::mean( sign_ccov * delta )
	weighted.indizes[ii] <- stats::weighted.mean( sign_ccov * delta , sqrt_N )
	parnames[ii] <- "ASSI"	
	#--- RATIO
	ii <- 3
	indizes[ii] <- base::sum( ccov * delta ) / base::sum( abs_ccov )
	weighted.indizes[ii] <- base::sum( ccov * delta * sqrt_N ) / base::sum( abs_ccov * sqrt_N )
	parnames[ii] <- "RATIO"	
	#--- MADCOV
	ii <- 4
	indizes[ii] <- 100 * base::mean( abs_ccov )
	weighted.indizes[ii] <- 100* stats::weighted.mean( abs_ccov , sqrt_N ) 
	parnames[ii] <- "MADCOV100"		
	#--- MCOV
	ii <- 5
	indizes[ii] <- 100 * base::mean( ccov )
	weighted.indizes[ii] <- 100* stats::weighted.mean( ccov , sqrt_N ) 
	parnames[ii] <- "MCOV100"		
		
	#******************************************
    #--- calculate weighted indizes
    res <- base::data.frame( "unweighted" = indizes , "weighted" = weighted.indizes )
    rownames(res) <- parnames
    return(res)
}
#-----------------------------------------------------------------------------------------------------

