

##############################################################
# compute item probabilities
xxirt_compute_itemprobs <- function( item_list , items , Theta , ncat ,
		partable , partable_index , item_index = NULL ){
		TP <- base::nrow(Theta)
		maxK <- base::max(ncat)
		if ( is.null(item_index) ){ 
		    I <- base::length(items)
			item_index <- 1:I
		}			
		# I <- length(items)
		I <- base::length(item_index)		
		# compute item probabilities as a function of theta
		probs <- base::array( 0 , dim=c(I,maxK,TP) ) 		
		for (jj in 1:I){
			# ii <- 1
			ii <- item_index[jj]
			item_ii <- item_list[[ii]]
			par_ii <- partable[ partable_index[[ii]] , "value" ]
			arg_ii <- base::list( par = par_ii  , Theta = Theta , ncat = ncat[ii] )
			probs_ii <- base::do.call( item_ii$P , arg_ii )
			probs[ jj, 1:ncat[ii] ,] <- base::t(probs_ii)
		}
		base::return(probs)					
}
#############################################################################						