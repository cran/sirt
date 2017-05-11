

###########################################################################
xxirt_compute_likelihood <- function( probs_items , dat , resp_index,
		dat_resp=NULL )
{		
# zz0 <- Sys.time()
	N <- nrow(dat)
	TP <- dim(probs_items)[3]
	I <- dim(probs_items)[1]
	maxK <- dim(probs_items)[2]
	p.xi.aj <- matrix( 1 , nrow=N , ncol=TP )
# cat("   * bef1 ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		
	probs <- matrix( probs_items , nrow=I , ncol=maxK*TP )
	p.xi.aj <- xxirt_compute_likelihood_helper( dat = dat, dat_resp = dat_resp ,
					probs = probs, TP=TP , maxK = maxK )
	return(p.xi.aj)
}
###########################################################################			
