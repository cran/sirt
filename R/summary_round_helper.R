
summary_round_helper <- function( obji , digits , exclude = NULL, print=TRUE){
	NC <- base::ncol(obji)
	ind <- 1:NC
	if ( ! base::is.null(exclude) ){
		ind2 <- base::which( base::colnames(obji) %in% exclude )
		ind <- base::setdiff( ind , ind2 )
	}
	obji[,ind] <- base::round( obji[,ind] , digits )
	base::rownames(obji) <- NULL
	base::print(obji)
	base::invisible(obji)	
}