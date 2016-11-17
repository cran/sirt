
weighted_stats_extend_wgt <- function( wgt , mat )
{
	N1 <- base::nrow(mat)
	N2 <- base::ncol(mat)
	if ( base::is.null(wgt) ){
		wgt <- base::rep( 1 , N1 )
	}
	if ( base::is.vector(wgt) ){
		wgt <- base::matrix( wgt , nrow = N1 , ncol= N2 )	
	}	
	base::return(wgt)
}