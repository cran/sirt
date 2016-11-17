

weighted_rowSums <- function( mat , wgt=NULL){
	wgt <- weighted_stats_extend_wgt( wgt=wgt , mat=mat )
	mat1 <- base::rowSums( mat * wgt , na.rm=TRUE) 
	base::return(mat1)
}