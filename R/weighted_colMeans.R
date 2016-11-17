
weighted_colMeans <- function( mat , wgt=NULL){
	wgt <- weighted_stats_extend_wgt( wgt=wgt , mat=mat )
	mat1 <- base::colSums( mat * wgt , na.rm=TRUE) 
	mat2 <- base::colSums( wgt , na.rm=TRUE) 
	mat1 <- mat1 / mat2
	base::return(mat1)
}