

mlnormalMatrix2List <- function( G , mat , freq_id ){
	matlist <- base::as.list(1:G)
	for (gg in 1:G){
		ind_gg <- base::seq( freq_id[ gg , "start"] , freq_id[ gg , "end"] )
		matlist[[gg]] <- mat[ ind_gg , base::seq( 1 , freq_id[ gg, "dim_id"] ) ]
	}
	base::return(matlist)
}