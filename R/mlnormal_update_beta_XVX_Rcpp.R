
################################################################
# update beta
mlnormal_update_beta_XVX_Rcpp <- function( NB , Z_index , G ,
		V1_list , X_list , y_list , rcpp_args , X , y){
		# PACKAGE <- "R"
		PACKAGE <- "sirt"
		res <- CallSwitch( "mlnormal_update_beta_rcpp_helper" ,
					rcpp_args$dim_id , rcpp_args$start , rcpp_args$end ,
					G , X , y , rcpp_args$V1 ,
					PACKAGE = PACKAGE )
		#--- output
		res$XV_list <- NULL
		base::return(res)
}
