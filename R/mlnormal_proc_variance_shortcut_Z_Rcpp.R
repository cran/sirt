

mlnormal_proc_variance_shortcut_Z_Rcpp <- function( Z_list , Z_index , G , freq_id){
	#-- create list with arguments for Rcpp functions
	rcpp_args <- base::list( 
					"Z_index" = base::as.vector(Z_index) ,
					"dim_Z_index" = base::dim(Z_index) ,
					"update_dim" = base::as.vector(freq_id[,"update_dim"]) ,
					"start_orig" = base::as.vector(freq_id[,"start_orig"])-1 ,	
					"end_orig" = base::as.vector(freq_id[,"end_orig"])-1 ,
					"orig_id" = base::as.vector(freq_id[,"orig_id"]) ,
					"dim_id" = base::as.vector(freq_id[,"dim_id"])	,
					"start" = base::as.vector(freq_id[,"start"])	,
					"end" = base::as.vector(freq_id[,"end"]),
					"N" = base::as.integer( base::max( freq_id[,"end"] ) ) ,
					"max_dim" = base::as.integer( base::max( freq_id[,"dim_id"] ))
							)

	#--- load Rcpp function
	res <- CallSwitch( "mlnormal_proc_variance_shortcut_Z_restructure" ,
				Z_list , rcpp_args$update_dim , rcpp_args$start_orig , 
				rcpp_args$end_orig , rcpp_args$dim_Z_index , rcpp_args$Z_index , 
				rcpp_args$orig_id, rcpp_args$dim_id , PACKAGE="sirt")
	freq_id$update_dim <- res$update_dim[,1]						
	
	#--- output
	res <- base::list( "freq_id" = freq_id , "rcpp_args" = rcpp_args )
	base::return(res)
}
		