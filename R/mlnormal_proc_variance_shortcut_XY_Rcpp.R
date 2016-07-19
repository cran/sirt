

mlnormal_proc_variance_shortcut_XY_Rcpp <- function(y , X , G , freq_id){
		freq_id <- base::as.matrix( freq_id )
		X <- base::as.matrix(X)
		#--- Rcpp function
		res <- CallSwitch( "mlnormal_proc_variance_shortcut_XY_restructure" ,
					freq_id ,  y , X, G , PACKAGE = "sirt" ) 
		base::return(res)
}
		