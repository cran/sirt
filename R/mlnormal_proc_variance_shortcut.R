
##########################################################################
# process data for shortcuts in variance estimation
mlnormal_proc_variance_shortcut <- function( id , y , X , Z_list ,
		Z_index , use_Rcpp , G ){
zz0 <- Sys.time()		
		#*** group sizes
		freq_id <- base::rowsum( 1+0*id , id )
		freq_id <- base::data.frame( 
						base::as.numeric( base::rownames(freq_id ) ) , 
						freq_id[,1] )
		base::colnames(freq_id) <- base::c("orig_id","dim_id")		
		freq_id$start_orig <-  1 + base::c(0 , 
								base::cumsum( freq_id[1:(G-1) , "dim_id"] ) )
		freq_id$end_orig <-  base::cumsum( freq_id[1:G , "dim_id"] )		
		freq_id <- freq_id[ base::order( freq_id[,2] ) , ]
		
		G <- base::nrow(freq_id)
		freq_id$id <- 1:G
		freq_id$update_dim <- base::c( 1,1 * ( base::diff(freq_id$dim_id) > 0 ) )
		freq_id$start <-  1 + base::c(0 , base::cumsum( freq_id[1:(G-1) , "dim_id"] ) )
		freq_id$end <-  base::cumsum( freq_id[1:G , "dim_id"] )
				
# cat("##### create freq_id") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1
 
		#--------------------------------------
		#---- check equality of Z_index and Z_list		
		if ( ! use_Rcpp ){
			mlnormal_proc_vs_Z <- mlnormal_proc_variance_shortcut_Z_R
		} else {
			mlnormal_proc_vs_Z <- mlnormal_proc_variance_shortcut_Z_Rcpp
		}			
# cat("**** vor mlnormal_proc_vs_Z")	

		res <- mlnormal_proc_vs_Z( Z_list=Z_list , Z_index = Z_index , G = G , 
		            freq_id = freq_id )
		freq_id <- res$freq_id
		rcpp_args <- res$rcpp_args
		
# cat("##### check equalities") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1									
 
		#--- do compute vector
		do_compute <- freq_id$update_dim == 1		
		if ( use_Rcpp ){
			rcpp_args$do_compute <- base::as.integer( do_compute )
		}
		
		#---------------------------------------
		#---- rearrange Z_index and Z_list
		Z_index <- Z_index[ freq_id[,1] , , , drop=FALSE]
		Z_list0 <- Z_list
		for (gg in 1:G){
			Z_list[[gg]] <- Z_list0[[  freq_id[gg,1] ]]
		}
#  cat("##### rearrange Z_list and Z_index") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1		

		#---------------------------------------
		#---- rearrange y and X
		if ( use_Rcpp){		
			mlnormal_proc_vs_XY <- mlnormal_proc_variance_shortcut_XY_R
		} else {  # The Rcpp function is slower than the R function
			mlnormal_proc_vs_XY <- mlnormal_proc_variance_shortcut_XY_R
		}				
		res <- mlnormal_proc_vs_XY(y=y, X=X, G=G, freq_id=freq_id)
		y <- res$y
		X <- res$X		
# cat("##### rearrange X and y") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1											
 
		id <- base::rep( 1:G , freq_id$dim_id )				
 
		#---------------------------------------
		#---------- output
		res <- base::list( id = id , y = y , X = X , Z_list = Z_list ,
					Z_index = Z_index , freq_id = freq_id , do_compute = do_compute ,
					rcpp_args = rcpp_args )
		base::return(res)
}
#############################################################################