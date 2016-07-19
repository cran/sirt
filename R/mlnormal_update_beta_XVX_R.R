
################################################################
# update beta
mlnormal_update_beta_XVX_R <- function( NB , Z_index , G ,
		V1_list , X_list , y_list , rcpp_args , X , y){

		XVX <- base::matrix( 0 , nrow=NB , ncol=NB)
		XVY <- base::matrix( 0 , nrow=NB , ncol=1)
		dimZ <- base::dim( Z_index )
		XV_list <- base::as.list(1:G)
	# zz0 <- Sys.time()
		for (gg in 1:G){
			# gg <- 1
			# compute V for group gg
#			Z_index_gg <- Z_index[gg,,,drop=FALSE]
#			Z_list_gg <- Z_list[[gg]] 
			V_gg1 <- V1_list[[gg]]		
			X_gg <- X_list[[gg]]
			y_gg <- y_list[[gg]]
			# XV_gg <- t(X_gg) %*% V_gg1	
			XV_gg <- base::crossprod( X_gg , V_gg1	)
			XV_list[[gg]] <- XV_gg
			XVX <- XVX + XV_gg %*% X_gg
			XVY <- XVY + XV_gg %*% y_gg
		}
		#--- output
		res <- base::list( XVX = XVX , XVY = XVY , XV_list = XV_list )
		base::return(res)
}
