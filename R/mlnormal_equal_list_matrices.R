

mlnormal_equal_list_matrices <- function( list_mat1  , list_mat2 , dim_list , 
        eps=1E-30){
		if ( base::is.null(dim_list) ){
			dim_list <- base::length(list_mat1)
		}
		v1 <- 1
		for ( dd in 1:dim_list ){
			v1 <- v1 * mlnormal_equal_matrix( mat1 =list_mat1[[dd]] , 
						mat2 = list_mat2[[dd]] , eps=eps )
		}
		res <- if ( v1 == 1 ){ TRUE } else { FALSE }
		base::return(res)
}