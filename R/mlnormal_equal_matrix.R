
mlnormal_equal_matrix <- function( mat1  , mat2 , eps=1E-30){
		res <- base::all( mat1 == mat2 )
		base::return(res)
}
####################
#			maxval <- base::max( base::abs( mat1 - mat2 ))
#			res <- if ( maxval < eps ){ TRUE } else { FALSE }
#		res <- base::identical( mat1 , mat2 )		
#		res <- base::all.equal( mat1 , mat2 )