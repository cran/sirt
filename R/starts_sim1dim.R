
#######################################################
# STARTS unidimensional model: simulation 
starts_sim1dim <- function( N , W , var_trait , 
	   var_ar , var_state , a )
{

	# var_state <- var_disturbance / ( 1 - a^2 )	
	#  Var(S) = a^2 Var(S) + Var(U)
	#   => Var(U) = (1-a^2)*Var(S)
	
	#--- renaming input arguments
	var_error <- var_state
	var_state <- var_ar
	#---
	
	var_disturbance <- (1-a^2)*var_state
	
	matr <- base::matrix( NA , nrow=N , ncol=W )
	base::colnames(matr) <- paste0("W",1:W)

	states <- matr

	trait <- stats::rnorm( N , sd = base::sqrt( var_trait ))

	ww <- 1
	states[,ww] <- stats::rnorm( N , sd = base::sqrt( var_state ))
	for ( ww in 2:W){
		states[,ww] <- a*states[,ww-1] + stats::rnorm( N , sd = base::sqrt(var_disturbance))
	}
	matr <- trait + states + stats::rnorm( N*W , sd = base::sqrt( var_error ))
	base::return(matr)
}
#######################################################################		