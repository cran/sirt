

IRT.se.xxirt <- function( object , ...){
	c1 <- coef(object)
	v1 <- vcov(object)
	par1 <- xxirt_partable_extract_freeParameters( object$partable )
	par2 <- xxirt_parTheta_extract_freeParameters( object$customTheta )
	N1 <- base::length(par1)
	N2 <- base::length(par2)
	dfr <- data.frame("partype"= base::c( base::rep("item",N1), base::rep("Theta",N2) ) )
	dfr$parlabel <- names(c1)
	dfr$value <- c1
	dfr$se <- base::sqrt( base::diag(v1) )
	base::return(dfr)
}