
##########################################
# coef S3 method for xxirt objects
coef.xxirt <- function(object,...){
		par1 <- xxirt_partable_extract_freeParameters( object$partable )
		par2 <- xxirt_parTheta_extract_freeParameters( object$customTheta )
		par <- base::c(par1 , par2)
		base::return(par)
}
##########################################		