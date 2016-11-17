
smirt_squeeze <- function( val , lower , upper, est)
{
	D <- 1
	is_matrix <- FALSE
	if ( base::is.matrix(est) ){
		D <- base::ncol(est)
		is_matrix <- TRUE
	}
	est <- base::matrix( est , ncol=D)
	val0 <- val
	for (dd in 1:D){
		val[,dd] <- base::ifelse( val[,dd] < lower , lower , val[,dd] )
		val[,dd] <- base::ifelse( val[,dd] > upper , upper , val[,dd] )
		ind_dd <- base::which(est[,dd] == 0)
		if ( base::length(ind_dd) > 0 ){
			val[ ind_dd , dd] <- val0[ ind_dd,dd]
		}
	}
	if ( ! is_matrix){
		val <- val[,1]
	}
	base::return(val)
}