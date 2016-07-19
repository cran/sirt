
############################################
# compute EAP and its standard deviation
xxirt_EAP <- function(p.aj.xi , Theta ){
	D <- base::ncol(Theta)
	e1 <- p.aj.xi %*% Theta
	base::colnames(e1) <- base::paste0("EAP.Dim",1:D)
	e2 <- p.aj.xi %*% Theta^2
	base::colnames(e2) <- base::paste0("SD.EAP.Dim",1:D)	
	e2 <- e2 - e1^2 
	res <- base::cbind( e1 , e2)
	res <- res[ , base::rep( c(0,D) , D ) + 1:D ]
	base::return(res)
}
################################################		