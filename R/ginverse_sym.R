
#######################################################
# code from Eugene Demidenko: book mixed effects models
ginverse_sym <- function(A, eps= 1E-8){
    # Generalized inverse of a symmetric matrix A
    PV <- base::eigen(A,symmetric=TRUE)
    V0 <- IV <- PV$values
	av0 <- base::abs(V0)
    IV[ av0 > eps] <- 1/V0[ av0 > eps]
    IV[ av0 <= eps] <- 0
    Ainv <- PV$vectors %*% ( IV*( base::t(PV$vectors) ) )
    base::return(Ainv)
}
#######################################################
