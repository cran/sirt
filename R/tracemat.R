
########################################
# trace of a matrix
tracemat <- function(A)
{ 
    res <- base::sum( base::diag(A) ) 
	base::return(res)
} 
########################################