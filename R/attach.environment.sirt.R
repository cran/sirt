#######################################
# attach all elements in an environment   
.attach.environment.sirt <- function( res , envir ){
#	e1 <- environment()
	CC <- length(res)
	for (cc in 1:CC){
		assign( names(res)[cc] , res[[cc]] , envir=envir )		
					}
			}
