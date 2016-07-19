 
#******************************************************************************************************
# Yen's Q3 statistic (1984)
##NS export(yen.q3)
Q3 <- function( dat , theta , b , progress = TRUE ){
        # INPUT:
        # dat       ... data frame
        # theta     ... theta estimate
        # b         ... item difficulty estimate
        base::cat("Yen's Q3 Statistic based on an estimated theta score \n*** " )
        I <- base::ncol(dat)
        base::cat(I , "Items | " )
        # expected probability
        expected <- .prob.rasch( theta , b )
        # residual 
        residual <- dat - expected  
        # initialize matrix of Q3 values
        I <- base::ncol(dat)
        q3.matr <- base::matrix( NA , nrow=I ,ncol= I )
        q3.long <- base::matrix( NA , nrow=(I-1)*I/2 , ncol=4 )
        base::colnames(q3.matr) <- base::rownames(q3.matr) <- base::colnames(dat)
		q3.matr <- stats::cor( residual , use = "pairwise.complete.obs")		
		nares <- 1 - base::is.na(residual)
#		NIP <- t(nares) %*% nares
		NIP <- base::crossprod(nares)
		itempairs <- base::t( utils::combn( I , 2 ) )
		q3.long[,3] <- q3.matr[ itempairs ]
		q3.long[,4] <- NIP[ itempairs ]
		q3.long[,1] <- base::colnames(q3.matr)[ itempairs[,1] ]
		q3.long[,2] <- base::colnames(q3.matr)[ itempairs[,2] ]			
        q3.long <- base::as.data.frame( q3.long )
        q3.long[,3] <- base::as.numeric( base::paste( q3.long[,3] ))
        q3.long <- q3.long[ base::order( q3.long[,3] ) , ]
        base::colnames(q3.long) <- base::c("Item1" , "Item2" , "Q3" , "N" )
        q3.long <-   q3.long[ !is.na( q3.long[,3] ) , ]
		base::cat( base::paste( base::nrow(q3.long) , "item pairs\n" ) )
		MQ3 <- base::mean( q3.long[,3] )
		SDQ3 <- stats::sd( q3.long[,3] )
		Q3.stat <- stats::quantile( q3.long[,3] , prob = c(  .10 , .25 , .50 , .75 , .90  ) )
		Q3.stat <- base::c("M" = MQ3 , "SD" = SDQ3 , "Min" = base::min(q3.long[,3] ) , 
						Q3.stat , "Max" = base::max(q3.long[,3] ) )
		base::cat("*** Q3 Descriptives\n")
		base::print( base::round(Q3.stat,3))
        res <- base::list( "q3.matrix" = q3.matr , "q3.long" = q3.long ,
				"expected" = expected , "residual" = residual , "Q3.stat" = Q3.stat )
        base::return(res)
}
#*****************************************************************************************************************

# Q3 <- yen.q3
