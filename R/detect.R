 

#-----------------------------------------------------------------------------------------------------
detect.index <- function( ccovtable , itemcluster ){
    # INPUT:
    # result from ccov.np 
    # itemcluster ... identifies an item cluster for each item
    #.............................
    # calculate delta
    ccovtable <- ccovtable$ccov.table
    ccovtable$delta <- ifelse( itemcluster[ ccovtable$item1ID ] == itemcluster[ ccovtable$item2ID ] , 1 , -1 )
    # calculate indizes
    indizes <- c( 100*mean( ccovtable$ccov * ccovtable$delta ) , 
                mean( sign( ccovtable$ccov ) * ccovtable$delta ) , 
                    sum( ccovtable$ccov  * ccovtable$delta ) / sum( abs( ccovtable$ccov ) ) )
	# calculate weighted indizes 
    weighted.indizes <- c( 100* stats::weighted.mean( ccovtable$ccov * ccovtable$delta , sqrt(ccovtable$N) ) ,
        stats::weighted.mean( sign( ccovtable$ccov ) * ccovtable$delta , sqrt(ccovtable$N) ) , 
        sum( ccovtable$ccov  * ccovtable$delta * sqrt(ccovtable$N) ) / sum( abs( ccovtable$ccov ) * sqrt(ccovtable$N) ) )
    res <- data.frame( "unweighted" = indizes , "weighted" = weighted.indizes )
    rownames(res) <- c("DETECT" , "ASSI" , "RATIO" )
    res
    }
#-----------------------------------------------------------------------------------------------------


#*********************************************************
# auxiliary function for creating a covariance matrix
.create.ccov <- function( cc , data ){
    ccc <- cc$ccov.table
    I <- max( ccc$item1ID , ccc$item2ID )
    ccov.matrix <- matrix( 0 , I , I)
    rownames(ccov.matrix) <- colnames(ccov.matrix) <- colnames(data)
    LL <- nrow(ccc)
    for (ll in 1:LL){ 
            ccov.matrix[ ccc$item1ID[ll] , ccc$item2ID[ll] ] <- ccc$ccov[ll]
            ccov.matrix[ ccc$item2ID[ll] , ccc$item1ID[ll] ] <- ccov.matrix[ ccc$item1ID[ll] , ccc$item2ID[ll] ]
                        }
    return( ccov.matrix) 
        }
#*********************************************************
