
##############################################################
# write elements from mcmcmlist into code file
mcmclist2coda <- function( mcmclist , name , coda.digits=5 ){
    m1 <- mcmclist[[1]]
    vars <- base::colnames(m1)
    #--- create codaIndex file
    BB <- base::nrow(m1)
    VV <- base::length(vars)
    c1 <- base::paste( vars , base::seq( 1 , BB*VV , BB ) , 
	                    base::seq( BB , BB*VV , BB ) )
    base::writeLines( c1 , base::paste0( name , "_codaIndex.txt" ) )
    #--- create coda file
    m2 <- base::matrix( m1 , ncol=1 )
    m2 <- base::paste( base::rep(1:BB , VV ) , base::round( m2[,1] , coda.digits ) )
    base::writeLines( m2 , base::paste0( name , "_coda1.txt" ) )
            }
#########################################################
