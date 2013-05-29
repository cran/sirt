 
# 0.01  2012-xx-yy


# 0.01  2012-06-23  o initial release


#-------------------------------------------------------




########################################################################################
# stratified Cronbach's Alpha
##NS export(stratified.cronbach.alpha)
stratified.cronbach.alpha <- function( data , itemstrata ){
	# function .cronbach.alpha
    .cronbach.alpha <- function( data ){ 
        # covariance
        c1 <- cov( data , use = "pairwise.complete.obs" )
        # mean covariance
        c1a <- c1 ; diag(c1a) <- 0
        I <- ncol(data)
        mc <- sum(c1a) / ( I^2 - I )
        # mean variance
        mv <- mean( diag(c1) )
        alpha <- I * mc / ( mv + (I-1) *mc )
        mean.tot <- mean( rowSums(data) )
        var.tot <- var( rowSums( data ) )
        res <- list( "I" = I , "alpha" = alpha , "mean.tot" = mean.tot ,  "var.tot" = var.tot )
        return(res)
            }
	#---------------------------------------
    dfr <- data.frame( "scale" = "total" , .cronbach.alpha(data) )
    # calculation of stratified alpha
    itemstrata.u <- sort(unique( itemstrata[,2] ))
    for (gg in itemstrata.u){
        #gg <- itemstrata.u[1]
        dfr1 <- data.frame( "scale" = gg , .cronbach.alpha( data[ , itemstrata[ itemstrata[,2] == gg , 1] ] ) )
        dfr <- rbind( dfr , dfr1 )
            }
    # stratified alpha
    dfr$alpha.stratified <- NA
    dfr$alpha.stratified[1] <- 1 - sum ( ( 1 -  dfr[ -1 , "alpha" ] ) * dfr[ -1 , "var.tot" ] ) / dfr[1,"var.tot" ]
	obji <- dfr
	obji[ , -c(1:2)] <- round( obji[,-c(1:2) ] , 3 )	
	print( obji )
    invisible(dfr)
    }
########################################################################################

