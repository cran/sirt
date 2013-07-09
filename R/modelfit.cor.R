
##########################################
# Modelfit in sirt
modelfit.sirt <- function( object ){
	#*****
	# rasch.mml
	if (class(object)=="rasch.mml"){
		mod <- object
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
			probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
			probs[ , 2 , ] <- t(prob1)
			probs[ , 1 , ] <- 1 - t(prob1)
		dat <- mod$dat
					}
	#*****
	# rasch.mirtlc
	if (class(object)=="rasch.mirtlc"){
		mod <- object$estep.res
		posterior <- mod$f.qk.yi
		prob1 <- mod$pjk
			probs <- array( NA , dim=c( ncol(prob1) , 2 , nrow(prob1)) )
			probs[ , 2 , ] <- t(prob1)
			probs[ , 1 , ] <- 1 - t(prob1)
		dat <- object$dat
					}	
	#******
	# rasch.pml
	if ( class(object) !="rasch.pml"){ pmlobject <- NULL } else {
		data <- NULL ; posterior <- NULL ; probs <- NULL ; pmlobject <- object }
	#*******
	# smirt	
	if (class(object) == "smirt"){
		# note that for polytomous response data some adaptations are
		# necessary: see modelfit in the CDM package
		mod <- object
		probs <- mod$probs
		posterior <- mod$f.qk.yi
		dat <- mod$dat
					}	
	# calculate modelfit.cor	
	res <- modelfit.cor( data = dat , posterior =posterior , probs = probs ,
			pmlobject=pmlobject)
	return(res)
	}
################################################################################

#############################################################################
modelfit.cor <-
function( data=NULL , posterior=NULL , probs=NULL , pmlobject=NULL ){
	if ( is.null(pmlobject)){
		data.resp <- 1 - is.na(data)
		data[ is.na(data) ] <- 9
		data1 <- data*data.resp
		I <- ncol(data)
		# calculate counts (ignore weights here!!) 
		n11 <- t(  ( data==1) * data.resp ) %*% ( ( data==1) * data.resp )
		n10 <- t(  ( data==1) * data.resp ) %*% ( ( data==0) * data.resp )
		n01 <- t(  ( data==0) * data.resp ) %*% ( ( data==1) * data.resp )
		n00 <- t(  ( data==0) * data.resp ) %*% ( ( data==0) * data.resp )
		
		p1 <- colMeans(  ( data==1) * data.resp ) 
		# p0 <- colMeans(  ( data==0) * data.resp ) 
		
		# expected counts
		exp1 <- rep(NA, I )
		for (ii in 1:I){
			# ii <- 1		
#			pr.ii1 <- matrix( probs[ii,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#			p3ii <-  pr.ii1 * posterior
#			exp1[ii] <- sum( rowSums( p3ii ) * data.resp[,ii ] ) / sum( data.resp[,ii] )
			exp1[ii] <- sum( colSums( posterior * data.resp[,ii] ) * probs[ii,2,] ) / sum( data.resp[,ii] )
					}
		#********************************
		# covariances 
		
		ip <- itempairs <- t( combn(I,2 ) )
		colnames(itempairs) <- c("item1" , "item2" )
		itempairs <- as.data.frame( itempairs )            
		itempairs$n11 <- n11[ ip ]
		itempairs$n10 <- n10[ ip ]
		itempairs$n01 <- n01[ ip ]
		itempairs$n00 <- n00[ ip ]
		itempairs$n <- rowSums( itempairs[ , c("n11","n10", "n01","n00") ] )		
        itempairs$Exp00 <- itempairs$Exp01 <- itempairs$Exp10 <- itempairs$Exp11 <- NA		
						}
	##################################################						
    if ( ! is.null( pmlobject ) ){

		ip0 <- pmlobject$itempairs
		I <- max(pmlobject$itempairs[,c("item1","item2")])
		itempairs <- ip0[ , c("item1","item2") ]
		itempairs$n11 <- ip0$f11
		itempairs$n10 <- ip0$f10		
		itempairs$n01 <- ip0$f01		
		itempairs$n00 <- ip0$f00		
		itempairs$n <- rowSums( itempairs[ , c("n11","n10", "n01","n00") ] )				
		itempairs$Exp00 <- ip0$p00 * itempairs$n
		itempairs$Exp10 <- ip0$p10 * itempairs$n
		itempairs$Exp01 <- ip0$p01 * itempairs$n
		itempairs$Exp11 <- ip0$p11 * itempairs$n	
	    p1 <- NULL		
					}
    itempairs$corExp <- itempairs$corObs <- NA
    
    m1 <- matrix( c(1,1,1,0,0,1,0,0) , 4 , 2 , byrow=T )
    
	# define further quantities
	itempairs$X2 <- NA
#	itempairs$G2 <- NA	
	itempairs$RESIDCOV <- NA		
	itempairs$Q3 <- NA			
	
	#***
	# calculate expected score for every person and every item
	exp.ii.jj <- posterior %*% t( probs[,2,] )
	#***
	
    for (ii in 1:(I-1) ){
        for (jj in (ii+1):I){
    # ii <- 1
    # jj <- 2
        ii1 <- which ( itempairs$item1 == ii &  itempairs$item2 == jj )
		ps.iijj <- colSums( posterior[ data.resp[,ii]*data.resp[,jj]>0 , ] )
		
	if ( is.null( pmlobject)){
        diijj <- data.resp[,ii ]*data.resp[,jj ]
#        pr.ii1 <- matrix( probs[ii,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#        pr.jj1 <- matrix( probs[jj,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )    
#        p3ii <-  pr.ii1 * pr.jj1 * posterior
#        itempairs[ii1,"Exp11"] <- sum( rowSums( p3ii ) * diijj )    
		itempairs[ii1,"Exp11"] <- sum( probs[ii,2,]*probs[jj,2,] * ps.iijj )
	
#        pr.ii1 <- matrix( probs[ii,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#        pr.jj1 <- matrix( probs[jj,1,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )    
#        p3ii <-  pr.ii1 * pr.jj1 * posterior
#        itempairs[ii1,"Exp10"] <- sum( rowSums( p3ii ) * diijj )
		itempairs[ii1,"Exp10"] <- sum( probs[ii,2,]*probs[jj,1,] * ps.iijj )
    
#        pr.ii1 <- matrix( probs[ii,1,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#        pr.jj1 <- matrix( probs[jj,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )    
#        p3ii <-  pr.ii1 * pr.jj1 * posterior
#        itempairs[ii1,"Exp01"] <- sum( rowSums( p3ii ) * diijj )
		itempairs[ii1,"Exp01"] <- sum( probs[ii,1,]*probs[jj,2,] * ps.iijj )		
    
#        pr.ii1 <- matrix( probs[ii,1,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#        pr.jj1 <- matrix( probs[jj,1,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )    
#        p3ii <-  pr.ii1 * pr.jj1 * posterior
#        itempairs[ii1,"Exp00"] <- sum( rowSums( p3ii ) * diijj )
		itempairs[ii1,"Exp00"] <- sum( probs[ii,1,]*probs[jj,1,] * ps.iijj )				
		
							}
							
        itempairs[ii1, "corObs"]  <-   .corr.wt( x = m1[,1,drop=FALSE] ,  y = m1[,2,drop=FALSE] , 
            w = as.numeric( itempairs[ii1,c("n11","n10","n01","n00") ] ) )
    
        itempairs[ii1, "corExp"]  <-   .corr.wt( x = m1[,1,drop=FALSE] ,  y = m1[,2,drop=FALSE] , 
            w = as.numeric( itempairs[ii1,c("Exp11","Exp10","Exp01","Exp00") ] ) )
		#***
		# Q3 statistic
	if ( is.null( pmlobject)){		
#        pr.ii1 <- matrix( probs[ii,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )
#        pr.jj1 <- matrix( probs[jj,2,] , nrow= nrow(data) , ncol= dim(probs)[3] , byrow=T )    
#        p3ii <-  pr.ii1 * posterior
#		expii <- rowSums( p3ii )
#        p3jj <-  pr.jj1 * posterior
#		expjj <- rowSums( p3jj )
		# calculate residuals
#		data.res <- data[, c(ii,jj) ] - cbind( expii , expjj )
#		data.res <- data.res[ diijj == 1 , ]
#		itempairs[ii1,"Q3"] <- cor(data.res)[1,2]
		data.res <- data[, c(ii,jj) ] - exp.ii.jj[ , c(ii,jj) ]
		data.res <- data.res[ diijj == 1 , ]
		itempairs[ii1,"Q3"] <- cor(data.res)[1,2]		
		
								}
                            }
# print(ii) ; flush.console()							
                }
    
	##############################
	itempairs$X2 <- ( itempairs$n00 - itempairs$Exp00 )^2 / itempairs$Exp00 +
					( itempairs$n10 - itempairs$Exp10 )^2 / itempairs$Exp10	+
					( itempairs$n01 - itempairs$Exp01 )^2 / itempairs$Exp01	+
					( itempairs$n11 - itempairs$Exp11 )^2 / itempairs$Exp11						
	# G2
#	itempairs$G2 <- itempairs$n00 * log( itempairs$Exp00 / max( itempairs$n00 , .01 ) ) +
#					itempairs$n01 * log( itempairs$Exp01 / max( itempairs$n01 , .01 ) ) +
#					itempairs$n10 * log( itempairs$Exp10 / max( itempairs$n10 , .01 ) ) +
#					itempairs$n11 * log( itempairs$Exp11 / max( itempairs$n11 , .01 ) ) 
#    itempairs$G2 <- -2*itempairs$G2					
	itempairs$RESIDCOV <- ( itempairs$n11 * itempairs$n00 - itempairs$n10 * itempairs$n01 ) / itempairs$n^2 -
				( itempairs$Exp11 * itempairs$Exp00 - itempairs$Exp10 * itempairs$Exp01 ) / itempairs$n^2	
	##############################
	# labels
    itempairs$item1 <- colnames(data)[ itempairs$item1 ]
    itempairs$item2 <- colnames(data)[ itempairs$item2 ]

	modelfit <- data.frame( "est" = c( 
			mean( abs( itempairs$corObs - itempairs$corExp ) ) ,
			mean( itempairs$X2 ) , # mean( itempairs$G2) ,
			mean( 100*abs(itempairs$RESIDCOV ) ) ,
			mean( abs( itempairs$Q3 ) )
						) )
	rownames(modelfit) <- c("MADcor" , "MX2" , # "MG2",
				"100*MADRESIDCOV" , "MADQ3" )
    
#     "pfit" <- data.frame( "item" = colnames(data) , "pObs" = p1 , "pExp" = exp1 )
	print( round(modelfit,5) , digits=5 )   
#    cat("MAD Correlation (Observed minus Expected)" , round( MADcor , 4 ) , "\n" )    
    res <- list( "modelfit" = modelfit , "itempairs" = itempairs ) # , "pfit" = pfit )
    return(res)
    }

#######################################################################	
# auxiliary function: weighted correlation	
.corr.wt <- function( x, y, w = rep(1,length(x))) {
#  stopifnot(length(x) == dim(y)[2] )
  w <- w / sum(w)
  # Center x and y, using the weighted means
  x <- x - sum(x * w)
  ty <- y - sum( y * w)
  # Compute the variance
  vx <- sum(w * x * x)
  vy <- sum(w * ty * ty)
  # Compute the covariance
  vxy <- sum(ty * x * w)
  # Compute the correlation
  vxy / sqrt(vx * vy)
}
##################################################################################