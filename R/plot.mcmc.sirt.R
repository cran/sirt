
######################################################
# plot results of objects of class mcmc.sirt
plot.mcmc.sirt <- function( x , layout=1 , conflevel=.90 , 
	round.summ=3 , lag.max=100 , col.smooth="red" , lwd.smooth=1 , 
	col.ci="orange" , cex.summ=1 , ask=FALSE , ... ){
	
	object <- x	# rename x into object
	
	# layout type
	# layout=1 : standard output from coda package	
	# library(coda)
    if (layout==1){  plot(object$mcmcobj , ask=ask , ...) }         

	#***************
	# layout =2
	if (layout==2){
		mcmcobj <- (object$mcmcobj)[[1]]
		lag.max <- min( nrow(mcmcobj) , lag.max )
		# index vector
		a1 <- attr(mcmcobj,"mcpar")
		iterindex <- seq(a1[1] , a1[2] )
		smcmcobj <- object$summary.mcmcobj
		VV <- ncol(mcmcobj)		
		ci.quant <- -qnorm( (1-conflevel)/2 )
		par( mfrow=c(2,2))
		for (vv in 1:VV){
#	vv <- 15
			x.vv <- as.vector( mcmcobj[,vv] )
			parm.vv <- colnames(mcmcobj)[vv]
			sparm.vv <- smcmcobj[ smcmcobj$parameter == parm.vv , ]
			#***
			# traceplot
			plot( iterindex , x.vv , type="l" ,  main= paste0( "Traceplot of " , parm.vv ) ,
				xlab="Iterations" , ylab="" , ... )
			x1 <- as.numeric( x.vv )
			xmin <- min(x1)
			xmax <- max(x1)
			l1 <- loess( x1 ~ iterindex ) 
			lines( iterindex ,l1$fitted  , col= col.smooth , lwd= lwd.smooth )
			#***
			# density estimate
			plot( density( x.vv ) , main= paste0( "Density of " , parm.vv ) )
			
			c1 <- quantile( x1 , ( 1 - conflevel  ) / 2 )			
			c2 <- quantile( x1 , 1 - ( 1 - conflevel  ) / 2 )			
#			lines( sparm.vv$Mean + c(-1,1)*ci.quant * sparm.vv$SD , c(0,0) , col=col.ci , lwd=3 )
			lines( c(c1,c2) , c(0,0) , col=col.ci , lwd=3 )
			points( sparm.vv$Mean , 0 , pch=17 , col= col.ci , cex=1.5)
			#***
			# autocorrelation function
			acf( x.vv , lag.max=lag.max ,
				main= paste0( "Autocorrelation of " , parm.vv ) )
			#***
			# numerical summary
			plot( c(0,1) , c(0,1) , axes=FALSE , xlab="" , ylab="", 
					main= paste0( "Summary of " , parm.vv ) , type="n" , ...)
			x0 <- 0 ; y0 <- 0
			heights.summ = c( .05 ,  .20 , .35 ,  .5 , .65 , .8 , .95)
			text( x0 + .0015 , y0 + heights.summ[7] , "Posterior Mean =" , cex= cex.summ  , pos=4)
			text( x0 + .5 , y0 + heights.summ[7] , 
				paste0( format.numb( x = mean( x1 ) , digits = round.summ)  )  , pos=4 )
			hvv <- heights.summ[6]	
			text( x0 + .0015 , y0 + hvv , "Posterior Mode =" , cex= cex.summ  , pos=4)				
			text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$MAP , digits = round.summ)  )  , pos=4 )			
				
			text( x0 + .0015 , y0 + heights.summ[5] , "Posterior SD   =" , cex= cex.summ  , pos=4)
			text( x0 + .5 , y0 + heights.summ[5] , 
				paste0( format.numb( x = sd( x1 ) , digits = round.summ)  )  , pos=4 )

			hvv <- heights.summ[4]
			text( x0 + .0015 , y0 + hvv , 
							paste( round(100*conflevel ) , "% Credibility Interval = " ,sep="") ,
							cex= cex.summ , pos=4 )

			hvv <- heights.summ[3]
				ci.lower <- format.numb( quantile( x1 , ( 1 - conflevel  ) / 2 ) , digits = round.summ )
				ci.upper <- format.numb( quantile( x1 , 1- ( 1 - conflevel  ) / 2 ) , digits = round.summ )            
			text( x0 + .25 , y0 + hvv , 
							paste( "[" , ci.lower ,    "," , ci.upper , "]" ,  sep="") ,
							cex= cex.summ  , pos=4)
			hvv <- heights.summ[2]
			text( x0 + .0015 , y0 + hvv , "Rhat =" , cex= cex.summ  , pos=4)
			text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$Rhat , digits = 2)  )  , pos=4 )
			hvv <- heights.summ[1]
			text( x0 + .0015 , y0 + hvv , "PercSEratio =" , cex= cex.summ  , pos=4)
			text( x0 + .5 , y0 + hvv , 
				paste0( format.numb( x = sparm.vv$PercSEratio , digits = 1)  )  , pos=4 )
			par(ask=ask)				
					}				
			par(mfrow=c(1,1))
			 }

                    }
#######################################################
# format numbers
format.numb <- function( x , digits ){ 
    a1 <- round( x , digits ) + 10^{-(digits +1 ) } 
    a1 <- substring( a1 , 1 , nchar(a1) - 1 )
    return(a1)
    }
#######################################################