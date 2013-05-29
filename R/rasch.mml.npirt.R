


# 1.01  2012-11-yy	O collect all nonparametric IRT functions into this file



# 1.0x  2012-11-yy
################################################################




#########################################################################################				 
.mstep.mml.npirt <- function( pjk , r.jk , n.jk , theta.k , npformula , npmodel, G , I){				 
					rjk0 <- r.jk[,,1]
					njk0 <- n.jk[,,1]
					if (G > 1){
						for (gg in 2:G){
							rjk0 <- rjk0 + r.jk[,,gg]
							njk0 <- njk0 + n.jk[,,gg]						
									}
								}
									
		if ( is.null( npformula ) ){ pjk <- t(rjk0 / njk0 )  }
				else {
						cat("ICC estimation |")	
						prbar <- floor( 10 * ( 1:I )	/ I )
						prbar <- c( 1 , diff(prbar))
						for (ii in 1:I){		
								#ii <- 3
							dfr1 <- data.frame( "theta" = theta.k , "y" = 1 , "wgt" = rjk0[ii,] )
							dfr0 <- data.frame( "theta" = theta.k , "y" = 0 , "wgt" = njk0[ii,] - rjk0[ii,] )
							dafr <- data.frame( rbind( dfr0 , dfr1 ) )
							if ( prbar[ii] == 1){ cat("~"); flush.console() }
							theta <- dafr$theta.k
							wgt <- dafr$wgt
							y <- dafr$y
							ICC_ <- model.matrix( npformula[[ii]] , dafr )						
							npmodel[[ii]] <- glm( y ~ 0 + ICC_ , weights = wgt , family="binomial" )						
#							npmodel[[ii]] <- glm( npformula[[ii]] , 
#										data = dafr , weights = dafr$wgt , family="binomial")					
							pjk[,ii] <- fitted( npmodel[[ii]] )[ seq( 1 , length(theta.k) ) + length(theta.k) ]
								}
							cat("\n")		
						}
					res <- list( "pjk" = pjk , "npmodel"=npmodel)				
							}
##########################################################################################

