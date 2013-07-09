

#------------------------------------------------------------------------
# --------------------------------------------------------------------#
# Semiparametric Maximum Likelihood Estimation in the Raschtype
#  Model                                                              
# item discrimination and guessing parameter can be fixed
rasch.mml2 <- function( dat , theta.k = seq(-6,6,len=21) , group = NULL , weights = NULL , 
						constraints = NULL , 
                        glob.conv = 10^(-5) , parm.conv = 10^(-4) , mitermax = 1 , 
                        mmliter = 1000 , progress = TRUE ,  
                        fixed.a = rep(1,ncol(dat)) , 
                        fixed.c = rep(0,ncol(dat)) , 
                        fixed.d = rep(1,ncol(dat)) , 
						fixed.K = rep(3,ncol(dat)) , 
						b.init = NULL , 
						est.a = NULL , 
						est.b = NULL , 
						est.c = NULL , est.d = NULL , 
						min.c = 0 , max.c = 1 , 
						min.d = 0 , max.d = 1 ,
						est.K = NULL , min.K = 1 , max.K = 20 , 
                        pid = 1:( nrow(dat) ) , trait.weights = NULL ,
                        center.trait = TRUE , alpha1 = 0 , alpha2 = 0 ,
                        est.alpha = FALSE , equal.alpha = FALSE , 
                        designmatrix = NULL , alpha.conv = parm.conv , 
						numdiff.parm = 0.001 , numdiff.alpha.parm= numdiff.parm , 
						distribution.trait = "normal" , 
						Qmatrix = NULL , 
						variance.fixed = NULL , 
						mu.fixed = cbind( seq(1,ncol(Qmatrix)) , rep(0,ncol(Qmatrix)) ) ,
						irtmodel = "raschtype" , npformula = NULL ,  
						use.freqpatt = is.null(group) , 
						... ){
    #******************************************************************************************##
    # INPUT:                                                                                ***##
    # dat           ... data frame with item responses
    # theta.k       ... grid of theta values where the trait density is evaluated
    # group         ... vector of group entries (numbered from 1 to G)
    # weights       ... sample weights or absolute frequency of item response patterns
    # glob.conv     ... global convergence criterion
    # conv1         ... convergence of parameters
    # mitermax      ... maximum number of iterations within the M step
    # WLE           ... should WLEs be estimated?
    # mmliter       ... maximum number of iterations during Maximum Likelihood Estimation
    # progress      ... should the estimation progress being displayed?
    # constraints   ... matrix with two columns: first column is item label (name or number), 
	#						second column is fixed item difficulty value (b parameter)
    # fixed.a       ... vector of fixed item discriminations
    # fixed.c       ... vector of fixed guessing parameters
    # b.init        ... initial estimates of item difficulties
	# est.a			... estimation of discrimination parameter
	# est.c, est.d	... estimated groups of parameters for lower and upper asymptote
	# max.c, min.d  ... maximal c and minimal d parameters to be estimated
    # pid           ... labels for subject ID's (-> pid ... person IDs)
    # trait.weights ... vector of fixed trait distribution
    # center.trait  ... set the mean of trait distribution equal to zero? (default = TRUE)
    # nplausible    ... number of plausible values to be drawn for item fit estimation
    # est.alpha     ... should alpha parameters be estimated?
    # equal.alpha   ... should equal alpha's be estimated?
    # designmatrix  ... Q matrix of item parameter restrictions    
    # numdiff.parm  ... step parameter for numerical differentiation
    # normal.trait  ... normal distribution assumption of the trait
	# ramsay.qm		... estimate quotient model of ramsay?
    #******************************************************************************************##
    #****
    # specifications
    conv1 <- parm.conv ; nplausible = 5 
	dat <- as.matrix(dat)
# a0 <- Sys.time()	
	# models
	npirt <- ramsay.qm <- FALSE
	I <- ncol(dat)
	if ( irtmodel == "ramsay.qm" ){ 
		ramsay.qm <- TRUE
		kG <- NULL
		}
	if ( irtmodel == "npirt" ){ 
			npirt <- TRUE 
			I <- ncol(dat)
			if ( ! is.null(npformula) ){						
				if ( length( npformula) == 1 ){
					npformula <- rep( npformula , I )
								}
				npformula0 <- npformula
				npformula <- list( 1:I) 							
				for (ii in 1:I){ 
					npformula[[ii]] <- as.formula( npformula0[ii] ) 
								}
								}
				npmodel <- list(1:I)
				}	
	# multidimensional model
	D <- 1
	if ( ! is.null( Qmatrix ) ){
		D <- ncol(Qmatrix)
		if ( D==2){ theta.k <- expand.grid( theta.k , theta.k ) }
		if ( D==3){ theta.k <- expand.grid( theta.k , theta.k , theta.k) }
		if ( D==4){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k) }
		if ( D==5){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k, theta.k) }
		if ( D==6){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k, theta.k, theta.k) }
		if ( D==7){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k, theta.k, theta.k, theta.k) }
		if ( D==8){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k, theta.k, theta.k, theta.k, theta.k) }
		if ( D==9){ theta.k <- expand.grid( theta.k , theta.k , theta.k, theta.k, theta.k, theta.k, theta.k, theta.k, theta.k) }		
#		Qmatrix <- sapply( 1:max(item.dim) , FUN = function(dd){
#						1*(item.dim == dd) } )		
					}
	if (D > 1){
		if ( is.null(variance.fixed) & ( sum(est.a) > 0) ){ 
			variance.fixed <- as.matrix( cbind( 1:D , 1:D , 1 ) )
			library(mvtnorm)
				}
			}				
	Sigma.cov <- diag(D)			
	mu <- rep(0,D)
# cat("114") ; a1 <- Sys.time(); print(a1-a0) ; a0 <- a1		
#	ramsay.qm <- FALSE
    if ( ! ramsay.qm) { raschtype <- TRUE  	}
	if (ramsay.qm | npirt ){
		raschtype <- FALSE
		# no alpha, a, c or d parameters can be estimated
		est.alpha <- FALSE	
		est.a <- est.c <- est.d <- NULL	
		pow.qm <- 1	# This parameter is ignored in analyses
					}
    # computation time
    s1 <- Sys.time()
    if (est.alpha){ 
            if (is.null(alpha1) ){ alpha1 <- 0 }
            if (is.null(alpha2) ){ alpha2 <- 0 }
                }
    #**************************************************************************************
    # some data checks
	ag1 <- NULL
	if( max( colMeans( is.na( dat ) ) ) == 1 ){
		stop("Remove items which have no observations!")
					}
    if ( ! is.null(group) ){ 
            t1 <- table(sort(group) )
			group.orig <- group
			group <- match( group.orig , sort(unique( group.orig)) )	
			ag1 <- aggregate( group , list( group.orig) , mean )
			colnames(ag1) <- c("group" , "groupindex" )
                            }
    # center trait: if there exists constraints, then do not center
	
	if ( is.null( colnames(dat) ) ){ 
			colnames(dat) <- paste( "I" , 1:ncol(dat) , sep="")
						}
	  if ( ! is.null( constraints ) ){ 
            center.trait <- F 
          if( ! is.numeric( constraints[,1] ) ){
            constraints[,1] <- match( paste(constraints[,1]) , colnames(dat) )
                                       }           
             constraints <- na.omit(constraints)
             constraints <- constraints[ constraints[,1] <= ncol(dat) , ]
                                }
    if ( ! is.null( designmatrix) ){
            if ( ncol(dat) != nrow(designmatrix) ){ 
                    stop( "Row dimension of designmatrix should be equal to number of items")
                                }
                        }
	# est.b parameters
	if (! is.null(est.b) ){
#		bG <- unique( est.b ) 
		bG <- unique( setdiff( est.b ,0 )) 
		if ( is.null( b.init) ){ b <- rep(0 , I ) }

		designmatrix <- matrix( 0 , ncol(dat) , length(bG) )	
		for (bb in bG){
			# bb <- bG[1]
#			designmatrix[ which( est.b == bb ) , bb ] <- 1
			designmatrix[ which( est.b == bb ) , match(bb,bG) ] <- 1
						}
				}
	# set starting values for estimated c and d parameters
	if ( sum(est.c) > 0 ){	fixed.c[ est.c > 0 ] <- .10 }
	if ( sum(est.d) > 0 ){	fixed.d[ est.d > 0 ] <- .95 }


    #****************************************************************************************
     WLE <- FALSE 
     pure.rasch <- -9   # this parameter is only included for historical reasons of this program.
    # specify weights
    if ( is.null(weights) ){ weights <- rep( 1 , nrow(dat) ) }
    # display
    if ( progress & ( npirt ) ){
        	cat("------------------------------------------------------------\n")
        cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
        cat("Nonparametric IRT Model (Rossi, Wang & Ramsay, 2002) \n") 
        	cat("------------------------------------------------------------\n")
        flush.console()
      }
    if ( progress & ( ramsay.qm ) ){
        	cat("------------------------------------------------------------\n")
        cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
        cat("Quotient Model (Ramsay, 1989) \n") 
#        if (normal.trait){ cat("Normal trait distribution \n") } else { cat("Nonparametric trait distribution \n") }
#		if (ramsay.qm){ cat("Log Normal Distribution of Theta with Power of" , pow.qm , "\n") }
        	cat("------------------------------------------------------------\n")
        flush.console()
      }
    if ( progress & (raschtype) ){
       	cat("------------------------------------------------------------\n")
        cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
        if ( est.alpha ){ 
            cat(paste( "Raschtype Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
                    } else {
            cat(paste( "Raschtype Model with generalized logistic link function: alpha1=",alpha1 , " , alpha2=", alpha2 , " \n") )
                        }
		if ( sum(est.c) > 0){ cat(paste( "Estimated guessing parameter groups \n") )}  ## estimated guessing parameters
		if ( sum(est.d) > 0){ cat(paste( "Estimated slipping parameter groups \n") )}  ## estimated slipping parameters
        	cat("------------------------------------------------------------\n")
        flush.console()
      }
    # revise guessing parameter (if necessary)
    if ( !is.null(fixed.c) ){
        # calculate itemmeans
        itemmean <- colMeans( dat ,  na.rm = T )
        if ( any( itemmean < fixed.c) ){  
                cat ( "revise fixed guessing estimates\n")
                fixed.c[ itemmean < fixed.c] <- 0
                }
            }
        # data preparations 
		if ( ! is.null( group ) ){ use.freqpatt <- FALSE }
        dp <- .data.prep( dat , weights = weights , use.freqpatt = use.freqpatt)
        dat1 <- dp$dat1
        dat2 <- dp$dat2
        dat2.resp <- dp$dat2.resp
        freq.patt <- dp$freq.patt
        n <- dp$n
        I <- dp$I
        # probability weights at theta.k
        if (D==1){
			pi.k <- dnorm( theta.k ) 
			pi.k <- pi.k / sum( pi.k )
				}
		if (D > 1){
			pi.k <- dmvnorm( theta.k , mean = rep(0,D) , sigma = Sigma.cov )	
			pi.k <- pi.k / sum( pi.k )
				}	
		G <- 1
		pi.k <- matrix( pi.k , nrow=length(pi.k) , ncol=G )
        # group calculations
        if ( !is.null( group )){ 
            G <- length( unique( group ) )
            pi.k0 <- pi.k
            pi.k <- matrix( 0 , nrow=length(pi.k0) , ncol=G)
            for (gg in 1:G){
                pi.k[,gg] <- pi.k0
                              }
                    }
		sd.trait <- mean.trait <- rep(0,G)					
        # initial estimates for item difficulties
        if ( is.null(b.init) & is.null(est.b) ){   
			b <- - qlogis( colMeans( dat , na.rm=T ) )
				if ( FALSE ){ 			
#				if ( ramsay.qm ){ 			
						b <-   - log( ( fixed.K * colSums( dat , na.rm=T ) ) / 
									( colSums( 1 - dat , na.rm=T ) ) ) 
								}
						} 
        if ( (!is.null(b.init) ) & is.null(est.b) ){   						
					b <- b.init 
					}
		if ( G == 1 ){ 	group <- rep(1, nrow(dat1)) }					
		# missing data indicators
		ind.ii.list <- list(1:I)
		for (ii in 1:I){
			ind.ii.list[[ii]] <- which( dat2.resp[,ii] == 1 )
						}
        mean.trait <- rep(0,G)					
        sd.trait <- rep(1,G)							
        # initial iteration index
        iter <- 0
        par.change <- dev.change <- 3
        dev <- 99
		apmax <- 0
        maxalphachange <- 1
		  # display
		  disp <- "...........................................................\n"	

#		old_increment.d <- old_increment.c <- rep( .2 , I )
	if( sum( est.d ) > 0 ){ 
		old_increment.d <- rep( .2 , length( unique( est.d[ est.d > 0 ] ) ) )
			}
	if( sum( est.c ) > 0 ){ 
		old_increment.c <- rep( .2 , length( unique( est.c[ est.c > 0 ] ) ) )
			}
	old_increment_b <- rep( 2 , I )

	# initialize standard errors
	se.alpha <- se.K <- se.b <- se.a <- se.c <- se.d <- NULL
	
	# Ramsay QM
#	if ( irtmodel == "ramsay.qm" ){ normal.trait <- TRUE }

#################################################	
        #--------------------------------#
        # MML Iteration Algorithm        #
        while ( ( dev.change > glob.conv | par.change > conv1 | maxalphachange > alpha.conv ) & iter < mmliter ){
		if (progress){ 
		  cat(disp)	
		  cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )
          flush.console()
					}
#zz0 <- Sys.time()					
				b0 <- b
                dev0 <- dev    		
                # perform E Step

               if ( ramsay.qm ){ 
					e1 <- .e.step.ramsay( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
									fixed.K , group , pow.qm = pow.qm , ind.ii.list )
                                } 
				if (raschtype & D==1){ 
                    e1 <- .e.step.raschtype( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
                                    fixed.a , fixed.c , fixed.d ,  alpha1 , alpha2 , group )        
									 }
				if (raschtype & D>1){ 
                    e1 <- .e.step.raschtype.mirt( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
                                    fixed.a , fixed.c , fixed.d ,  alpha1 , alpha2 , group ,
									 mu , Sigma.cov , Qmatrix)        									
									 }									 								 
									 
				if (npirt ){
					if (iter == 0){	
							pjk <- plogis( outer( theta.k , b , "-" ) ) 
								}
					e1 <- .e.step.ramsay( dat1 , dat2 , dat2.resp , theta.k , pi.k , I , n , b ,
									fixed.K , group , pow.qm = pow.qm , ind.ii.list ,
									pjk=pjk )				
							}
                n.k <- e1$n.k
                n.jk <- e1$n.jk
                r.jk <- e1$r.jk
                pjk <- e1$pjk
                f.qk.yi <- e1$f.qk.yi
				f.yi.qk <- e1$f.yi.qk
                dev <- -2*e1$ll		
#cat("e step") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
		
                # perform M Step
				#****
				# Ramsay QM
                 if ( ramsay.qm ){ 
                    m1 <- .m.step.ramsay( theta.k , b , n.k , n , n.jk , r.jk , I , 
							conv1 , constraints , 
							mitermax , pure.rasch ,  trait.weights , fixed.K , 
							designmatrix = designmatrix , group = group ,      
							numdiff.parm=numdiff.parm , pow.qm = pow.qm )
					se.b <- m1$se.b
											} 
				 # generalized Rasch type model
				 if (raschtype ){					 
                    m1 <- .m.step.raschtype( theta.k , b , n.k , n , n.jk , r.jk , pi.k , 
							I , conv1 , constraints , mitermax , pure.rasch ,        
							trait.weights , fixed.a , fixed.c , fixed.d ,  alpha1 , 
							alpha2 , designmatrix = designmatrix ,
							group = group , numdiff.parm=numdiff.parm ,
							Qmatrix = Qmatrix , old_increment=old_increment_b , 
							est.b=est.b )
					se.b <- m1$se.b							
                                }
#cat("m step raschtype") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	

								
				# nonparametric IRT model
				 if (npirt ){
				    pjk0 <- pjk
					res <- .mstep.mml.npirt( pjk , r.jk , n.jk , theta.k , 
									npformula , npmodel, G , I)					
					pjk <- res$pjk
					npmodel <- res$npmodel												
#					pjk <- t( rjk0 / njk0 )
					apmax <- max( pi.k[,1]*abs( pjk - pjk0)/.40 )
					m1 <- list( "b"=b , "G" = G , 
						"pi.k" = pi.k , "center" = FALSE )
						}

		# update mean and covariance in multidimensional models
		if ( D > 1){
			theta.k <- as.matrix(theta.k)
#			delta.theta <- (theta.k[2,1] - theta.k[1,1])^D
			delta.theta <- 1
				hwt <- e1$f.qk.yi		
				hwt <- hwt / rowSums(hwt)
				thetabar <- hwt%*%theta.k
				# calculation of mu
				mu <- colSums( thetabar * dat1$Freq ) / sum( dat1$Freq )

				if ( ! is.null(mu.fixed ) ){
				  if (is.matrix(mu.fixed) ){	
				    mu0 <- mu
					mu[ mu.fixed[,1] ] <- mu.fixed[,2]
					if ( sum( as.vector(mu.fixed[1,1:2]) - c(1,0))==0 ){ 
						mu[-1] <- -mu0[1] + mu[-1]
											}
										}
#				  if ( mu.fixed == "center"){
#						mu <- mu - mean(mu)
#										}
									}				
				# calculation of the covariance matrix
				theta.k.adj <- theta.k - matrix( mu , nrow=nrow(theta.k) , 
									ncol=ncol(theta.k) , byrow=TRUE)
				for (dd1 in 1:D){
					for (dd2 in dd1:D){
						tk <- theta.k.adj[,dd1]*theta.k.adj[,dd2]
						h1 <- dat1$Freq * ( hwt %*% tk ) * delta.theta
						Sigma.cov[dd1,dd2] <- sum( h1 ) / sum( dat1$Freq )
						if (dd1 < dd2 ){ Sigma.cov[dd2,dd1] <- Sigma.cov[dd1,dd2] }
										}
									}
				if ( ! is.null(variance.fixed ) ){
					Sigma.cov[ variance.fixed[,1:2] ] <- variance.fixed[,3]
					Sigma.cov[ variance.fixed[,c(2,1)] ] <- variance.fixed[,3]		
									}
				diag(Sigma.cov) <- diag(Sigma.cov) + 10^(-10)
				pi.k <- matrix( dmvnorm( theta.k , mean = mu , sigma = Sigma.cov )	, ncol=1 )		
				m1$pi.k <- pi.k <- pi.k / sum( pi.k )	

							}
			# end MIRT
			#*****
                b <- m1$b
                # distribution
                G <- m1$G
                pi.k <- m1$pi.k	
				if (!is.null( trait.weights) ){	
						pi.k <- matrix( trait.weights , ncol=1 )
#						print("bin ich hier??")
							}	
				#****************************************************
                # latent ability distribution
                if (distribution.trait=="normal" & D == 1){ 
					delta.theta <- 1
#					delta.theta <- theta.k[2] - theta.k[1]
#                    sd.trait <- mean.trait <- rep(0,G)			
				  h <- .0005
                    for (gg in 1:G){ 
						pi.k0 <- pi.k
						f.yi.qk.gg <- e1$f.yi.qk[group==gg,]
						dat1.gg <- dat1[group==gg,2]
						X1 <- rep(1,nrow(f.yi.qk.gg) )
					if ( gg > 1 | ( ! center.trait ) ){
						#*********************************	
						# mean estimation	
						d.change <- .est.mean( dat1.gg , f.yi.qk.gg , X1 , pi.k , pi.k0 , gg ,
								mean.trait , sd.trait , theta.k , h)
						mean.trait[gg] <- mean.trait[gg] + d.change		
                        pi.k[,gg] <- dnorm( theta.k , mean = mean.trait[gg] , sd = sd.trait[gg] )
                        pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
													}
						if (center.trait){ mean.trait[1] <- 0 }			
						#*********************************	
						# SD estimation	
						if ( ( gg > 1 ) | ( sum(est.a) == 0 )  ){ 
							d.change <- .est.sd( dat1.gg , f.yi.qk.gg , X1 , pi.k , pi.k0 , gg ,
								mean.trait , sd.trait , theta.k , h )						
							sd.trait[gg] <- sd.trait[gg] + d.change	
											}
					if ( ( ! is.null(est.a) ) | ( irtmodel == "npirt" )  ){ 
							sd.trait[1] <- 1 
									}
                        pi.k[,gg] <- dnorm( theta.k , mean = mean.trait[gg] , sd = sd.trait[gg] )
                        pi.k[,gg] <- pi.k[,gg] / sum( pi.k[,gg] )
								}
						}  # end normal distribution
		#######################################
               if (distribution.trait!="normal" & D == 1){ 
			   for (gg in 1:G){ 
					pik1 <-	n.k[,gg] / sum(n.k[,gg] )
					pik1 <- pik1 + 10^(-10)
					lpik1 <- log( pik1 )
					tk <- theta.k
					if ( distribution.trait=="smooth2"){ 
							formula1 <- lpik1 ~ tk + I(tk^2)
											}				
					if ( distribution.trait=="smooth3"){ 
							formula1 <- lpik1 ~ tk + I(tk^2) + I(tk^3)
											}
					if ( distribution.trait=="smooth4"){ 
							formula1 <- lpik1 ~ tk + I(tk^2) + I(tk^3)+I(tk^4)
											}
					mod <- lm( formula1 , weights = pik1 )
					pik2 <- exp( fitted(mod))
					pi.k[,gg] <- pik2 / sum(pik2)
					if (center.trait & gg==1){
						mmm1 <- weighted.mean( theta.k , pik2 )	
						theta.k <- theta.k - mmm1
											}
					if ( ( ! is.null(est.a) ) | ( irtmodel == "npirt" )  ){ 
					   if (gg== 1){
						 sd1 <- sqrt( sum( theta.k^2 * pi.k[,1] ) - sum( theta.k * pi.k[,1] )^2 )
						 theta.k <- theta.k / sd1 
									}
									}
							}
 						}  # end non-normal distribution		
						
        ##############################
        # estimation of alpha, c and d parameters
        alpha.change <- 0
        maxalphachange <- 0
		a1a <- a1b <- 0
		a1K <- a1c <- 0
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of a parameters
		if ( sum( est.a ) > 0 & raschtype ){
                h <- numdiff.parm 
				fixed.a0 <- fixed.a
				# identify different a parameter groups
				aG <- setdiff(unique( est.a ) , 0 )
				# a estimation
				res <- .mml.raschtype.est.a( theta.k , b , fixed.a , fixed.c , fixed.d ,
					pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.a , Qmatrix)		
				fixed.a <- res$fixed.a
				se.a <- res$se.a					
				a1a <- max( abs( fixed.a - fixed.a0 ) )
						}	
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of c parameters
		if ( sum( est.c ) > 0 & raschtype ){
                h <- numdiff.parm 
				fixed.c0 <- fixed.c
				# identify different c parameter groups
				cG <- setdiff( unique( est.c ) , 0 )			
				res <- .mml.raschtype.est.c( theta.k , b , fixed.a , fixed.c , fixed.d ,
					pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.c ,
					min.c , max.c , iter , old_increment.c , Qmatrix)
				fixed.c <- res$fixed.c
				se.c <- res$se.c					
                a1b <- max( abs( fixed.c - fixed.c0 ) )
						}
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of d parameters
		if ( sum( est.d ) > 0 & raschtype ){
                h <- numdiff.parm 
				fixed.d0 <- fixed.d
				# identify different c parameter groups
				dG <- setdiff( unique( est.d ) , 0 )
				res <- .mml.raschtype.est.d( theta.k , b , fixed.a , fixed.c , fixed.d ,
					pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.d ,
					min.d , max.d , iter , old_increment.d,Qmatrix)
				fixed.d <- res$fixed.d
				se.d <- res$se.d					
                a1c <- max( abs( fixed.d - fixed.d0 ) )
						}
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# estimation of K parameters in Ramsay's quotient model
		if ( sum( est.K ) > 0 & ramsay.qm ){
                h <- numdiff.parm 
				fixed.K0 <- fixed.K
				# identify different c parameter groups
				kG <- setdiff( unique( est.K ) , 0 )

				res <- .mml.ramsay.est.K( theta.k , b , fixed.a , fixed.c , fixed.d ,
					fixed.K , pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.K ,
					min.K , max.K , iter , pow.qm )
				fixed.K <- res$fixed.K
				se.K <- res$se.K									
				# convergence is indicated in metric guess.K = 1 / ( fixed.K + 1 )
                a1K <- max( abs( 1/(1+fixed.K) - 1/(1+fixed.K0) ) )
						}	
		#***************************
		# estimation of alpha
        if ( est.alpha ){ 
                alpha1.old <- alpha1
                h <- numdiff.alpha.parm 
                fixed.c.M <- outer( rep(1,nrow(pjk)) , fixed.c )
                fixed.d.M <- outer( rep(1,nrow(pjk)) , fixed.d )            
                #****
                # alpha1
                pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a ,Qmatrix)
                    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
                pjk <- ( pjk + .000000005 ) / 1.00000001 
                pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
                pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1+h , alpha2 , fixed.a,Qmatrix)
                    pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
                pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
                pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
                pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1-h , alpha2 , fixed.a,Qmatrix)
                    pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
                pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
                pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
                # first order derivative
                # f(x+h) - f(x-h) = 2* f'(x) * h
                ll0 <- ll1 <- ll2 <- rep(0,G)
                for (gg in 1:G){ 
                    ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
                    ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
                    ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
                            }                
                ll0a1 <- ll0 <- sum(ll0)
                ll1a1 <- ll1 <- sum(ll1)
                ll2a1 <- ll2 <- sum(ll2)
                d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
#				d1 <- ( ll1 - ll0 ) / h
                # second order derivative
                # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
                d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
                # change in item difficulty
                alpha.change <- - d1 / d2
                alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
                alpha1 <- alpha1 + alpha.change
                a1 <- abs(alpha.change )
				se.alpha <- sqrt( 1 / abs(d2) )
                #****
                # alpha2
                pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a,Qmatrix)
                    pjk <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk
                pjk <- ( pjk + .000000005 ) / 1.00000001 
                pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
                pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 +h , fixed.a,Qmatrix )
                    pjk1 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk1
                pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
                pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
                pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 -h , fixed.a,Qmatrix)
                    pjk2 <- fixed.c.M + ( fixed.d.M - fixed.c.M) * pjk2
                pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
                pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
                # first order derivative
                # f(x+h) - f(x-h) = 2* f'(x) * h
                ll0 <- ll1 <- ll2 <- rep(0,G)
                for (gg in 1:G){ 
                    ll0[gg] <- sum( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
                    ll1[gg] <- sum( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
                    ll2[gg] <- sum( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) )
                            }                
                ll0 <- sum(ll0)
                ll1 <- sum(ll1)
                ll2 <- sum(ll2)
                d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
#				d1 <- ( ll1 - ll0 ) / h

                # second order derivative
                # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
                d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
                # change in item difficulty
                alpha.change <- - d1 / d2
                alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )
                alpha2 <- alpha2 + alpha.change
                a2 <- abs(alpha.change)
                maxalphachange <- max(a1, a2)
				se.alpha <- c( se.alpha , sqrt( 1 / abs(d2) ) )
                if (equal.alpha){
                        ll0 <- ll0a1 + ll0 
                        ll1 <- ll1a1 + ll1                 
                        ll2 <- ll2a1 + ll2                 
                        d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
                        d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2                
                        alpha.change <- - d1 / d2
                        alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )
                        alpha2 <- alpha1 <- alpha1.old + alpha.change
                        a2 <- abs(alpha.change)
                        maxalphachange <- max(a2)
						se.alpha <- sqrt( 1 / abs(d2) )
                            }
						}
# cat("distribution / rest") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1	
						
            ############################
                # possibly incorrect deviance calculated at the M step
                #                dev <- -2*m1$ll
                # iteration index
                dev.change <- abs( ( dev - dev0)/ dev0 )
                par.change <- max( c( abs(b - b0 ) , abs(alpha.change ) , a1a , a1b , a1c , a1K ,
								apmax) )
                # display convergence
                if ( progress  ){   
						cat( paste( "   Deviance = "  , 
                             round( dev , 4 ) , 
							 if (iter > 0 ){ " | Deviance change = " } else {""} ,
							 if( iter>0){round( - dev + dev0 , 6 )} else { ""}	,"\n",sep=""))
					if ( ! npirt ){							 
						cat( paste( "    Maximum b parameter change = " , 
                             round( max(abs(b - b0 )) , 6 ) ,  " \n"   )  )  
									}
                                    if ( est.alpha ){           
                                        cat( paste( "    alpha1=" , round(alpha1,3) , " | alpha2=", round( alpha2,3) , 
                                                    " | max alpha change " , round( maxalphachange ,7 ) , "\n" , sep=""))
												}
                                    if ( sum(est.a) > 0  ){           
                                        cat( paste( "    Maximum a parameter change = " , 
												paste( round(a1a ,6) , collapse=" " ) , "\n" , sep=""))
												}
									if ( sum(est.c) > 0  ){           
                                        cat( paste( "    Maximum c parameter change = " , 
												paste( round(a1b ,6) , collapse=" " ) , "\n" , sep=""))												
												}
                                    if ( sum(est.d) > 0  ){           
                                        cat( paste( "    Maximum d parameter change = " , 
												paste( round(a1c ,6) , collapse=" " ) , "\n" , sep=""))
												}
                                    if ( npirt  ){           
                                        cat( paste( "    Maximum weighted ICC change = " , 
												paste( round(apmax ,6) , collapse=" " ) , "\n" , sep=""))
												}
									if ( sum(est.K) > 0  ){           
                                        cat( paste( "    Maximum K parameter change = " , 
												paste( round(a1K ,6) , collapse=" " ) , "\n" , sep=""))
												}
								if ( D > 1 ){
								  cat("    Mean              | " )
								  cat( round(as.vector(mu),3))
								  cat("\n    Covariance Matrix | " )
								  cat( round(Sigma.cov[!upper.tri(Sigma.cov)],3))
								  cat("\n")
															}							
									flush.console() 
                                    }
                iter <- iter + 1
			
                     }
		####################################### end iterations #####################
		############################################################################

			if (npirt & ( ! is.null(npformula ) ) ){
				item <- NULL
				for (ii in 1:I){
					item.ii <- data.frame( "item" = colnames(dat)[ii] )
					smod.ii <-  summary(npmodel[[ii]])
					item.ii <- data.frame( cbind( item.ii , rownames(smod.ii$coef) , 
										smod.ii$coef[,1:2] ) )
					colnames(item.ii)[-1] <- c("par" , "est" , "se" )
					rownames(item.ii) <- NULL
					item <- rbind( item , item.ii )
							}			
						}
					  
        #**********************************************
        # standard error for item parameter
		# ...
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
#		ic$traitpars <- ic$itempars <- NA
        if ( distribution.trait=="normal"){
				ic[[ "np" ]] <-  ( G - 1 ) + ncol(dat) + ( G - 0 )
							}
        if ( distribution.trait=="smooth2"){
				ic[[ "np" ]] <-  ( G - 1 ) + ncol(dat) + ( G - 0 )
							}
        if ( distribution.trait=="smooth3"){
				ic[[ "np" ]] <- ( G - 1 ) + ncol(dat) + ( G - 0 ) + G
							}
        if ( distribution.trait=="smooth4"){
				ic[[ "np" ]]<-  ( G - 1 ) + ncol(dat) + ( G - 0 ) + 2*G
							}
#		ic$itempars <- ic$traitpars - ncol(dat)
        # subtract fixed constraints
        if ( ! is.null( constraints) ){ 
			ic$np <- ic$np - nrow(constraints) 
#			ic$itempars <- ic$itempars - nrow(constraints)
				}
        # subtract constraints due to designmatrix
        if ( ! is.null( designmatrix ) ){ 
				ic$np <- ic$np - ncol(dat) + ncol(designmatrix) 
#				ic$itempars <- ic$itempars - ncol(dat) + ncol(designmatrix) 			
						}
        # alpha estimation
        ic$np <- ic$np + est.alpha * 2 - equal.alpha *1
#        ic$itempars <- ic$itempars + est.alpha * 2 - equal.alpha *1		
		# guessing, slipping and discrimination parameter estimation
		if ( sum(est.c) > 0 ){ 
					ic$np <- ic$np + length(cG) 
						}
		if ( sum(est.d) > 0 ){ 
				ic$np <- ic$np + length(dG) 
							}		
		if ( sum(est.a) > 0 ){ 
				ic$np <- ic$np + length(aG) - 1 
						}		
		if ( sum(est.K) > 0 ){ 
				ic$np <- ic$np + length(kG)
							}
		# parameters for multiple dimensions
		if (D>1){ 
			# mean vector
			MM <- nrow(mu.fixed )
#			if ( mu.fixed == "center" ){ MM <- 1 }
			if ( is.null(mu.fixed) ){ MM <- 0 }
			ic$np <- ic$np + length(mu) - MM
			# covariance matrix
			ic$np <- ic$np - 1*(sum(est.a)==0) + D*(D+1)/2		# SD's
			if ( ! is.null(variance.fixed)){ ic$np <- ic$np - nrow( variance.fixed ) }
				}
		
		# item parameter for nonparametric models
		if (npirt & ( ! is.null(npformula ) ) ){
				ic$np <- nrow(item) }
		if (npirt & (  is.null(npformula ) ) ){
				ic$np <- prod( dim(pjk)) }
    	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC (conistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
        # item statistics
		if ( npirt & ( ! is.null(npformula ) ) ){ item0 <- item }
        item <- data.frame( "item" = colnames(dat) , "N" = colSums( weights*(1 - is.na(dat)) ) , 
                            "p" = colSums( weights*dat   , na.rm=T) / colSums( weights*(1 - is.na(dat)) ) , 
                            "b" = b  )
		if ( ! is.null( constraints) ){   
				est.b <- 1:I
				est.b[ constraints[,1] ] <- 0
				item$est.b <- est.b
					}							
		if ( npirt & ( ! is.null(npformula ) ) ){ item <- merge( x = item[,1:3] , y = item0 , by = "item" ) }							
		if ( ! npirt ){		
			if (is.null(est.b)){ item$est.b = seq(1,I) } else { item$est.b <- est.b }
			# fixed parameters
			item$a <- fixed.a
			if ( ! is.null( est.a) ){  item$est.a <- est.a } else { item$est.a <- rep(0,I) }
			# include threshold
			item$thresh <- item$a*item$b		
			# guessing parameter
			item$c <- fixed.c
			if ( ! is.null( est.c) ){ item$est.c <- est.c } else { item$est.c <- rep(0,I) }
			item$d <- fixed.d
			if ( ! is.null( est.d) ){  item$est.d <- est.d } else { item$est.d <- rep(0,I) }
			if (m1$center){  if ( is.null(constraints) ){ #    item[I,4] <- NA 
																		} 
						else { item[ constraints[,1] ,4] <- NA } }
			rownames(item) <- colnames(dat)
				}
        # latent ability distribution
        skewness.trait <- sd.trait <- mean.trait <- rep(0,G)
		if ( D==1){
			for (gg in 1:G){ 
				mean.trait[gg] <- weighted.mean( theta.k , pi.k[,gg] ) 
				sd.trait[gg] <- sqrt( weighted.mean( ( theta.k - mean.trait[gg] )^2 , pi.k[,gg] ) ) 
				skewness.trait[gg] <- sum( ( theta.k - mean.trait[gg]  )^3 * pi.k[,gg] ) / sd.trait[gg]^3
				if (gg == 1 & npirt ){ sd.trait[gg] <- 1 }
						}
					}
			# center trait distribution
#			if ( center.trait & G < 1 ){
#					theta.k <- theta.k - mean.trait
#					b <- b - mean.trait
#					item$itemdiff <- b
#					mean.trait <- 0
#							}      
			trait.distr <- data.frame( "theta.k" = theta.k , "pi.k" = pi.k )
        # item response pattern
		if ( D==1 ){
			ability.est <- data.frame( dat1 , theta.k[ whichrowMaxs( f.qk.yi )$arg ] )
			colnames(ability.est) <- c("pattern" , "AbsFreq" , "mean" , "MAP" )
					}
		if (D>1){					
			ability.est <- data.frame( dat1 , theta.k[ whichrowMaxs( f.qk.yi )$arg ,] )
			colnames(ability.est) <- c("pattern" , "AbsFreq" , "mean" , 
					paste("MAP.Dim",1:D,sep="") )		
				}
		if (D==1){
			ability.est$EAP <- rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , theta.k  )  )
			ability.est$SE.EAP <- sqrt( rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , 
							theta.k^2  )  ) - ability.est$EAP^2 )
				 }
		if (D>1){
		   for (dd in 1:D){
			ability.est[ , paste("EAP.Dim",dd,sep="")] <- 
					rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , theta.k[,dd]  )  )
			ability.est[ , paste("SE.EAP.Dim",dd,sep="")] <- 
				sqrt( rowSums( f.qk.yi * outer( rep(1,nrow(ability.est)) , theta.k[,dd]^2  )  ) - 
							ability.est[,paste("EAP.Dim",dd,sep="")]^2 )
								}
						}				 
        # posterior distribution
        rownames(f.qk.yi) <- dat1[,1]		
        # merging ability estimates
#        if ( ! is.null(group)){  
		if ( G > 1 ){
                    ability.est2 <- cbind( freq.patt , ability.est[,-1] ) 
                            } else {
                ability.est2 <- merge( freq.patt , ability.est , 1 , 1 )
                            }
        ability.est2 <- ability.est2[ order(ability.est2$index) , -c(3:5) ]   
			
		# EAP reliability estimate
		reliability <- NULL
		if (D==1){ 		
			reliability$eap.reliability <- 
					1 - mean(ability.est2$SE.EAP^2) / ( mean(ability.est2$SE.EAP^2) + var(ability.est2$EAP) )
						}
		if (D>1){ 		
			r1 <- rep(0,D)
			for (dd in 1:D){
			r1[dd] <- 1 - mean(ability.est2[,paste("SE.EAP.Dim",dd,sep="")]^2) / 
								( mean(ability.est2[,paste("SE.EAP.Dim",dd,sep="")]^2) + 
									var(ability.est2[,paste("EAP.Dim",dd,sep="")]) )
							}
		  if ( is.null( colnames(Qmatrix) ) ){
			dimnamesPars <- paste( "Dim",1:D , sep="")
							} else { dimnamesPars <- colnames(Qmatrix) }							
          names(r1) <- dimnamesPars
          reliability$eap.reliability <- r1
		  names(mu) <- dimnamesPars
		  rownames(Sigma.cov) <- colnames(Sigma.cov) <- dimnamesPars
						}
        # include person ID
        ability.est2$pid <- pid
		# match ability patterns

		ind1 <- match( ability.est2$freq.patt , ability.est$pattern  )
		ability.est <- ability.est[ ind1  , ]
		f.qk.yi <- f.qk.yi[ind1,]
		f.yi.qk <- f.yi.qk[ind1,]		
		
        # output fixed.a and fixed.c         
        if ( is.null(fixed.a ) & is.null(fixed.c) ){  fixed.a <- rep(1,I) ; fixed.c <- rep(0,I) }
        # include item discrimination
		if (D==1){
			i1 <- item$emp.discrim <- round( item.discrim( dat ,  ability.est2$MAP ) , 3 )
				}
		if (npirt){
			i1 <- data.frame( "item" = colnames(dat) , "emp.discrim" = i1 )
			item$emp.discrim <- NULL
			item <- merge( x = item , y = i1 , by = "item" )
				}
		if ( ! npirt ){
			item$alpha1 <- alpha1
			item$alpha2 <- alpha2
					}
		#---------------------------------------------------------
		# item summary Ramsay QM
		item2 <- NULL
		if ( ramsay.qm){
		if ( is.null(est.K) ){ est.K <- rep(0,I) }
			item2 <- data.frame( "item" = item$item  , "N" = item$N	 , "p" = item$p , 
							"K" = fixed.K , "est.K" = est.K , 
							"b" = exp(b) , "log_b" = b , "est.b" = item$est.b , 
							"guess.K" = 1/(fixed.K+1) , 
							"emp.discrim" = item$emp.discrim )
							}
        # result
        res <- list( "dat" = dat , "item" = item , "item2" = item2 , 
					"trait.distr" = trait.distr , "mean.trait" = mean.trait , 
					"sd.trait" = sd.trait , "skewness.trait" = skewness.trait , 
                    "deviance" = dev ,  "pjk" = pjk ,   "person" = ability.est2 , "pid" = pid , 
                    "ability.est.pattern" = ability.est , "f.qk.yi" =  f.qk.yi , "f.yi.qk" =  f.yi.qk,
                    "pure.rasch" = pure.rasch  , "fixed.a" = fixed.a , "fixed.c" = fixed.c , 
					"G" = G ,"alpha1"=alpha1 , "alpha2" = alpha2 , 
					"se.b" = se.b , "se.a" = se.a , "se.c" = se.c , "se.d" = se.d ,
					"se.alpha" = se.alpha , "se.K" = se.K , 
					"iter" = iter ,
					"reliability" = reliability , "ramsay.qm" = ramsay.qm ,
					"irtmodel" = irtmodel , "D" = D , "mu" = mu , 
					"Sigma.cov"=Sigma.cov , "theta.k" = theta.k , 
					"trait.weights" = trait.weights  # , "pi.k" = pi.k 
# collect results of npmodel
							) 
        class(res) <- "rasch.mml"
        res$ic <- ic
		res$est.c <- est.c
		res$groupindex <- ag1
        # computation time
        s2 <- Sys.time()
		res$s1 <- s1
		res$s2 <- s2
		res$Rfcttype <- "rasch.mml2"
        if (progress){ 
				cat("------------------------------------------------------------\n")
                cat("Start:" , paste( s1) , "\n")
                cat("End:" , paste(s2) , "\n")
                cat("Difference:" , print(s2 -s1), "\n")
                cat("------------------------------------------------------------\n")
                    }       
        #..................
        return( res )
        }
#---------------------------------------------------------------------------









#*******************************************************
# Summary for rasch.mml object                         *
##NS S3method(summary,rasch.mml)
summary.rasch.mml <- function( object , ... ){
    # object      ... object from rasch.mml                #
	
	npirt <- object$irtmodel == "npirt"	
	D <- object$D
	
	a5 <- 1*( npirt & ncol( object$item) == 5 )

	cat("------------------------------------------------------------\n")
		d1 <- packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Semiparametric Marginal Maximum Likelihood Estimation \n")
	if ( object$Rfcttype == "rasch.mml" ){ cat("Function 'rasch.mml' \n\n") }
	if ( object$Rfcttype == "rasch.mml2" ){ cat("Function 'rasch.mml2' \n\n") }	
	if ( ! object$ramsay.qm ){
		cat("Rasch Type Model with Fixed Discrimination, Guessing and Slipping Parameters \n") 
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		moments <- genlogis.moments( alpha1=object$alpha1 , alpha2=object$alpha2)
		cat("Moments: \n" ); print(round(moments,2)) ; cat("\n")
							}
	    if ( object$ramsay.qm ){  cat("Quotient Model (Ramsay, 1989) \n") 	}
	    if ( object$irtmodel == "npirt" ){  cat("Nonparametric IRT \n") 	}		
        flush.console()
	if ( sum(object$est.c) > 0){ cat(paste( "Estimated guessing parameter groups \n") )}  
				## estimated guessing parameters
    if ( object$G > 1 ){ cat("\nMultiple Group Estmation with",object$G , "Groups \n") 
			print(object$groupindex) ; cat("\n")
			}
	cat("------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , 
				round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , 
					round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , 
				round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

	#------------------------------------
	if ( object$D == 1){	
	if ( is.null( object$trait.weights) ){
		cat( "Trait Distribution (" , length(object$trait.distr[,1]) , " Knots )\n" , 
				  "Mean=" , round( object$mean.trait,3) , " SD=" , round( object$sd.trait , 3) ,
				  " Skewness=" , round( object$skewness.trait , 3) 
				  ) 
						}
	if ( ! is.null( object$trait.weights) ){
		M1 <- weighted.mean( object$theta.k ,object$trait.weights )
		S1 <- sqrt( weighted.mean( object$theta.k^2 ,object$trait.weights ) - M1^2	)
		cat( "Fixed Trait Distribution (" , length(object$trait.distr[,1]) , " Knots )\n" , 
				  "Mean=" , round( M1 ,3 ) , 
					" SD=" , round( S1 , 3) 
								)
						}						
						}
	if ( object$ramsay.qm ){ cat( "      Note: log theta distribution is parametrized!") }
	cat("\n")
	if ( D > 1){
		cat("Mean Vector\n") ; print( round( object$mu , 3 ) )
		cat("\nCovariance Matrix\n") ; print( round( object$Sigma.cov , 3 ) )	
		cat("\n")
		covmat <- object$Sigma.cov 
		covmat2 <- cov2cor( covmat )
		diag(covmat2) <- sqrt( diag(covmat) )
		cat("\nStandard Deviations / Correlation Matrix\n") ; print( round( covmat2 , 3 ) )	
		cat("\n")
				}
	
	if ( object$irtmodel != "npirt" ){	
		cat( "Item Difficulty Distribution (" , nrow(object$item) , " Items )\n" , 
				  "Mean=" , round( mean(object$item$b) ,3) , " SD=" , 
							round( sd(object$item$b) , 3) , "\n") 
							}
    cat( "Distribution of Items Administered (" , nrow(object$item) , " Items )\n" , 
              "Mean=" , round( mean(rowSums( 1 - is.na(object$dat) )) ,3) , " SD=" , 
                        round( sd(rowSums( 1 - is.na(object$dat) )) ,3) , "\n\n") 
	cat( "EAP Reliability: ") 
	cat(round( object$reliability$eap.reliability,3 ) )
	cat( "\n")
	if ( a5 == 0 ){
	cat("------------------------------------------------------------\n")
		cat("Item Parameter \n")
		if ( ! object$ramsay.qm ){ obji <- object$item } else { obji <- object$item2 }
		rvars <- seq( 2 , ncol(obji ) )
		ind <- which( colMeans( is.na( obji )) == 1 )
		roundvars <- setdiff( rvars , ind )
		if (npirt ){ 
		  roundvars <- c("p","est","se" )
				}		
		for (vv in roundvars ){ 
			obji[,vv] <- round( obji[,vv] , 3 ) 
					}
		rownames(obji) <- NULL
		print( obji )                
			}
                }
#*******************************************************



