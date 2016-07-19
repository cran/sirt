
###############################################
# Bradley-Terry model in sirt
btm <- function( data ,	ignore.ties = FALSE  , fix.eta=NULL , fix.delta=NULL ,
			fix.theta = NULL , 
			maxiter=100 , conv = .0001 , eps=.3){

		s1 <- base::Sys.time()
		CALL <- base::match.call()		
		admiss <- c(0,1,.5)
		est.delta <- TRUE
		est.eta <- TRUE
		delta <- -1
		eta <- 0
		center.theta <- TRUE
		
		if ( ! is.null( fix.eta ) ){
			eta <- fix.eta
			est.eta <- FALSE
						}
		
		
		if ( ignore.ties ){ 
			admiss <- admiss[1:2] 
			delta <- -99
			est.delta <- FALSE
				}
		
		data <- data[ ! is.na( data[,3] ) , ]
		data <- data[ data[,3] %in% admiss , ]		
		dat0 <- dat <- data
					
		# teams
		teams <- base::unique( base::sort( c( base::paste(dat[,1] ) , base::paste(dat[,2] ) ) ) )
		dat0[,1] <- base::match( base::paste( dat0[,1] ) , teams )
		dat0[,2] <- base::match( base::paste( dat0[,2] ) , teams )

		dfr <- base::data.frame( "res1" = 1 * ( dat0[,3] == 1 ),
					"res0" = 1 * ( dat0[,3] == 0 ) ,
					"resD" = 1 * ( dat0[,3] == 1/2 ) )            

		# calculate scores for each team
		TP <- base::length(teams)
		r1 <- base::rowsum( dfr , dat0[,1] )
		ind1 <- base::as.numeric(rownames(r1))
		r2 <- base::rowsum( dfr , dat0[,2] )
		ind2 <- base::as.numeric(rownames(r2))
		r2 <- r2[, c(2,1,3) ]
		r0 <- base::as.data.frame( base::matrix( 0 , nrow=TP , ncol=3 ))
		r0[ ind1 , ] <- r0[ ind1,] + r1
		r0[ ind2 , ] <- r0[ ind2,] + r2
		r3 <- r0
		colnames(r3) <- colnames(dfr)
		score <- r3[,1]*1 + r3[,3]*1/2
		maxscore <- base::rowSums(r3)
		# epsilon adjustment
		score <- eps + ( maxscore-2*eps)*score / maxscore
		propscore <- score / maxscore

		# initial ability for each team
		theta <- stats::qlogis( propscore )
		# eliminate individuals with extreme scores
		elim_persons <- FALSE
		if ( base::sum( propscore %in% c(0,1) ) > 0 ){
			elim_persons <- TRUE
			elim_persons_index <- base::which( propscore %in% c(0,1) )
			theta_elim <- theta[ elim_persons_index ]
			# define matrix with sets probabilities to zero for 
			# comparisons which are excluded
			indicator_elim <- dat0
			indicator_elim[,3] <- 0
			indicator_elim[ indicator_elim[,1] %in% elim_persons_index , 3 ] <- 1
			indicator_elim[ indicator_elim[,2] %in% elim_persons_index , 3 ] <- 1
										}

		
		some.fixed.theta <- FALSE
		if ( ! is.null( fix.theta) ){
			fix.theta.index <- base::match( names(fix.theta) , teams )
			if ( base::sum( is.na( fix.theta.index ) ) > 0 ){
				base::stop( paste0( "Cannot find all individuals with fixed values\n" ,
						"  in 'fixed.theta'\n") )
								}			
			some.fixed.theta <- TRUE
			center.theta <- FALSE
									}
		
		# number of dyads
		ND <-  nrow(dat0)
		max.change <- 1E5
		iter <- 0			
		incrfac <- .98
		maxincr <- 1
		se.delta <- NA
		se.eta <- NA
		
		while( ( iter < maxiter ) & ( max.change > conv) ) {

			theta0 <- theta
			delta0 <- delta
			eta0 <- eta
			
			M1 <- base::matrix(0 , nrow=ND , ncol=3)
			M1[,1] <- theta[ dat0[,1] ] + eta
			M1[,2] <- theta[ dat0[,2] ] 
			M1[,3] <- delta + ( theta[ dat0[,1] ] + theta[ dat0[,2] ] + eta ) / 2
			
			M1 <- base::exp(M1)
			M1 <- M1 / base::rowSums(M1)

			if ( elim_persons ){
				M1[ indicator_elim[,3] == 1 , ] <- 0
								}
			
			maxincr <- maxincr * incrfac
			
			#***********************************
			# derivatives with respect to delta	
			if ( est.delta ){			
				d1 <- base::sum(r1[,3]) - base::sum( M1[,3] )						
				d2 <- base::sum( M1[,3] * ( 1 - M1[,3] ) )
				incr <- d1 / d2
				incr <- base::ifelse( base::abs(incr) > maxincr , maxincr*base::sign(incr) , incr )
				delta <- delta + incr
				se.delta <- base::sqrt( 1 / d2 )			
							}
			#***********************************
			# derivatives with respect to eta		
			if ( est.eta ){
				d1 <- base::sum(r1[,1]+r1[,2]/2) - base::sum( M1[,1] + M1[,3]/2 )						
				d2 <- base::sum( M1[,1] * ( 1 - M1[,1] - M1[,3]/2 ) +
								M1[,3]/2 * ( 1/2 - M1[,1] - M1[,3]/2 ) )
				incr <-  d1 / d2
				incr <- base::ifelse( base::abs(incr) > maxincr , 
								maxincr*base::sign(incr) , incr )
				eta <- eta + incr
				se.eta <- base::sqrt( 1 / d2 )			
						}
			
			#******************
			# derivatives with respect to theta
			# first derivative
			fac1 <-  M1[,1] + M1[,3]/2
			deriv_theta_i1 <- fac1
			fac2 <-  M1[,2] + M1[,3]/2
			deriv_theta_i2 <- fac2
			# d1a <- ( rowsum( deriv_theta_i1 , dat0[,1] )[,1] + 
			# 			rowsum( deriv_theta_i2 , dat0[,2] )[,1] )
			h1 <- base::rowsum( deriv_theta_i1 , dat0[,1] )
			h2 <- base::rowsum( deriv_theta_i2 , dat0[,2] )
			d1a <- base::rep(0,TP)
			d1a[ ind1 ] <- d1a[ ind1 ] + h1
			d1a[ ind2 ] <- d1a[ ind2 ] + h2
			d1 <- score -  d1a
			# second derivative
			d2a <- M1[,1]  * (1 - M1[,1] - M1[,3] / 2 ) + 
							M1[,3]/2 * (1/2 - M1[,1] - M1[,3] / 2 )
			d2b <- M1[,2]  * (1 - M1[,2] - M1[,3] / 2 ) + 
							M1[,3]/2 * (1/2 - M1[,2] - M1[,3] / 2 )
			# d2 <- ( rowsum( d2a , dat0[,1] ) + rowsum( d2b , dat0[,2] ) )[,1]

			h1 <- base::rowsum( d2a , dat0[,1] )
			h2 <- base::rowsum( d2b , dat0[,2] )
			d2 <- base::rep(1E-20,TP)
			d2[ ind1 ] <- d2[ ind1 ] + h1
			d2[ ind2 ] <- d2[ ind2 ] + h2			
			incr <- d1/d2
			incr <- base::ifelse( base::abs(incr) > maxincr , maxincr*base::sign(incr) , incr )
			theta <- theta + incr

			theta2 <- theta								
			se.theta <- base::sqrt( 1 / d2 )
			if ( elim_persons ){
				theta[ elim_persons_index ] <- theta_elim 
				se.theta[ elim_persons_index ] <- NA
				theta2[ base::abs(theta) %in% Inf ] <- NA				
								}			
			
			if (center.theta){
				theta <- theta - base::mean(theta2 , na.rm=TRUE)	
							}
							
			if (some.fixed.theta){
				theta[ fix.theta.index ] <- fix.theta
				se.theta[ fix.theta.index ] <- NA
								}
													
			# theta <- ifelse( abs(theta) > theta.max , sign(theta)*theta.max , theta )
			# assess convergence			
			iter <- iter + 1	
			
			theta_ch <- base::abs( theta - theta0 )

			if (elim_persons){
				theta_ch <- base::ifelse( theta_ch == Inf , NA , theta_ch )
							}
			
			theta.change <- base::max( theta_ch , na.rm=TRUE)
			delta.change <- base::max( base::abs( delta - delta0 ))
			eta.change <- base::max( base::abs( eta - eta0 ))
			max.change <- base::max( c(theta.change,delta.change, eta.change) )
			
			base::cat( base::paste0("**** Iteration " , iter , 
						" | Maximum parameter change = " , base::round(max.change , 7) , "\n")
								)
			utils::flush.console()
			
						}
						
			# arrange output			
			pars <- base::data.frame("parlabel" = c("Ties" , "Home") ,
						"par" = c("delta" , "eta") )
			pars$est <- c( delta , eta )
			pars$se <- c( se.delta , se.eta )	

			# estimated individual effect
			effects <- base::data.frame( "individual" = teams , 
							"id" = base::seq( 1 , base::length(teams) ) )
			effects$Ntot <- base::rowSums(r3)
			effects$N1 <- r3[,1]
			effects$ND <- r3[,3]
			effects$N0 <- r3[,2]
			effects$score <- score
			effects$propscore <- propscore
			effects$theta <- theta
			effects$se.theta <- se.theta
			effects <- effects[ base::order( effects$theta , decreasing=TRUE) , ]
			rownames(effects) <- NULL
			# summary of effects parameters
			theta2 <- theta
			theta2[ base::abs(theta) %in% Inf ] <- NA	
			summary.effects <- base::data.frame( 
			        "M" = base::mean(theta2,na.rm=TRUE) ,
					"median" = stats::median(theta) ,
					"SD" = stats::sd(theta2,na.rm=TRUE) , 
					"min" = base::min(theta) ,
					"max" = base::max(theta) )
			
			# probabilities
			probs <- M1
			colnames(probs) <- c("p1" , "p0" , "pD")
						
			# fit statistics
			res0 <- btm_fit( probs , dat0 , ind1 , ind2 , TP)
			effects$outfit <- res0$outfit
			effects$infit <- res0$infit
			
			# MLE reliability
			v2 <- base::mean(effects$se.theta^2)
			v0 <- stats::var(effects$theta) 
			mle.rel <- 1 - v2 / v0
			sep.rel <- base::sqrt( v0 / v2 )
			# output list
			res <- base::list( effects = effects , pars=pars , summary.effects=summary.effects,
						mle.rel = mle.rel , sepG = sep.rel , 
						probs=probs , data = dat0 )			
			res$CALL <- CALL
			res$iter <- iter			
			ic <- base::list( "n" = length(teams) , "D" = nrow(dat0) )
			res$ic <- ic							
			s2 <- base::Sys.time()
			res$s1 <- s1
			res$s2 <- s2			
			base::class(res) <- "btm"
			base::return(res)
			}
#########################################################################			