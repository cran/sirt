



#*************************************************************************************
# E Step Multidimensional Latent Class Rasch Model                                 #
.e.step.mirtlc.mlc1 <- function( dat1 , dat2 , dat2.resp , pi.k , pjk , I , 
                   b , a  , group , G , theta.k ,  D , dimensions , Qmatrix , 
				   f.qk.yi=NULL  ){
    #...................................                    
		#***
		# array notation of probabilities
		if ( D== 1){
			pjk <- .prob.raschtype.genlogis( theta.k , b ,
					alpha1=0 , alpha2=0 , fixed.a=a )
					}
		if (D>1){
			thetaPred <- theta.k %*% t(Qmatrix )
			bPred <- matrix( b , nrow=nrow(theta.k) , ncol= I , byrow=TRUE)
			aPred <- matrix( a , nrow=nrow(theta.k) , ncol= I , byrow=TRUE)
			pjk <- plogis( aPred*(thetaPred - bPred ) )
				}
		pjkL <- array( NA , dim=c(2 , nrow(pjk) , ncol(pjk) ) )
		pjkL[1,,] <- 1 - pjk
		pjkL[2,,] <- pjk	
		if (D==1){ NT <- length(theta.k)  } else {NT <- nrow(theta.k) }
		f.yi.qk <- matrix( 1 , nrow(dat2) , NT )
		for (ii in 1:ncol(dat2)){
		#	ii <- 1
			ind.ii <- which( dat2.resp[,ii] == 1 )
			f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * pjkL[ dat2[ind.ii,ii]+1 , ,ii]
						}
		#******
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k , ncol=1 ) }
    for ( gg in 1:G){ 
        f.qk.yi[ group == gg , ] <- f.yi.qk[ group == gg , ] * 
				outer( rep( 1 , nrow(dat2[ group==gg,]) ) , pi.k[,gg] )
                    }				
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )		
	# expected counts at theta.k
	if (D==1){ NT <- length(theta.k) } else { NT <- nrow(theta.k ) }
    n.k <- matrix( 0 , NT , G )
    r.jk <- n.jk <- array( 0 , dim=c( ncol(dat2) , NT , G) )
    ll <- rep(0,G)
    for (gg in 1:G){ 
        n.k[,gg] <- colSums( dat1[group==gg,2] * f.qk.yi[group==gg,,drop=FALSE]  )
        # expected counts at theta.k and item j
        n.jk[,,gg] <- ( t(dat2.resp[group==gg,]) * outer( rep(1,I) , dat1[group==gg,2] ) ) %*% 
					f.qk.yi[group==gg, ]
        # compute r.jk (expected counts for correct item responses at theta.k for item j
        r.jk[,,gg] <- ( t(dat2[group==gg,]) * t( dat2.resp[group==gg,]) * 
						outer( rep(1,I) , dat1[group==gg,2] ) ) %*% f.qk.yi[ group==gg,]
        # compute log-Likelihood
        ll[gg] <- sum( dat1[group==gg,2] * log( rowSums( f.yi.qk[group==gg,] * 
					outer( rep( 1,nrow(f.yi.qk[group==gg,,drop=FALSE]) ) , pi.k[,gg] ) ) ) )
				}
    res <- list( "n.k" = n.k , "n.jk" = n.jk , "r.jk" = r.jk , "f.qk.yi" = f.qk.yi , "pjk" = pjk  ,
            "f.yi.qk" = f.yi.qk , "ll" = sum(ll) )
    return(res)
    }
#*************************************************************************************





#*************************************************************************************
# E Step Raschtype Model                                                        #
.e.step.mirtlc.lc <- function( dat1 , dat2 , dat2.resp , pi.k , pjk , I , 
                   group , G , theta.k ,  f.qk.yi=NULL  ){
    #...................................                    
		#***
		# array notation of probabilities
		pjkL <- array( NA , dim=c(2 , nrow(pjk) , ncol(pjk) ) )
		pjkL[1,,] <- 1 - pjk
		pjkL[2,,] <- pjk	
		f.yi.qk <- matrix( 1 , nrow(dat2) , length(theta.k) )
		for (ii in 1:ncol(dat2)){
		#	ii <- 1
			ind.ii <- which( dat2.resp[,ii] == 1 )
			f.yi.qk[ind.ii,] <- f.yi.qk[ind.ii,] * pjkL[ dat2[ind.ii,ii]+1 , ,ii]
						}
		#******
    f.qk.yi <- 0 * f.yi.qk
    if ( G==1 ){ pi.k <- matrix( pi.k , ncol=1 ) }
    for ( gg in 1:G){ 
        f.qk.yi[ group == gg , ] <- f.yi.qk[ group == gg , ] * outer( rep( 1 , nrow(dat2[ group==gg,]) ) , pi.k[,gg] )
                    }				
    f.qk.yi <- f.qk.yi / rowSums( f.qk.yi )		
	# expected counts at theta.k
    n.k <- matrix( 0 , length(theta.k) , G )
    r.jk <- n.jk <- array( 0 , dim=c( ncol(dat2) , length(theta.k) , G) )
    ll <- rep(0,G)
    for (gg in 1:G){ 
        n.k[,gg] <- colSums( dat1[group==gg,2] * f.qk.yi[group==gg,,drop=FALSE]  )
        # expected counts at theta.k and item j
        n.jk[,,gg] <- ( t(dat2.resp[group==gg,]) * outer( rep(1,I) , dat1[group==gg,2] ) ) %*% 
					f.qk.yi[group==gg, ]
        # compute r.jk (expected counts for correct item responses at theta.k for item j
        r.jk[,,gg] <- ( t(dat2[group==gg,]) * t( dat2.resp[group==gg,]) * 
						outer( rep(1,I) , dat1[group==gg,2] ) ) %*% f.qk.yi[ group==gg,]
        # compute log-Likelihood
        ll[gg] <- sum( dat1[group==gg,2] * log( rowSums( f.yi.qk[group==gg,] * 
					outer( rep( 1,nrow(f.yi.qk[group==gg,,drop=FALSE]) ) , pi.k[,gg] ) ) ) )
				}
    res <- list( "n.k" = n.k , "n.jk" = n.jk , "r.jk" = r.jk , "f.qk.yi" = f.qk.yi , "pjk" = pjk  ,
            "f.yi.qk" = f.yi.qk , "ll" = sum(ll) )
    return(res)
    }
#*************************************************************************************



########################################################
# calculate class probabilities
.m.step.mirtlc.mlc1 <- function( pjk , n.k , r.jk , n.jk , G , Nclasses ,
			theta.k , b , a , I , ref.item , mstep.maxit ,
			des.theta , des.b , theta.fixed , theta.normal , f.qk.yi,D ,
			distribution.trait , est.a , Qmatrix ,modeltype , range.b , range.a){
    if (G==1){ 
		pi.k <- n.k / sum( n.k )
			}
	if ( G> 1){
		pi.k <- n.k / matrix( colSums(n.k ) , nrow=Nclasses , ncol=G , byrow=TRUE)
			}							
	#*****
	# perform logistic regression			
	if (G== 1){
		c1 <- r.jk[,,1]
		i1 <- n.jk[,,1] - c1
		 }
	if (G > 1){
		c1 <- apply( r.jk , c(1,2) , sum )
		i1 <- apply( n.jk - r.jk , c(1,2) , sum )
				}

	# class weights correct response
	wc1 <- matrix( c1 , nrow=I*Nclasses , 1 )
	wi1 <- matrix( i1 , nrow=I*Nclasses , 1 )
	# outcome
	y <- c( rep(1,I*Nclasses) , rep(0,I*Nclasses) )
	wgt <- c( wc1 , wi1 )
	if ( theta.fixed ){
		if (D==1){ 
			theta.offset <- rep( theta.k , each=I ) 
     		theta.offset <- c( theta.offset , theta.offset )
				}
		if (D>1){
			tk1 <- matrix( theta.k , nrow=1 , ncol=ncol(des.theta) )
			tk1 <- matrix( tk1 , nrow=nrow(des.theta) , ncol=ncol(tk1) , byrow=TRUE)
			theta.offset <- des.theta * tk1
			theta.offset <- rowSums( theta.offset )
				}
					}
	if ( ( ! theta.fixed ) ){
		# structure theta
		# theta_1: class_1 , ... , class_D
		# theta_2: class_1 , ... , class_D
		# ... theta_D: class_1 , ... , class_D
		# unconstrained theta optimization
		if (modeltype == "MLC2"){
			des1 <- matrix( a , nrow=nrow(des.theta) , ncol = 1)
			des.theta <- des.theta * outer( des1[,1] , rep(1,ncol(des.theta) ) )
			des.b <- des.b * outer( des1[,1] , rep(1,ncol(des.b) ) )
								}
    		mod2 <- glm( y ~ 0 + des.theta + des.b  , weights = wgt , family ="binomial" ,
				control=list(maxit=mstep.maxit ) )
		if (D==1){ 
			theta.k <- coef(mod2)[ 1:Nclasses ]	# theta
			b0 <- coef(mod2)[ -c( 1:Nclasses ) ]
			    }
		if (D>1){
			theta.k0 <- theta.k
			theta.k <- coef(mod2)[ 1:(D*Nclasses) ]	# theta
			theta.k <- matrix( theta.k , nrow=Nclasses , ncol=D )
			b0 <- coef(mod2)[ -c( 1:(D*Nclasses) ) ]
			    }
		b[ setdiff( 1:I , ref.item ) ] <- b0
						}
	if ( theta.fixed ){
		# constrained theta optimization
		if (modeltype == "MLC2"){
			des1 <- matrix( a , nrow=nrow(des.theta) , ncol = 1)
			theta.offset <- des1[,1] * theta.offset
			des.b <- des.b * outer( des1[,1] , rep(1,ncol(des.b) ) )
								}
		mod2 <- glm( y ~ 0 + offset(theta.offset) + des.b  , weights = wgt , family ="binomial" ,
				control=list(maxit=mstep.maxit ) )
		b0 <- coef(mod2)
		b[ setdiff( 1:I , ref.item ) ] <- b0	
		############################################
		# normal distribution assumption
		if (distribution.trait == "normal" & ( D==1) ){ # D=1	
			for (gg in 1:G){
				pik1 <-	pi.k[,gg]
				m1 <- sum( theta.k * pik1 )
				sd1 <- sqrt( sum( theta.k^2 * pik1 ) - m1^2 )
				pik2 <- dnorm( theta.k , mean=m1 , sd = sd1 )
				pi.k[,gg] <- pik2 / sum(pik2)
							}
						}
		############################################
		# log linear smoothing: allow for skewness
		if (distribution.trait %in% c("smooth2" , "smooth3" , "smooth4") & 
						( D==1) ){ # D=1	
			for (gg in 1:G){
				pik1 <-	pi.k[,gg]
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
						formula1 <- lpik1 ~ tk + I(tk^2) + I(tk^3) + I(tk^4)
										}
				mod <- lm( formula1 , weights = pik1 )
				pik2 <- exp( fitted(mod))
				pi.k[,gg] <- pik2 / sum(pik2)
							}
						}
		##################################################		
						}
	# range restrictions
	b[ b < range.b[1] ] <- range.b[1]
	b[ b > range.b[2] ] <- range.b[2]
	## 2PL estimation
	if (modeltype == "MLC2" ){
		a <- .mirtlc.est.a( theta.k=theta.k , b=b , fixed.a=a , 
				pjk=pjk , alpha1=0 , alpha2=0 , 
				h=.0001 , G=G , I=I , r.jk=r.jk , n.jk=n.jk , est.a=est.a , Qmatrix=Qmatrix )	
 		a[ a < range.a[1] ] <- range.a[1]
	    a[ a > range.a[2] ] <- range.a[2]
						}
	res <- list( "pi.k" = pi.k , "pjk" = pjk , "theta.k" = theta.k , 
			"b" = b , "a" = a )
    return(res)
            }
########################################################




########################################################
# calculate class probabilities
.m.step.mirtlc.lc <- function( pjk , n.k , r.jk , n.jk , G , Nclasses ){
    if (G==1){ 
		pi.k <- n.k / sum( n.k )
		pi.k <- matrix( pi.k , nrow=ncol(pjk) , ncol=nrow(pjk) )
			}
	if ( G> 1){
		pi.k <- n.k / matrix( colSums(n.k ) , nrow=Nclasses , ncol=G , byrow=TRUE)
			}
    for (cc in 1:Nclasses ){
    # cc <- 1
    if (G==1){ 
			pjk[ cc , ] <- r.jk[  , cc , 1] / n.jk[ , cc , 1] 
				}
    if (G>1){ 
			pjk[ cc , ] <- rowSums( r.jk[  , cc , ] ) / rowSums( n.jk[ , cc , ]  )
				}
                    }
    res <- list( "pi.k" = pi.k , "pjk" = pjk )
    return(res)
            }
########################################################


#*********************************************************************		
# Estimation of a parameter (discrimination parameter)		
.mirtlc.est.a <- function( theta.k , b , fixed.a , 
					pjk , alpha1 , alpha2 , h , G , I , r.jk , n.jk , est.a , Qmatrix ){
				# cc <- cG[1]	
#					est.aa <- 1 * (est.a == aa )
					#****
					# a
					pjk <- .prob.raschtype.genlogis( theta.k , b , alpha1 , alpha2 , fixed.a ,
								Qmatrix)
					pjk <- ( pjk + .000000005 ) / 1.00000001 
					pjk.M <- t(pjk) ;  qjk.M <- 1 - pjk.M
					pjk1 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a + h,Qmatrix)			
					pjk1 <- ( pjk1 + .000000005 ) / 1.00000001 
					pjk1.M <- t(pjk1) ; qjk1.M <- 1 - pjk1.M
					pjk2 <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , fixed.a - h,Qmatrix)
					pjk2 <- ( pjk2 + .000000005 ) / 1.00000001 
					pjk2.M <- t(pjk2) ; qjk2.M <- 1 - pjk2.M
					# first order derivative
					# f(x+h) - f(x-h) = 2* f'(x) * h
					ll0 <- ll1 <- ll2 <- matrix(0, I , G)
				# was ist hier das G? => G ist hier Anzahl der Gruppen
					for (gg in 1:G){ 
						ll0[,gg] <- rowSums( r.jk[,,gg] * log( pjk.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk.M ) )
						ll1[,gg] <- rowSums( r.jk[,,gg] * log( pjk1.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk1.M ) )
						ll2[,gg] <- rowSums( r.jk[,,gg] * log( pjk2.M ) + ( n.jk[,,gg] - r.jk[,,gg]  ) * log( qjk2.M ) 	)		
								   }  
					ll0 <- rowSums(ll0)
					ll1 <- rowSums(ll1)
					ll2 <- rowSums(ll2)
					# aggregate with respect to estimation of a
					a1 <- aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.a) , sum , na.rm=T)	
					a1 <- a1[ a1[,1] > 0 , ]					
					ll0 <- a1[,2]
					ll1 <- a1[,3]
					ll2 <- a1[,4]			
					d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
					# second order derivative
					# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
					d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
					d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
					# change in item difficulty
					a.change <- - d1 / d2
#					a.change <- ifelse( abs( a.change ) > .2 , .2*sign(a.change) , a.change )              
					# dampening parameter as in tam
					old_increment <- .2
					ci <- ceiling( abs(a.change) / ( abs( old_increment) + 10^(-10) ) )
					a.change <- ifelse( abs( a.change) > abs(old_increment)  , 
										a.change/(2*ci) , a.change )					
	#					a.change <- a.change
					a.change <- a.change[ match( est.a , a1[,1] ) ]
					if ( any( est.a == 0 ) ){
						a.change[ est.a == 0 ] <- 0
									}									
					fixed.a <- fixed.a + a.change
					fixed.a[ fixed.a < 0 ] <- 0
				return(fixed.a)
					}
#*********************************************************************