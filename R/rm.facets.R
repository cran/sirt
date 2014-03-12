
#################################################################
# Facets Model for Raters:
# MML estimation
rm.facets <- function( dat , pid=NULL , rater=NULL ,
	Qmatrix=NULL , theta.k=seq(-9,9,len=30) , 
	est.b.rater=TRUE , est.a.item=FALSE , est.a.rater=FALSE ,
	max.b.increment=1 , numdiff.parm=.00001 , maxdevchange=.10 ,
	globconv=.001 , maxiter=1000 , msteps=4 , mstepconv=.001){
	s1 <- Sys.time()
	dat <- as.matrix(dat)
	if ( is.null(rater)){	
		rater <- rep(1,nrow(dat)) 
		est.b.rater <- FALSE
		est.a.rater <- FALSE		
		    }
	if ( is.null(pid)){  pid <- seq(1,nrow(dat) ) }
	pcm.param <- FALSE	
	theta.k0 <- theta.k
	pi.k <- dnorm( theta.k )
	pi.k <- pi.k / sum( pi.k )
	# process data
	procdata <- res <- .prep.data.rm( dat=dat , rater=rater , pid=pid )
	dat2 <- as.matrix(res$dat2)
	dat2.resp <- as.matrix(res$dat2.resp)
	rater.index1 <- res$rater.index
	dataproc.vars <- res$dataproc.vars
	VV <- res$VV
	RR <- res$RR
	item.index <- res$dataproc.vars$item.index 
	rater.index <- res$dataproc.vars$rater.index 
	
	# maximum categories
	maxK <- apply( dat , 2 , max , na.rm=T )
	K <- max( maxK )
	if ( is.null(Qmatrix) ){
		Qmatrix <- matrix( 1:K , nrow=VV , ncol=K , byrow=T)
					}
	TP <- length(theta.k)
	I <- VV*RR

	# define constraints on tau.item parameters
	# if not all categories are observed
	tau.item.fixed <- NULL
	if ( min(maxK) < K ){
		for (vv in 1:VV){
		#vv <- 1
			K.vv <- maxK[vv]
			if ( K.vv < K ){
				for (zz in (K.vv+1):K ){
					d1 <- data.frame( "item"= vv, 
							"categ"=zz , "val"=99 )	
					tau.item.fixed <- rbind( tau.item.fixed , d1 ) 				
									}
							}
					}
		tau.item.fixed <- as.matrix(tau.item.fixed )
				}
	# starting values for item difficulties
	b.item <- - qlogis( colMeans( dat , na.rm=T ) / maxK  )
	if ( ! pcm.param ){ b.item <- 0*b.item	}
	
	tau.item <- matrix( 0 , nrow=VV , ncol=K )
	rownames(tau.item) <- colnames(dat)
	tau.item <- matrix( seq( -2 , 2 , len=K ) , nrow=VV , ncol=K , byrow=T )

	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
	N <- aggregate( N , list( rater.index ) , sum )[,2]
	M1 <- aggregate( M1 , list( rater.index ) , mean )[,2]		
	b.rater <- - qlogis( M1 / K )
	b.rater <- b.rater - mean( b.rater )
	a.item <- rep(1,VV)
	a.rater <- rep(1,RR)
	if ( ! est.b.rater ){ b.rater <- rep(0,RR) }

	# init standard errors
	se.b.rater <- NA*b.rater
	se.a.rater <- NA*a.rater
	se.a.item <- NA*a.item
		
	# inits
	iter <- 0
	dev0 <- dev <- 0
	conv <- devchange <- 1000
	sigma <- 1
	disp <- "...........................................................\n"	
	
	#****************************************************
	# start EM algorithm
    while( ( ( maxdevchange < devchange ) | (globconv < conv) ) &
			( iter < maxiter )
						){
		cat(disp)	
		cat("Iteration" , iter+1 , "   " , paste( Sys.time() ) , "\n" )	
		
		# previous values
		b.item0 <- b.item
		b.rater0 <- b.rater
		tau.item0 <- tau.item
		dev0 <- dev
		sigma0 <- sigma
		a.item0 <- a.item
		a.rater0 <- a.rater
# zz0 <- Sys.time()
		# calculate probabilities
		probs <- .rm.facets.calcprobs2( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
				theta.k ,RR )				
# cat("facets.calcprob   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1						
		# calculate posterior
		res <- .rm.posterior( dat2 , dat2.resp , TP , pi.k , K, I , probs )
# cat("posterior   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1								
		f.yi.qk <- res$f.yi.qk
		f.qk.yi <- res$f.qk.yi
		n.ik <- res$n.ik
		N.ik <- res$N.ik
		pi.k <- res$pi.k
		ll <- res$ll
		# estimate b.rater parameter
		if( est.b.rater ){
			if (iter ==0){	max.b.increment -> b.rater.incr }
			res <- .rm.facets.est.b.rater( b.item , b.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=b.rater.incr,theta.k ,
					msteps , mstepconv  )
			b.rater <- res$b.rater
			se.b.rater <- res$se.b.rater
			b.rater.incr <- abs( b.rater0 - b.rater )
						}
# cat("est b rater   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1														
		# estimate tau.item parameters
		if (iter ==0){	max.b.increment -> tau.item.incr }
		res <- .rm.facets.est.tau.item( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=tau.item.incr  , theta.k ,
				msteps, mstepconv , tau.item.fixed )
		tau.item <- res$tau.item
		se.tau.item <- res$se.tau.item
		tau.item.incr  <- abs( tau.item0 - tau.item )
# cat("est tau item   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1															
		
		# estimate a.item parameter
		if (est.a.item){
			res <- .rm.facets.est.a.item( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1,theta.k ,
				msteps, mstepconv )		
			a.item <- res$a.item
			se.a.item <- res$se.a.item
				}
# cat("est a item   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1																
		# estimate a.rater parameter
		if (est.a.rater){
			res <- .rm.facets.est.a.rater( b.item , b.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , a.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1,theta.k ,
				msteps, mstepconv )		
			a.rater <- res$a.rater
			se.a.rater <- res$se.a.rater
				}
# cat("est a rater   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1																			
		flush.console()		
		# update distribution
		w2 <- sum( theta.k^2 * pi.k )
		sigma <- sqrt(w2)
		pi.k <- dnorm( theta.k , mean=0 , sd=sigma )
		pi.k <- pi.k / sum( pi.k )
		dev <- -2*ll
		# convergence criteria
		conv <- max( abs(b.rater-b.rater0) , abs( a.rater-a.rater0) , 
					abs( tau.item0-tau.item) , abs( a.item - a.item0 ) )
		iter <- iter+1
		devchange <- abs( ( dev - dev0 ) / dev0  )
# cat("calc ll   ") ; zz1 <- Sys.time(); print(zz1-zz0) ; zz0 <- zz1																			

		
		#****
		# print progress			
		cat( paste( "   Deviance = "  , round( dev , 4 ) , 
			if (iter > 1 ){ " | Deviance change = " } else {""} ,
			if( iter>1){round( - dev + dev0 , 6 )} else { ""}	,"\n",sep="") )
		cat( paste( "    Maximum b.rater parameter change = " , 
				paste( round(max(abs(b.rater0-b.rater)) ,6) , collapse=" " ) , "\n" , sep=""))				
		cat( paste( "    Maximum a.rater parameter change = " , 
				paste( round(max(abs(a.rater0-a.rater)) ,6) , collapse=" " ) , "\n" , sep=""))															
		cat( paste( "    Maximum tau.item parameter change = " , 
				paste( round(max(abs(tau.item0-tau.item)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste( "    Maximum a.item parameter change = " , 
				paste( round(max(abs(a.item0-a.item)) ,6) , collapse=" " ) , "\n" , sep=""))
		cat( paste(" Trait SD = " , round( sigma , 3 ) , sep="") , "\n")
		# flush.console()			
				}
				
	# *********
	# arrange OUTPUT
	#---
	# Information criteria
	ic <- list( "deviance" = dev , "n" = nrow(dat2) )
	ic$VV <- VV
	ic$RR <- RR
	ic$np.item <- sum( maxK) + est.a.item*(VV-1)
	ic$np.rater <- est.b.rater*(RR-1) + est.a.rater*(RR-1)
	ic$np <- 1 + ic$np.item + ic$np.rater
    # AIC
    ic$AIC <- dev + 2*ic$np
    # BIC
    ic$BIC <- dev + ( log(ic$n) )*ic$np
    # CAIC (conistent AIC)
    ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
	# corrected AIC
    ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )		
	
	#---
	# person
	person <- procdata$person.index
	NP <- nrow(person)
	person$score <- rowSums( dat2 * dat2.resp )
	mkrr <- rep( maxK , RR )
	person$maxscore <- rowSums( dat2.resp * outer( rep(1,NP) , mkrr ) )
	person$EAP <- rowSums( f.qk.yi * outer( rep(1,NP) , theta.k) )
	person$SE.EAP <- sqrt( rowSums( f.qk.yi * outer( rep(1,NP) , theta.k^2) ) - 
			( person$EAP) ^2 )
	EAP.rel <- 1 - mean( person$SE.EAP^2 ) / 
				( mean( person$SE.EAP^2 ) + var( person$EAP ) )

	
	#---
	# item
	if (!is.null(tau.item.fixed)){
		tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA
		se.tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- NA		
							}
	
    item <- data.frame( "item" = colnames(dat) , 
			"N" = colSums( 1-is.na(dat)) , 
			"M" = colMeans( dat , na.rm=T ) )
	for (kk in 1:K){ item[ , paste0("tau.Cat",kk) ] <- tau.item[,kk] }
    item$a <- a.item

	obji <- item
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
	cat("*********************************\n")
	cat("Item Parameters\n")
    print( obji )		
	
	#---
	# rater
	M1 <- colSums( dat2 ) / colSums( dat2.resp )
	N <- colSums( dat2.resp )
	N <- aggregate( N , list( rater.index ) , sum )[,2]
	M1 <- aggregate( M1 , list( rater.index ) , mean )[,2]	
    rater <- data.frame( "rater" = rater.index1[,1] , 
			"N" = N , 
			"M" = M1 , 
			"b" = b.rater ,
			"a" = a.rater )
	rater$thresh <- rater$a * rater$b
	obji <- rater
	for (vv in seq(2,ncol(obji) )){
		obji[,vv] <- round( obji[,vv],3 ) }
	cat("*********************************\n")
	cat("Rater Parameters\n")		
    print( obji )		
	cat("*********************************\n")
	cat("EAP Reliability = " , round(EAP.rel,3) , "\n")		
	
	s2 <- Sys.time()
	
    res <- list("deviance" = dev , "ic"=ic , "item"=item , "rater"=rater ,
		"person" = person , "EAP.rel"=EAP.rel , 
		"sigma"=sigma , 
		"tau.item"=tau.item , "se.tau.item"=se.tau.item ,
		"a.item"=a.item , "se.a.item"=se.a.item ,
		"b.rater"=b.rater , "se.b.rater"=se.b.rater , 
		"a.rater"=a.rater , "se.a.rater"=se.a.rater , 
		"f.yi.qk"=f.yi.qk , "f.qk.yi"=f.qk.yi , "probs"=probs ,
		"n.ik"=n.ik , "maxK"=maxK , "procdata" =procdata , "iter"=iter , 
		"s1"=s1 , "s2"=s2 , "tau.item.fixed"=tau.item.fixed)
	class(res) <- "rm.facets"
	return(res)

		}
