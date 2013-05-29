


#------------------------------------------------------------------
# Copula estimation in a raschtype model
rasch.copula2 <- function( dat , itemcluster , 
						copula.type = "bound.mixt" ,
						progress = TRUE , mmliter = 1000 , delta = NULL ,
						theta.k = seq(-4,4,len=21) , alpha1=0 , alpha2=0 , numdiff.parm = .000001 ,
						est.b = seq(1,ncol(dat)) , est.a = rep(1,ncol(dat)) , est.delta = NULL , 
						b.init = NULL , a.init = NULL , 
						est.alpha = FALSE , 
						glob.conv = .001 , alpha.conv = .0001 , conv1 = .001 ,
						dev.crit = .2
#						pattern.off = FALSE
										){
	###############################################################
	# INPUT:
	# dat	... data frame
	# itemcluster ... vector of integers
	# progress	... display progress?
	# mmliter 	... see rasch.mml
	# delta ... initial delta estimate (is relevant for fixing delta parameters)
	# theta.k... number of grid theta points
	# alpha1, alpha2 ... rasch type parameter
	# numdiff.parm ... numerical differentiation parameter
	# est.b	... which b parameters shall be estimated
	###############################################################
	s1 <- Sys.time()
	group = NULL	
	# arrange item clusters item clusters
    t1 <- table(itemcluster)
	itemcluster[ which( itemcluster %in% names(t1)[ t1 == 1	 ] ) ] <- 0
	itemcluster <- match( itemcluster , unique( sort(itemcluster[itemcluster!=0]) ) )
	itemcluster[ is.na( itemcluster ) ] <- 0
	t1 <- table(itemcluster)	
	t1c <- t1[ names(t1) != 0 ]
	t1b <- as.numeric( names(t1))
	t1c <- t1[ names(t1) != 0 ]
	
    if ( any( t1c == 1) ){
           stop( "There should be at least two items in an item cluster\n" )
                                    }			
            x1 <- seq( 1 , max(t1b) ) 
            x2 <- sort(setdiff( as.numeric(sort(names(t1))) , 0 ))
            if ( sum( abs( x1-x2) ) > 10^(-10) ){
               stop( "Item cluster identifiers must be recoded to 1, ..., C\n" )
                                 }
	CC <- length(x1) # number of clusters
    # calculation of number of itemclusters	
    if ( progress  ){
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function rasch.copula2\n")
        cat("---------------------------------------------------------------------------------------------------------- \n")
        flush.console()
      }
	 # arrange copula types
	if ( length( copula.type ) == 1 ){ copula.type <- rep( copula.type , CC ) }
	  
	I <- ncol(dat)
	if ( is.null( colnames(dat))){
		colnames(dat) <- paste("Item" , 1:I , sep="")
				}
	# remove cases where all responses are missings
	ind <- which( rowMeans( is.na(dat) ) < 1 )
	N0 <- nrow(dat)
	dat <- dat[ ind , ]
	Nmiss <- N0 - nrow(dat)
	if (Nmiss > 0){ 
		cat(paste("Removed ", Nmiss , " cases with only missing reponses\n",sep=""))
				}
	# groups
	if ( ! is.null( group) ){
		cat("Multiple groups are ignored!\n")
		cat("This option is only implemented in 'rasch.copula'\n")
#		groups <- unique( group )		
#		G <- length( groups )	
		G <- 1 ; group <- NULL
			} else { G <- 1 }
	GG <- G
	# data preparation (frequencies)
	dat10 <- dat00 <- dat
	dat10[ is.na( dat10) ] <- 9
	patt <- paste("P",dat10[,1] ,sep="")
	for (ii in 2:I){ patt <- paste( patt , dat10[,ii],sep="") }
	pattern <- data.frame( table(patt) )
	colnames(pattern) <- c("pattern" , "freqwgt")
	# patttern in data
#	pattern.in.data <- data.frame(match( patt , pattern$pattern )
	pattern.in.data <- patt
	# calculate frequencies in multiple group case
	if (G > 1 ){ 
		for (gg in 1:G){
		# gg <- 1
		t1 <- table( patt[ group == gg  ] )
		pattern <- merge( pattern , t1 , by.x=1 , by.y = 1 , all=T )
					}
		pattern[ is.na(pattern) ] <- 0	
		colnames(pattern)[-c(1:2)] 	<- paste("freqwgt" , 1:G , sep="")		
				}
	dat0 <- matrix( 0 , nrow=nrow(pattern) , I )
	for (ii in 1:I){ 
			dat0[,ii] <- as.numeric( substring( paste( pattern[,1] ) , ii+1 , ii+1 ) )
				}
	dat0[ dat0 == 9 ] <- NA
	colnames(dat0) <- colnames(dat)
	# define mu and sigma
	mu <- rep(0,G)
	sigma <- rep(1,G)
	if ( G > 1){ 
#			mu[2] <- 1 ; sigma[2] <- 1.2 
				}
	#------------------------
	dat <- dat0
	dat2 <- dat
	dat2.resp <- 1 * ( !is.na( dat2) )			
	dat2[ is.na(dat2) ] <- 0
		# data preparation in case of item clusters
		dat2.ld <- NULL		# dat2.ld is NULL if there are no item clusters
			if ( is.null(delta)){ 
					delta <- runif( CC , .3 , .7 ) 
					delta <- ifelse( copula.type == "frank" , 1.3 , delta )
#					delta <- ifelse( copula.type == "cook.johnson" , 1 , delta )										
							}	# initial estimate of delta
			if ( is.null(est.delta)){ est.delta <- 1:CC }
			dp.ld <- as.list( 1:CC )
			# item pattern
			for (cc in 1:CC){	
	#			cc <- 1
				icl.cc <- which( itemcluster == cc )
				dp.ld.cc <- .calc.copula.itemcluster( D = length(icl.cc) )
				dp.ld.cc$items <- icl.cc
				dp.ld.cc$N.items <- NCC <- length(icl.cc)
				dp.ld.cc$itemlabels <- colnames(dat)[icl.cc]
				# item selection for independent items
				m1 <- outer( rep(1,2^NCC) , icl.cc)
				m2 <- dp.ld.cc$patt * m1   + ( 1 - dp.ld.cc$patt ) * ( m1 + I )
				m2 <- matrix( t(m2) ,  nrow=1 , byrow=T )[1,]
				res <- list( "items" = m2 )
				res1 <- rep(NCC , 2^NCC )
				names(res1) <- rownames( dp.ld.cc$patt )
				res$N.Index1 <- res1
				dp.ld.cc$independent <- res
				# item selection for dependent items
				m2 <- ( 1 - dp.ld.cc$patt ) * ( m1 + I )	
				m2 <- matrix( t(m2) ,  nrow=1 , byrow=T )[1,]				
				m2 <- c( m2[ m2> 0 ] , 2*I + 1 )
				res <- list( "items" = m2 )				
				res$N.Index1 <- rowSums( dp.ld.cc$patt  == 0 )
				res$N.Index1[ length(res$N.Index1) ] <- 1
				dp.ld.cc$dependent <- res
				dp.ld[[cc]] <- dp.ld.cc
							}
			# create data frame with item response pattern
			dat2.ld <- matrix(0 , nrow(dat2) , CC )
			dat3.ld <- as.list( 1:CC )
			for (cc in 1:CC){
				# cc <- 1	
				dp.cc <- dp.ld[[cc]]
				dat2.cc <- dat2[ , dp.cc$items ]
				l1 <- "P"
				for ( vv in seq(1,ncol(dat2.cc))){
					l1 <- paste( l1 , dat2.cc[,vv] , sep="")
						}
				dat2.ld[ , cc ] <- match( l1 , rownames( dp.cc$patt ) )
				dat2.ld[ rowSums( dat2.resp[ , dp.cc$items ]	) < length(dp.cc$items) , cc ] <- NA			
				NRR <- nrow( dp.cc$patt )
				dat3.ld.cc <- sapply( seq( 1 , NRR) , FUN = function(aa){ 1*(dat2.ld[,cc] == aa ) } )
				dat3.ld.cc[ is.na(dat3.ld.cc) ] <- 0
				dat3.ld[[cc]] <- dat3.ld.cc
							}					
				# response indicator
				dat2.ld.resp <- 1 - is.na( dat2.ld )
				# set missings in dat2 to some arbitrary category
				dat2.ld[ is.na( dat2.ld ) ] <- 1
				# create dat2 data sets for local independence items
				itemcluster0 <- ind2 <- which( itemcluster == 0 )
				bdat2.li.resp <- dat2.li.resp <- bdat2.li <- dat2.li <- NULL
				if ( length(ind2) > 0 ){ 
						dat2.li <- dat2[ , ind2 , drop=FALSE]
						dat2.li.resp <- dat2.resp[ , ind2 , drop=FALSE]
						bdat2.li <- cbind( dat2.li , 1 - dat2.li )
						bdat2.li.resp <- cbind( dat2.li.resp , dat2.li.resp )
										}
		# descriptives itemcluster
		Ncat.ld <- max( unlist( lapply( dp.ld , FUN = function(ll){ nrow(ll$patt) } ) ))
				
		# design table for estimating item difficulties
		b.design <- NULL
		if (length(itemcluster0) > 0){
			b.design <- data.frame( "itemcluster"=0 , 
					"item" = itemcluster0 , "b.indexgroup" = 1 )
									}
		for (cc in 1:CC){
		#	cc <- 1
			g.cc <- ( dp.ld[[cc]] )$items
			bb <- data.frame("itemcluster"= cc , 
						"item" = g.cc , "b.indexgroup" = seq( 1 , length(g.cc) ) )
			b.design <- rbind( b.design , bb )
					}
		b.design <- b.design[ order( b.design$item ) , ]
		b.design$est.b <- est.b
		t1 <- table( setdiff( est.b , 0 ) )
		if ( max(t1) > 1 ){	b.design$b.indexgroup <- 1:I }
		if ( sum( est.b == 0 ) > 0 ){ 
				b.design[ est.b == 0 , "b.indexgroup"  ] <- b.design
							}
		b1 <- setdiff( 1:max(b.design$b.indexgroup) , 
						setdiff( b.design$b.indexgroup , 0 ) )
		if ( length(b1) > 0 ){ 	b.design$b.indexgroup <- 1:I }						
		
	#########################################################################
	#--------------------------------------------------
	# initial estimate of item difficulty
#	b <- rasch.pairwise( dat , progress = FALSE)$item$itemdiff
	I <- ncol(dat2)
	if ( is.null( b.init) ){
		b <- - qlogis( ( colMeans( dat00 , na.rm=T ) + .005 ) / 1.01 )
				} else { b <- b.init }
	# initial estimate of (mean) item discrimination
	if ( is.null(a.init) ){ a <- rep( 1 , I ) } else { a <- a.init }
	# density weights
	wgt.theta <- dnorm(theta.k)
	wgt.theta <- wgt.theta / sum( wgt.theta )
    if ( G > 1){
		wgt.theta <- matrix(0 , length(theta.k) , G )
		for ( gg in 1:G){
			wgt.theta[,gg] <- dnorm( theta.k , mean = mu[gg] , sd = sigma[gg] )
			wgt.theta[,gg] <- wgt.theta[,gg] / sum( wgt.theta[,gg] )			
						}
				}
	iter <- 0

	#**********************************
	# BEGIN MARGINAL MAXIMUM LIKELIHOOD ESTIMATION
	dev <- 1 ; absdev.change <- par.change <- dev.change <- 1000 
	res.posterior <- NULL
	maxalphachange <- 0
	while ( ( ( absdev.change > dev.crit | dev.change > glob.conv | par.change > conv1 | maxalphachange > alpha.conv ) & iter < mmliter ) ){
		cat( paste(rep("-" , 70), collapse="") , "\n")
		k1 <- floor( log10(iter+1) )
		x1 <- "        |" 
		x1 <- substring( x1 , k1+1 )
		s1c <- Sys.time()
		cat( paste( paste( "MML EM Iter." , iter + 1 ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
						s1c  , "  " ,
						if ( iter > 0 ){ paste( round(difftime(s1c ,s1b , units='secs' ),4) , "secs" ) } , 
						"\n" ,sep="") ) # 
		s1b <- Sys.time()
		h <- numdiff.parm 
		dev0 <- dev
		#************************************
		# estimation of b parameters
		b0 <- b
		# identify different b parameter groups
#		bG <- setdiff( unique( est.b ) , 0 )
		bG <- setdiff( unique( b.design$b.indexgroup ) , 0 )		
		prbar <- seq( 1 , 10 , len = length(bG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of b:     |")		
		for (bb in bG){
			est.bb <- 1 * ( b.design$b.indexgroup == bb )
			b.design.bb <- b.design[ b.design$b.indexgroup == bb , ]
			if (bb == 1 ){ 
				rescop <- .ll.rasch.copula20( theta.k , b0 , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG , copula.type , 
							Ncat.ld	)
				res.posterior <- rescop
							}
			# is this really necessary?
#			wgt.theta <- rescop$pik
						
			rest1 <- .update.ll.rasch.copula21( theta.k , b0 + h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type)
					
			rest2 <- .update.ll.rasch.copula21( theta.k , b0 - h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type)
			ll0 <- ll1 <- ll2 <- rep(0,I)
			# numerical derivatives independent items
			if ( rescop$calc.ind ){
			for (jj in 1:(length(itemcluster0) ) ){
				ll0[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rescop$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rescop$pjk.theta.k0[,jj] ) )
				ll1[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rest1$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rest1$pjk.theta.k0[,jj] ) )
				ll2[ itemcluster0[jj] ] <- sum( rescop$rjk0.1[,jj] * log( rest2$pjk.theta.k0[,jj] ) +
								rescop$rjk0.0[,jj] * log( 1-  rest2$pjk.theta.k0[,jj] ) )								
										} }
			itemcluster1 <- b.design.bb[ b.design.bb$itemcluster > 0 , "item"]
			rjkCC <- rescop$rjkCC
			for (jj in 1:(length(itemcluster1) ) ){
					cc <- b.design.bb[ b.design.bb$item == itemcluster1[jj] , "itemcluster" ]
					ll0[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rescop$pjk.theta.kCC[[cc]] ) )
					ll1[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rest1$pjk.theta.kCC[[cc]] ) )
					ll2[itemcluster1[jj]] <- sum( rjkCC[[cc]] * log( rest2$pjk.theta.kCC[[cc]] ) )			
								}
			a1 <- aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.b) , sum , na.rm=T)					
			ll0 <- a1[,2]
			ll1 <- a1[,3]
			ll2 <- a1[,4]			
			b.change <- nr.numdiff( ll0=ll0 , ll1=ll1 , ll2=ll2 , h=h )	
			hstep <- .5^( log(iter) )
			b.change <- ifelse( abs( b.change ) > hstep , hstep*sign(b.change) , b.change )              			
			b.change <- b.change[ match( est.b , a1[,1] ) ]		
			b <- b + b.change
			cat( paste( rep( "-" , prbar[bb]), collapse="") )
			flush.console()				
					}
        a1b <- max( abs( b - b0 ) )
		cat("|     max. parm. change" , round( a1b , 5),"\n")
		
		wm1 <- sum( theta.k * rescop$pik )
		wsd <- sqrt( sum( ( theta.k - wm1 )^2 * rescop$pik ) )	
		#******************************************************************************
		# estimation of a parameters
		a0 <- a
		# identify different a parameter groups
		aG <- setdiff( unique( est.a ) , 0 )
		prbar <- seq( 1 , 10 , len = length(aG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of a:     |")
		for (aa in aG){
			est.aa <- 1 * (est.a == aa )
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type , Ncat.ld)
			ll0 <- rescop$ll
			ll1 <- .update.ll.rasch.copula20( theta.k , b, alpha1 , alpha2 , a + h*est.aa , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type )$ll
			ll2 <- .update.ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a - h*est.aa , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2					
 			if ( abs(d2) < 10^(-20) ){ d2 <- 10^20 }
			a.change <- - d1 / d2
			a.change <- ifelse( abs( a.change ) > .3 , .3*sign(a.change) , a.change )              
			a.change <- a.change * est.aa
			a <- a + a.change
			a[ a < 0 ] <- .01			
#			cat( aa , " ") ; 
			cat( paste( rep( "-" , prbar[aa]), collapse="") )
			flush.console()
							}
		if ( length(aG) < 2 ){ cat( paste( rep( "-" , 10 - length(aG) ), collapse="") ) }
		a1a <- max( abs( a - a0 ) )
		cat("|     max. parm. change" , round( a1a , 5),"\n")
		#******************************************************************************
		# estimation of delta parameters
		delta0 <- delta
		dG <- setdiff( unique( est.delta ) , 0 )
		prbar <- seq( 1 , 10 , len = length(dG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of delta: |")
		if ( length(dG) == 0 ){ cat( paste(rep("-",10),collapse="") ) }
		if ( length(dG) > 0 ){
		# identify different a parameter groups
			est.cc <- 1 # * ( est.delta == cc )
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
								CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , 
								delta , wgt.theta , I , 
								bdat2.li , bdat2.li.resp  , pattern, GG , copula.type, Ncat.ld )
#			ll0 <- rescop$ll 
			rest1 <- .update.ll.rasch.copula21( theta.k , b, alpha1 , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta + h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type)
			rest2 <- .update.ll.rasch.copula21( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta - h*est.cc , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG , copula.type)
			ll0 <- ll1 <- ll2 <- rep(0,length(dG))
			rjkCC <- rescop$rjkCC
			for (cc in 1:CC){
					ll0[cc] <- sum( rjkCC[[cc]] * log( rescop$pjk.theta.kCC[[cc]] ) )
					ll1[cc] <- sum( rjkCC[[cc]] * log( rest1$pjk.theta.kCC[[cc]] ) )
					ll2[cc] <- sum( rjkCC[[cc]] * log( rest2$pjk.theta.kCC[[cc]] ) )				
								}
			a1 <- aggregate( cbind( ll0 , ll1 , ll2 ) , list(est.delta) , sum , na.rm=T)				
			ll0 <- a1[,2]
			ll1 <- a1[,3]
			ll2 <- a1[,4]		
			delta.change <- nr.numdiff( ll0=ll0 , ll1=ll1 , ll2=ll2 , h=h )			

			ct <- sapply( dG , FUN = function(dd){
					( copula.type[ est.delta == dd ] )[1] } )
			maxstep <- ifelse( copula.type=="bound.mixt" , .2 , .9 )
			hstep <- maxstep^( log( 2+iter))
			delta.change <- ifelse( abs( delta.change ) > hstep , 
							hstep*sign(delta.change) , delta.change )              														
			delta.change <- delta.change[ match( est.delta , a1[,1] ) ]
			delta <- delta + delta.change
			delta[ delta <= 0 ] <- 2*numdiff.parm		
			delta <- ifelse( copula.type == "bound.mixt" & ( delta > 1 ) ,
						1 - 2 * numdiff.parm , delta )
			cat( paste( rep( "-" , 10 ), collapse="") )
			flush.console()
							}
        a1d <- max( abs( delta - delta0 ) )
		cat("|     max. parm. change" , round( a1d , 5),"\n")
		#******************************************************************************
		# estimation of alpha parameters
		alpha10 <- alpha1
		alpha20 <- alpha2
		prbar <- 5
		cat(" Estimation of alpha: |")		
		# alpha1
		if (est.alpha){
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld)
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula20( theta.k , b, alpha1 + h , alpha2 , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster  , pattern , GG , copula.type)$ll
			ll2 <- .update.ll.rasch.copula20( theta.k , b , alpha1 - h , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster , pattern , GG, copula.type)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k1 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha1 <- alpha1 + alpha.change
			}
			cat( paste( rep( "-" , prbar), collapse="") )
			flush.console()		
		# alpha2
		if (est.alpha){
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , pattern , GG , copula.type, Ncat.ld)
			ll0 <- rescop$ll 
			ll1 <- .update.ll.rasch.copula20( theta.k , b, alpha1 , alpha2+h , a  , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG , copula.type)$ll
			ll2 <- .update.ll.rasch.copula20( theta.k , b , alpha1 , alpha2 -h, a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta  , wgt.theta , I , 
							bdat2.li , bdat2.li.resp , rescop , itemcluster ,pattern , GG , copula.type)$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2		
			alpha.change <- - d1 / d2
			a1k2 <- alpha.change <- ifelse( abs( alpha.change ) > .1 , .1*sign(alpha.change) , alpha.change )              
			alpha2 <- alpha2 + alpha.change
			}			
			cat( paste( rep( "-" , prbar), collapse="") )
			flush.console()		
        a1k <- max( abs( c( alpha1 - alpha10, alpha2 - alpha20 )) )
		cat("|     max. parm. change" , round( a1k , 5),"\n")
		#******************************************************************************
		# estimation of mu parameters
		a1m <- 0
		if (G>1){
		mu0 <- mu
		# identify different a parameter groups
		muG <- 1:(GG-1)
		prbar <- seq( 1 , 10 , len = length(muG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of mu:    |")
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type , Ncat.ld)
			ll0 <- rescop$ll
			# mu + h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg] + h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld)
			ll1 <- rescop$ll
			# mu - h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg] - h , sd = sigma[gg] )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type, Ncat.ld )
			ll2 <- rescop$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			mu.change <- - d1 / d2
			mu.change <- ifelse( abs( mu.change ) > .3 , .3*sign(mu.change) , mu.change )      		
			mu.change <- mu.change * ( ( 1:G ) == gg )
			mu <- mu + mu.change
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] )
			wgt.theta[,gg] <- w2 / sum(w2)
			cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
			flush.console()
							}
		if ( length(muG) < 2 ){ cat( paste( rep( "-" , 10 - length(muG) ), collapse="") ) }
		a1m <- max( abs( mu - mu0 ) )
		cat("|     max. parm. change" , round( a1m , 5),"\n")
				}			# end mu
		######################################################################
		#******************************************************************************
		# estimation of sigma parameters
		a1s <- 0
		if (G>1){
		sigma0 <- sigma
		# identify different a parameter groups
#		h <- 10 * numdiff.parm
		sigmaG <- seq(1,GG-1)
		prbar <- seq( 1 , 10 , len = length(sigmaG) )
		prbar <- floor( prbar )
		prbar <- c( prbar[1] , diff(prbar) )
		cat(" Estimation of sigma: |")
		for (gg in 2:G){
#			est.aa <- est.a * (est.a == aa )		
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG,copula.type,Ncat.ld )
			ll0 <- rescop$ll
			# sigma + h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] +h)
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type,Ncat.ld)
			ll1 <- rescop$ll
			# sigma - h		
			w1 <- wgt.theta
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg]-h )
			w1[,gg] <- w2 / sum(w2)
			rescop <- .ll.rasch.copula20( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , 
							wgt.theta = w1 , I , 
							bdat2.li , bdat2.li.resp  , pattern , GG , copula.type,Ncat.ld)
			ll2 <- rescop$ll			
			d1 <- ( ll1 - ll2  ) / ( 2 * h )    
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
			sigma.change <- - d1 / d2
			sigma.change <- ifelse( abs( sigma.change ) > .3 , .3*sign(sigma.change) , sigma.change )      		
			sigma.change <- sigma.change * ( ( 1:G ) == gg )
			sigma <- sigma + sigma.change
			w2 <- dnorm( theta.k , mean = mu[gg]  , sd = sigma[gg] )
			wgt.theta[,gg] <- w2 / sum(w2)
#			cat( aa , " ") ; 
			cat( paste( rep( "-" , prbar[gg-1]), collapse="") )
			flush.console()
							}
		if ( length(sigmaG) < 2 ){ cat( paste( rep( "-" , 10 - length(sigmaG) ), collapse="") ) }
		a1s <- max( abs( sigma - sigma0 ) )
		cat("|     max. parm. change" , round( a1s , 5),"\n")
				}			# end sigma
		######################################################################		
		
		#**********************************************************************************
		# convergence display 
#		a1 <- aggregate( b , list( est.b) , mean )
#		cat("   b parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )
#		a1 <- aggregate( a , list( est.a) , mean )
#		cat("   a parameters: " , paste( round( a1[,2] , 3 ) , collapse= " " ) , "\n" )
#		cat("   delta parameters: " , paste( round( delta , 3 ) , collapse= " " ) , "\n" )
#		cat("   alpha parameters: " , paste( round( c(alpha1 , alpha2) , 3 ) , collapse= " " ) , "\n" )		
#		cat("   mu parameters: " , paste( round( mu , 3 ) , collapse= " " ) , "\n" )		
#		cat("   sigma parameters: " , paste( round( sigma , 3 ) , collapse= " " ) , "\n" )
		#******************************************************************************
		iter <- iter + 1 
		thetawidth <- diff( theta.k )[1]
		M2 <- outer( rep(1,nrow(pattern)), wgt.theta )
		#**********************
		# deviance rasch.mml
        # ll[gg] <- sum( dat1[,2] * log( rowSums( f.yi.qk[,] * 
		#			outer( rep(1,nrow(f.yi.qk[,])) , pi.k ) ) ) )		
		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm * M2 ) ) )
		#**********************
		# deviance tam
		# deviance <- - 2 * sum( pweights * log( res.hwt$rfx * thetawidth ) )		
#		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm) * thetawidth ) )
#		dev <- - 2 * sum( pattern$freqwgt * log( rowSums(rescop$post.unnorm * M2 ) ))
        dev.change <- abs( ( dev - dev0)/ dev0 )
		absdev.change <- abs( dev- dev0 )
        par.change <- max( a1a , a1b , a1d , a1k , a1m , a1s)
		cat( "Deviance = "  ,   round( dev , 5 ) , 
				" | Deviance change = " , round( absdev.change , 4 ) , 
				"| max. parm. change = " ,  round( par.change , 6 ) ,  " \n"   )  
		if ( ( dev > dev0 ) & ( iter > 4 ) ){ cat("   Deviance has increased! Convergence Problems?\n") }
		flush.console()
			}
	# end MML iterations
	#**********************************************************************************
	# Standard error estimation (This is a TO DO!)
	iterend <- iter
#	iter <- 1
#	cat( paste(rep("-" , 70), collapse="") , "\n")
#	k1 <- floor( log10(iter+1) )
#	x1 <- " |" 
#	x1 <- substring( x1 , k1+1 )	
#	cat( paste( paste( "Standard errors (SE's)"  ) , x1 , paste( rep( "*" , 10  ) , collapse="") , "|  " ,
#					Sys.time() , "\n" ,sep="") ) # 
#	h <- numdiff.parm 
#	dfr <- NULL
		#************************************
		# standard errors b
		# identify different b parameter groups
#		bG <- setdiff( unique( est.b ) , 0 )
#		prbar <- seq( 1 , 10 , len = length(bG) )
#		prbar <- floor( prbar )
#		prbar <- c( prbar[1] , diff(prbar) )
#		cat(" SE's of b:            |")	
#		for (bb in bG){
#			est.bb <- est.b * (est.b == bb )
#			if (bb == 1 ){ 
#				rescop <- .ll.rasch.copula2( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp )
#				ll0 <- rescop$ll
#							}							
#
#			ll1 <- .update.ll.rasch.copula( theta.k , b + h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp , rescop , itemcluster )
							
							
#			ll2 <- .update.ll.rasch.copula( theta.k , b - h*est.bb , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
#							CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
#							bdat2.li , bdat2.li.resp , rescop , itemcluster )
			# second order derivative
			# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
#			d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2				
#			d2 <- sqrt( 1 / sum( - d2 ) )
#			dfr.pp <- data.frame( "parameter" = "b" , "pargropup" = bb , 
#							"est" = b[ bb ] , "se" = d2 )
#			dfr <- rbind( dfr , dfr.pp )
#			cat( paste( rep( "-" , prbar[bb]), collapse="") )
#			flush.console()
#							}
#		cat("|\n")
#    print(dfr)

	################################################
	# evaluation of posterior distribution
	post <- rescop$post
	M3 <- outer( rep(1,nrow(post)) , theta.k ) 
	pattern$EAP <- rowSums( M3 * post )
	pattern$PostVar <- rowSums( M3^2 * post ) - pattern$EAP^2 
	M.EAP <- weighted.mean( pattern$EAP , pattern$freqwgt )
	Var.EAP <- sum(( ( pattern$EAP - M.EAP )^2 * pattern$freqwgt )) / ( sum( pattern$freqwgt ) )
	MVar.EAP <- weighted.mean( pattern$PostVar , pattern$freqwgt )
	EAP.Rel <- Var.EAP / ( Var.EAP + MVar.EAP )
	

	#********************************************************
	# information criteria
        # calculations for information criteria
        ic <- list( "deviance" = dev , "n" = nrow(dat00) )
        # number of parameters to be estimated
        # these formulas hold when assuming normal distributions
		bG <- setdiff( unique( est.b ) , 0 )
		aG <- setdiff( unique( est.a ) , 0 )		
		dG <- setdiff( unique( est.delta ) , 0 )				
        ic[[ "np" ]] <- length(bG) + length(aG) + length(dG) + 2*est.alpha
        # AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np	
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )			
	#**********************************************************************************
	# results item parameters			
	item <- data.frame( "item" = colnames(dat) , 
				"N" = colSums(!is.na(dat00)) , 
				"p" = colMeans( dat00 , na.rm=T ), 
				"b" = b , "est.b"= est.b , "a" = a , "est.a" = est.a )
	item$thresh <- item$a * item$b
	# add results dependency parameter for item clusters
	item$itemcluster <- itemcluster
	item$delta <- 0
    cat("---------------------------------------------------------------------------------------------------------- \n")
	for (cc in 1:CC){
		# cc <- 1
		dcc <- dp.ld[[cc]]
		item[ dcc$items , "delta"] <- delta[cc]
				}
	item$est.delta <- 0
	i1 <- which( item$itemcluster > 0 )
	icl <- item$itemcluster[ i1 ]		
	item$est.delta[ i1 ] <- est.delta[ icl ]	
	# print item summary
	cat("Parameter summary\n")
	.pr( item , digits=3 )		# print item statistics
	# dependency parameter
	cat("\nDependency parameters\n")
	summary.delta <- data.frame( "cluster" = 1:CC , "delta" = delta , 
				"est.delta" = est.delta , "copula.type" = copula.type	 )
	summary.delta$items <- sapply( 1:CC , FUN = function(cc){ 
		paste( colnames(dat)[ itemcluster == cc ] , collapse="-" )
				} )			
	.pr(summary.delta , digits = 3)
	cat(paste("\nEAP Reliability:" , round( EAP.Rel,3)),"\n\n")
	cat("Generalzed logistic link function\n")
	cat("alpha1=",round(alpha1,3)," alpha2=" , round(alpha2,3) , " \n\n")	
        # computational time
        s2 <- Sys.time()
       cat("---------------------------------------------------------------------------------------------------------- \n")
       cat("Start:" , paste( s1) , "\n")
       cat("End:" , paste(s2) , "\n")
       cat("Difference:" , print(s2 -s1), "\n")
       cat("---------------------------------------------------------------------------------------------------------- \n")

	
	datalist <- list( pattern.in.data = pattern.in.data , dat0 = dat0 ,
					dat2 = dat2 , dat2.resp = dat2.resp , dat2.li = dat2.li , 
					dat2.ld = dat2.ld , dat2.li.resp=dat2.li.resp , 
					dat2.ld.resp = dat2.ld.resp , dp.ld = dp.ld , CC = CC ,
					bdat2.li = bdat2.li , bdat2.li.resp = bdat2.li.resp , 
					itemcluster0 = itemcluster0 , dat3.ld = dat3.ld		
							)						
    # collect results
	v1 <- datalist$pattern.in.data	
	patternindex <- match( v1 , pattern$pattern )				
	res <- list( "N.itemclusters" = CC , "item" = item , "iter" = iterend , "dev" = dev ,
					"delta" = delta , "b" = b , "a" = a , "mu" = mu , "sigma" = sigma , 
					"alpha1"=alpha1 , "alpha2"=alpha2 , "ic" = ic , "theta.k" = theta.k , "deviance" = dev ,
					"pattern" = pattern	 , "datalist" = datalist	, "EAP.Rel" = EAP.Rel	,
					"copula.type" = copula.type	, "summary.delta" = summary.delta ,
					"f.qk.yi" =(res.posterior$post)[ patternindex ,] , 
					"f.yi.qk" =(res.posterior$post.unnorm)[ patternindex ,] , 					
					"s2" = s2 
								)	
	class(res) <- "rasch.copula2"
	return(res)
		}
#----------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------
# Function calculates necessary patterns for copula IRT models (Braeken, 2011)
.calc.copula.itemcluster <- function(D){
    require(gregmisc)
    res <- permutations(n=2, r=D, v=0:1, repeats.allowed=TRUE)
    rownames(res) <- apply( res , 1 , FUN = function(ll){ paste("P" , paste( ll ,collapse="") ,sep="") } )
    RR <- nrow(res) 
    matr <- matrix( 0 , RR , RR )
    rownames(matr) <- colnames(matr) <- rownames(res) 
    colnames(matr) <- gsub( "P" , "F" , colnames(matr) )
    vec <- 1:RR
    # calculation of formulas
    for (rr in vec){
        # rr <- 2
        res.rr <- outer( rep(1,nrow(res)) , res[rr,] ) - res
        a1.rr <- apply( res.rr , 1 , FUN = function(ll){ paste("F" , paste( ll ,collapse="") ,sep="") } )
        g1.rr <- ( (-1)^rowSums( res ))
        ind.rr <- which( apply( res.rr , 1 , min ) > -1 )
        a1.rr <- a1.rr[ind.rr]
        g1.rr <- g1.rr[ind.rr]
        matr[ rr , a1.rr ]  <- g1.rr
            }
    res1 <- list( "patt" = res , "calc" = matr )
    return(res1)
    }
#--------------------------------------------------------------------------------------



#*******************************************************
# Summary for rasch.copula object                         *
summary.rasch.copula2 <- function( object , ... ){
    # object      ... object from rasch.copula                #
        cat("---------------------------------------------------------------------------------------------------------- \n")
		d1 <- packageDescription("sirt")
#		cat( paste( d1$Package , d1$Version ,d1$Date ) , "\n" )
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
		cat( "Date " , paste( object$s2 ) , "\n" )
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function rasch.copula2\n")
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		cat("---------------------------------------------------------------------------------------------------------- \n")
		cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
		cat( "Number of persons = " , object$ic$n , " (" , nrow(object$pattern) , " Response Patterns)\n" )    
		cat( "Number of estimated parameters = " , object$ic$np , "\n" )   
		cat( "Number of iterations = " , object$iter , "\n" )   	
		cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
		cat( "AICc = " , round( object$ic$AICc , 2 ) , " | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )    
		cat(" (bias corrected AIC)\n" )   	
		cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , "\n" )  
		cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat( " (consistent AIC) \n\n" )         
		#******		
		cat( "Trait Distribution (" , length(object$theta.k) , " Knots )\n" , 
				  "Mean=" , 0 , " SD=" , 1 , "\n") 
		cat(paste("\nEAP Reliability:" , round( object$EAP.Rel,3)),"\n\n")			  
		cat("---------------------------------------------------------------------------------------------------------- \n")
		cat("Item Parameter \n")
		.pr( object$item , 3 )   
		cat("\nDependency parameters\n")
		.pr(object$summary.delta , digits = 3)
                }
#*******************************************************




#***************************************************
# This is an auxiliary function which helps for 
#   printing some progress
.pr <-function( object,digits){	
			ow <- options()$warn
			if ( length(dim(object)) == 2 ){
					options( warn = -1 )
					if ( nrow(object) >= 1 ){ 
						g1a <- apply( object , 2 , as.numeric )
								} else { g1a <- object }				
					g1a <- matrix(g1a , nrow= nrow(object) , ncol= ncol(object))				
					colnames(g1a) <- colnames(object)
					g1 <- colMeans( g1a )
					g1 <- which( ! is.na( g1 ) )
					options( warn = ow )
					object1 <- object
					object1[ , g1 ] <- round( object1[ , g1 ] , digits )
					# print( object1)					
					print( object1 )  } 
					else 
					{ print( round( object , digits ) ) }
				}
#***************************************************


#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula21 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0	
		eps1 <- 10^(-14)		
		# calculate necessary updates
		ind.b <- which( b != rescopula$b )
		ind.a <- which( a != rescopula$a )
		ind.delta <- which( delta != rescopula$delta )
		ind.alpha1 <- ( alpha1 != rescopula$alpha1 )	+ ( alpha2 != rescopula$alpha2 )
		if (ind.alpha1 > 0){ ind.alpha <- seq(1 , ncol(dat2.ld) ) } else { ind.alpha <- NULL }
		itemset <- union( ind.b , ind.a )
		itemset <- union( itemset , ind.alpha )	
		# update term local independence
		li.update <- 1 * ( sum( itemcluster0 %in% itemset ) > 0 )
		# update terms item dependence parameters
		ld.update <- sapply( 1:CC , FUN = function(cc){ 
				g1 <- intersect( which( itemcluster == cc )  , itemset )
				if ( length(g1)){ v1 <- cc } else { v1 <- NULL }
				v1
					} )
		ld.update <- unique( union( ind.delta , unlist( ld.update) ) )
		###########################################################################
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- rescopula$pjk.theta.kCC
		for (cc in ld.update){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){			
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
										}
			#-----------------------------------------------
			# Cook-Johnson Copula
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
					
								}
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
							} # end Frank copula
						}							
		#---------------------------------------------------------------

		
		#####################################
		# rearrange output
				res <- list( "pjk.theta.kCC"=pjk.theta.kCC , "pjk.theta.k0" = pjk.theta.k0  )
				return(res)
				}
#----------------------------------------------------------------------------------------------				



#----------------------------------------------------------------------------------------------
.update.ll.rasch.copula20 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , rescopula , itemcluster , pattern , GG , copula.type , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# look for items which change parameters for necessary update
		G <- GG
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0	
		eps1 <- 10^(-14)		
		# calculate necessary updates
		ind.b <- which( b != rescopula$b )
		ind.a <- which( a != rescopula$a )
		ind.delta <- which( delta != rescopula$delta )
		ind.alpha1 <- ( alpha1 != rescopula$alpha1 )	+ ( alpha2 != rescopula$alpha2 )
		if (ind.alpha1 > 0){ ind.alpha <- seq(1 , ncol(dat2.ld) ) } else { ind.alpha <- NULL }
		itemset <- union( ind.b , ind.a )
		itemset <- union( itemset , ind.alpha )	
		# update term local independence
		li.update <- 1 * ( sum( itemcluster0 %in% itemset ) > 0 )
		# update terms item dependence parameters
		ld.update <- sapply( 1:CC , FUN = function(cc){ 
				g1 <- intersect( which( itemcluster == cc )  , itemset )
				if ( length(g1)){ v1 <- cc } else { v1 <- NULL }
				v1
					} )
		ld.update <- unique( union( ind.delta , unlist( ld.update) ) )
		###########################################################################
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- rescopula$pjk.theta.kCC
		for (cc in ld.update){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){			
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
										}
			#-----------------------------------------------
			# Cook-Johnson Copula
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
					
								}
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1
							} # end Frank copula
						}							
		#---------------------------------------------------------------
		# calculate log likelihood
		rjk0.1 <- rescopula$rjk0.1
		rjk0.0 <- rescopula$rjk0.0
		rjkCC <- rescopula$rjkCC
		#.......................
		# one group
		if (G == 1){ 
			ll0 <- 	sum( rescopula$nk * log(rescopula$pik) )			
			 # likelihood part from independent items
			 if ( calc.ind ){	
				if ( nrow( rjk0.1) == 1 ){ 
					rjk.temp <- cbind( t( rjk0.1) , t(rjk0.0) )
							} else {
					rjk.temp <- cbind( rjk0.1 , rjk0.0 )
								}			 
					ll0 <- ll0 + sum( log(pjk.theta.k0) * rjk.temp )
							}
			 # likelihood part from dependent items
			for (cc in 1:CC){
				ll0 <- ll0 + sum( rjkCC[[cc]] * log( pjk.theta.kCC[[cc]] + 10^(-15) ) )
							}

								}
#					}
		#.......................
		# multiple groups 
		#	... to do	...		
			lli <- ll0		
		
		#####################################
		# rearrange output
				res <- list( "ll" = ll0 , "lli" = lli )
				return(res)
				}
#----------------------------------------------------------------------------------------------				




#----------------------------------------------------------------------------------------------
.ll.rasch.copula20 <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
		CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , wgt.theta , I , 
		bdat2.li , bdat2.li.resp , pattern , GG , copula.type , Ncat.ld , ... ){
		#-------------------------------------------------------
		# use this function for log likelihood calculation
		# calculation of terms of independent itemclusters?
		calc.ind <- length(itemcluster0) > 0
		G <- GG			# number of groups
		eps1 <- 10^(-14)
		ndat2 <- nrow(dat2.ld)
		M1 <- rep(1,ndat2)
		ntheta <- length(theta.k)
		M2 <- rep( 1, ntheta)
 		pjk.theta.k <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2 , a)
		pjk.theta.k01 <- cbind( pjk.theta.k , 1 - pjk.theta.k , 1  )
		
		################################################
		# E step
		################################################
		#.............................................		
		# probabilities for indepedent items
		if ( calc.ind ){
			pjk.theta.k0 <- pjk.theta.k01[ , c( itemcluster0 , itemcluster0 + I ) ]
								} else  {
			pjk.theta.k0 <- NULL
								}
		# probabilities for dependent items
		pjk.theta.kCC <- as.list( 1:CC )
	
		for (cc in 1:CC){
			# cc <- 2	# itemcluster cc
			dp.ld.cc <- dp.ld[[cc]]
			m1.cc <- pjk.theta.k01[ , dp.ld.cc$independent$items ]		
			v1.cc <- dp.ld.cc$independent$N.Index1
			#--------------------------------------------
			# Boundary Mixture Copula (Braeken, 2011)
			if (copula.type[cc] == "bound.mixt" ){
				# likelihood under independence				
				F0pjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
				# likelihood under dependence
				m1.cc <- pjk.theta.k01[ , dp.ld.cc$dependent$items ]		
				v1.cc <- dp.ld.cc$dependent$N.Index1
				pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
				F1pjk.cc <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
				pjk.theta.kCC[[cc]] <- ( 1 - delta[cc] ) * F0pjk.cc + delta[cc] * F1pjk.cc
							}
			#-----------------------------------------------
			# include other Copula models here
			# Cook-Johnson copula (Braeken et al., 2007)
			if (copula.type[cc] == "cook.johnson" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					pjk.cc[ , pp ] <- ( rowSums( ( F.Xr^(-delta.cc))^( outer( rep(1,ntheta) , ppcc )) ) 
													- R + 1 )^(-1/delta.cc)											
								}
#					pjk.theta.kCC[[cc]] <- t( dp.ld.cc$calc %*% t( pjk.cc ) )
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1					
							}
			#******************************************
			# Frank copula
			if (copula.type[cc] == "frank" ){			
				F.Xr <-  1- pjk.theta.k[ , dp.ld.cc$items ]
				R <- ncol(F.Xr)
				delta.cc <- delta[cc]
				prod.delta <- ( 1 - exp( - delta.cc ) )^(R-1)
				patt.cc <- dp.ld.cc$patt
				pjk.cc  <- matrix( 0 , nrow=ntheta , ncol= nrow(patt.cc) )
				for (pp in 1:( nrow(patt.cc) ) ){
					ppcc <- 1 - patt.cc[pp,]
					g1 <- rowProds2( ( 1 - exp( - delta.cc * F.Xr^( outer( rep(1,ntheta) , ppcc ))  ) ) )
					pjk.cc[,pp] <- - log( 1 - g1 / prod.delta ) / delta.cc										
								}
#					pjk.theta.kCC[[cc]] <- t( dp.ld.cc$calc %*% t( pjk.cc ) )	
					temp1 <- t( dp.ld.cc$calc %*% t( pjk.cc ) )												
					temp1[ temp1 < 0 ] <- eps1					
					pjk.theta.kCC[[cc]]	<- temp1					
							}  # end Frank copula
						}
								
#print( str(pjk.theta.kCC))				
# cat( "\n probabilities \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1	
						
		#.............................................
		# Calculate posterior distribution	
#		post0 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		# posterior distribution independent items
#		if ( length(itemcluster0) > 0){ 
#			post0 <- sapply( 1:ntheta , FUN = function(tt){ 
#					# tt <- 11
#					g1 <- outer( M1 , pjk.theta.k0[ tt , ] )
#					rowProds2( g1^( bdat2.li * bdat2.li.resp ) )
#						} )
#							}
							
		#****
		post0 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		I0 <- length(itemcluster0)
		if ( calc.ind ){
			pjkL <- array( NA , dim=c(2 , length(theta.k) , I0 ) )
			# pjkL: [ #categories , #thetagrid , #items ]
			pjkL[1,,] <- 1 - pjk.theta.k0[,1:I0]
			pjkL[2,,] <- pjk.theta.k0[,1:I0]
			for (ii in 1:I0 ){
				ind.ii <- which( bdat2.li.resp[,ii] == 1 )
				post0[ind.ii,] <- post0[ind.ii,] * pjkL[ bdat2.li[ind.ii,ii]+1 , ,ii]
						}
						}	
					
		# posterior distribution dependent items
#		post2 <- sapply( 1:ntheta , FUN = function(tt){ 		
				# tt <- 11
#				h1 <- 1
#				for (cc in 1:CC){
					# cc <- 1
#					pcc <- pjk.theta.kCC[[cc]]
#					h1 <- h1 * ( (pcc[tt,])[ dat2.ld[,cc] ] )^( dat2.ld.resp[,cc] )
#								}
#					h1
#							} )								
		#****
		post2 <- matrix( 1 , nrow=ndat2 , ncol=ntheta )
		pjkL <- array( NA , dim=c(Ncat.ld , length(theta.k) , CC ) )
	    for (cc in 1:CC){
			#	cc <- 1
			p1.cc <- t( pjk.theta.kCC[[cc]] )	
			pjkL[ seq( 1 , nrow(p1.cc) )  ,, cc ] <- p1.cc
						}
		for (cc in 1:CC ){
			ind.ii <- which( dat2.ld.resp[,cc] == 1 )
			post2[ind.ii,] <- post2[ind.ii,] * pjkL[ dat2.ld[ind.ii,cc] , ,cc]
					}						
						
				
		post <- post0 * post2		# product of independent and dependent parts
		post.unnorm <-  post 		
		post <- post * outer( M1 , wgt.theta )		
		post <- post / rowSums( post)	# standardization of posterior distribution

# cat( "\n posterior\n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1			

	
		#....................................................
		# Calculate expected counts
		# expected counts independent item responses
		njk0 <- rjk0.0 <- rjk0.1 <- NULL
		# pattern[,gg+1] is the frequency weight of a response pattern in group gg
		gg <- 1
		if ( calc.ind ){
#			FW1 <- outer( pattern[,gg+1] , rep( 1 , ncol(dat2.li.resp) ) )
			# students at items
#			njk0 <- t( sapply( 1:ntheta , FUN = function(tt){
#						colSums( FW1 * dat2.li.resp * outer( post[,tt] , rep( 1 , ncol(dat2.li.resp) ) ) )
#							} )	)
			njk0 <- t( post ) %*% ( pattern[,gg+1] * dat2.li.resp )
			# how many students solved correctly items
#			rjk0.1 <- t( sapply( 1:ntheta , FUN = function(tt){
#						colSums( dat2.li * FW1 * dat2.li.resp * outer( post[,tt] , rep( 1 , ncol(dat2.li.resp) ) ) )
#							} )	)
			rjk0.1 <- t( post ) %*% ( pattern[,gg+1] * dat2.li.resp * dat2.li )
			rjk0.0 <- njk0 - rjk0.1
									}
		# use loops in case of multiple groups !!!
				
		
		#*****************
		# expected counts dependent items
		rjkCC <- as.list( 1:CC )
		gg <- 1		
		for ( cc in 1:CC){	
			#	cc <- 1		
			rjkCC[[cc]] <- t(post) %*% ( pattern[,gg+1] * dat3.ld[[cc]] * dat2.ld.resp[,cc] )
					}
					
		# total counts
		gg <- 1
		tc <- outer( pattern[,gg+1] , rep( 1 , ntheta ) )
		nk <- colSums( tc * post)
		pik <- nk / sum(nk)

# cat( "\n expected counts \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1			
		
		#---------------------------------------------------------------
		# calculate log likelihood
		#.......................
		# one group
		if (G == 1){ 
			 ll0 <- sum( nk * log(pik) )					
#			 ll0 <- sum( nk * log(wgt.theta) )					
			 # likelihood part from independent items
			 if ( calc.ind ){	
				if ( nrow( rjk0.1) == 1 ){ 
					rjk.temp <- cbind( t( rjk0.1) , t(rjk0.0) )
							} else {
					rjk.temp <- cbind( rjk0.1 , rjk0.0 )
								}
					ll0 <- ll0 + sum( log(pjk.theta.k0) * rjk.temp )
							}					
			 # likelihood part from dependent items
			for (cc in 1:CC){
				ll0 <- ll0 + sum( rjkCC[[cc]] * log( pjk.theta.kCC[[cc]] + 10^(-15) ) )
								}
								
					}
# cat( "\n Likelihood \n" ) ; aa1 <- Sys.time() ; print(aa1-aa0) ; aa0 <- aa1	
				
		#.......................
		# multiple groups 
		#	... to do	...
			
			lli <- ll0
		################################
		# arrange output
		res <- list( "ll"=ll0 ,  "b" = b , "a" = a , "delta"=delta ,
						"alpha1" = alpha1 , "alpha2" = alpha2 , "lli" = lli ,
						"post" = post , "post.unnorm" = post.unnorm ,
						"rjk0.1" = rjk0.1 ,
						"rjk0.0" = rjk0.0 , "rjkCC" = rjkCC ,
						"pjk.theta.kCC" = pjk.theta.kCC , 
						"pjk.theta.k0" = pjk.theta.k0 ,
						"nk" = nk , "pik" = pik ,"calc.ind" = calc.ind
													)
		return(res)
		#~~~~~~~~~~~~~~
		# Output Version < 1.0
		#				res <- list( "ll"=ll1 , "g1"=res , "b" = b , "a" = a , "delta"=delta ,
		#									"alpha1" = alpha1 , "alpha2" = alpha2 , "lli" = lli ,
		#									"ll.theta.post" = ll.theta.post ,
		#									"exp.theta.post" = exp.theta.post
		#											)
				}
#----------------------------------------------------------------------------------------------				



###############################################################################################
# auxiliary function: person parameter estimation
.likelihood.rasch.copula <- function( theta.k , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
        bdat2.li , bdat2.li.resp , eps=10^(-20) , ... ){
        pjk.theta.k.tt <- .prob.raschtype.genlogis( theta.k , b  , alpha1 , alpha2   , a)
        g1 <- matrix( 0 , nrow(pjk.theta.k.tt) , CC + 1 )
		ndat2 <- nrow(dat2.ld)
        M1 <- rep(1,ndat2)
        pqjk.theta.k.tt <- cbind( pjk.theta.k.tt , 1 - pjk.theta.k.tt , 1 )
                    if ( length(itemcluster0) > 0 ){
                pqjk.theta.k.tt0 <- pqjk.theta.k.tt[ , c( itemcluster0 , itemcluster0+I) ]
                # likelihood for independent items at theta tt
                ll.tt <- ( pqjk.theta.k.tt0 ^ bdat2.li )^bdat2.li.resp
                g1[,1] <- rowProds2( ll.tt )        } else { g1[,1] <- 1 }
 
            g1i <- g1    
            # likelihood for dependent items
            for (cc in 1:CC){
                #                cc <- 2
                dat3.ld.cc <- dat3.ld[[cc]] 
                dp.ld.cc <- dp.ld[[cc]]
                m1.cc <- pqjk.theta.k.tt[ , dp.ld.cc$independent$items ]
                v1.cc <- dp.ld.cc$independent$N.Index1
                # product under independence                
                Fpjk.cc <- .rowProds2.bundle( m1 = m1.cc , v1 = v1.cc)
                # evaluate likelihood
                g1i[,cc+1] <- g1.tt <- ( rowSums(Fpjk.cc * dat3.ld.cc ) )^dat2.ld.resp[  ,cc]
                # product under dependence  
                m1.cc <-  pqjk.theta.k.tt[ , dp.ld.cc$dependent$items ] 
                v1.cc <- dp.ld.cc$dependent$N.Index1
                F0pjk.cc <- .rowMins2.bundle( m1 = m1.cc , v1 = v1.cc)
                F1pjk.cc <- F0pjk.cc %*% t(dp.ld.cc$calc)
                g2.tt <- ( rowSums(F1pjk.cc * dat3.ld.cc) )^dat2.ld.resp[ ,cc]
                g3.tt <- ( 1 - delta[cc] ) * g1.tt + delta[cc] * g2.tt               
                g1[,cc+1] <- g3.tt
                            }
                res <- g1   
                g1 <- rowProds2( g1 )
                g1i <- rowProds2( g1i ) 
				# calculate log likelihood
				g1[ g1 < eps] <- eps
				g1i[ g1i < eps] <- eps				
				g1 <- log( g1 )
				g1i <- log( g1i )
                res <- list( "loglike.dep" = g1 , "loglike.ind" = g1i )
                return(res)
                    }                
###############################################################################################



##################################################################
# function for person parameter estimation in rasch copula models
person.parameter.rasch.copula <- function( raschcopula.object , numdiff.parm = .001 , 
					conv.parm = .001 , maxiter = 20 , stepwidth = 1 , 
					print.summary = TRUE , 					... ){
        dat2.li <- NULL					
        dat2 <- raschcopula.object$datalist$dat2
        dat2.resp <- raschcopula.object$datalist$dat2.resp
        dat2.ld <- raschcopula.object$datalist$dat2.ld
        dat2.ld.resp <- raschcopula.object$datalist$dat2.ld.resp
        dat2.li.resp <- raschcopula.object$datalist$dat2.li.resp
        dat3.ld <- raschcopula.object$datalist$dat3.ld
        bdat2.li <- raschcopula.object$datalist$bdat2.li
        bdat2.li.resp <- raschcopula.object$datalist$bdat2.li.resp
        CC <- raschcopula.object$datalist$CC
        dp.ld <- raschcopula.object$datalist$dp.ld
        itemcluster0 <- raschcopula.object$datalist$itemcluster0
        b <- raschcopula.object$b
        a <- raschcopula.object$a
        alpha1 <- raschcopula.object$alpha1
        alpha2 <- raschcopula.object$alpha2
        delta <- raschcopula.object$delta
		pattern <- raschcopula.object$pattern
        ndat2 <- nrow(dat2)
        ntheta <- nrow(dat2)
        I <- ncol(dat2)
		######################
		# missing response pattern
		mp <- paste("R" , dat2.resp[,1] , sep="")
		for (vv in seq( 2 , ncol(dat2) ) ){
			mp <- paste( mp , dat2.resp[,vv] , sep="")
					}
        # initial estimate of theta
#        theta0 <- rowMeans( dat2 * dat2.resp  ) 
        dat20 <- dat2 * dat2.resp
		dat20[ dat2.resp == 0 ] <- NA
		theta0 <- rowMeans( dat20 , na.rm=T )
        theta0 <- qlogis( theta0 )
        # theta0[1:2] <- c(0,1)
        theta0[ theta0 == - Inf] <- -9999
        theta0[ theta0 ==  Inf] <- 9999
        theta.init <- theta0i <- theta0
        ii <- 0
        h <- numdiff.parm
        a1m <- 990
        while( a1m > conv.parm & ii < maxiter ){
            # evaluation of likelihood at theta0
            rescop0 <- .likelihood.rasch.copula( theta.k = theta0 , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop1 <- .likelihood.rasch.copula( theta.k = theta0 + h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop2 <- .likelihood.rasch.copula( theta.k = theta0 - h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop0i <- .likelihood.rasch.copula( theta.k = theta0i , b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop1i <- .likelihood.rasch.copula( theta.k = theta0i + h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )
            rescop2i <- .likelihood.rasch.copula( theta.k = theta0i - h, b , alpha1 , alpha2 , a , dat2.li , itemcluster0 , 
                        CC , dp.ld , dat2.ld , dat3.ld , dat2.ld.resp , dat2.li.resp , delta , I , 
                        bdat2.li , bdat2.li.resp  )                        
            #******
            # estimation assuming dependence
            ll1 <- rescop1$loglike.dep
            ll2 <- rescop2$loglike.dep
            ll0 <- rescop0$loglike.dep            
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2d <- d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2       
            theta.change <- - d1 / d2
            theta.change[ abs( theta.init ) == 9999 ] <- 0
			theta.change[ is.na( theta.change ) ] <- 0
            a1t1 <- theta.change <- ifelse( abs( theta.change ) > stepwidth , stepwidth*sign(theta.change) , theta.change )                        
            theta0 <- theta0 + theta.change
			ind1 <- ( abs( theta.change ) < conv.parm )

            #******
            # estimation assuming independence
            ll1 <- rescop1i$loglike.ind
            ll2 <- rescop2i$loglike.ind
            ll0 <- rescop0i$loglike.ind            
            d1 <- ( ll1 - ll2  ) / ( 2 * h )    
            # second order derivative
            # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
            d2i <- d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2       
            theta.change <- - d1 / d2
            theta.change[ abs( theta.init ) == 9999 ] <- 0
			theta.change[ is.na( theta.change ) ] <- 0
            a1t2 <- theta.change <- ifelse( abs( theta.change ) > stepwidth , stepwidth*sign(theta.change) , theta.change )                        
            theta0i <- theta0i + theta.change
            a1m <- max( abs(a1t1) , abs(a1t2) )
            ii <- ii+1
            cat("Iteration" , ii , ":   max. parm. change" , round( a1m , 5))
			cat("   |" , sum(ind1) , "out of" , ndat2 , "cases converged")
			cat(" (", round(100*sum(ind1)/length(ind1),1) , "%)\n")
			flush.console()
                    }
		theta0[ abs( theta.init  ) == 9999 ] <- NA
		theta0i[ abs( theta.init  ) == 9999 ] <- NA
        res <- data.frame( "pattern" = pattern[,1] , 
					"missing.pattern" = match( mp , unique(mp) ) , 
					"freqwgt" = pattern$freqwgt , 
					"converged" = 1*(ind1 == 1) , 
				"score" = rowSums(dat2) , "max" = rowSums( dat2.resp) , 
                 "theta.dep" = theta0 , "theta.ind" = theta0i )		
		res$setheta.dep <- sqrt( - 1 / d2d )
		res$setheta.ind <- sqrt( - 1 / d2i )
		res$setheta.dep[ is.na(theta0i) ] <- NA
		res$setheta.ind[ is.na(theta0i) ] <- NA
		res$seinflat <- res$setheta.dep / res$setheta.ind
		res[ is.na(res$theta.dep) , "converged" ] <- NA
		x1 <- seq( grep( "theta" , colnames(res) )[1] , ncol(res) ) 
		for (vv in x1){ 
				res[ paste(res$converged) == 0, vv] <- NA 
						}	
		res0 <- res <- res[ order( paste( 10000+ res$missing.pattern , 10000+ res$score)) , ]
		# calculate a summary
		index <- rep( seq(1,nrow(res0)) , res$freqwgt )
		res <- res[ index , ]
		a1 <- aggregate( res[ , c("theta.dep" , "theta.ind")] , list( res$missing.pattern ,  res$score , res$max) , mean , na.rm=T )
		colnames(a1) <- c("missing.pattern" , "score" , "max" , "M.theta.dep" , "M.theta.ind" )
		a1$N <- aggregate( 1+0*res[ , c("theta.dep")] , list( res$missing.pattern ,res$score, res$max) , sum  , na.rm=T )[,4]	
		a1 <- a1[ , c(1:3,6,4,5) ]
		a1$SD.theta.dep <- aggregate( res[ , c("theta.dep")] , list( res$missing.pattern ,res$score, res$max) , sd  , na.rm=T)[,4]		
		a1$Min.theta.dep <- aggregate( res[ , c("theta.dep")] , list( res$missing.pattern ,res$score, res$max) , min  , na.rm=T)[,4]
		a1$Max.theta.dep <- aggregate( res[ , c("theta.dep" )] , list( res$missing.pattern ,res$score, res$max) , max , na.rm=T )[,4]
		if ( abs(alpha1) + abs( alpha2) > 0 ){
			a1$SD.theta.ind <- aggregate( res[ , c("theta.ind")] , list( res$missing.pattern ,res$score, res$max) , sd  , na.rm=T)[,4]		
			a1$Min.theta.ind <- aggregate( res[ , c("theta.ind")] , list( res$missing.pattern ,res$score, res$max) , min  , na.rm=T )[,4]
			a1$Max.theta.ind <- aggregate( res[ , c("theta.ind" )] , list( res$missing.pattern ,res$score, res$max) , max  , na.rm=T)[,4]
					}		
		a1$M.seinflat <- aggregate( res[ , c("seinflat")] , list( res$missing.pattern ,res$score, res$max) , mean  , na.rm=T)[,4]
		a1$M.setheta.dep <- aggregate( res[ , c("setheta.dep" )] , list( res$missing.pattern ,res$score, res$max) , mean , na.rm=T )[,4]		
		a1$M.setheta.ind <- aggregate( res[ , c("setheta.ind" )] , list( res$missing.pattern ,res$score, res$max) , mean  , na.rm=T)[,4]		
		if (print.summary){
			cat("\n..................................................................\n")
			cat("Mean percentage standard error inflation\n\n")
			a4 <- aggregate( res[,"seinflat"] , list( res$missing.pattern) , mean , na.rm=T)
			a4[,2] <- round( 100*(a4[,2] - 1) , 2 )
			colnames(a4) <- c("missing.pattern" , "Mperc.seinflat")
			print(a4)	
			cat("\n..................................................................\n")
			cat("Summary theta estimation\n\n")
			a1b <- a1
			a1b[,seq(5,ncol(a1))] <- round( a1[ , seq(5,ncol(a1)) ] , 4 )
			print(a1b)
				}
		# person parameter estimates		
				
		ind <- match( raschcopula.object$datalist$pattern.in.data , res0$pattern )		
 
		# results		
		res <- list( "person" = res0[ ind , ] , "se.inflat" = a4 , 
						"theta.table" = res0 , "pattern.in.data" = raschcopula.object$datalist$pattern.in.data ,
						"summary.theta.table" = a1)
        return(res)
           }
#########################################################################################



								
##################################################################################
# product of rows in a matrix m1 | bundlewise calculated by a vector v1
.rowProds.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1prod <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	for (ll in 1:L1){
		if (v1[ll] > 1 ){ 
				m1prod[,ll] <- rowProds( m1[ , v1min[ll]:v1max[ll] ] )
					} else { m1prod[ll] <- m1[ , v1min[ll] ] }
						}
	m1prod
		}
##################################################################################
.rowProds2.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1prod <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	for (ll in 1:L1){
		if (v1[ll] > 1 ){ 
				m1prod[,ll] <- rowProds2( m1[ , v1min[ll]:v1max[ll] ] )
					} else { m1prod[ll] <- m1[ , v1min[ll] ] }
						}
	m1prod
		}
#*********************************************************************************
##################################################################################
# product of rows in a matrix m1 | bundlewise calculated by a vector v1
#.rowMins.bundle <- function( m1 , v1 ){
#	L1 <- length(v1)
#	m1min <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
#	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
#	v1max <- cumsum(v1)
#	m1min[ , which(v1==1)] <- m1[ , v1min[ v1 == 1 ] ]
#	for (ll in (1:L1)[ v1 > 1] ){
#				m1min[,ll] <- rowMins( m1[ , v1min[ll]:v1max[ll] ] )
#						}
#	m1min
#		}
.rowMins2.bundle <- function( m1 , v1 ){
	L1 <- length(v1)
	m1min <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )
	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )
	v1max <- cumsum(v1)
	m1min[ , which(v1==1)] <- m1[ , v1min[ v1 == 1 ] ]
	for (ll in (1:L1)[ v1 > 1] ){
				m1min[,ll] <- rowMins2( m1[ , v1min[ll]:v1max[ll] ] )
						}
	m1min
		}
##################################################################################



#************************************************************************
# Function rowProds2 
rowProds2 <- function(matr){
	y <- matr[,1]
	for (ii in 2:dim(matr)[2]){
		y <- y * matr[,ii] }
	return(y)
		}
#...................................................................
