
################################################################
# C code
# calculation of probabilities
probraterfct1 <- function (crater,drater,dimA,B,dimB){ 
.Call("file21a045a71b7", crater,drater,dimA,B,dimB, PACKAGE = "sirt")
					}
# array multiplication					
arraymult1 <- function (A,dimA,B,dimB){ 
.Call("file211c19faab", A, dimA, B, dimB, PACKAGE = "sirt")
					}					

################################################################
# calculate probabilities
.rm.hrm.calcprobs  <- function(  c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				theta.k ,RR , prob.item=NULL , prob.rater = NULL ){
	# calculate probabilities for true ratings
	a <- a.item
	b <- tau.item
	if ( is.null( prob.item ) ){ 
		res <- .rm.pcm.calcprobs( a , b , Qmatrix=Qmatrix , theta.k , I=VV , K , TP )
				} else { res <- prob.item }
	# calculate probabilities for raters
	calc.rater <- FALSE
	if (is.null(prob.rater)){
		calc.rater <- TRUE
						}
# calc.rater <- TRUE						
#	prob.categ <- array( 0 , dim=c(I , K+1 , TP ) )	
#	if (calc.rater){
#		for (ii in 1:I){
	#		ii <- 1	
#			h1[ii,,] <- matrix( c.rater[ii,] , nrow=K+1 , ncol=K,byrow=T) - 
#				         matrix( (0:K) * d.rater[ii] , nrow=K+1 , ncol=K,byrow=F )
#		h1[ii,,] <- - outer( (0:K)*d.rater[ii] , as.vector(c.rater[ii,]) , "-" )
#						}
#		h1 <- plogis( h1 )
#		prob.rater[,1,] <- h1[,,1]
#		for (kk in 1:(K-1)){ prob.rater[,kk+1,] <- h1[,,kk+1] - h1[,,kk] }
#		prob.rater[,K+1,] <- 1-h1[,,K]
#				}
				
	dimA <- c(I , K+1, K+1 )
	res2 <- res[ item.index ,,]	
	dimB <- dim(res2)					
	BM <- matrix( res2 , dimA[1]*dimB[2] , dimB[3] )
	#****	
	# if prob.rater is calculated	
	if (calc.rater){
		res2 <- probraterfct1( crater=c.rater , drater=d.rater , 
				dimA=dimA,B=BM,dimB=dimB)
		prob.categ <- array( res2$probtotal , dim= c(dimA[c(1,2)],dimB[3]) )
		prob.rater <- array( res2$PRA , dim=dimA )		
					}
	#***
	# if prob.rater is not calculated
	AM <- matrix( prob.rater , dimA[1]*dimA[2] , dimA[3] )
    if ( ! calc.rater ){
		y <- arraymult1( AM , dimA , BM , dimB )
		prob.categ <- array( y , dim= c(dimA[c(1,2)],dimB[3]) )    
					}
	
	res <- list("prob.total"=prob.categ , "prob.rater"=prob.rater , 
		"prob.item"	= res )
	return(res)	
			}
#############################################################
###################################################			
# c.rater
.rm.hrm.est.c.rater <- function(  c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=1,theta.k ,
					msteps , mstepconv , est.c.rater , prob.item  ){
    h <- numdiff.parm
	if (est.c.rater=="r"){ diffindex <- rater.index }
	if (est.c.rater=="i"){ diffindex <- item.index }
	if (est.c.rater=="e"){ diffindex <- rep(1,I) }	
	if (est.c.rater=="a"){ diffindex <- 1:I }		
	RR <- I/VV
	cat("  M steps c.rater parameter    |")
	it <- 0 ;	conv1 <- 1000
	se.c.rater <- 0 * c.rater
#	max.b.increment <- max.b.increment + 0*c.rater	
	
	Q0 <- 0 * c.rater
	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){
		b0 <- c.rater
		for (kk in 1:K){	
	#		kk <- 1
			Q1 <- Q0
			Q1[,kk] <- 1			
			pjk <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=prob.item, prob.rater=NULL)$prob.total				
			pjk1 <- .rm.hrm.calcprobs( c.rater+h*Q1, Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=prob.item, prob.rater=NULL)$prob.total				
			pjk2 <- .rm.hrm.calcprobs( c.rater-h*Q1 , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=prob.item, prob.rater=NULL)$prob.total		
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.b.increment , numdiff.parm )				
					
			c.rater[,kk] <- c.rater[,kk] + res$increment[diffindex]
			se.c.rater[,kk] <- ( sqrt( abs(-1/res$d2) ) )[diffindex ]
			if (kk>1){ 
				ind <- which( c.rater[,kk] < c.rater[,kk-1] )
				if (length(ind)>0){
					l1 <- c.rater[ind,kk-1]
					c.rater[ind, kk-1] <- c.rater[ind,kk]				
					c.rater[ ind , kk ] <- l1
							}
						}
					}
	#		max.b.increment <- abs( b.rater - b0 )
		conv1 <- max( abs( c.rater - b0 ) )
		it <- it+1
		cat("-")  #; flush.console()
			}
	cat(" " , it , "Step(s) \n")	
    res <- list("c.rater" = c.rater , "se.c.rater" = se.c.rater , 
			"ll" = sum(res$ll0) 
				)
    return(res)
			}				

			
###################################################			
# d.rater
.rm.hrm.est.d.rater <- function(  c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
					n.ik , numdiff.parm, max.b.increment=1,theta.k ,
					msteps , mstepconv , d.min , d.max , est.d.rater , prob.item ){
    h <- numdiff.parm
	
	if (est.d.rater=="r"){ diffindex <- rater.index }
	if (est.d.rater=="i"){ diffindex <- item.index }
	if (est.d.rater=="e"){ diffindex <- rep(1,I) }	
	if (est.d.rater=="a"){ diffindex <- 1:I }			
	
	RR <- I/VV	
	cat("  M steps d.rater parameter    |")
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){
		b0 <- d.rater
		r1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
				prob.item=prob.item,prob.rater=NULL)				
		pjk <- r1$prob.total				
		pjk1 <- .rm.hrm.calcprobs( c.rater, Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater+h , item.index , rater.index , theta.k,RR,
				prob.item=prob.item,prob.rater=NULL)$prob.total				
		pjk2 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater-h , item.index , rater.index , theta.k,RR,
				prob.item=prob.item,prob.rater=NULL)$prob.total		
		# numerical differentiation			
		res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				max.increment=max.b.increment , numdiff.parm )		
		d.rater <- d.rater + res$increment[diffindex]
		d.rater[ d.rater < d.min ] <- d.min		
		d.rater[ d.rater > d.max ] <- d.max				
#		max.b.increment <- abs( b.rater - b0 )
		conv1 <- max( abs( d.rater - b0 ) )
		it <- it+1
		cat("-")  #; flush.console()
			}
	cat(" " , it , "Step(s) \n")	
    res <- list("d.rater" = d.rater , "se.d.rater" = sqrt( abs(-1/res$d2) ) , 
			"ll" = sum(res$ll0) 
				)
    return(res)
			}				

			
			
			
#####################################################################
.rm.hrm.est.tau.item <- function( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1  , theta.k ,
				msteps, mstepconv , tau.item.fixed , prob.rater ){
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(c.rater)	
	Q0 <- matrix(0,nrow=VV, ncol=K)
	se.tau.item <- Q0
	cat("  M steps tau.item parameter   |")
	it <- 0 ;	conv1 <- 1000
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		tau.item0 <- tau.item
		for (kk in 1:K){
	#		kk <- 1
			Q1 <- Q0
			Q1[,kk] <- 1
			r1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR, 
					prob.item=NULL , prob.rater=prob.rater )
			pjk <- r1$prob.total				
			pjk1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item+h*Q1 ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL , prob.rater=prob.rater)$prob.total				
			pjk2 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item-h*Q1 ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL , prob.rater=prob.rater)$prob.total		
			# numerical differentiation			
			res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
					max.increment=max.b.increment , numdiff.parm )					
			increment <- Q1*matrix( res$increment , nrow=VV , ncol=K)	
			tau.item <- tau.item + increment
			se.tau.item[,kk] <- sqrt(abs(-1/res$d2)	)
					}
		conv1 <- max( abs( tau.item - tau.item0 ) )
		it <- it+1
		cat("-") # ; flush.console()
		if (!is.null(tau.item.fixed)){
			tau.item[ tau.item.fixed[,1:2,drop=FALSE] ] <- tau.item.fixed[,3]
								}
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()
	res <- list("tau.item" = tau.item , "se.tau.item" = se.tau.item , 
			"ll" = sum(res$ll0) , "prob.item"=r1$prob.item )
    return(res)
					}
#########################################################################
.rm.hrm.est.a.item <- function( c.rater , Qmatrix , tau.item ,
				VV , K , I , TP , a.item , d.rater , item.index , rater.index ,
				n.ik , numdiff.parm , max.b.increment=1,theta.k ,
				msteps, mstepconv , prob.rater ){
    h <- numdiff.parm
	diffindex <- item.index
	RR <- length(c.rater)
	cat("  M steps a.item parameter     |")
	it <- 0 ;	conv1 <- 1000	
	while( ( it < msteps ) & ( conv1 > mstepconv ) ){	
		a.item0 <- a.item	
			r1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL,prob.rater=prob.rater)
			pjk <- r1$prob.total				
			pjk1 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item+h , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL,prob.rater=prob.rater)$prob.total				
			pjk2 <- .rm.hrm.calcprobs( c.rater , Qmatrix , tau.item ,
					VV , K , I , TP , a.item-h , d.rater , item.index , rater.index , theta.k,RR,
					prob.item=NULL,prob.rater=prob.rater)$prob.total			
		# numerical differentiation			
		res <- .rm.numdiff.index( pjk , pjk1 , pjk2 , n.ik , diffindex , 
				max.increment=max.b.increment , numdiff.parm )					
		a.item <- a.item + res$increment
		a.item[ a.item < .05 ] <- .05
#		a.item <- a.item - mean(a.item ) + 1
		b1 <- mean( log( a.item ) )
		a.item <- a.item / exp( b1 )
		conv1 <- max( abs( a.item - a.item0 ) )
		it <- it+1
		cat("-") # ; flush.console()	
			}
	cat(" " , it , "Step(s) \n")	#; flush.console()	
    res <- list("a.item" = a.item , "se.a.item" = sqrt( abs(-1/res$d2 )) , 
			"ll" = sum(res$ll0) , "prob.item" = r1$prob.item )
    return(res)
			}			
###############################################################				