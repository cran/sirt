
tetrachoric2 <- function( dat , delta=.007 , maxit = 1000000 ,
	cor.smooth=TRUE ){
	# process data
	dat <- as.matrix(dat)
	# calculate tau
	tau <- - qnorm( colMeans(dat,na.rm=TRUE ) )
	dat.resp <- 1 - is.na(dat)
	dat[ dat.resp==0] <- 9
	I <- ncol(dat)
	# calculate frequencies
	dfr <- data.frame( "item1" = rep(1:I,I) , "item2" = rep(1:I, each=I ) )
	h1 <- t( ( dat==1 ) ) %*% ( dat==0 )
	dfr$f11 <- matrix( t( ( dat==1 ) ) %*% ( dat==1 ) , ncol=1 , byrow=T ) 
	dfr$f10 <- matrix( t( ( dat==1 ) ) %*% ( dat==0 ) , ncol=1 , byrow=T ) 
	dfr$f01 <- matrix( t( ( dat==0 ) ) %*% ( dat==1 ) , ncol=1 , byrow=T ) 
	dfr$f00 <- matrix( t( ( dat==0 ) ) %*% ( dat==0 ) , ncol=1 , byrow=T ) 
	dfr$ftot <- dfr$f11 + dfr$f10 + dfr$f01 + dfr$f00
	dfr$p11 <- dfr$f11 / dfr$ftot
	dfr$pi1 <- ( dfr$f11 + dfr$f10 ) / dfr$ftot
	dfr$pi2 <- ( dfr$f11 + dfr$f01 ) / dfr$ftot
	# subdata of dfr
	dfr <- dfr[ dfr$item1 > dfr$item2 , ]
	dfr <- dfr[ dfr$ftot > 0 , ]
	dfr$qi1 <- qnorm( dfr$pi1)
	dfr$qi2 <- qnorm( dfr$pi2)
	# functions defined by Cengiz Zopluoglu   
	L <- function(r,h,k) {(1/(2*pi*sqrt(1-r^2)))*exp(-((h^2-2*h*k*r+k^2)/(2*(1-r^2))))}
	S <- function(x,r,h,k) { x-(L(r,h,k)*delta) }
	A0 <- dfr$A0 <- dfr$p11 - dfr$pi1 * dfr$pi2 
	dfr$r0 <- delta / 2
	dfr$iter <- 0

	dfr$iter <- 0
	dfr$conv <- 0
	ii <- 0

	vars <-  c("A0","r0","iter","conv")
	cat("    0 : " )
	while( ( mean( dfr$conv) < 1 ) & ( ii < maxit ) ){
		# iterations
		dfr0 <- dfr
		#***
		ii <- ii+1
		dfr$A0 <- dfr$A0 - delta * L( dfr$r0 , dfr$qi1 , dfr$qi2 )
		dfr$r0 <- dfr$r0 + delta
		dfr$iter <- ii
		ind <- which( dfr$A0 < 0 )
		if ( length(ind) > 0 ){ 
			h1 <- dfr$r0 - delta/2 + dfr$A0 / L( dfr0$r0 , dfr$qi1 , dfr$qi2 )
			dfr$r0[ind] <- h1[ind]
			dfr$conv[ind] <- 1
							}
		i2 <- which( dfr$conv==1 & dfr0$conv==1 )
		if (length(i2) > 0){    dfr[ i2 , vars] <- dfr0[ i2 , vars] }
		if (ii %% 100 == 0 ){ cat(".") ; flush.console() }
		if (ii %% 1000 == 0 ){ cat("\n" ,ii, ": ") }    
		# print(dfr) ; flush.console()
			}
	cat("\n")
	TC <- matrix(NA , I , I )
	diag(TC) <- 1
	TC[ as.matrix(dfr[ , c("item1","item2") ] ) ] <- dfr$r0
	TC[ as.matrix(dfr[ , c("item2","item1") ] ) ] <- dfr$r0
	if (cor.smooth){ 
		library(psych)
		TC <- cor.smooth(TC) 
			}
	res <- list("tau"=tau , "rho" = TC )
	return(res)
	}