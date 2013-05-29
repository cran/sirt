 
# 0.01  2012-xx-yy


# 0.01  2012-06-23  o initial release
# 0.02  2012-06-24  o corrected a bug in summary table of rasch.pairwise 
# 0.03  2012-06-24  o included function rasch.pairwise.itemcluster
# 0.04  2012-07-14  o worked on a bug:
#						PL Iter. 1 : max. parm. change =  NA 
#						Fehler in while (max.change > conv) { : 
#						Fehlender Wert, wo TRUE/FALSE nötig ist



# 0.0x  2012-0x-yy
#-------------------------------------------------------




#---------------------------------------------------------------#
# Rasch estimation (Approximate Method)                         #
# also called MINCHI method                                     #
# Handbook of Statistics Vol. 26                                #
# Chapter of G. Fischer: p. 544 ff.                             #
# pairwise likelihood method
##NS export(rasch.pairwise)
rasch.pairwise <- function( dat , conv = 0.0001 , maxiter = 3000 , progress = TRUE , b.init = NULL ){
    # INPUT:                                               #
    # dat      ... data frame                              #
    # conv     ... convergence in epsilon parameters       #
    # maxiter  ... maximal number of iterations            #
    # progress ... display progress in MINCHI method       #
    # beta.init ... initial beta values                    #
    s1 <- Sys.time()
        # should items being excluded?
        item.elim <- which( colMeans( dat , na.rm=T ) %in% c(0,1))
        if (length(item.elim)>0){ dat <- dat[ , - item.elim ] }

        dat <- as.matrix(dat)
        # data preparation
        dat.resp <- 1 - is.na( dat )
        dat.9 <- dat
        dat.9[ is.na(dat) ] <- 9
        # calculate n_{ij}
        n.ij <- t( dat.9 * dat.resp ) %*%  ( ( 1 - dat.9 ) * dat.resp  )
        # which item pairs occur in estimation procedure
        delta.ij <- 1 * ( n.ij + t( n.ij ) > 0 ) 

        # initial values for beta
        if( is.null( b.init) ){ beta <- - qlogis( colMeans( dat , na.rm=T ) ) } else { beta <- b.init } 
        # calculate y_{ij} values
        y.ij <- n.ij / ( n.ij + t( n.ij) )
        y.ij[ delta.ij == 0 ] <- 0
        y.ji <- t( y.ij )

        eps <- exp( - beta )
        change <- 1
        iter <- 0
        while( change > conv & iter < maxiter ){
                eps0 <- eps
                eps <- sqrt( rowSums( y.ij * eps * delta.ij ) / colSums( y.ij * 1 / eps ) )
                change <- max( abs( eps0 - eps ) ) 
                iter <- iter + 1
                if ( progress ){
                    cat( "PL Iter." , iter , ": max. parm. change = " , 
                                round( max(abs( log(eps0) - log(eps))) , 6 ) , "\n")
                    flush.console()
                        }                
                }
        item <- data.frame( "N" = colSums(1 -is.na(dat)) , "p" = colMeans( dat , na.rm=T ) , 
#                                "itemdiff" = scale( - log(eps) , scale=F ) 
						"b" = - log(eps) ,
						"itemcluster"= rep(0,ncol(dat))
									)
        s2 <- Sys.time()									
        res <- list( "b" = - log( eps ) , "eps" = eps , "iter" = iter , "conv" = conv , "dat" = dat ,
                    "freq.ij" = n.ij  , "item" = item , "fct" = "rasch.pairwise",
					"s1"=s1 , "s2"=s2  ) 
        class(res) <- "rasch.pairwise"
        return(res)
    }    
#-------------------------------------------------------------------




#**********************************************************
# Summary for rasch.minchi object                         *
##NS S3method(summary,rasch.pairwise)
summary.rasch.pairwise <- function( object , ...){
    cat("------------------------------------------- \n")
    d1 <- packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat("  Function" , object$fct , "\n") 
    cat("------------------------------------------- \n")
    cat("Pairwise likelihood estimation \n")
    cat("Rasch Model \n")
    cat("------------------------------------------- \n")
    cat("Item Parameters \n")
    print( round( object$item , 3 ))                
                }
#*******************************************************





############################################################
# Pairwise estimation with itemclusters
##NS export(rasch.pairwise.itemcluster)
rasch.pairwise.itemcluster <- function( dat , itemcluster = NULL ,
			conv = .00001 , maxiter = 3000 , progress = TRUE , b.init = NULL){
    if ( is.null(b.init) ){ 
#			b.init <- - qlogis( colMeans( dat , na.rm=T ) ) 
			b.init <- - qlogis( colMeans( dat , na.rm=T ) ) 
				}
	s1 <- Sys.time()
    I <- ncol(dat)
	dat <- as.matrix(dat)
	dat0 <- dat
	dat[ is.na(dat) ] <- 9
    b <- b.init
    # create count tables
    Aij <- t( dat == 0 ) %*% ( dat == 1 )
    Aji <- t( dat == 1 ) %*% ( dat == 0 )
	# set some entries to zero for itemclusters
	clusters <- unique( itemcluster[ itemcluster != 0 ] )
	CC <- length(clusters)
	for (cc in clusters){
		icc <- which( itemcluster == cc )
		Aji[icc,icc] <- Aij[icc,icc] <- 0
				}
    nij <- Aij + Aji
    eps0 <- eps <- exp(  b )
    max.change <- 10
    iter <- 1
    while( max.change > conv ){
        b0 <- b
		eps0 <- eps
#		g1 <- sum( nij[ii,] * ( eps0[ii] + eps0 )^(-1) )	
		m1 <- matrix( eps0 , I , I , byrow=TRUE ) + matrix( eps0 , I , I ) 
		g1 <- rowSums( nij / m1 )
		eps <- rowSums( Aij ) / g1 
        b <-  log(  eps )
        max.change <- max(abs( b - b0 ))
        if ( progress ){
                cat( "PL Iter." , iter , ": max. parm. change = " , 
                        round( max.change , 6 ) , "\n")
                flush.console()
                    } 
        iter <- iter + 1               
                }				
        item <- data.frame( "N" = colSums(1 -is.na(dat0)) , "p" = colMeans( dat0 , na.rm=T ) , 
                        "b" =  log(eps) )
		if ( is.null(itemcluster) ){ itemcluster <- rep(0,I) }
		item$itemcluster <- itemcluster
		s2 <- Sys.time()									
        res <- list( "b" = b , "eps" = eps , "iter" = iter , "conv" = conv , "dat" = dat0 ,"item" = item ,
			"fct" = "rasch.pairwise.itemcluster" , "s1"=s1 , "s2"=s2  ) 
        class(res) <- "rasch.pairwise"
        return(res)
       }
#-------------------------------------------------------------------
