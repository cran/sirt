

########################################################
# mcmclist descriptives
mcmc_summary <- function( mcmcobj , quantiles=c(.025,.05,.50,.95,.975) ){
 	summary.mcmcobj <- summary(mcmcobj , quantile=quantiles)	
    dat.bugs <- mcmcobj  # [[1]]
    vars <- base::colnames(dat.bugs)
    n.iter <- base::nrow(dat.bugs)
        discret <- base::round( base::seq( 1 , n.iter , length= 4 ) )
        MAP <- Rhat <- base::rep(0,length(vars) )
        for (ii in base::seq( 1 , length(vars )) ){
            vv <- vars[ii]
            dat.vv <- base::as.vector( dat.bugs[ , vv ] )
            # extract 3 subchains
            l1 <- dat.vv[ discret[1]:discret[2] ]
            l2 <- dat.vv[ (discret[2]+1):discret[3] ]
            l3 <- dat.vv[ (discret[3]+1):discret[4] ]
            W <- ( stats::var(l1) + stats::var(l2) + stats::var(l3) ) / 3 
            S1 <- n.iter / 3
            est.chains <- c( base::mean(l1) , base::mean(l2) , base::mean(l3) )
            B <-  S1 / 2 * base::sum( ( est.chains - base::mean( est.chains ) )^2 ) 
            Rhat[ii] <- ( ( S1-1 ) / S1 * W + B / S1 ) / W
            # mode estimation   
            m1 <- stats::density( dat.vv , from = base::min(dat.vv) ,
              			   to = base::max(dat.vv) )
            MAP[ii] <- m1$x[ base::which( m1$y  == base::max( m1$y) ) ]         
                            }						
    res <- base::data.frame( "MAP" = MAP , "Rhat" = Rhat )
	base::rownames(res) <- vars
	smc3  <- res
	smc2 <- summary.mcmcobj$statistics
	base::colnames(summary.mcmcobj$quantiles) <- base::paste0( "Q" , 100*quantiles )
	# calculate effective sample size
	effSize <- coda::effectiveSize( mcmcobj )	
	statis <- summary.mcmcobj$statistics
	statis <- base::cbind( statis[ , c(1,2) ] , 
				base::apply( as.matrix(mcmcobj) , 2 , stats::mad ) , 
				base::apply( as.matrix(mcmcobj) , 2 , skewness.sirt ) ,
				statis[,c(3,4) ]	)
	base::colnames(statis)[3:4] <- c("MAD" , "skewness" )
	dfr <- base::data.frame( "parameter" = base::rownames(smc3) , 
				statis , smc3 , "SERatio" = smc2[,4] / smc2[,2] ,
				"sampSize" = base::nrow(as.matrix(mcmcobj)) , "effSize" = effSize , 
				summary.mcmcobj$quantiles )	
    base::rownames(dfr) <- NULL		
    base::return(dfr)
        }
###########################################		
