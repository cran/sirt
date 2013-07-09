
########################################################
# mcmc.list descriptives
mcmc.list.descriptives <- function( mcmcobj , quantiles=c(.025,.05,.1,.9,.95,.975) ){
    library(coda)
 	summary.mcmcobj <- summary(mcmcobj , quantile=quantiles)	
    dat.bugs <- mcmcobj[[1]]
    vars <- colnames(dat.bugs)
    n.iter <- nrow(dat.bugs)
        discret <- round( seq( 1 , n.iter , length= 4 ) )
        MAP <- Rhat <- rep(0,length(vars) )
        for (ii in seq( 1 , length(vars ))){
            vv <- vars[ii]
            dat.vv <- as.vector( dat.bugs[ , vv ] )
            # extract 3 subchains
            l1 <- dat.vv[ discret[1]:discret[2] ]
            l2 <- dat.vv[ (discret[2]+1):discret[3] ]
            l3 <- dat.vv[ (discret[3]+1):discret[4] ]
            W <- ( var(l1) + var(l2) + var(l3) ) / 3 
            S1 <- n.iter / 3
            est.chains <- c(mean(l1) , mean(l2) , mean(l3) )
            B <-  S1 / 2 * sum( ( est.chains - mean( est.chains ) )^2 ) 
            Rhat[ii] <- ( ( S1-1 ) / S1 * W + B / S1 ) / W
            # mode estimation   
            m1 <- density( dat.vv , from = min(dat.vv) , to = max(dat.vv) )
            MAP[ii] <- m1$x[ which( m1$y  == max( m1$y) ) ]         
                            }						
    res <- data.frame( "MAP" = MAP , "Rhat" = Rhat )
	rownames(res) <- vars
	smc3  <- res
	smc2 <- summary.mcmcobj$statistics
	colnames(summary.mcmcobj$quantiles) <- paste0( "Q" , 100*quantiles )
	dfr <- data.frame( "parameter" = rownames(smc3) , 
		summary.mcmcobj$statistics , smc3 , "PercSEratio" = 100*smc2[,4] / smc2[,2] ,
		summary.mcmcobj$quantiles )	
    return(dfr)
        }