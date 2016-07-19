
###########################################
# confidence interval
mcmc_confint <- function( mcmcobj , parm , level = .95 , 
			exclude="deviance" ){	
	mcmcobj <- mcmcobj[ , ! ( colnames(mcmcobj) %in% exclude ) ]
	if ( ! base::missing(parm) ){
		mcmcobj <- mcmcobj[,parm]
			}
	q1 <- ( 1 - level ) / 2		
	h1 <- base::apply( mcmcobj , 2 , stats::quantile , q1 )
	q2 <- 1 - ( 1 - level ) / 2		
	h2 <- base::apply( mcmcobj , 2 , stats::quantile , q2 )	
	res <- base::data.frame( h1 , h2)	
	colnames(res)[1] <- base::paste0( base::round( 100*q1 ,1 ) , " %")
	colnames(res)[2] <- base::paste0( base::round( 100*q2 ,1 ) , " %")
	rownames(res) <- colnames(mcmcobj)
 	base::return(res)
		}
###############################################		