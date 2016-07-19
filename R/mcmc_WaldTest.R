
##########################################################
# Wald Test for a set of hypotheses
mcmc_WaldTest <- function( mcmcobj , hypotheses ){
	NH <- base::length(hypotheses)
	n1 <- base::ncol(mcmcobj)
	mcmcobj <- mcmc_derivedPars( mcmcobj , hypotheses)
	n2 <- base::ncol(mcmcobj)
	mcmcobj <- mcmcobj[ , base::seq(n1+1,n2) ]
	v1 <- mcmc_vcov(mcmcobj)
	s1 <- mcmc_summary(mcmcobj)
	c1 <- s1$MAP
	# compute test statistic
	W <- base::t(c1) %*% base::solve(v1) %*% c1
	stat <- base::c( "chi2" = W , "df" = NH)	
	stat["p"] <- 1 - stats::pchisq( W , df = stat["df"])
	res <- base::list( "hypotheses_summary" = s1 ,
			  "chisq_stat" = stat
					)
	base::class(res) <- "mcmc_WaldTest"
	base::return(res)	
		}
##############################################################		
# summary of Wald Test based on MCMC output
summary.mcmc_WaldTest <- function( object , digits = 3 , ... ){
    base::cat("Wald Test\n")
	W1 <- base::sprintf( paste0("%." , digits , "f" ) , object$chisq_stat["chi2"] )

	v1 <- base::paste0("Chi^2 = " ,  W1 , ", df = " , object$chisq_stat["df"])
	v1 <- base::paste0( v1 , ", p = " , base::sprintf( paste0("%." , digits , "f" ) , 
					object$chisq_stat["p"] ) )
	base::cat(v1)

	base::cat("\n\nSummary Hypotheses\n")	
	obji <- object$hypotheses_summary
	vars <- base::c("parameter","MAP","SD", "Q2.5", "Q97.5" , "Rhat","SERatio",
					"effSize" )
	obji <- obji[,vars]
	NO <- base::ncol(obji)
	for (vv in 2:NO ){
		obji[,vv] <- base::round( obji[,vv] , digits )
		}
	obji[,NO] <- base::round( obji[,NO] )	
	base::print(obji)
			}
##################################################################			