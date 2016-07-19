
amh_compare_estimators <- function(res_MAP , res_Mean , res_pmle ){
	
	pars0 <- res_MAP$pars
	NP <- base::length(pars0)	
	dfr <- base::matrix( 0 , nrow=NP+3 , ncol=3)
	colnames(dfr) <- base::c("MAP" , "mMAP" , "Mean")
	
	
	res_pmle$ll <- res_pmle$loglik
	res_pmle$priorval <- res_pmle$logprior
	
	for (cc in 1:3){			
		# cc <- 1
		res0 <- base::switch( cc , 
					res_MAP ,
					res_pmle ,
					res_Mean )		
		dfr[(1:NP)+3,cc] <- res0$pars
		dfr[1,cc] <- res0$ll
		dfr[2,cc] <- res0$priorval
		dfr[3,cc] <- res0$posteriorval
	}
	
	dfr <- base::data.frame( "parm" = base::c( base::c("Log Likelihood" , 
				"Log Prior", "Log Posterior") , base::names(pars0) ) , dfr )		
	base::return(dfr)
}