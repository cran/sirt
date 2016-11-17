

##########################################################################
# alternating least squares for Haberman linking
linking_haberman_als <- function(logaM , wgtM , maxiter , conv,
	   progress , est.type , cutoff )
{
	#****************
    iter <- 0
	parchange <- 1000
	NS <- base::ncol(logaM)
	NI <- base::nrow(logaM)
	#-- initial values study parameters
	logaAt <- base::rep(0,NS)
	At_inits <- TRUE
	if ( At_inits ){
		logaAt <- base::colSums( logaM * wgtM , na.rm=TRUE ) / 
						base::colSums( wgtM , na.rm=TRUE )
		logaAt <- logaAt - logaAt[1]						
	}
	wgtM0 <- wgtM
	wgt_adj <- 1 + 0*wgtM
	eps <- 1E-5
	wgtM <- wgtM + eps
	#*** begin algorithm
	while( ( parchange > conv ) & (iter < maxiter) ){
		logaAt0 <- logaAt
		
		#---------------
        # calculate average item parameter
		logaAt_M <- base::matrix( logaAt, nrow=NI, ncol=NS, byrow=TRUE)
		# logaj <- base::rowSums( ( logaM - logaAt_M ) * wgtM , na.rm=TRUE)
		logaM_adj1 <- logaM - logaAt_M 		
		# logaj <- base::rowSums( logaM_adj1 * wgtM , na.rm=TRUE) /
		#			base::rowSums( wgtM , na.rm=TRUE)								
		logaj <- weighted_rowMeans( mat = logaM_adj1 , wgt= wgtM )								
		# calculate adjusted mean slope
		logaMadj <- logaM - logaj
		res <- linking_haberman_als_residual_weights( logaj = logaj , logaAt = logaAt,
					logaM=logaM , cutoff=cutoff , wgtM0=wgtM0 , eps=eps )
		wgtM <- res$wgtM
		# logaAt <- base::colSums( logaMadj * wgtM , na.rm=TRUE ) / 
		# 				base::colSums( wgtM , na.rm=TRUE )
		logaAt <- weighted_colMeans( mat = logaMadj , wgt = wgtM )					
		logaAt[1] <- 0
		# logaAt <- logaAt - logaAt[1]		
		# logaj <- logaj - base::mean(logaj)
		
		#*** calculate residual and weight
		res <- linking_haberman_als_residual_weights( logaj = logaj , logaAt = logaAt,
					logaM=logaM , cutoff=cutoff , wgtM0=wgtM0 , eps=eps )
		loga_resid <- res$loga_resid			
		wgtM <- res$wgtM
		wgt_adj <- res$wgt_adj
	
		parchange <- base::max( base::abs( logaAt0 - logaAt )  )
		if (progress){
			base::cat( paste0( "** " , est.type , " estimation | Iteration " , iter  , " | " , 
				"Max. parameter change = " , base::round( parchange , 6 ) ) , "\n")
			utils::flush.console()
		}
		iter <- iter + 1

	}
	if (progress){
		base::cat("\n")	
	}	

	#------- summary of regression
	
	# residual SD
	selitems <- base::which( base::rowSums( 1 - base::is.na( loga_resid ) ) > 1 )
	
	#--- calculation of standard errors of regression coefficients
	if ( stats::sd(logaAt) < 1E-10 ){
		res <- base::list( vcov = 0*base::diag(NS-1) , se = base::rep(0,NS-1)  )
	} else {
		res <- linking_haberman_als_vcov( regr_resid = loga_resid , 
				   regr_wgt = wgtM , selitems = selitems ,
				   transf_pars = logaAt )
	}
	
	#--- item statistics
	item_stat <- data.frame( "study" = colnames(wgtM0) )	
	item_stat$N_items <- colSums( wgtM0 > 0 , na.rm=TRUE)
	item_stat$sumwgt_items <- colSums( wgt_adj  , na.rm=TRUE )
	#-------
	#*** end algorithm
	res <- base::list( "logaAt"=logaAt , "logaj" = logaj ,
				"loga_resid" = loga_resid , "loga_wgt" = wgtM ,
				"loga_wgt_adj" = wgt_adj ,
				"vcov" = res$vcov , "se" = base::c( NA , res$se) ,
				item_stat = item_stat )
	base::return(res)
}
##########################################################################
