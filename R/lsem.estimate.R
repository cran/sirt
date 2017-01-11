
#################################################
# estimate LSEM model

lsem.estimate <- function( data , moderator , moderator.grid ,
		lavmodel , type="LSEM" , h = 1.1 , 
		residualize=TRUE, fit_measures = c("rmsea","cfi","tli","gfi","srmr"),
		standardized = FALSE , 
		standardized_type = "std.all" ,
		eps=1E-8 , verbose=TRUE , ... ){

		
	CALL <- base::match.call()
	s1 <- base::Sys.time()	
		
	lavaan.args <- base::list(...)	
	if (standardized){
		if ( type == "MGM"){
			base::stop("standardized=TRUE cannot be applied for type='MGM'")
		}
	}
	# group moderator if type="MGM"
	out <- lsem.group.moderator( data , type , moderator.grid , moderator ,
				residualize , h)
    data <- out$data
	moderator.grouped <- out$moderator.grouped
	h <- out$h
	residualize <- out$residualize
	moderator.grid <- out$moderator.grid
	
	# residualize input data	
	out <- lsem.residualize( data , moderator , moderator.grid ,
				lavmodel , h , residualize , eps , verbose )		
	G <- out$G
    data <- out$data	
	weights <- out$weights
    data$index <- base::seq(1,nrow(data))	
	residualized_interceps <- out$residualized_interceps
	# unweighted fit of lavaan model
	dat <- data
	lavfit <- lavaan::sem(model = lavmodel, data=dat,  ... )
	fM <- lavaan::fitMeasures( lavfit )
	fit_measures <- base::intersect( fit_measures , base::names(fM) )
	NF <- base::length(fit_measures)
	pars <- lavaan::parameterEstimates(lavfit)
 	if (standardized){			
		sol <- lavaan::standardizedSolution( lavfit , type=standardized_type)
		colnames(sol)[ base::which( colnames(sol) == "est.std" ) ] <- "est"
		sol$lhs <- base::paste0( "std__" , sol$lhs)
		pars <- plyr::rbind.fill( pars , sol )	
	} 
	pars <- base::apply( pars[ , c("lhs" , "op" , "rhs" ) ] , 1 , FUN = function(ll){
				base::paste0( ll[1] , ll[2] , ll[3] ) } )
	# fit LSEM for all moderator groups
	out2 <- lsem.fitsem( dat , weights , lavfit ,
			  fit_measures , NF , G , moderator.grid , verbose , pars ,
			  standardized )	
	parameters <- out2$parameters
	# fits <- out2$fits
			
	rownames(parameters) <- base::paste0( parameters$par ,
						"__" , parameters$grid_index )			
		
	#****************************
	# parameter and fit statistics summary
	parameters_summary <-  lsem.parameter.summary( parameters , 
		     moderator.density=out$moderator.density , verbose )
	out$moderator.density$Neff <- base::colSums(weights)
	
	obji0 <- obji <- out$moderator.density	
	obji$moderator <- obji$moderator 
	obji$wgt <- obji$wgt
	obji$Neff <- obji$Neff
	dfr <- base::data.frame( "M" = base::colMeans( obji0[,-1] ) , 
			"SD"= base::apply( obji0[,-1] , 2 , stats::sd ) ,
			"min" = base::apply( obji0[,-1] , 2 , base::min ) ,
			"max" = base::apply( obji0[,-1] , 2 , base::max ) 
					)		
	dfr0 <- base::data.frame("M"= base::mean( data[,moderator] , na.rm=TRUE ) , 
				"SD" = out$sd.moderator ,
				"min"= base::min( data[ , moderator ] , na.rm=TRUE ) ,
				"max"= base::max( data[ , moderator ] , na.rm=TRUE )  
							)										
	obji <- base::rbind( dfr0 , dfr )		
	rownames(obji) <- NULL
	moderator.stat <- base::data.frame("variable"=c("moderator" ,
				"wgt" , "Neff") , obji )
		
	# output
	s2 <- base::Sys.time()	
	res <- base::list( "parameters"=parameters , "weights"=weights , 				 
				 "parameters_summary" = parameters_summary , 
				 "bw"=out$bw , "h"=h , "N"=out$N , 
				 "moderator.density"=out$moderator.density , 
				 "moderator.stat" = moderator.stat , 
				 "moderator.grouped" = moderator.grouped , 
				 "m.moderator" = base::mean( data[,moderator] , na.rm=TRUE ) ,
				 "sd.moderator"=out$sd.moderator , "moderator"=moderator ,
				 "moderator.grid" = moderator.grid ,
				 "lavmodel"=lavmodel , "residualize"=residualize ,
				 "data"=data , "residualized.intercepts" = residualized_interceps , 
				 "lavaan.args"=lavaan.args ,
				 "fit_measures"=fit_measures , "s1"=s1 , "s2"=s2 ,
				 "standardized" = standardized , 
				 "standardized_type" = standardized_type , 
				 "type"=type , "CALL" = CALL )	
	base::class(res) <- "lsem"	
	base::return(res)	
		
}