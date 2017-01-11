
##############################################################
lsem.fitsem <- function( dat , weights , lavfit ,
			fit_measures , NF , G , moderator.grid , verbose ,
			pars , standardized ){

    parameters <- NULL
	fits <- NULL 	
	pars0 <- pars	
	
	if (verbose){
		base::cat( "** Fit lavaan model\n")
		G1 <- base::min(G,10)	
		pr <- base::round( base::seq(1,G , len=G1) )
		base::cat("|")
		base::cat( base::paste0( base::rep("*",G1) , collapse="") )
		base::cat("|\n")
		base::cat("|")
	}

	for (gg in 1:G){
		# gg <- 1			
		dat$weight <- weights[,gg]
		datsvy <- survey::svydesign(id=~index,   weights=~weight ,    data=dat)
		# fit the model using weighted data
		survey.fit <- lavaan.survey::lavaan.survey(lavaan.fit=lavfit, 
							survey.design=datsvy )
		dfr.gg <- pars <- lavaan::parameterEstimates(survey.fit) 		
		
		if (standardized){			
			sol <- lavaan::standardizedSolution( survey.fit )
			colnames(sol)[ base::which( colnames(sol) == "est.std" ) ] <- "est"
			sol$lhs <- paste0( "std__" , sol$lhs)
			pars <- plyr::rbind.fill( pars , sol )	
			dfr.gg <- pars
		} 							
		pars <- base::paste0( pars$lhs , pars$op , pars$rhs )					
		NP <- base::length(pars0)
		ind <- base::match( pars0 , pars )
		dfr.gg <- dfr.gg[ ind , ]
		dfr.gg <- base::data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
						  "par"= pars0 , "parindex" = 1:NP , dfr.gg	)
		dfr.gg0 <- base::data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
						  "par"= fit_measures , "parindex" = NP + 1:NF , 
						  "est"= lavaan::fitMeasures(survey.fit , fit.measures= fit_measures ) ,
						  "op"="fit" )
		vars <- base::setdiff( colnames(dfr.gg) , colnames(dfr.gg0) )
		for (vv in vars){ dfr.gg0[,vv] <- NA }
		dfr.gg <- base::rbind( dfr.gg , dfr.gg0[ , colnames(dfr.gg) ] )		
		parameters <- base::rbind( parameters , dfr.gg ) 
		# fits <- rbind( fits , dfr.gg ) 
		if (verbose){
			if ( gg %in% pr ){
				base::cat("-")
				utils::flush.console()
			}
		}
	}
	if (verbose){
		base::cat("|\n")
		utils::flush.console()
	}

	parameters <- parameters[ base::order(parameters$parindex) , ]	
#	fits <- fits[ order(fits$fitindex) , ]	
#	rownames(fits) <- NULL
				
	res <- base::list( "parameters" = parameters ) #  , "fits" = fits )
	base::return(res)	
			}
#######################################################################			