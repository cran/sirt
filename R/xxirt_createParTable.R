
##################################################
# create parameter table
xxirt_createParTable <- function( dat , itemtype , customItems=NULL ){
	I <- base::ncol(dat)
	ncat1 <- base::apply( dat , 2 , base::max , na.rm=TRUE ) + 1
	items <- base::colnames(dat)
	if ( base::length(itemtype) == 1 ){
		itemtype <- base::rep( itemtype , I )
	}
	dfr <- NULL
	CI <- base::length(customItems)	
	for (ii in 1:I){	
		# ii <- 1
		type_ii <- itemtype[ii]
		item_ii <- NULL
		for ( vv in 1:CI){
			ci_ii <- customItems[[vv]]
			if ( ci_ii$name == type_ii ){
				item_ii <- ci_ii
			}
		}
		if ( base::is.null(item_ii) ){
				base::stop( paste0( "Item type " , type_ii , " not found!") )
		}
		NP <- base::length( item_ii$par )
		dfr1 <- base::data.frame( "item" = base::rep( items[ii] , NP ) )		
		dfr1$itemnr <- ii
		dfr1$ncat <- ncat1[ii]		
		dfr1$class <- base::class(item_ii) 
		dfr1$type <- item_ii$name
		dfr1$parname <- base::names(item_ii$par)
		dfr1$value <- item_ii$par
		dfr1$est <- item_ii$est
		dfr1$lower <- item_ii$lower
		dfr1$upper <- item_ii$upper	
		dfr1$prior <- NA
		dfr1$prior_par1 <- NA
		dfr1$prior_par2 <- NA
		if ( ! base::is.null( item_ii$prior ) ){
			item_ii_prior <- base::names(item_ii$prior)
			ind_ii <- base::match( item_ii_prior , dfr1$parname )
			dfr1[ ind_ii , "prior" ] <- item_ii$prior
			dfr1[ ind_ii , "prior_par1" ] <- item_ii$prior_par1
			dfr1[ ind_ii , "prior_par2" ] <- 	item_ii$prior_par2
		}			
		dfr <- base::rbind( dfr , dfr1 )
	}
	#**** create parameter indices
	NP <- base::nrow(dfr)
	dfr$rowindex <- 1:NP
	# parameter index
	dfr$parindex <- base::cumsum( dfr$est )
	#*** parameter label
	dfr$parlabel <- base::paste0( dfr$item , "_" , dfr$parname )	
	base::attr(dfr , "ncat" ) <- ncat1
	base::attr(dfr , "items" ) <- items
	base::return(dfr)	
}
##################################################		

	#**** create parameter indices
#	NP <- nrow(dfr)
#	dfr$rowindex <- 1:NP
#	# parameter index
#	dfr$parindex <- cumsum( dfr$est )
#	#*** parameter label
#	dfr$parlabel <- paste0( dfr$item , "_" , dfr$parname )	
#	attr(dfr , "ncat" ) <- ncat1
#	attr(dfr , "items" ) <- items