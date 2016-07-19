
#######################################################
xxirt_postproc_parameters <- function( partable , customTheta , 
		items , probs_items ){	
	#**** item parameters
	p1 <- partable[ partable$parfree == 1 , ]
	par_items <- p1$value
	base::names(par_items) <- p1$parlabel
	#*** theta distribution parameters
	cs <- customTheta 
	par_Theta <- cs$par[ cs$est ]
	#*** structured form of parameters
	I <- base::length(items)
	parnames <- base::sort( base::unique( base::paste( partable$parname) ) )
    PN <- base::length(parnames) 
	m1 <- base::matrix(NA , nrow=I , ncol=PN)
	base::rownames(m1) <- items
	base::colnames(m1) <- parnames
	for (pp in 1:PN){
		p1 <- partable[ partable$parname == parnames[pp]  , ]
		m1[ p1$item , parnames[pp] ] <- p1$value
	}	
	p1 <- partable[ ! base::duplicated(partable$item ) , ]
	dfr <- base::data.frame( "item" = items , "type" = base::paste(p1$type) , m1  )
	base::rownames(dfr) <- NULL

	#*** probabilities item parameters
	pi_dim <- base::dim(probs_items)
	base::dimnames(probs_items)[[1]] <- items
	base::dimnames(probs_items)[[2]] <- base::paste0("Cat", base::seq(0,pi_dim[2]-1) )	
	#*** lower and upper bounds
	p1 <- partable[ partable$parfree == 1 , c("item" , "type" , "parname" ,
				"value" , "lower" , "upper" , "parindex" , "parlabel" ) ]
	p1$active <- 1 * ( p1$value > p1$lower )
	p1$active <- p1$active * ( p1$value < p1$upper )
	par_items_bounds <- p1	
	
	#*** output
	res <- base::list( par_items = par_items , par_Theta = par_Theta ,
	              probs_items = probs_items , par_items_summary = dfr ,
				  par_items_bounds = par_items_bounds )
    base::return(res)
}
#######################################################	