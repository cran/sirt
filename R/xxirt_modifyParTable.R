
################################################################
# modify parameter table
xxirt_modifyParTable <- function( partable , parname , item = NULL , value=NULL ,
		est = NULL , parlabel = NULL , parindex = NULL , lower=NULL ,
		upper = NULL , prior=NULL ,
		prior_par1 = NULL , prior_par2 = NULL
		){
	if ( is.null(item) ){
		item <- base::paste( base::unique(partable$item ))
	}
	ind <- base::which( (partable$parname %in% parname) & (partable$item %in% item ) )	
	vv <- value ; vv_name <- "value"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}
	vv <- est ; vv_name <- "est"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}
	vv <- lower ; vv_name <- "lower"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	vv <- upper ; vv_name <- "upper"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}			  	  	  
	vv <- prior ; vv_name <- "prior"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	vv <- prior_par1 ; vv_name <- "prior_par1"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	vv <- prior_par2 ; vv_name <- "prior_par2"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	vv <- parlabel ; vv_name <- "parlabel"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	vv <- parindex ; vv_name <- "parindex"	
	  if ( ! base::is.null(vv) ){ partable[ ind , vv_name ] <- vv	}		
	#*** output
	base::return(partable)	
}
##################################################		