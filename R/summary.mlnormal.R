#*******************************************************
# Summary for mlnormal object
summary.mlnormal <- function( object , digits = 4 , file=NULL , ...){

    # open sink
    CDM::osink( file = file , suffix = paste0( file, "__SUMMARY.Rout") )

	base::cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	base::cat( base::paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	base::cat( "Date of Analysis:" , base::paste( object$s2 ) , "\n" )
	base::cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
	
	base::cat("Call:\n", base::paste( base::deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
	
	base::cat( object$descriptions["des_method"] , "\n\n")
	
#	modeltype <- object$irtmodel
#		base::cat( "   " , object$ic$n , "Cases, " , object$ic$I , "Items, " , 
#		        object$G , "Group(s)", # "," ,
#				"\n")  
	
	if ( object$REML & ! ( object$vcov_REML ) ){
		base::cat(" * Standard errors for random effects are not properly calculated! \n")
		base::cat(" * Use vcov = TRUE as an argument in 'mlnormal' if valid standard errors\n")
		base::cat("   for random effects parameters are requested.\n\n")
	}
	
    base::cat("-----------------------------------------------------------------\n")
	
	base::cat( "Number of observations =" , object$ic$N , "\n" )	
	base::cat( "Number of clusters =" , object$ic$G , "\n\n" )	
	
	base::cat( "Number of iterations =" , object$iter , "\n\n" )
    base::cat( "Deviance = " , base::round( object$deviance , 2 ) , "\n" )
    base::cat( base::paste0( object$descriptions["log_like_verbose"] 
					, " =" ) , base::round( object$ic$loglike , 2 ) , "\n" )	
	
	#---- print posterior
	if ( object$posterior_obj$display_posterior ){
		base::cat( "Log prior =" , 
				base::round( object$posterior_obj$log_prior , 2 ) , "\n" )	
		base::cat( "Log posterior =" , 
				base::round( object$posterior_obj$log_posterior , 2 ) , "\n" )		
	}
	
	base::cat("\n")	
    base::cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    base::cat( "  Number of estimated beta parameters = " , object$ic$np.beta , 
				"\n" )    	
    base::cat( "  Number of estimated theta parameters = " , object$ic$np.theta , 
				"\n\n" )    
	
	if ( ! object$REML ){	
		base::cat( "AIC  = " , base::round( object$ic$AIC , 1 ) , " | penalty =" , 
				base::round( object$ic$AIC - object$ic$deviance ,2 ) , 
				"   | AIC = -2*LL + 2*p  \n" )    
	}
				
#    base::cat( "AICc = " , round( object$ic$AICc , 0 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
#		base::cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
#    base::cat( "BIC  = " , round( object$ic$BIC , 0 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
#			"   | BIC = -2*LL + log(n)*p  \n" )  
#    base::cat( "CAIC = " , round( object$ic$CAIC , 0 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
#		base::cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

    base::cat("-----------------------------------------------------------------\n")
	base::cat("Beta Parameters\n")
	excl <- base::c("parm","prior")
	
	obji <- object$beta_summary
	a <- summary_round_helper(obji, digits=digits, exclude = excl, print=TRUE)

    base::cat("-----------------------------------------------------------------\n")
	base::cat("Theta Parameters\n")
	obji <- object$theta_summary
	a <- summary_round_helper(obji, digits=digits, exclude = excl , print=TRUE)
	
	# close sink
    CDM::csink( file = file )		
}
#*******************************************************
