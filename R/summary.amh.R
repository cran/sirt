
#############################################################
summary.amh <- function( object , digits=3 , file=NULL ,... ){

    # open sink
    CDM::osink( file = file , suffix = paste0( file, "__SUMMARY.Rout") )

	cat("-----------------------------------------------------------------\n")
	d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	
	cat( object$description , "\n\n")
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")			
	
	cat( "Date of Analysis:" , "\n" )
	cat( "   Start:" , paste( object$time$start ) , "\n" )	
	cat( "   End  :" , paste( object$time$end ) , "\n" )		
	cat("Computation time:" , print(object$time$end - object$time$start), "\n\n")    

	cat( "Number of iterations =" , object$n.iter , "\n" )
	cat( "Number of burnin iterations =" , object$n.burnin , "\n" )			
	cat( "Number of saved iterations =" , object$n.saved , "\n\n" )			

    cat("-----------------------------------------------------------------\n")
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n\n" )    
												
    cat( "AIC  = " , round( object$ic$AIC , 0 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 0 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 0 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "DIC  = " , round( object$ic$DIC , 0 ) , " | penalty =" , 
				round( 2*object$ic$pD ,2 ) , 
			"   | DIC = -2*LL + 2*pD  \n" )    			
    cat( "CAIC = " , round( object$ic$CAIC , 0 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   
	
	
	cat("-----------------------------------------------------------------\n")
	cat("Parameter Summary \n")	

	obji <- object$amh_summary
	vars <- c("parameter","MAP","SD", "Q2.5", "Q97.5" , "Rhat","SERatio",
					"effSize" , "accrate")
	obji <- obji[-1,vars]
	NO <- ncol(obji) - 1
	for (vv in c( 2:(NO-1) , NO+1) ){
		obji[,vv] <- round( obji[,vv] , digits )
		}
	obji[,NO] <- round( obji[,NO] )	
	print(obji)
	invisible(obji)

	# close sink
    CDM::csink( file = file )		
	
	}
#############################################################	