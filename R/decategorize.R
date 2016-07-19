
###########################################################
decategorize <- function( dat , categ_design=NULL ){

	# preliminarites
	dat4 <- dat3 <- dat
	dfr <- categ_design 
	
	#****************************
	# handle categories
	if ( ! is.null( dfr ) ){	
		vars <- base::sort( base::unique( base::paste( dfr$variable )))
		VV <- base::length(vars)
		for (vv in 1:VV){
			# vv <- 3
			dfr.vv <- dfr[ base::paste(dfr$variable) == vars[vv] , ]
			dat4[, vars[vv] ] <- dfr.vv[ base::match( dat3[,vars[vv]] , dfr.vv$recode ) , "orig"] 
					}				
			}	
	#***************************				
	base::return(dat4)
		}	
#################################################################