
##################################################
# Linking Haberman ETS Research Report
linking.haberman <- function( itempars , conv = .00001 , maxiter=1000 ,
	progress=TRUE ){
	#****
    # convert itempars to data frame
	itempars <- as.data.frame( itempars )
	# include wgt if there does not exist a fifth columm
	if ( ncol(itempars) == 4){
		itempars$wgt <- 1
				}
	# extract studies
	studies <- sort( paste( unique( itempars[,1] ) ) )
	NS <- length(studies)
	# extract items
	items <- sort( paste( unique( itempars[,2] ) ) )
	NI <- length(items)
	# define a and b matrices
	wgtM <- bM <- aM <- matrix(NA , nrow=NI , ncol=NS)
	rownames(wgtM) <- rownames(bM) <- rownames(aM) <- items
	colnames(wgtM) <- colnames(bM) <- colnames(aM) <- studies
	# define item parameters
	for (ss in studies){
		# ss <- studies[1]
		itempars.ss <- itempars[ itempars[,1] == ss , ]
		aM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,3] 
		bM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,4] 
		wgtM[ paste(itempars.ss[,2]) , ss ] <- itempars.ss[,5] 
					}
	wgtM <- wgtM / matrix( rowSums( wgtM , na.rm=TRUE ) , nrow=NI , ncol=NS )
	#*****
	# estimation of A
	logaM <- log( aM )	
	resA <- .linking.haberman.als(logaM=logaM , wgtM=wgtM , maxiter=maxiter , 
				 conv=conv , progress = progress , est.type="A")
	aj <- exp( resA$logaj )
	At <- exp( resA$logaAt )
	#******
	# estimation of B
	bMadj <- bM / matrix( At , NI , NS , byrow=TRUE )
	resB <- .linking.haberman.als(logaM=bMadj , wgtM=wgtM , maxiter=maxiter , 
				 conv=conv , progress = progress , est.type="B")
	Bj <- resB$logaj
	Bt <- resB$logaAt
	#*****
	# transformations
	transf.itempars <- data.frame( "study" = studies , "At" = At , "Bt" = Bt )
	rownames(transf.itempars) <- NULL
	# This is the transformation for item parameters.
	transf.personpars <- transf.itempars
	transf.personpars$At <- 1/transf.itempars$At
	transf.personpars$Bt <-  - transf.itempars$Bt / transf.itempars$At	
	# new item parameters
	joint.itempars <- data.frame("item" = items , "aj"=aj , "bj"=Bj )
    # print output
	if (progress){
	    cat("Linear transformation for item parameters\n")
		obji <- transf.itempars
		obji[,-1] <- round( obji[,-1] , 3)
		print( obji )  	
	    cat("\nLinear transformation for person parameters\n")
		obji <- transf.personpars
		obji[,-1] <- round( obji[,-1] , 3)
		print( obji )  			
				}
	res <- list( "transf.itempars" = transf.itempars , 
	  "transf.personpars" = transf.personpars , 
		"joint.itempars" = joint.itempars )
	return(res)
		}
#######################################################################		
	
	