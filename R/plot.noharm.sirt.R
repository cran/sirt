

#####################################################
plot.noharm.sirt <- function( x , what="est" , efa.load.min=.3 , ... ){

	object <- x
	
	#***********
	# redefine object in case of factor analysis
	if ( object$model.type == "EFA" ){
		Fval <- object$promax
	 	object$Fpatt <- Fpatt <- 1 * ( abs(Fval) > efa.load.min	 )
    	 Fval[ Fpatt == 0 ] <- 0
         Pval <- object$factor.cor
		 D <- ncol(object$Ppatt)
		 object$Ppatt <- Ppatt <- 1-diag(D)
		  what <- "est"
		  object$loadings.theta <- Fval
		  object$factor.cor <- Pval
			}
	#*********
	# lavaan model
	Fpatt <- object$Fpatt
	Ppatt <- object$Ppatt
	if ( what=="est" ){
		Fval <- object$loadings.theta
		Pval <- object$factor.cor
				} else {
		Fval <- object$standardized.solution$Fval
		Pval <- object$standardized.solution$Pval
					}
	varnames <- colnames( object$dat )
#	colnames(Pval) <- colnames(Ppatt)
	colnames(Pval) <- rownames(Pval)	
	Psipatt <- object$Psipatt
	Psival <- object$residcorr

	#*****
	# create lavaan model
	
	# look for manifest variables with are unconnected to all other variables
	indvar_unconnected <- which( rowSums( Fpatt != 0 ) == 0  )
	I <- nrow(Fval)
	indvar_connected <- setdiff( 1:I , indvar_unconnected )
	var_connected <- setdiff( varnames , varnames[indvar_unconnected] )
	
	# thresholds
    if (length(var_connected) > 0 ){
			lavmodel <- paste0( paste0( var_connected , collapse="+" ) , " | t1 \n" )
					} else {
			lavmodel <- NULL 
					}
	# loadings for all dimensions
	D <- ncol(Fval)
	# loading matrix
	desF <- NULL
	for (dd in 1:D){
		# dd <- 1
		ind <- which( ( Fpatt[,dd] > 0 ) | ( Fval[,dd] != 0 ) )
		if ( length(ind) > 0 ){
			dfr1 <- data.frame( "dim" = dd , "dimlabel" = colnames(Pval)[dd] , 
					"ind" = ind , "item" = varnames[ind] , 
					"est" = round(Fval[ind,dd],4) , "Fpatt" = Fpatt[ind,dd] ,
					"lhs" = colnames(Pval)[dd] , "op" = "=~" , "rhs" = varnames[ind]
					 )	 		 
			desF <- rbind( desF , dfr1 )
			l1 <- ifelse( dfr1$Fpatt == 0 , paste0( dfr1$est , "*" ) , "" )        
			lavstr <- paste0( dfr1$dimlabel[1] , " =~ " , 
						paste0( paste0( l1  , dfr1$item ) , collapse="+" ) )
			lavmodel <- paste0( lavmodel , lavstr , "\n" , collapse="")
							}
					}					
	# matrix of factor correlations
	desP <- NULL
	for (ii in 1:D){
	for (jj in ii:D){
		# ii <- 1
		# jj <- 1
		dfr1 <- data.frame( "dim1" = colnames(Pval)[ii] , "dim2" = colnames(Pval)[jj]  )
		ind <- ( Ppatt[ii,jj] > 0 ) | ( Pval[ii,jj] != 0 )
		dfr1$est <- round( Pval[ii,jj] , 4 )
		dfr1$Ppatt <- Ppatt[ii,jj]
		l1 <- ifelse( dfr1$Ppatt == 0  , paste0( dfr1$est , "*" ) , "" )
		lavstr <- paste0( dfr1$dim1 , " ~~ " , l1 , dfr1$dim2 )
		if ( ind ){
			desP <- rbind( desP , dfr1 )
			lavmodel <- paste0( lavmodel , lavstr , "\n" , collapse="")
				}
		   }
		}
	# matrix of residual correlations	
	desPsi <- NULL
if ( ! is.null( Psival) ){
	
	for (ii in 1:(I-1) ){
	  for (jj in (ii+1):I ){
		dfr1 <- data.frame( "var1" = varnames[ii] , "var2" = varnames[jj] )
		ind <- ( Psipatt[ii,jj] > 0 ) | ( Psival[ii,jj] != 0 )
		dfr1$est <- Psival[ii,jj]
		dfr1$Psipatt <- Psipatt[ii,jj]
		l1 <- ifelse( dfr1$Psipatt == 0  , paste0( dfr1$est , "*" ) , "" )
		lavstr <- paste0( dfr1$var1 , " ~~ " , l1 , dfr1$var2 )
		if ( ind ){
			desPsi <- rbind( desPsi , dfr1 )
			lavmodel <- paste0( lavmodel , lavstr , "\n" , collapse="")
				}
						}
					}
				}
		unconvars <- length(indvar_unconnected) > 0				
		if ( unconvars ){
			
			IV <- length(indvar_unconnected)
			for (ii in 1:IV){
				# ii <- 1
				varii <- varnames[ indvar_unconnected[ii] ]
				lavstr <- paste0( varii , " ~~ 1*" , varii )
				lavmodel <- paste0( lavmodel , lavstr , "\n" , collapse="")
							}						
						}
						
	# lavaanified model
	lavmodelM <- lavaan::lavaanify(lavmodel)
	
	# collect estimates
	est <- c( object$thresholds[ indvar_connected ]  , desF$est , desP$est )
	if ( ! is.null( desPsi) ){ est <- c( est , desPsi$est ) }
	if ( unconvars ){ est <- c( est , rep(1,IV) ) }
	# create pseudo lavaan object
	# data(..., envir = environment())	
    HolzingerSwineford1939 <- NULL	
	data(HolzingerSwineford1939, package="lavaan" , envir= environment() )
	HS.model <- ' visual  =~ x1 + x2 + x3
				  textual =~ x4 + x5 + x6
				  speed   =~ x7 + x8 + x9 '
	lavfit <- lavaan::cfa(HS.model, data=HolzingerSwineford1939)
	# modify resulting object
	lavfit1 <- lavfit
	lpt <- lavfit1@ParTable
	lavfit1@ParTable <- as.list( lavmodelM )
	lavfit1@Fit@est <- est
	lavfit1@Fit@se <- 0
	# sem paths
	plotres <- semPlot::semPaths(object= lavfit1 , what="est" , ... )
#	plotres$lavmodel
#	plotres$lavfit1
	invisible(plotres)
		}
##############################################################