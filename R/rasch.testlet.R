 

#***********************************************************
# rasch.testlet
#***********************************************************


# 0.01  2012-06-23  o initial release
# 0.02  2012-07-05  o restructure rasch.testlet.glmer function
# 0.03  2012-07-12	o fixed a bug (.., .prnum missing)
# 0.04  2012-07-13	o include EAP and posterior variance estimation
# 0.05  2012-07-13	o allow for cloglog and loglog link functions
# 0.06  2012-07-13	o changed arguments in .prnum function
# 0.07  2012-07-13	o fixed a bug: smod <- summary(mod) in lme4 only allows 
#					  smod@coefs and not coef(smod)	
# 0.08  2012-07-14	o worked on a bug (coef??)
# 0.09  2012-07-14	o ...

# 0.01  2012-xx-yy

#-------------------------------------------------------



########################################################################################
# Estimation of testlet models in GLMM (lme4)
##NS export(rasch.testlet.glmer)
rasch.testlet.glmer <- function( dat , testlet.matrix = NULL , 
			link = "logit" , verbose = TRUE , progress=TRUE  ){
    # load R packages
	library(lme4)
	#***
    # setup data
    # dat1 <- data.frame( "pid" = seq( 1 , nrow(dat)) , dat )
    I <- ncol(dat)
    n <- nrow(dat)
    dat1 <- matrix( t(dat) , ncol= 1 , byrow=T )
    dat1 <- data.frame( "pid" = rep( 1:n  , each=I ) , "resp" = dat1 , 
                "item" = rep( colnames(dat) , n  ) )
    dat1 <- na.omit(dat1)
	# loglog link is not implemented in lme4
	# but: loglog for correct responses corresponds to the model
	#      of cloglog link for incorrect responses
	#--------------------------------------------------------
	# Small function which helps for printing purposes
	.prnum <- function( matr , digits , columns = 1:( ncol(matr)) ){
#			VV <- ncol(matr)
			for (vv in columns){
			# vv <- 1
			if ( is.numeric( matr[,vv]) ){ matr[,vv] <- round( matr[,vv] , digits ) }
							}
			print(matr)
					   }
	#--------------------------------------------------------	
	##########################################################
    # define testlet matrix
	if ( is.null( testlet.matrix)){
					est.testlet <- FALSE 	# no testlet model
					testlets <- NULL
							} else {	
					est.testlet <- TRUE		# estimation of the testlet model			
			if ( is.vector(testlet.matrix) ){ 
				tm <- sort( unique( testlet.matrix[ testlet.matrix > 0 ] ) )
				dfr <- NULL
				for (tt in tm){ 
						# tt <- tm[1]
						dfr.tt <- data.frame( "testlet" = paste( "Testlet" , tt , sep="") , 
									"item" = paste(colnames(dat)[ which( testlet.matrix == tt ) ] ) )
						dfr <- rbind( dfr , dfr.tt )
							} 
				testlet.matrix <- dfr
						}
			# define testlets in data frame
			testlets <- sort( unique( paste( testlet.matrix[,1] ) ) )
			for (tt in testlets){ 
				# tt <- testlets[1]
				items.tt <- paste( testlet.matrix[ paste(testlet.matrix[,1]) == tt , 2] )
				dat1[ , tt ] <- 1 * ( paste( dat1$item ) %in% items.tt )
								}
						}								
	##########################################################
    # no testlet model
	if ( ( ! est.testlet ) ){
		cat( "------------------------------------------------\n")
		modr <- "resp ~ 0 + as.factor(item) + ( 1 | pid )"
		s1 <- Sys.time()
		# logit link (Rasch model)
		if ( link == "logit"){ 
			cat("Estimation of the Rasch model using lme4\n\n")
			cat( modr , "\n\n"); flush.console()
			mod <- lmer( as.formula(modr)  , data = dat1 , family = "binomial" ,
							verbose = verbose  )
							}
		# linear probability model
		if ( link == "linear"){ 
			cat("Estimation of the linear probability model using lme4\n\n")
			cat( modr , "\n\n"); flush.console()
			mod <- lmer( as.formula(modr)  , data = dat1 , 
							verbose = verbose )
							}
		# Rasch type model (Cloglog link)
		if ( link == "cloglog"){ 
			cat("Estimation of the Rasch type cloglog link model using lme4\n\n")
			cat( modr , "\n\n"); flush.console()
			mod <- lmer( as.formula(modr)  , data = dat1 , 
							verbose = verbose , family = binomial( link = "cloglog" ) )
							}
		# loglog link
		if ( link == "loglog"){ 
			cat("Estimation of the Rasch type loglog link model using lme4\n\n")
			cat( modr , "\n\n"); flush.console()
			dat2 <- dat1
			dat2$resp <- 1 - dat1$resp			
			mod <- lmer( as.formula(modr)  , data = dat2 , 
							verbose = verbose , family = binomial( link = "cloglog" ) )
							}							
							}
	#########################################################################
	# testlet model
	if ( est.testlet ){
	    cat( "\n\n------------------------------------------------\n")
		modr <- "resp ~ 0 + as.factor(item) + ( 1 | pid )"				
		for (tt in testlets ){ modr <- paste( modr , "+ ( - 1 +" , tt , "| pid )" )}
		s1 <- Sys.time()
		if ( link == "logit"){
		    cat("Estimation of the Rasch testlet model using lme4\n\n")
			cat( modr , "\n\n") ; flush.console()					
			mod <- lmer( as.formula(modr)  , data = dat1 , family = "binomial" ,
							verbose = verbose )
							}
		if ( link == "linear"){
		    cat("Estimation of the linear probability testlet model using lme4\n\n")
			cat( modr , "\n\n") ; flush.console()					
			mod <- lmer( as.formula(modr)  , data = dat1 , verbose = verbose )
							}
		if ( link == "cloglog"){
		    cat("Estimation of the Rasch type cloglog link testlet model using lme4\n\n")
			cat( modr , "\n\n") ; flush.console()					
			mod <- lmer( as.formula(modr)  , data = dat1 , verbose = verbose ,
							family = binomial( link = "cloglog") )
							}
		if ( link == "loglog"){
		    cat("Estimation of the Rasch type loglog link testlet model using lme4\n\n")
			cat( modr , "\n\n") ; flush.console()			
			dat2 <- dat1
			dat2$resp <- 1 - dat1$resp			
			mod <- lmer( as.formula(modr)  , data = dat2 , verbose = verbose ,
							family = binomial( link = "cloglog") )
							}
					}
	##############################################################################
    # person parameter estimation respective prediction
#	re <- ranef( mod , postVar = TRUE )
		# This function does not seem to work for uncorrelated effects
	re <- ranef( mod )
	vc <- VarCorr( mod )
	vc1 <- vc$pid[1,1]
	eap <- (re$pid)[,1]	
	if ( link == "loglog" ){ eap <- - eap }
	# calculate EAP reliability
	eap.rel <- var( eap ) / vc1 
	##############################################################################	
		cat("\n\n")
#		print( summary(mod));
		print(mod)
		flush.console()
		s2 <- Sys.time()	
		# print needed computational time
		if (verbose){ 
					cat("---------------------------------\n")
					cat("Start:" , paste( s1) , "\n")
					cat("End:" , paste(s2) , "\n")
					cat("Difference:" , print(s2 -s1), "\n")
					cat("---------------------------------\n")
     				} 
	###################################################################################
	# print summary
	#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    cat( "\n\n------------------------------------------------\n")
    cat("*******************  RESULTS  ******************\n")
#    smod <- summary(mod)
     f1 <- fixef(mod)
#	f1 <- smod@coefs
#    rownames(f1) <- gsub( "as.factor(item)" , "" , rownames(f1) , fixed=TRUE)
    names(f1) <- gsub( "as.factor(item)" , "" , names(f1) , fixed=TRUE)
#    f1 <- f1[ colnames(dat) , ]	
#    f1[,c(1,3)] <- - f1[,c(1,3)]
    item <- data.frame( "item" = names(f1) , "b" = -f1 )
	if ( link == "loglog"){ item$b <- - item$b }
	if ( est.testlet ){ 
			item$testlet <-  testlet.matrix[ match( names(f1) , testlet.matrix[,2]  ) , 1 ] 
					} 	
	item$link <- link					
    sd2 <- sqrt( unlist( VarCorr(mod)) )
    sdmatr <- data.frame( "effect" = c("trait" , testlets ) , "SD" = sd2 )
	# include residual SD
	if ( link == "linear" ){
		sdmatr1 <- data.frame( "Residual" , attr( VarCorr(mod) , "sc" ) )
							}
	if ( link == "logit" ){
		sdmatr1 <- data.frame( "Residual" , sqrt(pi^2/3)	)
					}
	if ( link %in% c("loglog","cloglog" ) ){
		sdmatr1 <- data.frame( "Residual" , NA	)
					}
	colnames(sdmatr1) <- colnames(sdmatr)
	sdmatr <- rbind( sdmatr , sdmatr1 )
	rownames(sdmatr) <- NULL
	#'''''''''''''''''''''''''''''''''''''''''
    cat("\nItem difficulties\n")
	cat("\n")
    .prnum( item , 3 , columns=2)
    cat("\n\nStandard deviations\n\n")
    .prnum( sdmatr , 3 , columns=2)
    cat("\nEAP reliability: ")
    cat( round( eap.rel, 3 ) , "\n")
    # collect all results into a result list
    res <- list( "mod" = mod ,"item" = item ,  "sd.testlet" = sdmatr ,
			"dat1" = dat1 , "link" = link , "testlet.matrix" = testlet.matrix ,
			"eap" = eap , "eap.rel" = eap.rel )
#	class(res) <- "rasch.testlet.glmer"
    return(res)
	}
#####################################################################	
	
