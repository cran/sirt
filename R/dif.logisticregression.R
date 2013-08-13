 
# 0.01  2012-xx-yy


# 0.01  2012-06-23  o initial release


#-------------------------------------------------------



##NS export(dif.logistic.regression)
#---------------------------------------------------------------------------------------##
# This function performs itemwise DIF analysis by using logistic regression methods     ##
# uniform and nonuniform DIF                                                            ##
dif.logistic.regression <- function( dat , group , score ){
    # INPUT:
    # dat       ... data frame (must only include item responses)
    # group     ... group identifier (this has to be a dummy variable)
    # score     ... matching criterion
    
    I <- ncol(dat)
    matr <- NULL
    for (ii in 1:I){
     # ii <- 6
        dat.ii <- na.omit(data.frame( "y" = dat[,ii] , "score" = score , "group" = group ))
        mod1 <- glm( y  ~ score , data = dat.ii , family="binomial")
        mod2 <- glm( y  ~ score + group , data = dat.ii , family="binomial")
        mod3 <- glm( y  ~ score + group + score*group , data = dat.ii , family="binomial")
        # calculate correlation matrix
#        dat.ii$score.group <- dat.ii$score * dat.ii$group
#        c.ii <- cor( dat.ii )
# print( c.ii )
        # calculate R² according to Jodoin and Gierl (2001)
#        r2.mod1 <- mod1$coef[-1] * c.ii[2,1]    
#        r2.mod2 <- sum( mod2$coef[-1] * c.ii[2:3,1]  )
#        r2.mod3 <- sum( mod3$coef[-1] * c.ii[2:4,1]  )    
		
		h1 <- data.frame( "item" = colnames(dat)[ii] , 
				"N" = sum( 1- is.na( dat[,ii] ) , na.rm=T) )
		h1$R <- min(group)
		h1$F <- max(group)
		h1$nR <- sum(  ( 1- is.na( dat[,ii] ) )* (1-group) , na.rm=T)
		h1$nF <- sum(  ( 1- is.na( dat[,ii] ) )* (group) , na.rm=T)		
		h1$p <- mean(  dat[,ii], na.rm=T) 
        a1 <- aggregate( dat[,ii] , list( group) , mean , na.rm=T )[,2]
		h1$pR <- a1[1]
		h1$pF <- a1[2]
		h1$pdiff.adj <- NA	
		h1$uniformDIF <- mod2$coef[3]
		h1$se.uniformDIF <- sqrt( diag( vcov(mod2)) )[3]
		h1$t.uniformDIF <- mod2$coef[3] / sqrt( diag( vcov(mod2) ) )[3] 
		h1$sig.uniformDIF <- ""
		if ( h1$t.uniformDIF > 1.96 ){ h1$sig.uniformDIF <- "+" }
		if ( h1$t.uniformDIF < - 1.96 ){ h1$sig.uniformDIF <- "-" }
		h1$nonuniformDIF <- mod3$coef[4]
		h1$se.nonuniformDIF <- sqrt( diag( vcov(mod3)) )[4]
		h1$t.nonuniformDIF <- mod3$coef[4] / sqrt( diag( vcov(mod3) ) )[4] 
		h1$sig.nonuniformDIF <- ""
		if ( h1$t.nonuniformDIF > 1.96 ){ h1$sig.nonuniformDIF <- "+" }
		if ( h1$t.nonuniformDIF < - 1.96 ){ h1$sig.nonuniformDIF <- "-" }
		matr <- rbind( matr , h1 )
        cat( ii , " " ) ; flush.console()
        if ( ii %% 15 == 0 ){ cat("\n") }
        }
    cat("\n")
    # include variable of adjusted p values
#    ind <- which( colnames(matr) == "pF" )       
    matr[ , "pdiff.adj" ] <- matr$pR - matr$pF - mean( matr$pR - matr$pF  )   
    # matr <- data.frame( "item" = colnames(dat) , matr )
    return(matr)
    }
#------------------------------------------------------------------------------

