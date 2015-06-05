

###############################################################
lsem.parameter.summary <- function( parameters , moderator.density , verbose ){

	mod.density <- moderator.density
	NP <- max(parameters$parindex)
	if (verbose){
		cat("** Parameter summary \n\n")
		flush.console()
				}	
	parameters_summary <- NULL
	for (pp in 1:NP){
		# pp <- 1
		par.pp <- parameters[ parameters$parindex == pp , ]
		pars1 <- data.frame( "par" = paste(par.pp$par[1]) , "parindex"=pp)
		x <- par.pp[,"est"]
		pars1$M <- weighted.mean( x , mod.density[,2] )
		pars1$SD <- sirt.wtdSD( x ,mod.density[,2] ) 
		pars1$MAD <- sum( mod.density[,2] * abs( x - pars1$M ) )
		pars1$Min <- min(x)
		pars1$Max <- max(x)		
		mod1 <- lm( x ~ mod.density[,1]  , weights = mod.density[,2] )
		pars1$lin_int <- coef(mod1)[1]
		pars1$lin_slo <- coef(mod1)[2]
		pars1$SD_nonlin <- sirt.wtdSD( resid(mod1) ,mod.density[,2] )  
		parameters_summary <- rbind( parameters_summary , pars1 )
					}		
					
	return( parameters_summary )
			}
###############################################################	