

###############################################################
lsem.parameter.summary <- function( parameters , moderator.density , verbose ){

	mod.density <- moderator.density
	NP <- base::max(parameters$parindex)
	if (verbose){
		base::cat("** Parameter summary \n\n")
		utils::flush.console()
				}	
	parameters_summary <- NULL
	for (pp in 1:NP){
		# pp <- 1
		par.pp <- parameters[ parameters$parindex == pp , ]
		pars1 <- base::data.frame( "par" = base::paste(par.pp$par[1]) , "parindex"=pp)
		x <- par.pp[,"est"]
		pars1$M <- stats::weighted.mean( x , mod.density[,2] )
		pars1$SD <- sirt.wtdSD( x ,mod.density[,2] ) 
		pars1$MAD <- base::sum( mod.density[,2] * base::abs( x - pars1$M ) )
		pars1$Min <- base::min(x)
		pars1$Max <- base::max(x)		
		mod1 <- stats::lm( x ~ mod.density[,1]  , weights = mod.density[,2] )
		pars1$lin_int <- stats::coef(mod1)[1]
		pars1$lin_slo <- stats::coef(mod1)[2]
		pars1$SD_nonlin <- sirt.wtdSD( stats::resid(mod1) ,mod.density[,2] )  
		parameters_summary <- base::rbind( parameters_summary , pars1 )
					}				
	base::return( parameters_summary )
			}
###############################################################	