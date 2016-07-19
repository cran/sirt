
##############################################################
lsem.MGM.stepfunctions <- function( object , moderator.grid ){

		# object <- lsem.object
		moderator.grouped <- object$moderator.grouped
		parameters <- object$parameters
		
		dfr <- NULL
		
		G <- base::length(moderator.grid)

		for (gg in 1:G){
			# gg <- 4	
			mod.gg <- moderator.grid[gg]
			ind.gg <- base::which( ( moderator.grouped$min <= mod.gg ) & 
								( moderator.grouped$max > mod.gg ) )

			grouped.gg <-  object$moderator.grid[ind.gg]  
			parameters.gg <- parameters[ parameters$moderator ==  grouped.gg  , ]
			parameters.gg$moderator <- mod.gg
			dfr <- base::rbind( dfr , parameters.gg )
					}
		dfr <- dfr[ base::order( 10E20*dfr$parindex + dfr$moderator ) , ]
		base::return(dfr)
		}
##############################################################