
####################################################
# grouping a moderator variable
lsem.group.moderator <- function( data , type , moderator.grid,
			moderator , residualize , h ){
	    moderator.grouped <- NULL	
		if (type=="MGM"){
			G1 <- base::length(moderator.grid)
			moderator.grouped <- data.frame( "min" = moderator.grid[-G1] ,
					"max" = moderator.grid[-1] )
			moderator.grouped$mid <- base::rowMeans( moderator.grouped)
			v1 <- data[ , moderator ]
			v2 <- moderator.grouped$mid[1]
			for (gg in 2:G1){
				v2 <- base::ifelse( v1 > moderator.grouped$max[gg-1] ,
							moderator.grouped$mid[gg] , v2 )
							}
			data[,moderator] <- v2							
			# residualize <- FALSE
			h <- 1E-5
			moderator.grid <- moderator.grouped$mid
						}
		res <- base::list( "data" = data , "moderator.grouped " = moderator.grouped ,
					"residualize"=residualize , "h"=h ,
					"moderator.grid" = moderator.grid )
		base::return(res)			
			}
#########################################################			