

###########################################################################
xxirt_compute_likelihood <- function( probs_items , dat , resp_index ){		
		N <- nrow(dat)
		TP <- dim(probs_items)[3]
		I <- dim(probs_items)[1]
		p.xi.aj <- matrix( 1 , nrow=N , ncol=TP )
		for (ii in 1:I){
			# ii <- 1
			resp_ii <- resp_index[[ii]]
			v1 <- probs_items[ ii , dat[ resp_ii , ii ] + 1 , ]
			p.xi.aj[ resp_ii , ] <- p.xi.aj[ resp_ii ,  ] * v1
						}
		return(p.xi.aj)
			}
###########################################################################			
