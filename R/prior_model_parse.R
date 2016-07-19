
#############################################################
prior_model_parse <- function( prior_model ){
	ps <- base::strsplit( prior_model , split="\n" , fixed=TRUE)[[1]]
	# clean string
	ps <- prior_model_pars_CleanString( ps )	
	NP <- base::length(ps)
	prior <- base::as.list( 1:NP )
	for (pp in 1:NP){
		# pp <- 1
		ps_pp <- ps[pp]
		ps_pp1 <- base::strsplit( ps_pp , split="~" , fixed=TRUE)[[1]]
		# extract name
		base::names(prior)[pp] <- ps_pp1[1]
		prior[[pp]] <- base::as.list(1:2)
		# extract distribution
		ps_pp2 <- ps_pp1[2]
		ps_pp2a <- base::strsplit( ps_pp2 , split="(" , fixed=TRUE)[[1]]
		prior[[pp]][[1]] <- ps_pp2a[1]
		ps_pp3 <- ps_pp2a[2]
		ps_pp3 <- base::strsplit( ps_pp3 , split=")" , fixed=TRUE)[[1]][1]	
		ps_pp3 <- base::strsplit( ps_pp3 , split="," , fixed=TRUE)[[1]]
		NV <- base::length(ps_pp3)
		prior_pp2 <- as.list( 1:NV )		
		for (vv in 1:NV){
			# vv <- 1
			ps_vv <- ps_pp3[vv]
			h_vv <- base::strsplit(ps_vv , split="=" , fixed=TRUE)[[1]]
			len_h_vv <- base::length(h_vv)
			if ( len_h_vv == 1){
				prior_pp2[[vv]] <- base::suppressWarnings( base::as.numeric(h_vv) )
			}
			if ( len_h_vv == 2){
				prior_pp2[[vv]] <- base::suppressWarnings( base::as.numeric(h_vv[2]) )
				base::names(prior_pp2)[vv] <- h_vv[1]
			}			
		}
		prior[[pp]][[2]] <- prior_pp2
	}	
	base::return(prior)
}
###############################################################		
	