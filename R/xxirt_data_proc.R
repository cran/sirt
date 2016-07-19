
##############################################################
xxirt_data_proc <- function(dat , group = NULL, weights=NULL ){ 
		ncat <- base::apply( dat , 2 , base::max , na.rm=TRUE ) + 1
		I <- base::ncol(dat)
		items <- base::colnames(dat)
		N <- base::nrow(dat)
		if ( is.null(group) ){	
			group <- base::rep(1,N)
		}
		if ( is.null(weights) ){	
			weights <- base::rep(1,N) 
		}					
		group0 <- group
		groups_unique <- base::sort( base::unique( group ) )
		G <- base::length(groups_unique)
		group <- base::match( group0 , groups_unique )	
		maxK <- base::max(ncat)
		#*** group_index
		group_index <- base::as.list( 1:G )
		for (gg in 1:G){
			group_index[[gg]] <- base::which( group == gg )
		}
		#*** data with response indices
		dat_na <- base::is.na(dat)
		dat_resp <- 1 - dat_na
		resp_index <- base::as.list( 1:I )				
		for ( ii in 1:I){
			resp_index[[ii]] <- base::which( dat_resp[,ii] == 1 )
		}
		dat1 <- base::as.matrix(dat)
		dat1[ dat_na ] <- 0		
								
		#*** output
		res <- base::list( N = N , I = I , group = group , items = items , 
					group0 = group0 , G = G , groups_unique = groups_unique,
					maxK = maxK , ncat = ncat, weights=weights,
					group_index = group_index , dat_resp = dat_resp ,
					resp_index = resp_index , dat1 = dat1 )
		base::return(res)		
}
##############################################################		