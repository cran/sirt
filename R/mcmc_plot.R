

######################################################
# mcmc plot
mcmc_plot <- function(mcmcobj , ...){	
	x <- base::list( "mcmcobj" = mcmcobj )
	x$amh_summary <- mcmc_summary(mcmcobj)
	base::class(x) <- "amh"
	plot.amh(x, ... )
		}
