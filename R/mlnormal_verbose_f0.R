
mlnormal_verbose_f0 <- function(verbose, disp)
{
	if (verbose){		
		base::cat(disp)
	    base::cat("Preprocess data   " , 
		          base::paste( base::Sys.time() ) , "\n" )
		utils::flush.console()						
	}	
}		
			