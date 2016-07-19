
mlnormal_verbose_f1 <- function(verbose, disp , iter)
{
    if ( verbose ){
	    base::cat(disp)	
	    base::cat("Iteration" , iter , "   " , 
		              base::paste( base::Sys.time() ) , "\n" )
	}		
}		
			