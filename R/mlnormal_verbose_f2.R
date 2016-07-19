
mlnormal_verbose_f2 <- function(verbose, disp , iter, descriptions ,
		objfun , objfun0 , beta_change , theta_change )
{
    if (verbose){			
		base::cat( base::paste0( "   " , descriptions["log_like_verbose2"] , 
			         " = "  , base::round( objfun , 4 ) , 
			if (iter > 1){ " | Change = " } else {""} ,
			if( iter > 1){ base::round(  objfun - objfun0 , 6 )} else { ""}	,"\n",sep="") )
		base::cat( base::paste0( "    Maximum beta parameter change = " , 
					base::paste0( base::round( beta_change  ,6) , collapse=" " ) , "\n" , sep=""))
		base::cat( base::paste0( "    Maximum theta parameter change = " , 
					base::paste0( base::round( theta_change  ,6) , collapse=" " ) , "\n" , sep=""))				
		utils::flush.console()						
	}	
}		