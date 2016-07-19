
mlnormal_adjust_numdiff_parameter <- function( h , pars ){
    #** select h parameters according to size of parameters
    abs_par <- base::abs(pars)       
    hvec <- h * base::ifelse( abs_par > 1 , abs_par , 1 )
	base::return(hvec)
}