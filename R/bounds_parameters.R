

bounds_parameters <- function( pars , lower = NULL , upper = NULL)
{
	if ( ! base::is.null(lower)){
		pars <- base::ifelse( pars < lower , lower , pars )
	}
	if ( ! base::is.null(upper)){
		pars <- base::ifelse( pars > upper , upper , pars )
	}
	base::return(pars)	
}