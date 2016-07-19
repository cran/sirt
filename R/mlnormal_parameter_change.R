
mlnormal_parameter_change <- function( pars , pars0 )
{
	pars_change <- base::max( base::abs( pars - pars0) )
	base::return(pars_change)
}