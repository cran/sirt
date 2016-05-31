

#########################################################
# evaulate log-likelihood for amh fit
amh_loglike <- function( model , amh_summary , data , pars){
	pars_MAP <- amh_summary$MAP[-1]
	names(pars_MAP) <- names(pars)
	ll <- do.call( model , list( pars = pars_MAP , data = data ) )
	return(ll)
		}