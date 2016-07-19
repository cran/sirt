

#########################################################
# evaulate log-likelihood for amh fit
amh_loglike <- function( model , data , pars){
	ll <- base::do.call( model , base::list( pars = pars , data = data ) )
	base::return(ll)
}