
###############################################################################
xxirt_compute_priorDistribution <- function( Theta , customTheta , G ){
		P_Theta <- customTheta$P	
		arg_Theta <- base::list( "Theta" = Theta , "par" = customTheta$par , "G" = G )
        prior_Theta <- base::do.call( P_Theta , arg_Theta )
		base::return(prior_Theta)
}
###############################################################################				