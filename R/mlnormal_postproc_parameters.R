

#############################################################
# postprocessing of parameters
mlnormal_postproc_parameters <- function( theta , theta_init ,
		beta , beta_init , theta_infomat , XVX , level , prior_args ){ 
	
	theta_names <- base::names(theta_init)
	beta_names <- base::names(beta_init)	
	
	theta <- mlnormal_as_vector_names(pars = theta , parnames = theta_names	)
	beta <- mlnormal_as_vector_names(pars = beta , parnames = beta_names )
	NB <- base::length(beta)
	NT <- base::length(theta)
	
	theta_vcov <- base::solve(theta_infomat)
	base::rownames(theta_vcov) <- base::colnames(theta_vcov) <- theta_names
													
	theta_summary <- base::data.frame(
						"parm" = theta_names ,
						"prior" = NA , 
						"penalty" = NA , 
						"est" = theta , 
						"se" = mlnormal_sqrt_diag(matr=theta_vcov) )
    theta_summary <- parmsummary_extend( dfr = theta_summary , level = level )
	
	beta_vcov <- base::solve(XVX)
	beta_summary <- base::data.frame(
						"parm" = beta_names ,
						"prior" = NA , 	
						"penalty" = NA , 						
						"est" = beta , 
						"se" = mlnormal_sqrt_diag(matr=beta_vcov) )
    beta_summary <- parmsummary_extend( dfr = beta_summary , level = level )

	if ( prior_args$use_prior ){
		dens <- base::paste0(prior_args$dens$prior)
		beta_summary$prior <- dens[ 1:NB ]
		theta_summary$prior <- dens[ NB + 1:NT ]
	} else {
		beta_summary$prior <- NULL
		theta_summary$prior <- NULL
	}				

	if ( prior_args$use_penalty ){
		pens <- prior_args$penalty_pars
		beta_summary$penalty <- pens$lambda_beta * pens$weights_beta
		theta_summary$penalty <- pens$lambda_theta * pens$weights_theta 
	} else {
		beta_summary$penalty <- NULL
		theta_summary$penalty <- NULL
	}	
	
	#--- collect output
	coefs <- base::c( beta , theta )
		
	NP <- NB + NT
	parnames <- base::c( beta_names , theta_names )
	base::names(coefs) <- parnames
	vcovs <- base::matrix( 0 , nrow=NP , ncol=NP)
	base::rownames(vcovs) <- base::colnames(vcovs) <- parnames 
	vcovs[ 1:NB , 1:NB ] <- beta_vcov
	vcovs[ NB + 1:NT , NB + 1:NT ] <- theta_vcov
		
	#---------------------------------
	# OUTPUT
	res <- base::list( beta = beta , theta = theta , beta_summary = beta_summary ,
				theta_summary = theta_summary , beta_vcov = beta_vcov ,
				theta_vcov = theta_vcov , coefs = coefs , vcovs = vcovs)
	base::return(res)	
}
####################################################################