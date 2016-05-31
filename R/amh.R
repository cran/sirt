
#######################################################
# adaptive Metropolis-Hastings sampler
amh <- function( data , nobs , pars , model ,  prior , proposal_sd ,
        pars_lower = NULL , pars_upper = NULL , derivedPars = NULL , 
		n.iter = 5000 , n.burnin = 1000 , n.sims=3000 ,
		acceptance_bounds = c(.45,.55) ,
		proposal_refresh = 50 , print_iter = 50 
					){

		time <- list( "start" = Sys.time() )	
		CALL <- match.call()
					

		NP <- length(pars)
		if ( is.null(pars_lower)){
			pars_lower <- rep(-Inf, NP)
				}
		if ( is.null(pars_upper)){
			pars_upper <- rep(Inf, NP)
				}		
		acceptance_target <- mean( acceptance_bounds )
				
		#*** define objects for acceptance rates
		acceptance_parameters <- matrix( 0 , nrow=NP , ncol=2 )
		rownames(acceptance_parameters) <- names(pars)
		colnames(acceptance_parameters) <- c("accepted" , "sampled" )

		#*** at which iteration refreshing should be made?
		iter_refresh <- seq( 0 , n.burnin, proposal_refresh )
		
		#*** create object for saving chains
		n.sims1 <- n.iter - n.burnin 
		n.sims <- min( n.sims , n.sims1 )
		thin0 <- n.sims1 / n.sims
		thin1 <- floor(thin0)
		iter_thinned <- seq(n.iter , n.burnin+1 , - thin1)
		iter_thinned <- iter_thinned[ iter_thinned > n.burnin ]
		NS <- length(iter_thinned)
		pars_chain <- matrix( NA , nrow=NS , ncol=NP+1)
		colnames(pars_chain) <- c("deviance",names(pars))	
		# iter_save <- seq( n.burnin+1 , n.iter )
		iter_save <- iter_thinned
		it <- 1
		ss <- 1
		while( it <= n.iter){
			res <- amh_sampling( pars , data , model , prior ,
                       proposal_sd , acceptance_parameters , pars_lower ,
					   pars_upper )
			pars <- res$pars
			acceptance_parameters <- res$acceptance_parameters
			dev <- res$deviance
									
			#*** save parameters in chain
			if (it %in% iter_save){
				pars_chain[ss,-1] <- pars
				pars_chain[ss,1] <- dev
				ss <- ss+1
					}
					
			#*** refresh proposal SD
			if (it %in% iter_refresh ){
				res0 <- amh_proposal_refresh( acceptance_parameters , 
				            proposal_sd ,  acceptance_bounds )	
				acceptance_parameters <- res0$acceptance_parameters
				proposal_sd <- res0$proposal_sd	
					}			
					
			#*** print progress
			if ( it %% print_iter == 0 ){
				cat( paste0(" ** Iteration " , it , "\n")  )
				utils::flush.console()
					}			
			it <- it + 1
						}
		############# end iterations ######################
	
		#**** create mcmc object
		mcmcobj <- coda::mcmc(data= pars_chain , 
		                start = min(iter_thinned) , end = max(iter_thinned), 
		                thin = thin1)
		#--- include derived parameters in wanted
		ND <- 0
		if ( ! is.null(derivedPars) ){
			mcmcobj <- mcmc_derivedPars( mcmcobj=mcmcobj , derivedPars=derivedPars )		
			ND <- length(derivedPars)
							}
						
		amh_summary <- mcmc_summary(mcmcobj)
		accrate <- acceptance_parameters[,1] / acceptance_parameters[,2]
		amh_summary$accrate <- c( NA , accrate , rep(NA,ND) )
		
		#*** evaluate log-likelihood
		ll <- amh_loglike( model=model, amh_summary=amh_summary, data=data , pars=pars)
		#*** compute information criteria
		ic <- amh_ic( dev = - 2*ll, N = nobs , pars = pars , amh_summary = amh_summary,
					model = model , data = data )
		#*** methods
		coef1 <- amh_summary$MAP[-1]
		names(coef1) <- amh_summary$parameter[-1]		
		vcov1 <- var(mcmcobj)[-1,-1]
		
		time$end <- Sys.time()
		
		#**** output list
		res <- list( 
		    pars_chain = pars_chain ,
			acceptance_parameters = acceptance_parameters ,
			amh_summary = amh_summary , 
			coef = coef1 , vcov = vcov1 , 
			mcmcobj = mcmcobj , loglik = ll , ic = ic , deviance = -2*ll ,
			n.iter = n.iter , n.burnin = n.burnin ,
			thin = thin1 , n.saved = NS ,  proposal_sd = proposal_sd , 
			time=time , CALL = CALL
					)
		res$description <- "Adaptive Metropolis Hastings Sampling"				
		class(res) <- "amh"
		return(res)
			}
#################################################################################			