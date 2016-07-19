
#######################################################
# adaptive Metropolis-Hastings sampler
amh <- function( data , nobs , pars , model ,  prior , proposal_sd ,
        pars_lower = NULL , pars_upper = NULL , derivedPars = NULL , 
		n.iter = 5000 , n.burnin = 1000 , n.sims=3000 ,
		acceptance_bounds = c(.45,.55) ,
		proposal_refresh = 50 , print_iter = 50 
					){

		time <- base::list( "start" = Sys.time() )	
		CALL <- base::match.call()
				
		#*** convert prior if needed
		res <- pmle_process_prior( prior , pars )	
		prior <- res$prior
		dens <- res$dens				
		
		#****							
		NP <- base::length(pars)
		if ( is.null(pars_lower)){
			pars_lower <- base::rep(-Inf, NP)
		}
		if ( is.null(pars_upper)){
			pars_upper <- base::rep(Inf, NP)
		}		
		acceptance_target <- base::mean( acceptance_bounds )
		
		#*** object pseudo likelihood estimation
		pmle_pars <- base::list( pars = pars , posteriorval = -Inf )
		
		#*** define objects for acceptance rates
		acceptance_parameters <- base::matrix( 0 , nrow=NP , ncol=2 )
		base::rownames(acceptance_parameters) <- base::names(pars)
		base::colnames(acceptance_parameters) <- base::c("accepted" , "sampled" )

		#*** at which iteration refreshing should be made?
		iter_refresh <- base::seq( 0 , n.burnin, proposal_refresh )
		LI <- base::length(iter_refresh)
		proposal_sd_history <- matrix( NA , nrow=LI , ncol=NP )
		proposal_sd_history[1,] <- proposal_sd
		base::colnames(proposal_sd_history) <- base::names(pars)
		base::rownames(proposal_sd_history) <- iter_refresh
		base::attr(proposal_sd_history,"include") <- 2
		
		#*** create object for saving chains
		n.sims1 <- n.iter - n.burnin 
		n.sims <- base::min( n.sims , n.sims1 )
		thin0 <- n.sims1 / n.sims
		thin1 <- base::floor(thin0)
		iter_thinned <- base::seq(n.iter , n.burnin+1 , - thin1)
		iter_thinned <- iter_thinned[ iter_thinned > n.burnin ]
		NS <- base::length(iter_thinned)
		pars_chain <- base::matrix( NA , nrow=NS , ncol=NP+1)
		base::colnames(pars_chain) <- base::c("deviance",names(pars))	
		# iter_save <- seq( n.burnin+1 , n.iter )
		iter_save <- iter_thinned
		it <- 1
		ss <- 1
		while( it <= n.iter){
			res <- amh_sampling( pars=pars , data=data , model=model , prior=prior ,
                       proposal_sd=proposal_sd , 
					   acceptance_parameters=acceptance_parameters , pars_lower=pars_lower ,
					   pars_upper=pars_upper , pmle_pars=pmle_pars )
			pars <- res$pars
			acceptance_parameters <- res$acceptance_parameters
			dev <- res$deviance
			pmle_pars <- res$pmle_pars

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
				hh <- base::attr(proposal_sd_history,"include")
				proposal_sd_history[ hh , ] <- proposal_sd
				base::attr(proposal_sd_history,"include") <- hh + 1
			}			
					
			#*** print progress
			if ( it %% print_iter == 0 ){
				base::cat( base::paste0(" ** Iteration " , it , "\n")  )
				utils::flush.console()
			}			
			it <- it + 1
		}
		############# end iterations ######################
	
		#**** create mcmc object
		mcmcobj <- coda::mcmc( data= pars_chain , 
		                start = base::min(iter_thinned) , 
						end = base::max(iter_thinned), 
		                thin = thin1 )
		#--- include derived parameters in wanted
		ND <- 0
		if ( ! is.null(derivedPars) ){
			mcmcobj <- mcmc_derivedPars( mcmcobj=mcmcobj , derivedPars=derivedPars )		
			ND <- base::length(derivedPars)
		}
						
		amh_summary <- mcmc_summary(mcmcobj)
		accrate <- acceptance_parameters[,1] / acceptance_parameters[,2]
		amh_summary$accrate <- base::c( NA , accrate , base::rep(NA,ND) )
		
		#*** evaluate log-likelihood for MAP estimator	
		parnames <-	base::names(pars)
		pars_MAP <- amh_summary[ 2:(NP+1) , "MAP"]
		base::names(pars_MAP) <- parnames	
		res_MAP <- amh_posterior( pars = pars_MAP , model , prior , data)
		ll <- res_MAP$ll

		#*** evaluate log-likelihood for mean
		pars_Mean <- amh_summary[ 2:(NP+1) , "Mean"]
		base::names(pars_Mean) <- parnames	
		res_Mean <- amh_posterior( pars = pars_Mean , model , prior , data)
			
		#*** compare different estimators
		dfr_estimators <- amh_compare_estimators(res_MAP=res_MAP , res_Mean=res_Mean , 
					res_pmle = pmle_pars )
	
		#*** compute information criteria
		ic <- amh_ic( dev = - 2 * res_MAP$ll , N = nobs , pars = pars , 
					amh_summary = amh_summary,
					model = model , data = data , priorval = res_MAP$priorval  )
										
		#*** methods
		coef1 <- amh_summary$MAP[-1]
		base::names(coef1) <- amh_summary$parameter[-1]		
		vcov1 <- stats::var(mcmcobj)[-1,-1]
		
		time$end <- base::Sys.time()
		
		#**** output list
		res <- base::list( 
		    pars_chain = pars_chain ,
			acceptance_parameters = acceptance_parameters ,
			amh_summary = amh_summary , 
			coef = coef1 , vcov = vcov1 , 
			pmle_pars = pmle_pars , 
			comp_estimators = dfr_estimators , 
			mcmcobj = mcmcobj , loglik = ll , ic = ic , 
			deviance = -2*ll ,
			model = model , prior = prior , data = data , nobs = nobs , 
			prior_summary = dens , 
			n.iter = n.iter , n.burnin = n.burnin ,
			thin = thin1 , n.saved = NS ,  proposal_sd = proposal_sd , 
			proposal_sd_history = proposal_sd_history , 
			time=time , CALL = CALL
					)
		res$description <- "Adaptive Metropolis Hastings Sampling"				
		base::class(res) <- "amh"
		base::return(res)
}
#################################################################################			