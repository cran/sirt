

#################################################################
# refreshing the proposal SD
amh_proposal_refresh <- function( acceptance_parameters , proposal_sd ,
		acceptance_bounds ){
	target <- base::mean(acceptance_bounds)
	acc <- acceptance_parameters[,1] / acceptance_parameters[,2]
	SD.pp <- proposal_sd
	#*** compute new proposal SD
	SD.pp <- base::ifelse( acc < acceptance_bounds[1] ,
					SD.pp / ( 2 - acc / target ) , SD.pp )	
	SD.pp <- base::ifelse( acc > acceptance_bounds[2] ,
					SD.pp * ( 2 - (1-acc)/(1-target) ) , SD.pp )
	res0 <- base::list( proposal_sd = SD.pp ,
					acceptance_parameters = 0*acceptance_parameters)				
	base::return(res0)									
}								
#################################################################								