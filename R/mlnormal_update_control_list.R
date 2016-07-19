
mlnormal_update_control_list <- function(control, control0)
{
	# control0 is the output list
	# the elements in control are written into control0
	if ( ! base::is.null(control) ){
		N <- base::length(control)
		for (nn in 1:N){
			control0[[ base::names(control)[nn] ]] <- control[[nn]]
		}
	}
	base::return(control0)
}	