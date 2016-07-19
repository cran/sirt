

mlnormal_eval_penalty <- function( beta , theta , penalty_pars )
{
	penalty_beta <- penalty_pars$lambda_beta * penalty_pars$weights_beta * base::abs(beta)
	penalty_theta <- penalty_pars$lambda_theta * penalty_pars$weights_theta * base::abs(theta)
	res <- base::list( penalty_beta = base::sum(penalty_beta) , 
						penalty_theta = base::sum(penalty_theta))						
	base::return(res)
}