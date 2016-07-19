

mlnormal_as_vector_names <- function(pars , parnames)
{
	pars <- base::as.vector(pars)
	base::names(pars) <- parnames
	base::return(pars)
}