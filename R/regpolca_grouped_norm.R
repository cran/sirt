## File Name: regpolca_grouped_norm.R
## File Version: 0.01

regpolca_grouped_norm <- function(x)
{
    nx <- length(x)
    norm <- sqrt(nx)*sqrt(sum(x^2))
    return(norm)
}
