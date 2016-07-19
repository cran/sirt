

hard_thresholding <- function( x , lambda )
{
    x_abs <- base::abs(x)
    x <- base::ifelse( x_abs > lambda , x , 0 )
    base::return(x)
}
