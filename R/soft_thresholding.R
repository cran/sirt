

soft_thresholding <- function( x , lambda )
{
    x_abs <- base::abs(x)
    x <- base::ifelse( x_abs > lambda , x - base::sign(x) * lambda , 0 )
    base::return(x)
}


mlnormal_soft_thresholding <- soft_thresholding