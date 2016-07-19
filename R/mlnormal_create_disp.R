
mlnormal_create_disp <- function(symbol="." , length=30 , line_break = TRUE )
{
	v1 <- base::paste0( base::rep(symbol , length= length) , collapse="")
	v1 <- base::paste0( v1 , "\n")
	base::return(v1)
}