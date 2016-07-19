
###################################################
# improper density which is constant to 1
dimproper <- function(x){
	N <- base::length(x)
	dx <- base::rep(1,N)
	base::return(dx)
}
###################################################