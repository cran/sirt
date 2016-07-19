
###########################################
# Function for analyzing response patterns
md.pattern.sirt <- function(dat){
    dat <- base::as.matrix(dat)
	if ( base::ncol(dat)>1000 ){
	   base::stop("Function only works for datasets with fewer than 1000 variables!\n")
	}
    res <- md_pattern_rcpp( dat_=dat )
	rp_unique <- base::unique(res$unique_resp_patt)
	res$unique_resp_patt <- base::match( res$unique_resp_patt , rp_unique )	
	res$resp_patt <- base::match( res$resp_patt , rp_unique )
    res$dat.ordered <- res$dat[ base::order( res$resp_patt ) , ]
    base::return(res)
}
#*******************************************			
# calling the Rcpp function
md_pattern_rcpp <- function (dat_){ 
	base::.Call("md_pattern_csource", dat_ , PACKAGE = "sirt")
}			
#*******************************************					