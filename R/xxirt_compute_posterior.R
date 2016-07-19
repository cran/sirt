

###########################################################################
xxirt_compute_posterior <- function( prior_Theta , p.xi.aj , group ,
                 G , weights , dat1 , dat_resp , maxK , group_index ){
		N <- base::nrow(dat_resp)
		TP <- base::ncol(p.xi.aj)	
		I <- base::ncol(dat1)		
		# posterior distribution
		prior1 <- base::t( prior_Theta[ , group ] )
		p1 <- p.aj.xi <- prior1 * p.xi.aj
		p.aj.xi <- p.aj.xi / base::rowSums( p.aj.xi )
		# expected counts
		n.ik <- base::array( 0 , dim=c(I,maxK , TP,G) )		
		N.ik <- base::array( 0 , dim=c(I,maxK , TP) )		
		pi.k <- base::matrix( 0 , nrow=TP , ncol=G )		
		for (gg in 1:G){
			# gg <- 1
			ind_gg <- group_index[[gg]]	
			for (ii in 1:I){
				for (kk in 1:maxK){
					# ii <- 1
					# kk <- 1
					v1 <- weights[ind_gg] * p.aj.xi[ind_gg , ] * 
						   ( dat1[ind_gg , ii ] == (kk-1) ) * ( dat_resp[ind_gg , ii ] )
					n.ik[ ii , kk , , gg ] <- base::colSums( v1 )
				}  # end kk
			}  # end ii
			N.ik <- N.ik + n.ik[,,,gg]
			pi.k[,gg] <- base::colSums( p.aj.xi[ ind_gg , ] * weights[ ind_gg ] )
		}  # end gg
		res <- base::list( p.aj.xi = p.aj.xi , n.ik = n.ik , N.ik = N.ik, N.k = pi.k ,
						post_unnorm = p1 )
		base::return(res)
}
###########################################################################			