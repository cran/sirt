

############################################################
# Fit ISOP Model
isop.alg <- function( freq.correct , wgt , eps = .001 , 
		maxit=30 , progress=TRUE ){
    #****
    # initializations
    # monotone regression
    M1 <- as.matrix(freq.correct)
    wgt <- as.matrix(wgt)
    # isotonic row regression
    res <- .fit.isotonic.rows(X=M1 , wgt=wgt )
    RX <- res$RX
    IR <- res$IR
    # isotonic column regression
    res <- .fit.isotonic.cols(X=RX , wgt=wgt )
    CX <- res$CX
    IC <- res$IC    
    X <- M1
    deviation <- 1
    iter <- 0
    # ISOP algorithm
    while( ( iter < maxit) & ( deviation > eps ) ){    
        Xold <- X
        # isotonic rows
        res <- .fit.isotonic.rows(X=X-IC , wgt=wgt )
        RX <- res$RX
        IR <- res$IR    
        # isotonic columns
        res <- .fit.isotonic.cols(X=X-IR , wgt=wgt )
        CX <- res$CX
        IC <- res$IC
        # updated X
        X <- Xold - IR - IC
        # calculate deviation
        deviation <- sum( ( X - Xold )^2*wgt )
		iter <- iter + 1
        if (progress){ 
			cat( "ISOP Model | Iteration" , iter , "- Deviation =" ,  round( deviation , 5 ) , "\n")
			flush.console()			
					}
                    }        
    # deviation criterion
    wgt1 <- ( wgt / colSums( wgt ) ) / ncol(wgt)
    fit <- sqrt( sum( ( M1-CX  )^2 * wgt1  ) )
    # output
    res <- list( "fX" = CX , "ResX" = M1 - CX , "fit" = fit    )
    return(res)
        }
############################################################
#**********************
# fit isotonic rows
.fit.isotonic.rows <- function( X , wgt ){
    M2 <- X
    RR <- nrow(M2)
    for (rr in 1:RR){
        M2[rr,] <- monoreg(x=X[rr,], w=wgt[rr,] )$yf
                }
    RX <- M2
    IX <- X - RX
    res <- list( "X" = X , "RX"  = RX , "IR" = IX )
    return(res)
        }
		
#**********************
# fit isotonic columns
.fit.isotonic.cols <- function( X , wgt ){
    M2 <- X
    RR <- ncol(M2)
    for (rr in 1:RR){
        M2[,rr] <- monoreg(x=X[,rr], w=wgt[,rr] )$yf
                }
    RX <- M2
    IX <- X - RX
    res <- list( "X" = X , "CX"  = RX , "IC" = IX )
    return(res)
        }		