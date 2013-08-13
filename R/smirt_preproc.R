
################################################################
	# define estimation functions depending on condensation type
	#-#-#  calculation of probabilities
	.smirt.calcprob <- function( a , b, Q, thetak , c , d , irtmodel ){
			if ( irtmodel=="noncomp"){
				res <- calcprob.noncomp( a , b, Q, thetak , c , d )
									}
			if ( irtmodel=="comp"){
				res <- calcprob.comp( a , b, Q, thetak , c , d )
									}																																			
			return(res)
								}
################################################################								
	#-#-#	estimation of b parameters							 
	.smirt.est.b <- function(   b , a , c , d , Qmatrix , est.b , theta.k , 
			n.ik , I , K , TP , D ,  numdiff.parm, 
			max.increment ,	msteps ,  mstepconv , irtmodel , increment.factor ){
			if ( irtmodel=="noncomp"){			
			res <-	.smirt.est.b.noncomp(   b , a , c , d , Qmatrix , est.b , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment ,	msteps ,  mstepconv , increment.factor)	
								}
			if ( irtmodel=="comp"){			
			res <-	.smirt.est.b.comp(   b , a , c , d , Qmatrix , est.b , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment ,	msteps ,  mstepconv , increment.factor)	
								}														
			return(res)
									}
################################################################
	#-#-#   estimation of a parameters								
	.smirt.est.a <- function(  b , a , c , d , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , irtmodel  , increment.factor){		
			if ( irtmodel=="noncomp"){				
			 res <- .smirt.est.a.noncomp(  b , a , c , d , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , increment.factor)		
									}
			if ( irtmodel=="comp"){				
			 res <- .smirt.est.a.comp(  b , a , c , d , Qmatrix , est.a , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.a.increment, msteps ,  mstepconv , increment.factor)		
									}															
			return(res)
						}
						
################################################################						
	#-#-#   estimation of c parameters										
	.smirt.est.c <- function(  b , a , c , d , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , irtmodel  , increment.factor){		
			if ( irtmodel=="noncomp"){				
			 res <- .smirt.est.c.noncomp(  b , a , c , d , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
								}
			if ( irtmodel=="comp"){				
			 res <- .smirt.est.c.comp(  b , a , c , d , Qmatrix , est.c , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
								}															
			return(res)
						}
						
################################################################						
	#-#-#   estimation of d parameters										
	.smirt.est.d <- function(  b , a , c , d , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , irtmodel , increment.factor ){		
			if ( irtmodel=="noncomp"){
			 res <- .smirt.est.d.noncomp(  b , a , c , d , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
									}
			if ( irtmodel=="comp"){
			 res <- .smirt.est.d.comp(  b , a , c , d , Qmatrix , est.d , theta.k , 
				n.ik , I , K , TP , D ,  numdiff.parm, 
				max.increment, msteps ,  mstepconv , increment.factor)		
									}																												
			return(res)
						}