//  Code created: 2013-04-15 14:35:45


// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
}

// definition

SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix A(a);  
     Rcpp::NumericMatrix B(b);  
     Rcpp::NumericMatrix QQ(Q);  
     Rcpp::NumericMatrix THETA(thetak);  
     Rcpp::NumericVector CC(cc);  
     Rcpp::NumericVector DD(dd);  
       
     int I=A.nrow();  
     int D=A.ncol();  
     int TP=THETA.nrow();  
       
     // create matrix of probabilities  
     NumericMatrix prob (I,TP) ;  
     prob.fill(1);  
       
     NumericVector yy (TP);  
       
     for (int ii=0;ii<I;++ii){  
     for (int dd=0;dd<D;++dd){  
         if ( QQ(ii,dd)>0 ){  
             for (int tt=0; tt<TP; ++tt){  
                    yy[tt] =  A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)-B(ii,dd)   ;  
                         }   // end tt  
             NumericVector p1=plogis(yy) ;  
             for (int tt=0; tt<TP; ++tt){              
                 prob(ii,tt)=prob(ii,tt)*p1[tt] ;  
                   }  // end tt  
                 }       // end if Q(ii,dd)>0  
             }          // end dd   
     //***  
     // include guessing and slipping parameters  
      if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){   
         for (int tt=0;tt<TP;++tt){      
             prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;  
                     }   // end tt  
                 } // end if condition for guessing or slipping                  
         }       // end ii  
           
     ///////////////////////////////////////  
     /// OUTPUT      
     return(wrap(prob)) ;               
     // return( wrap(prob) );  
     
END_RCPP
}



