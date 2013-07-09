//  Code created: 2013-07-08 14:28:51


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
SEXP RM_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP K) ;
}

// definition

SEXP RM_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP K ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
     Rcpp::NumericVector KK(K);  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     int KKK = KK[0] + 1 ;   
       
     //*****  
     // calculate individual likelihood  
     NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( KKK*ii + DAT2(nn,ii) , tt ) ;  
                         }  
                     }  
                 }  
             }  
       
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return List::create(_["fyiqk"] = fyiqk );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}



