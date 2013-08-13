//  Code created: 2013-07-27 19:08:35


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
SEXP rowMaxsCPP_source( SEXP matr) ;
}

// definition

SEXP rowMaxsCPP_source( SEXP matr ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
       
     // create output vectors  
     NumericVector MAXVAL (NR) ;  
     NumericVector MAXIND (NR) ;  
     MAXIND.fill(1);  
       
     for (int nn=0;nn<NR;++nn){  
          MAXVAL[nn] = MATR( nn , 0 ) ;  
     	for (int cc=1 ; cc < NC ; ++cc ){  
     	    if ( MATR(nn,cc) > MAXVAL[nn] ){  
     	    	    MAXVAL[nn] = MATR(nn,cc) ;  
     	    	    MAXIND[nn] = cc + 1 ;  
     	    	           }  
     		}  
     	}  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     // return( wrap(prob) );  
     return List::create(_["maxval"] = MAXVAL , _["maxind"]=MAXIND ) ;     
     
END_RCPP
}



