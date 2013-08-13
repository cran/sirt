//  Code created: 2013-08-13 10:46:15


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
SEXP monoreg_rowwise_Cpp( SEXP yM, SEXP wM) ;
}

// definition

SEXP monoreg_rowwise_Cpp( SEXP yM, SEXP wM ){
BEGIN_RCPP
  
     //********************************************  
     // Adapted code from monoreg function  
     // contained in the fdrtool package  
     // original author: Korbinian Strimmer  
     //********************************************  
       
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix YM(yM);  
     Rcpp::NumericMatrix WM(wM) ;   
       
     // define row and column numbers  
     int NR=YM.nrow();  
     int NC=YM.ncol();  
     int nn=NC ;  
       
     // create output ghat matrix  
     NumericMatrix ghat (NR,NC) ;  
     NumericMatrix k (NR,NC) ;  
     NumericMatrix gew (NR,NC) ;  
              
     double neu ;  
     int zz ;   
     int c=0 ;  
     int j=0 ;  
       
     for (zz=0; zz < NR ; zz++){ // begin for zz  
     	c=0 ;  
     	j=0 ;	  
     	nn = NC ;  
     	k(zz,c) = 0;  
     	gew(zz,c) = WM(zz,0);  
     	ghat(zz,c) = YM(zz,0);  
     	  
     	//######  
     	for (j=1; j < nn; j++){   // begin for j ...		  
     	    c = c+1;  
     	    k(zz,c) = j;  
     	    gew(zz,c) = WM(zz,j);  
     	    ghat(zz,c) = YM(zz,j);  
     	    /* c is at least 1 as nn is > 1 */  
     	    while (ghat(zz,c-1) >= ghat(zz,c))  
     	    {  
     	      neu = gew(zz,c)+gew(zz,c-1);  
     	      ghat(zz,c-1) = ghat(zz,c-1)+(gew(zz,c)/neu)*(ghat(zz,c)-ghat(zz,c-1) );  
     	      gew(zz,c-1) = neu;  
     	      c = c-1;  
     	      if (c==0) break;  
     	    }  
     	  }     // end for j ...  
     	//##########    
     	while (nn >= 1){  
     	    for (j=k(zz,c); j < nn; j++){  
     	      ghat(zz,j) = ghat(zz,c);  
     	    		}  
     	    nn = k(zz,c);  
     	    c = c-1;  
     	  }	  	    
     } // end for zz  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     return( wrap(ghat) );  
               
END_RCPP
}



