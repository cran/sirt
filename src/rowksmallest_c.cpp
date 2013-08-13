//  Code created: 2013-07-28 12:17:03


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
SEXP rowKSmallest_C( SEXP matr, SEXP K, SEXP indexmatr, SEXP rnmatr) ;
}

// definition

SEXP rowKSmallest_C( SEXP matr, SEXP K, SEXP indexmatr, SEXP rnmatr ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     Rcpp::IntegerVector KK1(K) ;   
     Rcpp::NumericMatrix INDEXMATR(indexmatr);  
     Rcpp::NumericMatrix RNMATR(rnmatr);  
       
     // define row and column numbers  
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
     int KK = KK1[0] ;   
       
     // copy original matrix ;  
     Rcpp::NumericMatrix MATRK = Rcpp::clone(MATR);  
       
     // create output vectors  
     NumericMatrix SMALLVAL (NR,KK) ;  
     NumericMatrix SMALLIND (NR,KK) ;  
     // SMALLIND.fill(1);  
       
     // int kk=0 ;  
     double tmp1 ;  // temp value  
       
     for (int kk=0;kk<KK;++kk){ // begin for kk  
     for (int nn=0;nn<NR;++nn){ // begin for nn  
          SMALLIND(nn,kk) = kk+1 ;	  
          SMALLVAL(nn,kk) = MATRK( nn , kk ) ;  
     	for (int cc=kk+1 ; cc < NC ; ++cc ){  // begin for cc  
     	    if ( ( MATRK(nn,cc) < SMALLVAL(nn,kk) ) ||  
     	         (  ( MATRK(nn,cc) == SMALLVAL(nn,kk) ) && ( RNMATR(nn,cc)==1 ) )   
     	    		){ // begin comparison  
     	    	    tmp1 = SMALLVAL(nn,kk) ;  
     	    	    SMALLVAL(nn,kk) = MATRK(nn,cc) ;  
     	    	    MATRK(nn,kk) = MATRK( nn , cc );  
     	    	    MATRK(nn,cc) = tmp1 ;  
     	    	    tmp1 = SMALLIND(nn,kk) ;  
     	    	    SMALLIND(nn,kk) = INDEXMATR(nn,cc) ;  
     	    	    INDEXMATR(nn,kk) = INDEXMATR(nn,cc) ;  
     	    	    INDEXMATR(nn,cc)=tmp1 ;  
     	    	           } // end comparison  
     		}  // end for cc  
     	}   // end for nn  
       }  // end for kk  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     // return( wrap(prob) );  
     // return List::create(_["smallval"]=SMALLVAL , _["smallind"]=SMALLIND ,  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
      return List::create(_["smallval"]=SMALLVAL , _["smallind"]=SMALLIND ) ;     
     
END_RCPP
}



