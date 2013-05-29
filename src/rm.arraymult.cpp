
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
SEXP file211c19faab( SEXP A, SEXP dimA, SEXP B, SEXP dimB) ;
}

// definition

SEXP file211c19faab( SEXP A, SEXP dimA, SEXP B, SEXP dimB ){
BEGIN_RCPP
// define input and arrays
      // array A
      Rcpp::NumericMatrix AA(A);
      Rcpp::IntegerVector dimAA(dimA);
      // array B
      Rcpp::NumericMatrix BB(B);
      Rcpp::IntegerVector dimBB(dimB);
      Rcpp::NumericMatrix CC( dimAA[0]*dimAA[1] , dimBB[2] ) ;
      
      int a1 = dimAA[0];
      int a2 = dimAA[1];
      int a3 = dimAA[2];
      int b3 = dimBB[2]; 
      
      // ii -> loop within a matrix
      //for (int ii=0 ; ii < a1 ; ++ii ){
      for (int zz=0 ; zz < a2 ; ++zz){
          for (int ii=0 ; ii < a1 ; ++ii ){        
              for (int hh=0 ; hh < b3 ; ++hh){ // loop over columns
                  for (int kk=0 ; kk < a3 ; ++kk){
                      CC(ii+zz*a1,hh) += AA(ii+zz*a1,kk)*BB(ii+a1*kk,hh)  ; // *BB(kk,hh) ;
                              }
                          }
                      }
                  }
      
      // output
      return( wrap(CC) );
END_RCPP
}



