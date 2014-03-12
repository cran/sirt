//***  Code created: 2014-02-10 08:53:56
//***  gooijer_csn_table__2.01.Cpp


// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gooijercsntableaux.h"

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;




// declarations
extern "C" {
SEXP gooijer_csn_table( SEXP dat_, SEXP dat_perm_, SEXP RR_, SEXP NS_, 
	SEXP progress_, SEXP progress_vec_, SEXP score_index_) ;
}

// definition

SEXP gooijer_csn_table( SEXP dat_, SEXP dat_perm_, SEXP RR_, SEXP NS_, 
	SEXP progress_, SEXP progress_vec_, SEXP score_index_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix dat(dat_);  
     Rcpp::NumericMatrix dat_perm(dat_perm_);  
     int RR=as<int>(RR_) ;  
     int NS=as<int>(NS_) ;  
     int progress=as<int>(progress_) ;  
     Rcpp::NumericVector progress_vec(progress_vec_)  ;  
     Rcpp::NumericMatrix score_index(score_index_);  
       
       
     int N=dat.nrow() ;  
     Rcpp::NumericVector stat_perm(RR) ;  
       
     RNGScope scope;  
       
       
     Rcpp::NumericMatrix sampleM = dat_perm ;  
       
     // compute Gooijer statistic for original data  
     Rcpp::NumericVector stat = gta(dat) ;  
       
     //******* ;  
       
     Rcpp::NumericVector s1(1);  
     int zz=0;  
     if ( progress==1){  
     	Rcpp::Rcout << "|" <<   std::flush ;  
     			}  
     for (int rr=0;rr<RR;rr++){  
     	// permutation of original data ;  
     	sampleM = gooijer_permutation(sampleM , NS , N , score_index ) ;  
     	// compute statistic for permuted data ;  
     	s1 = gta(sampleM ) ;  
     	stat_perm[rr] = s1[0] ;  
     	if ( (progress==1) & ( rr==progress_vec[zz] ) ){  
     		zz = zz+1 ;  
     		if ( zz==10){	zz = 9 ; }   
     		Rcpp::Rcout << "-" <<   std::flush ;  
     				}  
     		}  
     if ( progress==1){  
     	Rcpp::Rcout << "|" <<   std::flush << std::endl ;  
     			}			  
      		  
        		  
     //*************************************************      
     // OUTPUT              
                   
      return Rcpp::List::create(    
         _["stat"] = stat ,    
         _["stat_perm"]=stat_perm  
         ) ;  
END_RCPP
}



