

// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "first_eigenvalue_sirt.h"


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
// #include <c:\Users\robitzsch\Dropbox\Eigene_Projekte\R-Routinen\IRT-Functions\sirt_package\Entwicklung\Eigenvalues\Code\first_eigenvalue_sirt__2.05.h>

// declarations
extern "C" {
SEXP eigenvaluesDsirt( SEXP X_, SEXP D_, SEXP maxit_, SEXP conv_) ;
}

// definition

SEXP eigenvaluesDsirt( SEXP X_, SEXP D_, SEXP maxit_, SEXP conv_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix Xr(X_);          
     int D = as<int>(D_);   
     int maxit = as<int>(maxit_);   
     double conv = as<double>(conv_) ;  
     Rcpp::NumericVector d1(2) ;  
       
     int K=Xr.nrow() ;  
       
     // Rcpp::List res2 ;  
     Rcpp::NumericVector dvec(D) ;  
     arma::mat u(K,D) ;  
       
     Rcpp::List res2;  
       
     //**********  
     // set matrices  
     arma::mat X(Xr.begin(), K, K , false);   
       
     arma::mat X0 = X ;  
       
     for (int dd=0;dd<D;dd++){  
     // int dd = 0;  
     	// estimate first eigenvalue of reduced matrix  
     	res2 = firsteigenvalsirt(X0,maxit,conv,K);  
     	d1 = res2["lambda1"] ;  
      // Rcpp::Rcout << " d1 = " << d1[0] << std::endl ;	  
     	dvec[dd] = d1[0] ;  
     	arma::mat u1 = res2["u"] ;  
     	u.col(dd) = arma::mat( u1.col(0) ) ;  
          for (int ii1=0;ii1<K;ii1++){
          	 X0(ii1,ii1) = X0(ii1,ii1) - dvec[dd] * u(ii1,dd)*u(ii1,dd) ;                     
          for (int ii2=ii1+1;ii2<K;ii2++){  
          	 X0(ii1,ii2) = X0(ii1,ii2) - dvec[dd] * u(ii1,dd)*u(ii2,dd) ;       
             X0(ii2,ii1) = X0(ii1,ii2) ; 
          	     			}  
          			}  
     	}  
                            
     ////////////////////////////////////  
     // OUTPUT:  
     return Rcpp::List::create(  
         _["d"]=dvec  ,   
         _["u"]= u   
            		) ;  
       
       
                   
                   
END_RCPP
}



///////////////////////////////////////////////////

Rcpp::List firsteigenvalsirt2(arma::mat X, int maxit, double conv, int K){

    double lambda_temp ;
    double lambda =0 ;
    double lambdadiff = 1000;
    double lambda_old ;

    Rcpp::NumericVector lambda_est(2);

    int iter = 0 ;

    //**********
    // set matrices
    arma::mat Xz ;


    arma::colvec z(K) ;
    double temp1 = 1 / sqrt( K ) ;
    for (int ii=0;ii<K;ii++){
        z[ii] = temp1 ;
                }

    ///////////////////////////////
    /// algorithm

    // for (int iter=0;iter < maxit;iter ++){
    while ( ( iter < maxit ) & ( lambdadiff > conv) ){
        lambda_old = lambda ;
        Xz = arma::mat( X * z ) ;
        lambda_temp = 0 ;
        for (int ii=0;ii<K;ii++){
            lambda_temp += Xz[ii]*Xz[ii] ;
                    }
        lambda = sqrt( lambda_temp ) ;
        lambdadiff = lambda - lambda_old ;
        if ( lambdadiff < 0 ){ lambdadiff = - lambdadiff ; }
        z = Xz / lambda ;
    // Rcpp::Rcout << "Iteration =" << iter << " lambda = " << lambda << std::endl ; 
    // Rcpp::Rcout << " difference = " << lambdadiff << std::endl ; 
            iter ++ ;
        }

    lambda_est[0] = lambda ;

    ////////////////////////////////////
    // OUTPUT:
    return Rcpp::List::create(
        _["u"]=z , 
        _["lambda1"]=lambda_est 
                ) ;
}






