

// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "evm_eigenvals2.h"

using namespace Rcpp;


// declarations
//extern "C" {
// Rcpp::List firsteigenvalsirt2(arma::mat X, int maxit, double conv, double K) ;
// Rcpp::NumericVector choppin_rowaveraging( arma::mat B , int I , double priorweight) ;
// Rcpp::List evm_aux( arma::mat B , int I , int powD , int maxit, double conv, double K ) ;
// Rcpp::List parameters_jackknife( Rcpp::NumericMatrix PARS ) ; 
// } 

 
//**********************************************************************
// compute means and covariances of estimators obtained by Jackknife
Rcpp::List parameters_jackknife( Rcpp::NumericMatrix PARS ){  
	int VV=PARS.nrow() ;
	int JJ=PARS.ncol() ;
	Rcpp::NumericVector PARS_means(VV);
	Rcpp::NumericMatrix PARS_vcov(VV,VV);
	double tmp3=0;
	// compute row means
	for (int vv=0;vv<VV;vv++){
		tmp3=0;
		//int vv=0 ;
		for (int jj=0;jj<JJ;jj++){
		      tmp3+= PARS(vv,jj) ;
				}
		PARS_means[vv] = tmp3 / JJ ;  
		PARS(vv,_) = PARS(vv,_) - PARS_means[vv] ; 
			}
	
	// compute covariance
	for (int vv1=0;vv1<VV;vv1++){
	for (int vv2=vv1;vv2<VV;vv2++){	
	for (int jj=0;jj<JJ;jj++){
	   PARS_vcov(vv1,vv2) += PARS(vv1,jj)*PARS(vv2,jj) ;
//	   if (vv1!=vv2){
//		PARS_vcov(vv2,vv1) += PARS_vcov(vv1,vv2) ;
//			}
		} // end jj
	PARS_vcov(vv1,vv2) = PARS_vcov(vv1,vv2) * (JJ-1) / JJ ;
	if (vv1!=vv2){
//	   PARS_vcov(vv2,vv1) = PARS_vcov(vv2,vv1) * (JJ-1) / JJ ;
           PARS_vcov(vv2,vv1) = PARS_vcov(vv1,vv2) ;
		}
		}
	 } 
	
	    return Rcpp::List::create(
		Rcpp::_["PARS_means"]= PARS_means,
		Rcpp::_["PARS_vcov"]= PARS_vcov
			) ;
	}

//************************************************************
// Eigenvector method
Rcpp::List evm_aux( arma::mat B , int I , int powD ,
	  int maxit, double conv, double K ){
    arma::mat TMP = B ;    
    arma::mat D= arma::zeros(I,I);
    for ( int hh=0;hh<(powD-1) ; hh++){
		TMP = TMP * B ;
			}   			
    //  calculate ratios in D
    for (int ii=0;ii<I;ii++){
        D(ii,ii) = 1 ;
    for (int jj=0;jj<I;jj++){
        if ( ii!= jj){
        D(ii,jj) = TMP(jj,ii) / TMP(ii,jj)    ;
                    }
                        }
                    }   
    // compute first eigenvalue                
    Rcpp::List res11 = firsteigenvalsirt2( D, maxit, conv, K) ;
    
    Rcpp::NumericVector u = res11["u"] ;
    double tmp1= 0 ;       

//    Rcpp::Rcout << "I " << I << "  u[0]" << u[0] << "  " <<   std::endl ;     
    // center difficulty vector                
    for ( int ii = 0 ; ii<I ; ii++){
	u[ii] = log( u[ii] ) ;
	 tmp1 += u[ii] ;
	 }
    Rcpp::NumericVector b = u - tmp1 / I ;       

// Rcpp::Rcout << "a100  " <<   std::endl ;      
    // compute eigenvalue and consistency index
    Rcpp::NumericVector tmp2 = res11["lambda1"] ;
   double lambda = tmp2[0] ; 
   double cons_index = ( lambda - I )   / ( I -1 ) ;    
 
// Rcpp::Rcout << "a200  " <<   std::endl ;      

    return Rcpp::List::create(
        Rcpp::_["lambda"]= lambda ,
        Rcpp::_["D"] = D , 
        Rcpp::_["cons_index"] = cons_index ,
      	Rcpp::_["b"] = b 
                ) ;   
    	}
    	
//***************************************                    
// Choppin's row averaging approach                    
Rcpp::NumericVector choppin_rowaveraging( arma::mat B , int I , double priorweight){ 
    Rcpp::NumericVector b_ra(I) ;
    arma::mat TMP2= arma::zeros(I,I) ;
    B = B + priorweight ;    
    //  calculate ratios in D
    for (int ii=0;ii<I;ii++){
    for (int jj=0;jj<I;jj++){
        if (ii!=jj){ 
            TMP2(ii,jj) = log(  B(jj,ii)  /  B(ii,jj)  )  ;    
                        }         
                       }
                   }                               
    for (int ii=0;ii<I;ii++){
    for (int jj=0;jj<I;jj++){
       b_ra[ii] += TMP2(ii,jj ) ;
                }    
       b_ra[ii] = b_ra[ii] / I ;            
            }
    return( wrap( b_ra) ) ;
        }
//**********************************************








// Rcpp::List firsteigenvalsirt2(arma::mat X, int maxit, double conv, double K){
