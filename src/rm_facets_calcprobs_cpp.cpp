

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
SEXP rm_facets_calcprobs_cpp( SEXP b_item_, SEXP b_rater_, SEXP Qmatrix_, 
	SEXP tau_item_, SEXP K_, SEXP I_, SEXP TP_, SEXP a_item_, SEXP a_rater_, 
	SEXP item_index_, SEXP rater_index_, SEXP theta_k_ ) ;
}

// definition

SEXP rm_facets_calcprobs_cpp( SEXP b_item_, SEXP b_rater_, SEXP Qmatrix_, 
	SEXP tau_item_, SEXP K_, SEXP I_, SEXP TP_, SEXP a_item_, SEXP a_rater_, 
	SEXP item_index_, SEXP rater_index_, SEXP theta_k_ ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericVector b_item(b_item_);  
     Rcpp::NumericVector b_rater(b_rater_);   
     Rcpp::NumericMatrix Qmatrix(Qmatrix_);  
     Rcpp::NumericMatrix tau_item(tau_item_);  
     // int VV = as<int>(VV_);  
     int K = as<int>(K_);  
     int I = as<int>(I_);  
     int TP = as<int>(TP_);  
     Rcpp::NumericVector a_item(a_item_);  
     Rcpp::NumericVector a_rater(a_rater_);   
     Rcpp::NumericVector item_index(item_index_);   
     Rcpp::NumericVector rater_index(rater_index_);  
     Rcpp::NumericVector theta_k(theta_k_);  
     // int RR = as<int>(RR_);  
     //    probs <- .rm.facets.calcprobs( b.item , b.rater , Qmatrix , tau.item ,  
     //           VV , K , I , TP , a.item , a.rater , item.index , rater.index ,  
     //           theta.k ,RR )       
       
     //***** calculate b  
     // b <- tau.item[ item.index , ]  
     Rcpp::NumericMatrix b(I,K) ;  
     // b0 <- ( matrix( b.rater , nrow= RR , ncol=K) )[ rater.index , ] * 	Qmatrix[ item.index ,]	   
     // b <- b + b0  
     for (int ii=0; ii<I ; ii++){  
        b.row(ii) = tau_item.row( item_index[ii] ) ;  
        for (int kk=0;kk<K;kk++){   	     
         b(ii,kk) = b(ii,kk) + b_rater[ rater_index[ii] ] * Qmatrix( item_index[ii] , kk ) ;  
         			}  
         		}  
       
     //****** calculate a    		  
     // a <- a.item[ item.index ] * a.rater[ rater.index ]    		  
     Rcpp::NumericVector a(I) ;  
     for (int ii=0;ii<I;ii++){  
     	a[ii] = a_item[ item_index[ii] ] * a_rater[ rater_index[ii] ] ;  
     			}  
     //******* calculate modified Q-matrix  
     // Qmatrix=Qmatrix[item.index,]  
     Rcpp::NumericMatrix Q(I,K) ;  
     for (int ii=0;ii<I;ii++){  
     	Q.row(ii) = Qmatrix.row( item_index[ii] ) ;  
     			}  
       
     ///**************************  
     // compute response probabilities according to the generalized partial credit model  
       
     //     probs <- array( 0 , dim=c(I,K+1,TP) )   # categories 0 , ... , K  
     Rcpp::NumericMatrix probs(I,(K+1)*TP) ;  
     //    for (kk in 1:K){  
     //        l0 <- matrix( - b[,kk] , nrow=I,ncol=TP)  
     //        l0 <- l0 + outer( a * Qmatrix[ , kk] , theta.k )  
     //        probs[,kk+1,] <- l0  
     //                }  
       
     // probs(ii , kk , tt ) ~ probs( ii , kk + tt * (K+1) )  
       
       
     double tmp1 = 0 ;   
       
     for (int tt=0;tt<TP;tt++){  
     for (int ii=0;ii<I;ii++){  
     	tmp1 = 1 ;	  
     	probs( ii , tt*(K+1) ) = 1 ;  
     	for (int kk=0;kk<K;kk++){  
     		probs( ii , kk+1 + tt*(K+1) ) = exp( - b(ii,kk) + a[ii] * Q(ii,kk) * theta_k[tt] );  
     		tmp1 += probs( ii , kk+1 + tt*(K+1) ) ;  
     				}  
     	for (int kk=0;kk<K+1;kk++){  
     		probs( ii , kk + tt*(K+1) ) = probs( ii , kk + tt*(K+1) ) / tmp1 ;  
     				}								  
     		}   // end ii  
     	} // end tt  
       
     return wrap( probs ) ;   
     	  
     ///////////////////////////////////////  
     /// OUTPUT                  
       
     // return List::create(  
     //	        _["b"] = b ,  
     //		_["a"]= a ,  
     //		_["probs"] = probs  
     //			) ;  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
     // return List::create(_["yM"]=YM , _["wM"]=WM ) ;     
     
END_RCPP
}



