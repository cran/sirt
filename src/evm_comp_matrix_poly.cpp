

// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "evm_eigaux_fcts.h"
// #include "first_eigenvalue_sirt.h"


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes
// #include <P: ... \Entwicklung\Eigenvector_Method\Code\evm_aux_fcts__3.20.h>

// declarations
extern "C" {
SEXP evm_comp_matrix_poly( SEXP dat_, SEXP dat_resp_, SEXP weights_, SEXP JJ_, 
	SEXP jackunits_, SEXP powD_, SEXP progress__, SEXP row_index_, SEXP col_index_) ;
}

// definition

SEXP evm_comp_matrix_poly( SEXP dat_, SEXP dat_resp_, SEXP weights_, SEXP JJ_, 
	SEXP jackunits_, SEXP powD_, SEXP progress__, SEXP row_index_, SEXP col_index_ ){
BEGIN_RCPP
  
       
     Rcpp::NumericMatrix dat(dat_);          
     Rcpp::NumericMatrix dat_resp(dat_resp_);          
     Rcpp::NumericVector weights(weights_) ;  
     int JJ = as<int>(JJ_) ;  
     Rcpp::NumericVector jackunits(jackunits_) ;  
     int powD = as<int>(powD_) ;  
     int progress_ = as<int>(progress__) ;  
     Rcpp::NumericMatrix row_index(row_index_);  
     Rcpp::NumericMatrix col_index(col_index_);  
       
     int N = dat.nrow() ;  
     int I = row_index.nrow() ;  
       
     double K = I ;  
     double conv= 0.0001 ;  
     int maxit=100 ;  
       
     // matrix for dichotomous responses  
     arma::mat B = arma::zeros( I , I  );  
     arma::mat Bjack = arma::zeros( I , I*JJ ) ;  
     arma::mat B2 = arma::zeros( I , I ) ;  
       
     Rcpp::NumericMatrix b_ra_jack(I , JJ ) ;  
     Rcpp::NumericMatrix b_evm_jack(I , JJ ) ;  
     Rcpp::NumericVector lambda_jack(JJ) ;  
     Rcpp::NumericVector cons_index_jack(JJ) ;  
     Rcpp::NumericVector lambda_jack_na(JJ) ;  
       
     //*****************************************  
     // start counting pairwise comparisons  
       
     if ( progress_==1){  
     Rcpp::Rcout << "*** Create pairwise comparison matrix " << std::flush <<  std::endl ;     
     }  
       
     int ii=0;  
     int jj=0;  
       
     for (int rr=0;rr<I;rr++){  
     for (int cc=0;cc<I;cc++){  
         ii = row_index(rr,0) ;  
         jj = col_index(cc,0) ;  
         if ( ii != jj ){  
         for (int nn=0;nn<N;nn++){  
             if ( dat_resp(nn,ii) * dat_resp(nn,jj) == 1 ){  
                 // x_ni = 1 & x_nj = 0  
                 if ( ( dat(nn,ii) ==  row_index(rr,1) ) & ( ( dat(nn,jj) == col_index(cc,1) ) ) ){  
                 B(rr,cc) += weights[nn] ;  
                 Bjack(rr,cc + I*jackunits[nn] ) += weights[nn] ;          
                     }  
                     }  // end if dat_resp(nn,ii)==1 & dat_resp(nn,jj)==1  
             }    // end nn      
             } // end if ii != jj   
         }  // end cc  
      } // end rr  
       
       
       
     ///***************************************  
     // eigenvector method  
       
     if ( progress_==1){  
     Rcpp::Rcout << "*** Perform eigenvector method "  << std::flush <<  std::endl ;     
     }  
       
     //***  original data  
     Rcpp::List res0 = evm_aux(B , I , powD , maxit, conv, K) ;  
     Rcpp::NumericVector lambda = res0["lambda"] ;  
     Rcpp::NumericMatrix D = res0["D"] ;   
     Rcpp::NumericVector cons_index = res0["cons_index"] ;  
     Rcpp::NumericVector b = res0["b"] ;   
       
     Rcpp::NumericVector tmp1 , tmp2 ;  
       
     //*** jackknife  
     int JJadj = JJ ;  
     for (int jj=0;jj<JJ;jj++){  
     // jj=1 ;  
     	B2=0*B2 ;  
     	B2 = B - Bjack.submat( arma::span(0,I-1) , arma::span( jj*I ,  I-1+jj*I  ) ) ;  
     	if (JJ==1){  
     	       B2 = B ;  
     		}  
     	Rcpp::List res2 = evm_aux(B2 , I , powD , maxit, conv, K) ;            
     	Rcpp::NumericVector v1 = res2["b"] ;  b_evm_jack(_,jj) = v1 ;  
     	tmp1 = res2["lambda"] ;  lambda_jack[jj] = tmp1[0] ;  
     	if ( R_IsNaN( tmp1[0] ) ){  
     		lambda_jack_na[jj] = 1 ;  
     		JJadj = JJadj - 1 ;  
     				}  
     	tmp2 = res2["cons_index"] ; cons_index_jack[jj] = tmp2[0] ;  
     	}  
        
       
     //************************************************              
     // statistical inference  
       
     // create matrix with all parameters  
     int VV = B.n_rows + 2 ;  
     Rcpp::NumericMatrix PARS( VV , JJadj ) ;  
     int jj2=0 ;  
     for (int jj=0;jj<JJ;jj++){  
        if ( lambda_jack_na[jj] == 0 ){	  
        	   for (int pp=0;pp<VV-2;pp++){	  
        		PARS(pp,jj2) = b_evm_jack(pp,jj) ;  
			PARS(VV-2,jj2 ) = lambda_jack[jj] ;  
			PARS(VV-1,jj2 ) = cons_index_jack[jj] ;  
     				}  
     		jj2 ++ ;				  
     			}  
     			}  
     // inference based on Jackknife  
     Rcpp::List res1 =  parameters_jackknife(PARS) ;			  
     Rcpp::NumericMatrix PARS_vcov = res1["PARS_vcov"] ;  
     Rcpp::NumericVector PARS_means = res1["PARS_means"] ;  
     		  
     //*************************************************      
     // OUTPUT              
                   
     return Rcpp::List::create(    
         _["D"] = D ,    
         _["B"] = B ,  
     //    _["Bjack"] = Bjack ,     
         _["b_evm"] = b ,  
         _["b_evm_jack"] = b_evm_jack ,  
         _["lambda"] = lambda ,  
         _["lambda_jack"] = lambda_jack ,   
         _["lambda_jack_na"] = lambda_jack_na ,   
         _["cons_index"] = cons_index ,  
         _["cons_index_jack"] = cons_index_jack ,      
         _["JJ"] = JJ ,  
         _["JJadj"] = JJadj ,   
         _["PARS_means"] = PARS_means ,  
         _["PARS_vcov"] = PARS_vcov         
         ) ;  
END_RCPP
}



