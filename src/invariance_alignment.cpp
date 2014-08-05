

// includes from the plugin
#include <RcppArmadillo.h>
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
SEXP ia_optim_lambda( SEXP lambda_, SEXP psi0_, SEXP psi0b_, SEXP align_scale_, 
		SEXP align_pow_, SEXP wgt_, SEXP eps_, SEXP group_combis_) ;
}

// definition

SEXP ia_optim_lambda( SEXP lambda_, SEXP psi0_, SEXP psi0b_, SEXP align_scale_, 
	SEXP align_pow_, SEXP wgt_, SEXP eps_, SEXP group_combis_ ){
BEGIN_RCPP
   
       
     // align.optim.lambda( lambda=lambda , psi0=psi0 , psi0b=psi0b ,  
     //                    align.scale=align.scale[1] , align.pow=align.pow[1] ,  
     //  wgt, eps=eps,group.combis)      
        
     Rcpp::NumericMatrix lambda(lambda_);  
     Rcpp::NumericVector psi0(psi0_);  
     Rcpp::NumericVector psi0b(psi0b_);  
     double align_scale = as<double>(align_scale_) ;  
     double align_pow = as<double>(align_pow_) ;  
     Rcpp::NumericMatrix wgt(wgt_);  
     double eps = as<double>(eps_) ;  
     Rcpp::NumericMatrix group_combis(group_combis_);  
       
       
     //    # optimization with respect to country SDs  
     //    lambda1 <- lambda / psi0  
     //    lambda1b <- lambda / psi0b      
     //    # function  
     //    fopt <- 0  
     //	I <- ncol(lambda)  
       
     int I = lambda.ncol() ;  
     int G = lambda.nrow() ;  
     Rcpp::NumericMatrix lambda1(G,I);  
     Rcpp::NumericMatrix lambda1b(G,I);  
       
     for (int ii=0;ii<I;ii++){  
     for (int gg=0;gg<G;gg++){  
     	lambda1(gg,ii) = lambda(gg,ii) / psi0[gg] ;  
             lambda1b(gg,ii) = lambda(gg,ii) / psi0b[gg] ;	  
     			}  
     		}  
       
     		  
       
     int GC = group_combis.nrow();  
     Rcpp::NumericVector fopt(GC);  
     Rcpp::NumericVector fopt1(GC);  
       
     for (int ii=0;ii<I;ii++){  
     for (int cc=0;cc<GC;cc++){  
     	//    fopt1 <- ( lambda1[ group.combis[,1] , ii ] - lambda1b[ group.combis[,2] , ii ] )^2  
     	fopt1[cc] = pow( lambda1( group_combis(cc,0) , ii) - lambda1b( group_combis(cc,1) , ii) , 2.0 ) ;  
     	//   fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] *   
     	//                  ( fopt1 / align.scale^2 + eps )^align.pow  
     	fopt[cc] += wgt( group_combis(cc,0) , ii) * wgt( group_combis(cc,1) , ii) *  
     			   pow( fopt1[cc] / ( align_scale*align_scale ) + eps , align_pow ) ;  
     	}  
     }  
       
     // sum over the same indices  
     Rcpp::NumericVector res(G);  
     for (int cc=0;cc<GC;cc++){  
     	res[ group_combis(cc,0) ] += fopt[cc] ;		  
     	}  
       
     return wrap(res) ;	  
     	  
     //*************************************************      
     // OUTPUT              
               
     //return Rcpp::List::create(   
     //    _["lambda"] = lambda  ,       
      //   _["fopt"] = fopt ,  
      //   ) ;    
       
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}

// declarations
extern "C" {
SEXP ia_optim_nu( SEXP lambda_, SEXP nu_, SEXP psi0_, SEXP psi0b_, SEXP alpha0_, 
	SEXP alpha0b_, SEXP align_scale_, SEXP align_pow_, SEXP wgt_, SEXP eps_, 
	SEXP group_combis_) ;
}

// definition

SEXP ia_optim_nu( SEXP lambda_, SEXP nu_, SEXP psi0_, SEXP psi0b_, SEXP alpha0_, 
	SEXP alpha0b_, SEXP align_scale_, SEXP align_pow_, SEXP wgt_, SEXP eps_, 
	SEXP group_combis_ ){
BEGIN_RCPP
   
       
     // align.optim.lambda( lambda=lambda , psi0=psi0 , psi0b=psi0b ,  
     //                    align.scale=align.scale[1] , align.pow=align.pow[1] ,  
     //  wgt, eps=eps,group.combis)      
        
     Rcpp::NumericMatrix lambda(lambda_);  
     Rcpp::NumericMatrix nu(nu_);  
     Rcpp::NumericVector psi0(psi0_);  
     Rcpp::NumericVector psi0b(psi0b_);  
     Rcpp::NumericVector alpha0(alpha0_);  
     Rcpp::NumericVector alpha0b(alpha0b_);  
     double align_scale = as<double>(align_scale_) ;  
     double align_pow = as<double>(align_pow_) ;  
     Rcpp::NumericMatrix wgt(wgt_);  
     double eps = as<double>(eps_) ;  
     Rcpp::NumericMatrix group_combis(group_combis_);  
       
       
     //    # optimization with respect to country SDs  
     //    nu1 <- nu - alpha0 * lambda   
     //    nu1b <- nu - alpha0b * lambda  
           
     int I = lambda.ncol() ;  
     int G = lambda.nrow() ;  
     Rcpp::NumericMatrix nu1(G,I);  
     Rcpp::NumericMatrix nu1b(G,I);  
       
     for (int ii=0;ii<I;ii++){  
     for (int gg=0;gg<G;gg++){  
     	nu1(gg,ii) = nu(gg,ii) - lambda(gg,ii) * alpha0[gg] ;  
             nu1b(gg,ii) = nu(gg,ii) - lambda(gg,ii) * alpha0b[gg] ;	  
     			}  
     		}  
       
     		  
       
     int GC = group_combis.nrow();  
     Rcpp::NumericVector fopt(GC);  
     Rcpp::NumericVector fopt1(GC);  
       
     for (int ii=0;ii<I;ii++){  
     for (int cc=0;cc<GC;cc++){  
     	//    fopt1 <- ( lambda1[ group.combis[,1] , ii ] - lambda1b[ group.combis[,2] , ii ] )^2  
     	fopt1[cc] = pow( nu1( group_combis(cc,0) , ii) - nu1b( group_combis(cc,1) , ii) , 2.0 ) ;  
     	//   fopt <- fopt + wgt[ group.combis[,1],ii] * wgt[ group.combis[,2],ii] *   
     	//                  ( fopt1 / align.scale^2 + eps )^align.pow  
     	fopt[cc] += wgt( group_combis(cc,0) , ii) * wgt( group_combis(cc,1) , ii) *  
     			   pow( fopt1[cc] / ( align_scale*align_scale ) + eps , align_pow ) ;  
     	}  
     }  
       
     // sum over the same indices  
     Rcpp::NumericVector res(G);  
     for (int cc=0;cc<GC;cc++){  
     	res[ group_combis(cc,0) ] += fopt[cc] ;		  
     	}  
       
     return wrap(res) ;	  
     	  
     //*************************************************      
     // OUTPUT              
               
     //return Rcpp::List::create(   
     //    _["lambda"] = lambda  ,       
      //   _["fopt"] = fopt ,  
      //   ) ;    
       
     // maximal list length is 20!  
       
       
     // Rcpp::Rcout << "tmp1 " <<  tmp1 <<  std::flush << std::endl ;  
       
       
     
END_RCPP
}





