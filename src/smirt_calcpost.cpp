//  Code created: 2013-04-15 14:31:27


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
SEXP SMIRT_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP dat2ind, SEXP pik, SEXP K) ;
}

// definition

SEXP SMIRT_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP dat2ind, SEXP pik, SEXP K ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
     Rcpp::NumericMatrix DAT2IND(dat2ind);  
     Rcpp::NumericVector PIK(pik);  
     Rcpp::NumericVector KK1(K) ;  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     int KK=KK1[0] ;   
       
     //*****  
     // calculate individual likelihood  
     NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;  
                         }  
                     }  
                 }  
             }  
       
     //****  
     // calculate posterior	  
     NumericMatrix fqkyi (N,TP) ;  
     for (int nn=0;nn<N;++nn){  
     	double total = 0; 	  
     	for (int tt=0;tt<TP;++tt){  
     		fqkyi(nn,tt) = fyiqk(nn,tt)*PIK[tt] ;		  
     		total += fqkyi(nn,tt) ;  
     			}  
     	for (int tt=0;tt<TP;++tt){  
     		fqkyi(nn,tt) = fqkyi(nn,tt)/total ;  
     			}  
     		}  
     //*****  
     // calculate counts  
     for (int tt=0;tt<TP ;++tt){  
     	PIK[tt] = 0 ;  
     	for (int nn=0;nn<N;++nn){	  
     		PIK[tt] += fqkyi(nn,tt) ;  
     				}  
     	PIK[tt] = PIK[tt] / N ;   
     			}  
       
     NumericMatrix nik (TP, I*(KK+1)) ;  
     NumericMatrix NIK (TP, I) ;			  
     for (int tt=0;tt<TP;++tt){  
     for (int ii=0; ii < I ; ++ii){  
     for (int kk=0;kk<KK+1;++kk){			  
     	for (int nn=0;nn<N;++nn){  
     //			nik( tt , ii + kk*I  ) += DAT2RESP(nn,ii)*(DAT2(nn,ii)==kk)*fqkyi(nn,tt)  ;  
     			nik( tt , ii + kk*I  ) += DAT2IND(nn,ii+kk*I) *fqkyi(nn,tt)  ;  
     				}  // end nn  
     		NIK(tt,ii) += nik(tt,ii+kk*I ) ; 				  
     			} // end kk  
     		}  // end ii  
     	}  // end tt  
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return List::create(_["fyiqk"] = fyiqk , _["f.qk.yi"]=fqkyi ,   
     	_["pi.k"] = PIK , _["n.ik"] = nik , _["N.ik"]=NIK );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}



