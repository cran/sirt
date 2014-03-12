

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
SEXP isop_tests_C( SEXP dat_, SEXP dat_resp_, SEXP weights_, SEXP jackunits_, SEXP JJ_) ;
}

// definition

SEXP isop_tests_C( SEXP dat_, SEXP dat_resp_, SEXP weights_, SEXP jackunits_, SEXP JJ_ ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix dat(dat_);  
     Rcpp::NumericMatrix dat_resp(dat_resp_) ;   
     Rcpp::NumericVector weights(weights_) ;  
     Rcpp::NumericVector jackunits(jackunits_) ;  
     int JJ=as<int>(JJ_) ;  
       
     int N=dat.nrow();  
     int I=dat.ncol();  
       
     Rcpp::NumericMatrix Esi(I,JJ+1);  
     Rcpp::NumericMatrix Edi(I,JJ+1);  
     Rcpp::NumericMatrix W1i(I,JJ+1);  
     Rcpp::NumericVector W1test(JJ+1);  
       
       
     // int ii=0 ;  
     for (int ii=0;ii<I;ii++){  
     for (int nn=0;nn<N;nn++){  
     for (int mm=0;mm<N;mm++){  
     if (mm!=nn){  
     for (int jj=0;jj<I;jj++){  
         if ( jj != ii ){  
             if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){  
     //            Esi( ii , 0 ) += weights[nn] ;  
                 Esi.row(ii) = Esi.row(ii) + weights[nn] ;  
                 Esi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Esi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){  
                 Esi.row(ii) = Esi.row(ii) + weights[nn] ;  
                 Esi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Esi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) < dat(mm,ii) ) & ( dat(nn,jj) > dat(mm,jj) ) ){  
     //            Edi( ii , 0 ) += weights[nn] ;  
                 Edi.row(ii) = Edi.row(ii) + weights[nn] ;  
                 Edi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Edi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }  
             if ( ( dat(nn,ii) > dat(mm,ii) ) & ( dat(nn,jj) < dat(mm,jj) ) ){  
                 Edi.row(ii) = Edi.row(ii) + weights[nn] ;  
                 Edi(ii, jackunits[nn]+1 ) +=  - weights[nn] ;  
                 if ( jackunits[nn] != jackunits[mm] ){   
                        Edi( ii , jackunits[mm] + 1 ) += - weights[nn] ;                                   
                                            }  
                                     }                                  
                                 }  // end if jj != ii  
                          }  // end for jj  
                      }  // end if mm!= nn  
                  }   // end for mm  
              } // end for nn  
          } // end for ii  
       
     // compute W1i  
     for (int ii=0;ii<I;ii++){       
         //int ii = 0 ;       
         W1i.row(ii) = ( Esi.row(ii) - Edi.row(ii)  ) /  ( Esi.row(ii) + Edi.row(ii)  ) ;  
             }  
       
     // compute statistic for the whole test          
       
     double tmp1 ;  
     double tmp2 ;  
              
              
     for (int jj=0;jj<JJ+1;jj++){          
         tmp2 = 0 ;  
         for (int ii=0;ii<I;ii++){          
             tmp1 =  ( Esi(ii,jj) + Edi(ii,jj) )   ;          
             tmp2 += tmp1 ;  
             W1test[jj] += tmp1 * W1i(ii,jj) ;       
                 }  
         W1test[jj] = W1test[jj] / tmp2 ;          
             }  
               
               
     ///////////////////////////////////////  
     /// OUTPUT                  
       
     return List::create(  
             _["W1test"] = W1test ,   
     		_["W1i"] = W1i ,  
     		_["Esi"] = Esi ,  
             _["Edi"] = Edi  
     			) ;  
     
END_RCPP
}



