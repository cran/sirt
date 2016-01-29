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


//********************************************************************
//********************************************************************
//********************************************************************
// unidimensionality test of Gooijer
//********************************************************************
//********************************************************************
//********************************************************************

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
         Rcpp::_["stat"] = stat ,    
         Rcpp::_["stat_perm"]=stat_perm  
         ) ;  
END_RCPP
}

//********************************************************************
//********************************************************************
//********************************************************************




//********************************************************************
//********************************************************************
//********************************************************************
// ISOP test
//********************************************************************
//********************************************************************
//********************************************************************


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
             Rcpp::_["W1test"] = W1test ,   
             Rcpp::_["W1i"] = W1i ,  
             Rcpp::_["Esi"] = Esi ,  
             Rcpp::_["Edi"] = Edi  
     			) ;  
     
END_RCPP
}



//********************************************************************
//********************************************************************
//********************************************************************


