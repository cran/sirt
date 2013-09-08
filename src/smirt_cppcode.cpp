
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;



//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPOST
/// smirt - calculation of posterior probabilities
////////////////////////////////////////////////////////////////////////
//**********************************************************************


// declarations
extern "C" {
SEXP SMIRT_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP dat2ind, SEXP pik, SEXP K) ;
}

// definition

SEXP SMIRT_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP dat2ind, SEXP pik, SEXP K ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::IntegerMatrix DAT2(dat2);  
     Rcpp::IntegerMatrix DAT2RESP(dat2resp);  
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

//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_COMP
/// smirt - calculation of probabilities in the compensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************

// declarations
extern "C" {
SEXP SMIRT_CALCPROB_COMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
}

// definition

SEXP SMIRT_CALCPROB_COMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ){
BEGIN_RCPP
  

    /////////////////////////////////////
    // INPUT
    Rcpp::NumericMatrix A(a);
    Rcpp::NumericVector B(b);
    Rcpp::NumericMatrix QQ(Q);
    Rcpp::NumericMatrix THETA(thetak);
    Rcpp::NumericVector CC(cc);
    Rcpp::NumericVector DD(dd);

    int I=A.nrow();
    int D=A.ncol();
    int TP=THETA.nrow();

    // create matrix of probabilities
    NumericMatrix prob (I,TP) ;
    prob.fill(1);

    NumericVector yy (TP);

    for (int ii=0;ii<I;++ii){
            for (int tt=0; tt<TP; ++tt){
            yy[tt] = - B[ii] ; 	
           for (int dd=0;dd<D;++dd){
                   yy[tt] += A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)  ;
                    } // end dd 
                        }   // end tt
            NumericVector p1=plogis(yy) ;
            for (int tt=0; tt<TP; ++tt){            
                prob(ii,tt)= p1[tt] ;
                  }  // end tt
    //***
    // include guessing and slipping parameters
     if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){ 
        for (int tt=0;tt<TP;++tt){    
            prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;
                    }   // end tt
                } // end if condition for guessing or slipping                
        }       // end ii
        
    ///////////////////////////////////////
    /// OUTPUT    
    return(wrap(prob)) ;             
    // return( wrap(prob) );
         
END_RCPP
}


//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_NONCOMP
/// smirt - calculations in the noncompensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************

// declarations
extern "C" {
SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd) ;
}

// definition

SEXP SMIRT_CALCPROB_NONCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix A(a);  
     Rcpp::NumericMatrix B(b);  
     Rcpp::NumericMatrix QQ(Q);  
     Rcpp::NumericMatrix THETA(thetak);  
     Rcpp::NumericVector CC(cc);  
     Rcpp::NumericVector DD(dd);  
       
     int I=A.nrow();  
     int D=A.ncol();  
     int TP=THETA.nrow();  
       
     // create matrix of probabilities  
     NumericMatrix prob (I,TP) ;  
     prob.fill(1);  
       
     NumericVector yy (TP);  
       
     for (int ii=0;ii<I;++ii){  
     for (int dd=0;dd<D;++dd){  
         if ( QQ(ii,dd)>0 ){  
             for (int tt=0; tt<TP; ++tt){  
                    yy[tt] =  A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)-B(ii,dd)   ;  
                         }   // end tt  
             NumericVector p1=plogis(yy) ;  
             for (int tt=0; tt<TP; ++tt){              
                 prob(ii,tt)=prob(ii,tt)*p1[tt] ;  
                   }  // end tt  
                 }       // end if Q(ii,dd)>0  
             }          // end dd   
     //***  
     // include guessing and slipping parameters  
      if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){   
         for (int tt=0;tt<TP;++tt){      
             prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;  
                     }   // end tt  
                 } // end if condition for guessing or slipping                  
         }       // end ii  
           
     ///////////////////////////////////////  
     /// OUTPUT      
     return(wrap(prob)) ;               
     // return( wrap(prob) );  
     
END_RCPP
}







