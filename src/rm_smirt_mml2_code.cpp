
// includes from the plugin

// #include <Rcpp.h>
#include <RcppArmadillo.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// array multiplication
/// The following function is used in arraymult1 in file rm.hrm_alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************




// declarations
extern "C" {
SEXP rm_arraymult1( SEXP A, SEXP dimA, SEXP B, SEXP dimB) ;
}

// definition

SEXP rm_arraymult1( SEXP A, SEXP dimA, SEXP B, SEXP dimB ){
BEGIN_RCPP
// define input and arrays
      // array A
      Rcpp::NumericMatrix AA(A);
      Rcpp::IntegerVector dimAA(dimA);
      // array B
      Rcpp::NumericMatrix BB(B);
      Rcpp::IntegerVector dimBB(dimB);
      Rcpp::NumericMatrix CC( dimAA[0]*dimAA[1] , dimBB[2] ) ;
      
      int a1 = dimAA[0];
      int a2 = dimAA[1];
      int a3 = dimAA[2];
      int b3 = dimBB[2]; 
      
      // ii -> loop within a matrix
      //for (int ii=0 ; ii < a1 ; ++ii ){
      for (int zz=0 ; zz < a2 ; ++zz){
          for (int ii=0 ; ii < a1 ; ++ii ){        
              for (int hh=0 ; hh < b3 ; ++hh){ // loop over columns
                  for (int kk=0 ; kk < a3 ; ++kk){
                      CC(ii+zz*a1,hh) += AA(ii+zz*a1,kk)*BB(ii+a1*kk,hh)  ; // *BB(kk,hh) ;
                              }
                          }
                      }
                  }
      
      // output
      return( wrap(CC) );
END_RCPP
}


//**********************************************************************
////////////////////////////////////////////////////////////////////////
// calculation posterior distribution
// RM_CALCPOST used in R function rm_calclike in file rm.alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************



// declarations
extern "C" {
SEXP RM_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP K) ;
}

// definition

SEXP RM_CALCPOST( SEXP dat2, SEXP dat2resp, SEXP probs, SEXP K ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
     Rcpp::NumericVector KK(K);  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
     int KKK = KK[0] + 1 ;   
       
     //*****  
     // calculate individual likelihood  
     NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( KKK*ii + DAT2(nn,ii) , tt ) ;  
                         }  
                     }  
                 }  
             }  
       
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
            Rcpp::_["fyiqk"] = fyiqk 
                );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}




//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// calculation of probabilities
// probraterfct1 in file rm.hrm_alg
////////////////////////////////////////////////////////////////////////
//**********************************************************************


// declarations
extern "C" {
SEXP rm_probraterfct1( SEXP crater, SEXP drater, SEXP dimA, SEXP B, SEXP dimB) ;
}

// definition

SEXP rm_probraterfct1( SEXP crater, SEXP drater, SEXP dimA, SEXP B, SEXP dimB ){
BEGIN_RCPP

      //////////////////////////////////////////
      //////// I N P U T ///////////////////////
      // c.rater
      Rcpp::NumericMatrix CRA(crater);
      // d.rater
      Rcpp::NumericVector DRA(drater) ;
      int K = CRA.ncol();
      int I = CRA.nrow();
      
      // define input and arrays
      // dimA (I,K+1,K+1)
      Rcpp::IntegerVector dimAA(dimA);
      // array B (I,K+1,TP)
      Rcpp::NumericMatrix BB(B);
      Rcpp::IntegerVector dimBB(dimB);
      Rcpp::NumericMatrix CC( dimAA[0]*dimAA[1] , dimBB[2] ) ;
      
      ////////////////////////////////////////////
      
      //*********************
      // create h1 matrix
      NumericMatrix h1(I*(K+1),K) ;
      int nrh1 = h1.nrow();
      
      for (int kk=0;kk<(K+1);++kk){
          for (int ii=0;ii<I;++ii){
              for (int jj=0;jj<K;++jj){
                  h1(ii+I*kk,jj) = CRA( ii , jj )  - DRA[ii]*kk ;
                                     }
                                 }
                          }
                          
      //**********************************
      // compute logistic distribution
      // for (int kk=0; kk<K; ++kk){
      for (int kk=0;kk<K;++kk){
          NumericVector inpv =h1(_,kk) ;
          NumericVector res0=plogis(inpv) ;
          for (int ii=0;ii<nrh1;++ii){
              h1(ii,kk) = res0[ii] ;
                               }
                      }
      
      //***********************************
      // compute matrix with rater probabilities  
      
      // P(X=0) ;                
      NumericMatrix PRA(I*(K+1),K+1) ;         
      for (int jj=0;jj<(K+1);++jj){                
          for (int ii=0;ii<I;++ii){
              PRA(ii,jj) = h1(ii+I*jj, 0 ) ;
                      }    
                  }
      // other categories
      for (int cc=1;cc<K;++cc){            
      for (int jj=0;jj<(K+1);++jj){            
          for (int ii=0;ii<I;++ii){
              PRA(ii+I*cc,jj) = h1(ii+I*jj, cc ) - h1(ii+I*jj, cc-1 ) ;
                      }               
                  }
              }
      // last category
      int cc=K ;        
      for (int jj=0;jj<(K+1);++jj){            
          for (int ii=0;ii<I;++ii){
              PRA(ii+I*cc,jj) = 1 - h1(ii+I*jj, cc-1 ) ;
                      }               
                  }
      
      //**************************
      // multiplication 
      
      NumericMatrix AA=PRA ;             
                  
      int a1 = dimAA[0];
      int a2 = dimAA[1];
      int a3 = dimAA[2];
      int b3 = dimBB[2]; 
      
      // ii -> loop within a matrix
      //for (int ii=0 ; ii < a1 ; ++ii ){
      for (int zz=0 ; zz < a2 ; ++zz){
          for (int ii=0 ; ii < a1 ; ++ii ){        
              for (int hh=0 ; hh < b3 ; ++hh){ // loop over columns
                  for (int kk=0 ; kk < a3 ; ++kk){
                      CC(ii+zz*a1,hh) += AA(ii+zz*a1,kk)*BB(ii+a1*kk,hh)  ; // *BB(kk,hh) ;
                              }
                          }
                      }
                  }
                           
      ///////////////////////////////////////////////////////
      ///////////// O U T P U T   ///////////////////////////
      return Rcpp::List::create(
         Rcpp::_["PRA"] = PRA , 
         Rcpp::_["probtotal"] = CC
                    );            
                  
      //return( wrap(res0) );
END_RCPP
}




//********************************************************************
//********************************************************************
//********************************************************************
// rm_facets_calcprobs
//********************************************************************
//********************************************************************
//********************************************************************

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
     return Rcpp::List::create(
        Rcpp::_["fyiqk"] = fyiqk , 
        Rcpp::_["f.qk.yi"]=fqkyi ,   
     	Rcpp::_["pi.k"] = PIK , 
        Rcpp::_["n.ik"] = nik , 
        Rcpp::_["N.ik"]=NIK 
                 );   
     	  
       
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



//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// SMIRT_CALCPROB_PARTCOMP
/// smirt - partially noncompensatory model
////////////////////////////////////////////////////////////////////////
//**********************************************************************

// declarations
extern "C" {
SEXP SMIRT_CALCPROB_PARTCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ,
	SEXP mui ) ;
}

// definition

SEXP SMIRT_CALCPROB_PARTCOMP( SEXP a, SEXP b, SEXP Q, SEXP thetak, SEXP cc, SEXP dd ,
	SEXP mui ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix A(a);  
     Rcpp::NumericMatrix B(b);  
     Rcpp::NumericMatrix QQ(Q);  
     Rcpp::NumericMatrix THETA(thetak);  
     Rcpp::NumericVector CC(cc);  
     Rcpp::NumericVector DD(dd);  
     Rcpp::NumericVector MUI(mui) ;
       
     int I=A.nrow();  
     int D=A.ncol();  
     int TP=THETA.nrow();  
               
     // create matrix of probabilities  
     NumericMatrix prob (I,TP) ;  
     prob.fill(1);  
       
     NumericVector yy1 (1);  
     NumericVector yy2 (1); 
     NumericVector yy3 (1);
     double tmp1 ;
     
     
     for (int ii=0;ii<I;++ii){  
     for (int tt=0; tt<TP; ++tt){
     	    yy1[0] = 0 ; 
            yy2[0] = 1 ;    	    
     for (int dd=0;dd<D;++dd){  
         if ( QQ(ii,dd)>0 ){    
             tmp1 = A(ii,dd)*QQ(ii,dd)*THETA(tt,dd)-B(ii,dd) ;
             yy1[0] = yy1[0] + tmp1  ;
             yy2[0] = yy2[0] * ( 1 + exp( tmp1 ) ) ;
                 }       // end if Q(ii,dd)>0  
             yy3[0] = exp( yy1[0] ) ;                              
                 }          // end dd   
          prob(ii,tt) = yy3[0]/(  MUI[ii] * yy2[0] + (1-MUI[ii])*(1+yy3[0]) ) ;
             
     //***  
     // include guessing and slipping parameters  
      if ( ( CC[ii] > 0 ) || ( DD[ii] < 1 ) ){        
             prob(ii,tt) = CC[ii] + ( DD[ii]-CC[ii] )* prob(ii,tt) ;  
                 } // end if condition for guessing or slipping                   
     		} // end tt                  
         }       // end ii  
           
     ///////////////////////////////////////  
     /// OUTPUT      
     return(wrap(prob)) ;               
     // return( wrap(prob) );  
     
END_RCPP
}


/////////////////////////////////////////////////////////////
/// rasch.mml2 : calculate counts

// user includes

// declarations
extern "C" {
SEXP MML2_RASCHTYPE_COUNTS( SEXP dat2, SEXP dat2resp, SEXP dat1, SEXP fqkyi, 
	SEXP pik, SEXP fyiqk) ;
}

// definition

SEXP MML2_RASCHTYPE_COUNTS( SEXP dat2, SEXP dat2resp, SEXP dat1, SEXP fqkyi, 
	SEXP pik, SEXP fyiqk ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericVector DAT1(dat1) ;  
     Rcpp::NumericMatrix FQKYI(fqkyi) ;  
     Rcpp::NumericVector PIK (pik) ;   
     Rcpp::NumericMatrix FYIQK(fyiqk) ;  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=FQKYI.ncol();  
       
     //**********************************  
     // calculate total counts n.k  
     NumericVector NK (TP);  
     NK.fill(0);  
       
     for (int tt=0;tt<TP;++tt){      
         for (int nn=0;nn<N;++nn){  
     //        NK[tt] = NK[tt] + FQKYI(nn,tt)*DAT1[nn] ;   
               NK[tt] += FQKYI(nn,tt)*DAT1[nn] ;   
                         } // end nn  
                     }  // end tt  
       
     //**********************************  
     // calculate n.jk and r.jk  
     NumericMatrix NJK (I,TP) ;  
     NumericMatrix RJK (I,TP) ;  
     NJK.fill(0) ;                   
     RJK.fill(0) ;    
                       
     for (int ii=0;ii<I;++ii){  
     for (int tt=0;tt<TP;++tt){  
         for (int nn=0;nn<N;++nn){  
             if (DAT2RESP(nn,ii)>0){  
                 NJK(ii,tt) += DAT1[nn] * FQKYI(nn,tt) ;   
                 RJK(ii,tt) += DAT1[nn] * FQKYI(nn,tt) * DAT2(nn,ii) ;               
                             }       // end if ...  
             }           //     end nn  
         } // end tt  
     }    // end ii  
       
     //*****************  
     // calculate loglikelihood  
       
     double LL=0;  
       
     //        ll[gg] <- sum( dat1[group==gg,2] * log( rowSums( f.yi.qk[group==gg,] *   
     //					outer( rep(1,nrow(f.yi.qk[group==gg,])) , pi.k[,gg] ) ) ) )  
       
     for (int nn=0;nn<N;++nn){  
         double total = 0; 	  
         for (int tt=0;tt<TP;++tt){  
               total += FYIQK(nn,tt) * PIK[tt] ;  
                             } // end tt  
         LL += log( total )*DAT1[nn] ;  
                 }  // end nn  
                           
                           
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
         Rcpp::_["nk"] = NK , 
         Rcpp::_["njk"] = NJK ,   
         Rcpp::_["rjk"] = RJK , 
         Rcpp::_["ll"] = LL 
                 );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}




//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution



// declarations
extern "C" {
SEXP MML2_CALCPOST_V1( SEXP dat2, SEXP dat2resp, SEXP probs) ;
}

// definition

SEXP MML2_CALCPOST_V1( SEXP dat2, SEXP dat2resp, SEXP probs ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
       
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
       
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
       Rcpp::_["fyiqk"] = fyiqk 
                   );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}



//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution (V2)


// declarations
extern "C" {
SEXP MML2_CALCPOST_V2( SEXP dat2, SEXP dat2resp, SEXP probs) ;
}

// definition

SEXP MML2_CALCPOST_V2( SEXP dat2, SEXP dat2resp, SEXP probs ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
       
     //*****  
     // calculate individual likelihood  
     NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
//             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;
	if ( ( DAT2(nn,ii) < 1 ) & ( DAT2(nn,ii) > 0 ) ){
             fyiqk(nn,tt) = fyiqk(nn,tt) * pow( PROBS(2*ii+1,tt),DAT2(nn,ii) )*
                	pow( PROBS( 2*ii  , tt ) , 1 - DAT2(nn,ii) ) ;
                		} else {
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;                			
                		}
                         }  
                     }  
                 }  
             }  
       
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
        Rcpp::_["fyiqk"] = fyiqk 
           );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}




//////////////////////////////////////////////////////////////////
/// rasch.mml2 - calculation of posterior distribution 
/// fuzzy log likelihood in the belief function framework


// declarations
extern "C" {
SEXP MML2_CALCPOST_V3( SEXP dat2, SEXP dat2resp, SEXP probs) ;
}

// definition

SEXP MML2_CALCPOST_V3( SEXP dat2, SEXP dat2resp, SEXP probs ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix DAT2(dat2);  
     Rcpp::NumericMatrix DAT2RESP(dat2resp);  
     Rcpp::NumericMatrix PROBS(probs);  
       
     int N=DAT2.nrow();  
     int I=DAT2.ncol();  
     int TP=PROBS.ncol();  
       
     //*****  
     // calculate individual likelihood  
     NumericMatrix fyiqk (N,TP) ;  
     fyiqk.fill(1);  
     for (int ii=0;ii<I;++ii){      
     for (int nn=0;nn<N;++nn){  
         if ( DAT2RESP(nn,ii)>0){  
         for (int tt=0;tt<TP;++tt){  
//             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;
	if ( ( DAT2(nn,ii) < 1 ) & ( DAT2(nn,ii) > 0 ) ){
             fyiqk(nn,tt) = fyiqk(nn,tt) * 
              (PROBS(2*ii+1,tt) * DAT2(nn,ii) + (PROBS( 2*ii , tt ) *(1-DAT2(nn,ii))) );
                		} else {
             fyiqk(nn,tt) = fyiqk(nn,tt) * PROBS( 2*ii + DAT2(nn,ii) , tt ) ;                			
                		}
                         }  
                     }  
                 }  
             }  
       
     			  
     ///////////////////////////////////////////////////////  
     ///////////// O U T P U T   ///////////////////////////  
     return Rcpp::List::create(
         Rcpp::_["fyiqk"] = fyiqk 
                        );   
     	  
       
     /// print output on R console  
     //		Rcpp::Rcout << "hier:" << std::endl << nn << std::endl;					
END_RCPP
}






