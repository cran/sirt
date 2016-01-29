

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

//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// interval_index_C
////////////////////////////////////////////////////////////////////////
//**********************************************************************

// declarations
extern "C" {
SEXP interval_index_C( SEXP matr, SEXP rn) ;
}

// definition

SEXP interval_index_C( SEXP matr, SEXP rn ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     Rcpp::NumericVector RN(rn) ;  
     	// rn random number for plausible value imputation  
       
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
       
     // create output vectors  
     NumericVector IND (NR) ;  
     IND.fill(0);  
       
     for (int nn=0;nn<NR;++nn){  
      	for (int cc=0 ; cc < NC ; ++cc ){  
     	    if ( MATR(nn,cc) > RN[nn] ){  
     	    	    IND(nn) = cc + 1 ;  
     	    	    break ;   
     	    	           }  
     		}  
     	}  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     return( wrap(IND) );  
     // return List::create(_["maxval"] = MAXVAL , _["maxind"]=MAXIND ) ;     
     
END_RCPP
}




//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// rowCumsums2_source
////////////////////////////////////////////////////////////////////////
//**********************************************************************


// declarations
extern "C" {
SEXP rowCumsums2_source( SEXP matr) ;
}

// definition

//# The C code was posted by Romain Francois at
//# http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-October/001198.html

SEXP rowCumsums2_source( SEXP matr ){
BEGIN_RCPP
     NumericMatrix input( matr ) ;  
          NumericMatrix output  = clone<NumericMatrix>( input ) ;  
       
          int nr = input.nrow(), nc = input.ncol() ;  
          NumericVector tmp( nr );  
          for( int i=0; i<nc; i++){  
              tmp = tmp + input.column(i) ;  
              NumericMatrix::Column target( output, i ) ;  
              std::copy( tmp.begin(), tmp.end(), target.begin() ) ;  
          }  
          return output ;
END_RCPP
}




//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// rowKSmallest_C
////////////////////////////////////////////////////////////////////////
//**********************************************************************



// declarations
extern "C" {
SEXP rowKSmallest_C( SEXP matr, SEXP K, SEXP indexmatr, SEXP rnmatr) ;
}

// definition

SEXP rowKSmallest_C( SEXP matr, SEXP K, SEXP indexmatr, SEXP rnmatr ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     Rcpp::IntegerVector KK1(K) ;   
     Rcpp::NumericMatrix INDEXMATR(indexmatr);  
     Rcpp::NumericMatrix RNMATR(rnmatr);  
       
     // define row and column numbers  
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
     int KK = KK1[0] ;   
       
     // copy original matrix ;  
     Rcpp::NumericMatrix MATRK = Rcpp::clone(MATR);  
       
     // create output vectors  
     NumericMatrix SMALLVAL (NR,KK) ;  
     NumericMatrix SMALLIND (NR,KK) ;  
     // SMALLIND.fill(1);  
       
     // int kk=0 ;  
     double tmp1 ;  // temp value  
       
     for (int kk=0;kk<KK;++kk){ // begin for kk  [donors]
     for (int nn=0;nn<NR;++nn){ // begin for nn  [cases]
          SMALLIND(nn,kk) = kk+1 ;	  
          SMALLVAL(nn,kk) = MATRK( nn , kk ) ;  
     	for (int cc=kk+1 ; cc < NC ; ++cc ){  
     			// begin for cc [matrix columns]  
     	    if ( ( MATRK(nn,cc) < SMALLVAL(nn,kk) ) ||  
     	         (  ( MATRK(nn,cc) == SMALLVAL(nn,kk) ) && ( RNMATR(nn,cc)==1 ) )   
     	    		){ // begin comparison  
     	    	    tmp1 = SMALLVAL(nn,kk) ;  
     	    	    SMALLVAL(nn,kk) = MATRK(nn,cc) ;  
     	    	    MATRK(nn,kk) = MATRK( nn , cc );  
     	    	    MATRK(nn,cc) = tmp1 ;  
     	    	    tmp1 = SMALLIND(nn,kk) ;  
     	    	    SMALLIND(nn,kk) = INDEXMATR(nn,cc) ;  
     	    	    INDEXMATR(nn,kk) = INDEXMATR(nn,cc) ;  
     	    	    INDEXMATR(nn,cc)=tmp1 ;  
     	    	           } // end comparison  
     		}  // end for cc  
     	}   // end for nn  
       }  // end for kk  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     // return( wrap(prob) );  
     // return List::create(_["smallval"]=SMALLVAL , _["smallind"]=SMALLIND ,  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
      return List::create(_["smallval"]=SMALLVAL , _["smallind"]=SMALLIND ) ;     
     
END_RCPP
}


//**********************************************************************
////////////////////////////////////////////////////////////////////////
/// rowMaxsCPP_source
////////////////////////////////////////////////////////////////////////
//**********************************************************************
// declarations
extern "C" {
SEXP rowMaxsCPP_source( SEXP matr) ;
}

// definition

SEXP rowMaxsCPP_source( SEXP matr ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix MATR(matr);  
     int NR=MATR.nrow();  
     int NC=MATR.ncol();  
       
     // create output vectors  
     NumericVector MAXVAL (NR) ;  
     NumericVector MAXIND (NR) ;  
     MAXIND.fill(1);  
       
     for (int nn=0;nn<NR;++nn){  
          MAXVAL[nn] = MATR( nn , 0 ) ;  
     	for (int cc=1 ; cc < NC ; ++cc ){  
     	    if ( MATR(nn,cc) > MAXVAL[nn] ){  
     	    	    MAXVAL[nn] = MATR(nn,cc) ;  
     	    	    MAXIND[nn] = cc + 1 ;  
     	    	           }  
     		}  
     	}  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     // return( wrap(prob) );  
     return List::create(_["maxval"] = MAXVAL , _["maxind"]=MAXIND ) ;     
     
END_RCPP
}


//*******************************************
// auxiliary function: calculation of probabilities in Copula model

// declarations
extern "C" {
SEXP rowmins2_bundle_C( SEXP m1_, SEXP v1_) ;
}

// definition

SEXP rowmins2_bundle_C( SEXP m1_, SEXP v1_ ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix m1(m1_);  
     Rcpp::NumericVector v1(v1_) ;  
       
       
     //	L1 <- length(v1)  
     int L1= v1.size();  
     int N1=m1.nrow() ;  
       
     //	m1min <- matrix( 0 , nrow=nrow(m1) , ncol= L1 )  
     Rcpp::NumericMatrix m1min(N1, L1) ;  
     Rcpp::NumericMatrix v1min(L1, 2) ;  
       
     //	v1min <- c(1 , cumsum(v1)[ - L1 ]+1 )  
     //	v1max <- cumsum(v1)  
     double t1=0 ;  
     v1min(0,0)=1 ;  
     for (int ll=0;ll<L1;ll++){  
     	t1 += v1[ll] ;  
     	v1min(ll,1) = t1 ;  
     	if (ll < L1 - 1 ){  
     		v1min(ll+1,0) = t1 + 1 ;  
     			}  
     	}  
       
       
     //	m1min[ , which(v1==1)] <- m1[ , v1min[ v1 == 1 ] ]  
     //	for (ll in (1:L1)[ v1 > 1] ){  
     //				m1min[,ll] <- rowMins2( m1[ , v1min[ll]:v1max[ll] ] )  
     //						}  
     //	m1min  
       
     int cc1=0;  
     int cc2=0;  
     double m1_vv = 0;  
       
     for (int ll = 0 ; ll < L1 ; ll++ ){  
     cc1 = v1min(ll,0) - 1 ;  
       
     for (int nn = 0 ; nn < N1 ; nn++){ // begin nn   
     	m1_vv = m1(nn , cc1 ) ;  
     	cc2 = v1min(ll,1) - 1 ;  
     	if (cc1 < cc2){  
     		for (int cc=cc1+1;cc<cc2+1;cc++){	  
     			if ( m1(nn,cc) < m1_vv ){  
     				m1_vv = m1(nn,cc) ;  
     						}  
     				} // end for cc  
     			} // end cc1 < cc2  
     	  
     	m1min(nn,ll) = m1_vv ;  
     	}  // end nn  
     	  
     } // end ll  
               
     return( wrap( m1min ) ) ;  
       
     ///////////////////////////////////////  
     /// OUTPUT                  
       
     // return List::create(  
     // _["v1min"] = v1min , _["m1min"] = m1min  
     //	) ;  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
     // return List::create(_["yM"]=YM , _["wM"]=WM ) ;     
     
END_RCPP
}


///************************************
// auxiliary function Copula model

// declarations
extern "C" {
SEXP calc_copula_itemcluster_C( SEXP D_, SEXP res_) ;
}

// definition

SEXP calc_copula_itemcluster_C( SEXP D_, SEXP res_ ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericVector DD(D_) ;  
     Rcpp::NumericMatrix res(res_) ;  
       
     int RR = res.nrow() ;  
     Rcpp::NumericMatrix matr(RR,RR) ;  
       
     double D = DD[0] ;  
     double t1=0;  
     double t2=0;  
       
     Rcpp::NumericMatrix res_rr(RR,D) ;  
     Rcpp::NumericVector rowsums_res(RR) ;  
     Rcpp::NumericVector g1_rr(RR) ;  
       
     // pattern in res  
       
     //    0 0 0 1   numerical value 2^0  
     //    0 0 1 0   numerical value 2^1  
     //    0 1 0 0   numerical value 2^2  
     //    1 0 0 0   numerical value 2^3   
     Rcpp::NumericVector res_patt(RR) ;  
       
       
     // for (int rr1=0;rr1<RR;rr1++){  
     // res_patt[rr1] = 0 ;  
     // for (int cc=0;cc<D;cc++){  
     //	res_patt[rr1] = res_patt[rr1] + res(rr1 , cc) * pow( 2 , D - cc - 1) ;   
     //	}  
     // }  
     	  
     //        g1.rr <- ( (-1)^rowSums( res ))  
       
       
       
     for (int ii=0;ii<RR;ii++){  
        t1 = 0 ;  
        for ( int cc=0;cc<D ; cc++){  
        	   t1 += res(ii,cc) ;  
        	   	}  
        rowsums_res[ii] = t1 ;  
        g1_rr[ii] = pow( -1.0 , rowsums_res[ii] ) ;  
        		}  
       
       
        		  
     for (int rr = 0 ; rr<RR;rr++){  
       
     // int rr = 2 ;  
       
     // res.rr <- outer( rep(1,nrow(res)) , res[rr,] ) - res  
     //        ind.rr <- which( apply( res.rr , 1 , min ) > -1 )  
     //        g1.rr <- g1.rr[ind.rr]  
     //        matr[ rr , a1.rr ]  <- g1.rr  
       
     for ( int zz=0;zz<RR ; zz++){  
     t2=3;  
     res_patt[zz]=0;  
     for ( int vv=0;vv<D ; vv++){
//                  resp_patt[nn] += pow(2.0,(double)(I-ii) ) - 1;       	     
           res_rr( zz , vv ) = res(rr,vv) - res( zz , vv ) ;  
           res_patt[zz] = res_patt[zz] + res_rr(zz , vv) * pow( 2.0 , (double)(D - vv - 1) ) ;         
           if ( res_rr(zz,vv) < t2 ){     
           		t2 = res_rr(zz,vv) ;   
           				}   
     			}  
           if ( t2 > -1 ){  
           	      matr( rr , res_patt[zz] ) = g1_rr[zz] ;      	       
           	      		}			  
     		}  
       
     //        g1.rr <- ( (-1)^rowSums( res ))  
       
      }    
       
               
     // return( wrap( matr ) ) ;  
       
     ///////////////////////////////////////  
     /// OUTPUT                  
       
      return List::create(  
      Rcpp::_["D"] = D , 
      Rcpp::_["res"] = res , 
      Rcpp::_["matr"] = matr , 
      Rcpp::_["res_rr"] = res_rr , 
      Rcpp::_["rowsums_res"] = rowsums_res ,  
      Rcpp::_["g1_rr"] = g1_rr , 
      Rcpp::_["res_patt"] = res_patt	
              ) ;  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
     // return List::create(_["yM"]=YM , _["wM"]=WM ) ;     
     
END_RCPP
}
//********************************************************************
//********************************************************************
//********************************************************************


//********************************************************************
//********************************************************************
//********************************************************************
// Extract missing data patterns
//********************************************************************
//********************************************************************
//********************************************************************
// declarations
extern "C" {
SEXP md_pattern_csource( SEXP dat_) ;
}

// definition

SEXP md_pattern_csource( SEXP dat_ ){
BEGIN_RCPP
         
     Rcpp::NumericMatrix dat(dat_);                        
     int I = dat.ncol();   
     int N = dat.nrow();  
     Rcpp::NumericMatrix dat_ind1(N,I) ;  
     Rcpp::NumericMatrix dat_ind0(N,I) ;  
     Rcpp::NumericVector resp_patt(N) ;  
     Rcpp::NumericVector freq1(I) ;  
     Rcpp::NumericVector freq0(I) ;  
            
     // define identifier of missingness pattern  
     for (int nn=0;nn<N; nn++){  
         for (int ii=0;ii<I ; ii++){  
             if (dat(nn,ii)==0){  
                 resp_patt[nn] += pow(2.0,(double)(I-ii) ) - 1;  
                 dat_ind0( freq0[ii] , ii) = nn+1 ;  
                 freq0[ii] ++ ;  
                     }  else {  
                 dat_ind1( freq1[ii] , ii) = nn+1 ;                      
                 freq1[ii] ++ ;  
                     }                      
         }  
     }  
       
     // unique response pattern  
     Rcpp::NumericVector unique_resp_patt = Rcpp::unique( resp_patt ) ;   
     // sort unique response pattern  
     std::sort(unique_resp_patt.begin(), unique_resp_patt.end());  
       
     // number of response patterns  
     int NP=unique_resp_patt.size();  
       
     Rcpp::NumericVector unique_resp_patt_freq(NP);  
     Rcpp::NumericVector unique_resp_patt_firstobs(NP);  
       
     for (int pp=0;pp<NP;pp++){  // begin pp  
         for (int nn=0;nn<N;nn++){    // begin nn  
             if ( resp_patt[nn] == unique_resp_patt[pp] ){ // begin if 1  
                  unique_resp_patt_freq[pp] ++ ;  
                   if( unique_resp_patt_firstobs[pp]==0){ // begin if 2  
                         unique_resp_patt_firstobs[pp]=nn+1 ;  
                   }  // end if 2  
                         } //end if 1  
                     }  
                 }  
       
     ////////////////////////////////////  
     // OUTPUT:  
     return Rcpp::List::create(  
         Rcpp::_["dat"] = dat ,  
         Rcpp::_["dat.resp1"] = dat_ind1 ,  
         Rcpp::_["dat.resp0"] = dat_ind0 ,  
         Rcpp::_["resp_patt"] = resp_patt ,  
         Rcpp::_["unique_resp_patt"] = unique_resp_patt ,  
         Rcpp::_["unique_resp_patt_freq"]=unique_resp_patt_freq ,  
         Rcpp::_["unique_resp_patt_firstobs"]=unique_resp_patt_firstobs ,              
         Rcpp::_["freq1"] = freq1 ,  
         Rcpp::_["freq0"] = freq0  
            		) ;  
     
END_RCPP
}
//********************************************************************
//********************************************************************
//********************************************************************



//********************************************************************
//********************************************************************
//********************************************************************
// Monotone rowwise regression in a matrix
//********************************************************************
//********************************************************************
//********************************************************************
// declarations
extern "C" {
SEXP monoreg_rowwise_Cpp( SEXP yM, SEXP wM) ;
}

// definition

SEXP monoreg_rowwise_Cpp( SEXP yM, SEXP wM ){
BEGIN_RCPP
  
     //********************************************  
     // Adapted code from monoreg function  
     // contained in the fdrtool package  
     // original author: Korbinian Strimmer  
     //********************************************  
       
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix YM(yM);  
     Rcpp::NumericMatrix WM(wM) ;   
       
     // define row and column numbers  
     int NR=YM.nrow();  
     int NC=YM.ncol();  
     int nn=NC ;  
       
     // create output ghat matrix  
     NumericMatrix ghat (NR,NC) ;  
     NumericMatrix k (NR,NC) ;  
     NumericMatrix gew (NR,NC) ;  
              
     double neu ;  
     int zz ;   
     int c=0 ;  
     int j=0 ;  
       
     for (zz=0; zz < NR ; zz++){ // begin for zz  
     	c=0 ;  
     	j=0 ;	  
     	nn = NC ;  
     	k(zz,c) = 0;  
     	gew(zz,c) = WM(zz,0);  
     	ghat(zz,c) = YM(zz,0);  
     	  
     	//######  
     	for (j=1; j < nn; j++){   // begin for j ...		  
     	    c = c+1;  
     	    k(zz,c) = j;  
     	    gew(zz,c) = WM(zz,j);  
     	    ghat(zz,c) = YM(zz,j);  
     	    /* c is at least 1 as nn is > 1 */  
     	    while (ghat(zz,c-1) >= ghat(zz,c))  
     	    {  
     	      neu = gew(zz,c)+gew(zz,c-1);  
     	      ghat(zz,c-1) = ghat(zz,c-1)+(gew(zz,c)/neu)*(ghat(zz,c)-ghat(zz,c-1) );  
     	      gew(zz,c-1) = neu;  
     	      c = c-1;  
     	      if (c==0) break;  
     	    }  
     	  }     // end for j ...  
     	//##########    
     	while (nn >= 1){  
     	    for (j=k(zz,c); j < nn; j++){  
     	      ghat(zz,j) = ghat(zz,c);  
     	    		}  
     	    nn = k(zz,c);  
     	    c = c-1;  
     	  }	  	    
     } // end for zz  
         
     ///////////////////////////////////////  
     /// OUTPUT                  
     return( wrap(ghat) );  
               
END_RCPP
}
//********************************************************************
//********************************************************************
//********************************************************************




