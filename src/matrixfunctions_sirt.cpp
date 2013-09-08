

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







