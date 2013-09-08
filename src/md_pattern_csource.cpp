//  Code created: 2013-08-26 17:38:39


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
                 resp_patt[nn] += pow(2,(I-ii) ) - 1;  
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
         _["dat"] = dat ,  
         _["dat.resp1"] = dat_ind1 ,  
         _["dat.resp0"] = dat_ind0 ,  
         _["resp_patt"] = resp_patt ,  
         _["unique_resp_patt"] = unique_resp_patt ,  
         _["unique_resp_patt_freq"]=unique_resp_patt_freq ,  
         _["unique_resp_patt_firstobs"]=unique_resp_patt_firstobs ,              
         _["freq1"] = freq1 ,  
         _["freq0"] = freq0  
            		) ;  
     
END_RCPP
}



