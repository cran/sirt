

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


//**************************************************************************
// mlnormal_proc_variance_shortcut_Z_restructure
//**************************************************************************


// declarations
extern "C" {
SEXP mlnormal_proc_variance_shortcut_Z_restructure( SEXP Z_list_, SEXP update_dim_, 
	SEXP start_orig_, SEXP end_orig_, SEXP dim_Z_index_, SEXP Z_index_, SEXP orig_id_, 
	SEXP dim_id_) ;
}

// definition
SEXP mlnormal_proc_variance_shortcut_Z_restructure( SEXP Z_list_, SEXP update_dim_, 
	SEXP start_orig_, SEXP end_orig_, SEXP dim_Z_index_, SEXP Z_index_, SEXP orig_id_, 
	SEXP dim_id_ ){
BEGIN_RCPP

     //@@ INPUT BEGIN
     Rcpp::List Z_list(Z_list_);  
     Rcpp::IntegerVector update_dim(update_dim_);  
     Rcpp::NumericVector start_orig(start_orig_);  
     Rcpp::NumericVector end_orig(end_orig_);  
     Rcpp::NumericVector dim_Z_index(dim_Z_index_);  
     Rcpp::NumericVector Z_index(Z_index_);  
     Rcpp::NumericVector orig_id(orig_id_);  
     Rcpp::NumericVector dim_id(dim_id_);  
     //@@ INPUT END              
     
     int G = dim_Z_index[0];  
     int NM = dim_Z_index[1];  
     int NP = dim_Z_index[2];  
       
     // copy of update_dim  
     Rcpp::NumericMatrix update_dim2(G,1);  
     update_dim2(_,0) = update_dim ;   
       
     int orig_gg = 0;  
     int orig_gg1 = 0;  
     double Z1 = 0.0;  
     double Z2 = 0.0;         
       
     //----- compare Z_index         
     for (int gg = 0; gg < G; gg++){       	  
     	orig_gg = orig_id[gg];  
     	if ( update_dim2(gg,0) == 0){  
     		orig_gg1 = orig_id[gg-1];						  
     	}       	  
     	for( int mm = 0; mm < NM ; mm++){  
     		if ( update_dim2(gg,0) == 0){  
     			for( int pp = 0; pp < NP ; pp++){  
     				Z1 = Z_index[ orig_gg + mm*G + pp*G*NM - 1] ;  
     				Z2 = Z_index[ orig_gg1 + mm*G + pp*G*NM - 1] ;				  
     				if ( Z1 != Z2 ){  
     					update_dim2(gg,0) = 1;					  
     				}  
     			}  // end parameter pp  
     		}  
     	}  // end matrices mm  
     }  // end groups gg  
              
     int dim_gg = 0;  
     double eps = 1E-15;  
     double val = 0.0;  
       
     Rcpp::NumericMatrix Z_gg_mm ;  
     Rcpp::NumericMatrix Z_gg1_mm ;  
     //----- compare Z_list          
     for (int gg=0; gg < G ; gg++){	  
     	orig_gg = orig_id[gg];  
     	if ( update_dim2(gg,0) == 0){  
     		orig_gg1 = orig_id[gg-1];						  
     	}		  
     	if ( update_dim2(gg,0) == 0 ){  
     		Rcpp::List Z_gg1 = Rcpp::as<Rcpp::List>( Z_list[ orig_gg1 - 1] ) ;  
     		Rcpp::List Z_gg = Rcpp::as<Rcpp::List>( Z_list[ orig_gg - 1] ) ;		  
     		dim_gg = dim_id[gg] ;	  
     		if ( update_dim2(gg,0) == 0 ){  
     			for (int mm=0; mm < NM; mm++){  
     				Rcpp::NumericMatrix Z_gg1_mm = Rcpp::as<Rcpp::NumericMatrix>( Z_gg1[mm] ) ;  
     				Rcpp::NumericMatrix Z_gg_mm = Rcpp::as<Rcpp::NumericMatrix>( Z_gg[mm] ) ;				  
     				for (int rr=0; rr < dim_gg; rr++){  
     					for (int cc=rr; cc < dim_gg; cc++){  
     						val = Z_gg1_mm(rr,cc) - Z_gg_mm(rr,cc) ;  
     						if ( val < 0){  
     							val = - val ;  
     						}  
     						if ( val > eps ){  
     							update_dim2(gg,0) = 1 ;  
     						}  
     					}  // end cc  
     				}  // end rr  
     			}   // end mm  
     		}  
     	}  
     }  // end gg	  
     	  
     //*************************************************      
     // OUTPUT                                 
     return Rcpp::List::create(    
         Rcpp::Named("NM") = NM   ,  
         Rcpp::Named("NP") = NP   ,      
         Rcpp::Named("Z_gg") = Z_gg_mm   ,  
         Rcpp::Named("update_dim") = update_dim2      
         ) ;    
                     	
END_RCPP
}

     // Rcout << "gg=" << gg << " | orig_gg=" << orig_gg <<  " | orig_gg1=" << orig_gg1 <<  std::endl;  
     
     // Rcout << "term compare mm,pp " << mm+1 << " "  << pp+1 << " | terms for gg and gg1 " << Z1  
     //  << " " << Z2 << std::endl;		

//**************************************************************************
// mlnormal_proc_variance_shortcut_XY_restructure
//**************************************************************************

// declarations
extern "C" {
SEXP mlnormal_proc_variance_shortcut_XY_restructure( SEXP freq_id_, SEXP y_, SEXP X_, SEXP G_) ;
}

// definition
SEXP mlnormal_proc_variance_shortcut_XY_restructure( SEXP freq_id_, SEXP y_, SEXP X_, SEXP G_ ){
BEGIN_RCPP
     
     //@@ INPUT BEGIN	
     Rcpp::NumericMatrix freq_id(freq_id_);  
     Rcpp::NumericVector y(y_);  
     Rcpp::NumericMatrix X(X_);  
     int G = Rcpp::as<int>(G_);  
     //@@ INPUT END  
     
     int N = X.nrow() ;  
     int V = X.ncol() ;  
       
     Rcpp::NumericMatrix X1(N,V);  
     Rcpp::NumericVector y1(N);  
             
     int hh=0;  
     int min_gg=0;  
     int max_gg=0;  
       
     for (int gg = 0 ; gg < G ; gg++){  
     	min_gg = freq_id(gg,2) - 1;  
     	max_gg = freq_id(gg,3) ;     
     	for (int mm = min_gg ; mm < max_gg ; mm++){  
     		y1[hh] = y[mm] ;  
     		for (int vv=0; vv < V; vv++){  
     			X1(hh,vv) = X(mm,vv);  
     		}  
     		hh ++  ;  
     	}  
     }	                
       
     //*************************************************      
     // OUTPUT                                 
     return Rcpp::List::create(    
         Rcpp::Named("N") = N ,  
         Rcpp::Named("V") = V ,  
         Rcpp::Named("X") = X1 ,  
         Rcpp::Named("y") = y1      
         ) ;  
END_RCPP
}


//**************************************************************************
// mlnormal_update_V_rcpp_helper
//**************************************************************************


// declarations
extern "C" {
SEXP mlnormal_update_V_rcpp_helper( SEXP Z_list_, SEXP Z_index_, SEXP dim_id_, 
	SEXP dim_Z_index_, SEXP startIndex_, SEXP endIndex_, SEXP N_, 
	SEXP max_dim_, SEXP do_compute_, SEXP theta_, SEXP use_ginverse_) ;
}

// definition

SEXP mlnormal_update_V_rcpp_helper( SEXP Z_list_, SEXP Z_index_, SEXP dim_id_, 
	SEXP dim_Z_index_, SEXP startIndex_, SEXP endIndex_, SEXP N_, 
	SEXP max_dim_, SEXP do_compute_, SEXP theta_, SEXP use_ginverse_ ){
BEGIN_RCPP
  
     //@@ INPUT BEGIN
     Rcpp::List Z_list(Z_list_);  
     Rcpp::NumericVector Z_index(Z_index_);  
     Rcpp::NumericVector dim_id(dim_id_);  
     Rcpp::NumericVector dim_Z_index(dim_Z_index_);  
     Rcpp::NumericVector startIndex(startIndex_);  
     Rcpp::NumericVector endIndex(endIndex_);  
     int N = Rcpp::as<int>(N_);  
     int max_dim = Rcpp::as<int>(max_dim_);  
     Rcpp::NumericVector do_compute(do_compute_);  
     Rcpp::NumericVector theta(theta_);  
     int use_ginverse = Rcpp::as<int>(use_ginverse_);  
     //@@ INPUT END        
       
     int G = dim_id.size();  
     int NM = dim_Z_index[1];  
     int NP = dim_Z_index[2];  
       
     Rcpp::NumericMatrix V( N , max_dim );  
     Rcpp::NumericMatrix V1( N , max_dim );  
       
     Rcpp::List V_list(G);  
     Rcpp::List V1_list(G);  
       
     arma::mat V_gg;  
     arma::mat V1_gg;  
     Rcpp::List Z_gg;  
     Rcpp::NumericMatrix Z_gg_mm ;  
     int dim_gg = 0;  
     double val=0; 
     double eps = 1E-15;
     double zval = 0;
       
     for (int gg = 0 ; gg < G; gg++){  
     	if ( do_compute[gg] == 1 ){  
     		dim_gg = dim_id[gg];  
     		Z_gg = Rcpp::as<Rcpp::List>( Z_list[gg] ) ;  
     		V_gg = arma::zeros<arma::mat>(dim_gg,dim_gg);  
     		double pow_gg_mm_pp = 0 ;  
     		for (int mm = 0; mm <NM; mm++){  
     			Z_gg_mm = Rcpp::as<Rcpp::NumericMatrix>( Z_gg[mm] ) ;  
     			val = 1 ;  
     			for (int pp = 0 ; pp < NP ; pp ++){   
     				pow_gg_mm_pp = Z_index[ gg + mm*G + pp*G*NM  ] ;  
     				if ( pow_gg_mm_pp > 0 ){  
     					val =  pow( theta[pp] , pow_gg_mm_pp ) * val ;  
     				}  
     			}	  
     			for (int rr=0;rr<dim_gg;rr++){  
     				for (int cc=0;cc<dim_gg;cc++){  
     					zval = Z_gg_mm(rr,cc) ;
     					if (rr <= cc ){
     						if ( ( zval > eps ) | ( zval < - eps ) ){
     							V_gg(rr,cc) = V_gg(rr,cc) + val* zval ;
     						}
     					} else {  
     						V_gg(rr,cc) = V_gg(cc,rr);  
     					}  
     				}  
     			}  
     		}  
     		// calculate inverse  
     		if (use_ginverse){  
     			V1_gg = arma::pinv( V_gg );  
     		} else {  
     			V1_gg = arma::inv( V_gg );  
     		}  
     	}  
     	// fill matrices V and V1  
     	for( int rr=0; rr < dim_gg ; rr++){  
     		for( int cc=0; cc < dim_gg ; cc++){  
     			V( startIndex[gg] - 1 + rr , cc ) = V_gg( rr , cc ) ;  
     			V1( startIndex[gg] - 1 + rr , cc ) = V1_gg( rr , cc ) ;  
     		}  
     	}    
     	// fill list  
     	V_list[gg] = V_gg ;  
     	V1_list[gg] = V1_gg ;	  
     }  
              
     //*************************************************      
     // OUTPUT                          
     return Rcpp::List::create(    
         Rcpp::Named("V") = V ,  
         Rcpp::Named("V1") = V1 ,  
         Rcpp::Named("V_list") = V_list ,  
         Rcpp::Named("V1_list") = V1_list  
         ) ;    
       
END_RCPP
}


//**************************************************************************
// mlnormal_update_beta_rcpp_helper
//**************************************************************************


// declarations
extern "C" {
SEXP mlnormal_update_beta_rcpp_helper( SEXP dim_id_, SEXP startIndex_, SEXP endIndex_, 
	    SEXP G_, SEXP X_, SEXP y_, SEXP V1_) ;
}

// definition

SEXP mlnormal_update_beta_rcpp_helper( SEXP dim_id_, SEXP startIndex_, SEXP endIndex_, 
	   SEXP G_, SEXP X_, SEXP y_, SEXP V1_ ){
BEGIN_RCPP
  
     //@@ INPUT BEGIN  
     Rcpp::NumericVector dim_id(dim_id_);  
     Rcpp::NumericVector startIndex(startIndex_);  
     Rcpp::NumericVector endIndex(endIndex_);  
     int G = Rcpp::as<int>(G_);  
     Rcpp::NumericMatrix X(X_);  
     Rcpp::NumericVector y(y_);  
     Rcpp::NumericMatrix V1(V1_);  
     //@@ INPUT END  
       
     int NB = X.ncol();  
     Rcpp::NumericMatrix XVX(NB,NB);  
     Rcpp::NumericMatrix XVY(NB,1);  
       
     int dim_gg=0;  
     int rr1 = 0;  
     int cc1 = 0;  
       
     for (int pp=0; pp < NB; pp++){  
     	for (int qq=pp; qq < NB; qq++){  
     		XVX(pp,qq) = 0 ;  
     		if (pp==qq){  
     			XVY(qq,0) = 0 ;  
     		}  
     		for (int gg=0; gg <G; gg++){  
     			dim_gg = dim_id[gg];  
     			for (int rr=0; rr < dim_gg; rr++){  
     				for (int cc=0; cc < dim_gg; cc++){  
     					rr1 = startIndex[gg] + rr - 1 ;  
     					cc1 = startIndex[gg] + cc - 1 ;		  
     					XVX(pp,qq) += X(rr1,pp) * V1(rr1,cc) * X(cc1,qq);  
     					if (pp==qq){  
     						XVY(pp,0) += X(rr1,pp) * V1(rr1,cc) * y[cc1];  
     					}  
     				}  // end cc	  
     			}   // end rr  
     			if (pp < qq){  
     				XVX(qq,pp) = XVX(pp,qq);  
     			}  // end if ( pp < qq )  
     		}  // end gg  
     	}   // end qq  
     }   // end pp  
                     
     //*************************************************      
     // OUTPUT                                 
     return Rcpp::List::create(    
         Rcpp::Named("XVX") = XVX ,  
         Rcpp::Named("XVY") = XVY   
         ) ;    
       
END_RCPP
}

       
       
     // Rcout << "gg=" << gg << " | orig_gg=" << orig_gg <<  " | orig_gg1=" << orig_gg1 <<  std::endl;  
       
     // Rcout << "term compare mm,pp " << mm+1 << " "  << pp+1 << " | terms for gg and gg1 " << Z1  
     //  << " " << Z2 << std::endl;				




