
##   
##   ##################################################################
##   # fslca ... Formann's structured latent class analysis
##   
##   Formann, A. K. (1985). Constrained latent class models: Theory and applications. 
##   British Journal of Mathematical and Statistical Psychology, 38(1), 87-111.
##   
##   Let p_{ic} the probability of item i in class c and
##   w_c the probability of class c.
##   
##   According to Formann (1985), it is assumed that
##   p_{ic} = plogis( x_ic ) and w_c = exp( z_c ) / sum_c  exp(z_c ).
##   
##   Then, the restricting equation of item response functions is
##   x_ic = sum_r q_icr lambda_r .
##   
##   The design matrix is defined as an array.
##   For estimating lambda_r it can be calculated at first,
##   which elements of q_icr are used for computation.
##   Consequently, tedious calculations can be avoided.
##   
##   The class probabilities are decomposed as
##   z_c = sum_r v_cr eta_r .
##   
##   In case of v_c1 = z_c^2 = theta_c^2 , 
##   a normal distribution can be approximated
##   with a dispersion parameter eta_1 .
##   
##   For a derivation of the first derivative, see Formann (1985).
##   
##   ##################################################################
