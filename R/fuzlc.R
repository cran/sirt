##
##
## fuzlc: fuzzy latent class model
##
## See paper of Deneoux (2013). Maximum likelihood estimation from
## uncertain data in the belief function framework.
## IEEE Transactions on Knowledge and Data Engineering.
##
## Notation in the paper:
##  i ... subject (person)
##  k ... class
##  j ... item
##  h ... category per item
##
## STEPS in estimation:
##
## Inits:
## - define number of latent classes K
## - set class probabilities equal to 1/K or sample them 
## - sample class-item-specific probabilities
##
## - arrange data with 'plausibilities' in array
##   dat[ subject , items , categories ]
##
## -------
## E-Step:
## (1) calculate posterior distribution t_{ik}
## (2) calculate gamma_{ik}^{jh} and beta_{ik}^{jh}
##     - only \sum_i beta_{ik}^{jh} is needed
## -------
## M-Step:
## (3) maximize integrated likelihood and compute class
##     probabilities pi_k and class-item-specific probabilities
##     alpha_k^{jh}
##
## - extend estimation to multiple group case  
##   -> class probabilities are group-specific
## 