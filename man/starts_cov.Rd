\name{starts_cov}
\alias{starts_cov}
\alias{starts_sim1dim}
%- Also NEED an '\alias' for EACH other topic documented here.



\title{
Some Functions for the STARTS model
}


\description{
Functions for computing the covariance matrix and simulating data
from the univariate STARTS model (Kenny & Zautra, 1995, 2001).
}

\usage{
starts_cov(W, var_trait, var_ar, var_state, a )

starts_sim1dim(N, W, var_trait, var_ar, var_state, a )
}


%- maybe also 'usage' for other objects documented here.
\arguments{
\item{W}{
Number of measurement waves.
}
\item{var_trait}{
Variance of trait component.
}
\item{var_ar}{
Variance of autoregressive component.
}    
\item{var_state}{
Variance of state component.
}    
\item{N}{
Sample size of persons
}    
\item{a}{
Stability parameter
}      
}

%%\details{  
%% xxx    
%% }


\value{
Covariance matrix or a data frame of simulated values
}

%% include references
\references{
Kenny, D. A., & Zautra, A. (1995). The trait-state-error 
model for multiwave data. 
\emph{Journal of Consulting and Clinical Psychology}, 
\bold{63}, 52-59.

Kenny, D. A., & Zautra, A. (2001). Trait-state models for 
longitudinal data. In L. M. Collins & A. G. Sayer (Eds.), 
\emph{New methods for the analysis of change} (pp. 243-263).
Washington, DC, US: American Psychological Association.

Wu, P.-C. (2016). Longitudinal stability of the Beck 
Depression Inventory II: A latent traitstate-occasion model. 
\emph{Journal of Psychoeducational Assessment}, 
\bold{34}, 39-53.
}

\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%xxx
%%}


\examples{
\dontrun{    
#############################################################################
# EXAMPLE 1: STARTS model | Data from Wu (2016)
#############################################################################

## define list with input data
## S ... covariance matrix, M ... mean vector

# read covariance matrix of data in Wu (older cohort, positive affect)
S <- matrix( c( 12.745, 7.046, 6.906, 6.070, 5.047, 6.110
    7.046, 14.977, 8.334, 6.714, 6.91, 6.624
    6.906, 8.334, 13.323, 7.979, 8.418, 7.951
    6.070, 6.714, 7.979, 12.041, 7.874, 8.099
    5.047, 6.91, 8.418, 7.874, 13.838, 9.117
    6.110, 6.624, 7.951, 8.099, 9.117, 15.132 ) ,  
    nrow=6 , ncol=6 , byrow=TRUE )
W <- 6   # number of measurement waves
data <- list( "S" = S , "M" = rep(0,W) , "n" = 660 , "W" = W  )

#*** likelihood function for the STARTS model
model <- function( pars , data ){ 
    # mean vector
    mu <- data$M
    # covariance matrix
    W <- data$W 
    var_trait <- pars["vt"]
    var_ar <- pars["va"]
    var_state <- pars["vs"]
    a <- pars["b"]
    Sigma <- sirt::starts_cov( W=W , var_trait=vt , var_ar=va , 
                var_state=vs , a=b )
    # evaluate log-likelihood
    ll <- sirt::loglike_mvnorm( S = data$S , Sigma = Sigma , M = data$M , mu = mu ,
                n = data$n , lambda = 1E-3)
    return(ll)
}
#** Note: 
#   (1) The function starts_cov calculates the model implied covariance matrix
#       for the STARTS model.
#   (2) The function loglike_mvnorm evaluates the loglikelihood for a multivariate
#       normal distribution given sample and population means M and mu, and sample 
#       and population covariance matrix S and Sigma.

#*** starting values for parameters
pars <- c( .33 , .33 , .33 , .75) 
names(pars) <- c("vt","va","vs","b")
#*** bounds for acceptance rates
acceptance_bounds <- c( .45 , .55 )
#*** starting values for proposal standard deviations
proposal_sd <- c( .1 , .1 , .1 , .1 )
names(proposal_sd) <- names(pars)
#*** lower and upper bounds for parameter estimates
pars_lower <- c( .001 , .001 , .001 , .001 )
pars_upper <- c( 10 , 10 , 10 , .999 ) 
#*** define prior distributions | use prior sample size of 3
prior_model <- "
    vt ~ dinvgamma(NA , 3 / 2 , 3 * .33 / 2 )
    va ~ dinvgamma(NA , 3 , 3 * .33 / 2)    
    vs ~ dinvgamma(NA , 3 , 3 * .33 / 2)        
    b ~ dbeta(NA, 4 , 4 )
        "
#*** estimate model with 'sirt::amh' function		
mod <- sirt::amh( data=data, nobs=data$n, pars=pars, model=model,  
            prior=prior_model, proposal_sd=proposal_sd, n.iter=n.iter, 
            n.burnin=n.burnin, pars_lower=pars_lower, pars_upper=pars_upper)
#*** model summary
summary(mod)

#---------------------------
# fitting the STARTS model in lavaan

library(lavaan)

## define lavaan model
lavmodel <- "
     #*** stable trait
     T =~ 1*W1 + 1*W2 + 1*W3 + 1*W4 + 1*W5 + 1*W6
     T ~~ vt * T
     W1 ~~ 0*W1
     W2 ~~ 0*W2
     W3 ~~ 0*W3
     W4 ~~ 0*W4
     W5 ~~ 0*W5     
     W6 ~~ 0*W6     
     #*** autoregressive trait
     AR1 =~ 1*W1
     AR2 =~ 1*W2     
     AR3 =~ 1*W3     
     AR4 =~ 1*W4          
     AR5 =~ 1*W5               
     AR6 =~ 1*W6              
     #*** state component      
     S1 =~ 1*W1
     S2 =~ 1*W2     
     S3 =~ 1*W3     
     S4 =~ 1*W4        
     S5 =~ 1*W5          
     S6 =~ 1*W6               
     S1 ~~ vs * S1
     S2 ~~ vs * S2
     S3 ~~ vs * S3
     S4 ~~ vs * S4 
     S5 ~~ vs * S5      
     S6 ~~ vs * S6      
     AR2 ~ b * AR1
     AR3 ~ b * AR2     
     AR4 ~ b * AR3          
     AR5 ~ b * AR4               
     AR6 ~ b * AR5                    
     AR1 ~~ va * AR1              
     AR2 ~~ v1 * AR2                   
     AR3 ~~ v1 * AR3     
     AR4 ~~ v1 * AR4   
     AR5 ~~ v1 * AR5       
     AR6 ~~ v1 * AR6 
     #*** nonlinear constraint
     v1 == va * ( 1 - b^2 )  
     #*** force variances to be positive
     vt > 0.001
     va > 0.001
     vs > 0.001
     #*** variance proportions
     var_total := vt + vs + va
     propt := vt / var_total
     propa := va / var_total
     props := vs / var_total
	"	
# estimate lavaan model
mod <- lavaan::lavaan(model = lavmodel, sample.cov = S, sample.nobs = 660)
# summary and fit measures
summary(mod)
lavaan::fitMeasures(mod)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
