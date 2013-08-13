\name{mcmc.3pno.testlet}
\alias{mcmc.3pno.testlet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3PNO Testlet Model
}
\description{
This function estimates the 3PNO testlet model (Wang, Bradlow & Wainer, 2002)
by Markov Chain Monte Carlo methods (Glas, 2012).
}
\usage{
mcmc.3pno.testlet(dat, testlets = rep(NA, ncol(dat)), param = 1, 
   weights = NULL, est.slope = TRUE, est.guess = TRUE, guess.prior = NULL, 
   testlet.variance.prior = c(1, 0.2), burnin = 500, iter = 1000, 
   N.sampvalues = 1000, progress.iter = 50, save.theta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous item responses for \eqn{N} persons and \eqn{I} items
}
  \item{testlets}{
An integer or character vector which indicates the allocation of items to 
testlets. Same entries corresponds to same testlets. 
If an entry is \code{NA}, then this item does not belong to any testlet. 
}
  \item{param}{
The parametrization of the testlet model (see Details). Possible values are
\code{1}, \code{2} or \code{3}. 
}
  \item{weights}{
An optional vector with student sample weights
}
  \item{est.slope}{
Should item slopes be estimated? The default is \code{TRUE}.
}
  \item{est.guess}{
Should guessing parameters be estimated? The default is \code{TRUE}.
}
  \item{guess.prior}{
A vector of length two or a matrix with \eqn{I} items and two columns 
which defines the beta prior distribution of guessing
parameters. The default is a non-informative prior, i.e. the Beta(1,1)
distribution.
}
  \item{testlet.variance.prior}{
A vector of length two which defines the (joint) prior for testlet variances
assuming an inverse chi-squared distribution.
The first entry is the effective sample size of the prior while the second
entry defines the prior variance of the testlet. The default of \code{c(1,.2)}
means that the prior sample size is 1 and the prior testlet variance is .2.
}
  \item{burnin}{
Number of burnin iterations
}
  \item{iter}{
Number of iterations
}
  \item{N.sampvalues}{
Maximum number of sampled values to save
}
  \item{progress.iter}{
Display progress every \code{progress.iter}-th iteration. If no progress
display is wanted, then choose \code{progress.iter} larger than \code{iter}.
}
  \item{save.theta}{
Should theta values be saved?
}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{
The testlet response model for person \eqn{n} at item \eqn{i} 
is defined in parametrization 1 (\code{param=1})
\deqn{ P(X_{ni} = 1 ) = c_i + ( 1 - c_i ) 
		\Phi ( a_i \theta_n + \gamma_{n,t(i)} + b_i	) \quad , \quad
	\theta_n \sim N ( 0  ,1 ) , \gamma_{n,t(i)} \sim N( 0 , \sigma^2_t ) }

For parametrization 2 (\code{param=2}) the response model is defined as	
\deqn{ P(X_{ni} = 1 ) = c_i + ( 1 - c_i ) 
		\Phi ( a_i \theta_n + a_i \gamma_{n,t(i)} + b_i	) \quad , \quad
	\theta_n \sim N ( 0  ,1 ) , \gamma_{n,t(i)} \sim N( 0 , \sigma^2_t ) }

Parametrization 3 (\code{param=3}) is the bifactor model
\deqn{ P(X_{ni} = 1 ) = c_i + ( 1 - c_i ) 
		\Phi ( a_i \theta_n + a_{ti} \gamma_{n,t(i)} + b_i	) \quad , \quad
	\theta_n \sim N ( 0  ,1 ) , \gamma_{n,t(i)} \sim N( 0 , 1 )}
with item slopes \eqn{a_i} on the general factor \eqn{\theta_n} and
item slopes \eqn{a_{ti}} on the testlet factor \eqn{\gamma_{n,t(i)}}.
	
In case of \code{est.slope=FALSE}, all item slopes \eqn{a_i} are set to 1. Then
a variance \eqn{\sigma^2} of the \eqn{\theta_n} distribution is estimated.

In case of \code{est.guess=FALSE}, all guessing parameters \eqn{c_i} are
set to 0.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
A list of class \code{mcmc.sirt} with following entries:
  \item{mcmcobj}{Object of class \code{mcmc.list}}
  \item{summary.mcmcobj}{Summary of the \code{mcmcobj} object. In this
  summary the Rhat statistic and the mode estimate MAP is included.
  The variable \code{PercSEratio} indicates the proportion of the Monte Carlo
  standard error in relation to the total standard deviation of the
  posterior distribution.}
  \item{ic}{Information criteria (DIC)}
  \item{burnin}{Number of burnin iterations}
  \item{iter}{Total number of iterations}
  \item{theta.chain}{Sampled values of \eqn{\theta_p} parameters}
  \item{deviance.chain}{Sampled values of Deviance values}  
  \item{EAP.rel}{EAP reliability}
  \item{person}{Data frame with EAP person parameter estimates for 
  	  \eqn{\theta_p} and their corresponding posterior standard
  	deviations}  
  \item{dat}{Used data frame}
  \item{weights}{Used student weights}
  \item{\dots}{Further values}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\references{
Glas, C. A. W. (2012). \emph{Estimating and testing the extended testlet model.}
LSAC Research Report Series, RR 12-03.

Wang, X., Bradlow, E. T., & Wainer, H. (2002). A general Bayesian model 
for testlets: Theory and applications. 
\emph{Applied Psychological Measurement}, \bold{26}, 109-128.
	
Wainer, H., Bradlow, E. T., & Wang, X. (2007). 
\emph{Testlet response theory and its applications}. 
Cambridge: Cambridge University Press.
}

\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
S3 methods: \code{\link{summary.mcmc.sirt}}, \code{\link{plot.mcmc.sirt}}	
	
For estimating testlet models using the \pkg{lme4} package see
\code{\link{rasch.testlet.glmer}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data.read
I <- ncol(dat)

# set burnin and total number of iterations here (CHANGE THIS!)
burnin <- 200
iter <- 500

#***
# Model 1: 1PNO model
mod1 <- mcmc.3pno.testlet( dat ,  est.slope=FALSE , est.guess=FALSE , 
        burnin=burnin, iter=iter )
summary(mod1)
plot(mod1,ask=TRUE) # plot MCMC chains in coda style
plot(mod1,ask=TRUE , layout=2) # plot MCMC output in different layout

#***
# Model 2: 3PNO model with Beta(5,17) prior for guessing parameters
mod2 <- mcmc.3pno.testlet( dat ,  guess.prior=c(5,17) ,
            burnin=burnin, iter=iter )
summary(mod2)

#***
# Model 3: Rasch (1PNO) testlet model
testlets <- substring( colnames(dat) , 1 , 1 )
mod3 <- mcmc.3pno.testlet( dat ,  testlets=testlets ,  est.slope=FALSE , 
           est.guess=FALSE , burnin=burnin, iter=iter )
summary(mod3)

#***
# Model 4: 3PNO testlet model with (almost) fixed guessing parameters .25
mod4 <- mcmc.3pno.testlet( dat ,  guess.prior=1000*c(25,75) , testlets=testlets , 
              burnin=burnin, iter=iter )
summary(mod4)
plot(mod4, ask=TRUE, layout=2)

#***
# Model 5: 2PNO testlet model (param=2)
mod5 <- mcmc.3pno.testlet( dat ,  est.guess=FALSE , testlets=testlets , 
            param=2 , burnin=burnin, iter=iter )
summary(mod5)

#***
# Model 6: 2PNO testlet model (param=3: bifactor model)
mod6 <- mcmc.3pno.testlet( dat ,  est.guess=FALSE , testlets=testlets , 
            param=3 , burnin=burnin, iter=iter )
summary(mod6)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Testlet model}
\keyword{Testlets}
\keyword{Markov Chain Monte Carlo (MCMC)}
% __ONLY ONE__ keyword per line
