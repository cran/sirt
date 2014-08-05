\name{rm.facets}
\alias{rm.facets}
\alias{summary.rm.facets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rater Facets Models with Item/Rater Intercepts and Slopes
}
\description{
This function estimates the unidimensional rater facets model (Lincare, 1994)
and an extension to slopes (see Details). The estimation
is conducted by an EM algorithm employing marginal
maximum likelihood.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
rm.facets(dat, pid=NULL, rater=NULL, Qmatrix=NULL, theta.k=seq(-9, 9, len=30), 
    est.b.rater=TRUE, est.a.item=FALSE, est.a.rater=FALSE, max.b.increment=1, 
    numdiff.parm=0.00001, maxdevchange=0.1, globconv=0.001, maxiter=1000, 
    msteps=4, mstepconv=0.001)

\method{summary}{rm.facets}(object,...)	    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Original data frame. Ratings on variables must be in rows,
i.e. every row corresponds to a person-rater combination.
}
  \item{pid}{
Person identifier.
}
  \item{rater}{
Rater identifier
}
  \item{Qmatrix}{
An optional Q-matrix. If this matrix is not provided,
then by default the ordinary scoring of categories
(from 0 to the maximum score of \eqn{K}) is used.
}
  \item{theta.k}{
A grid of theta values for the ability distribution.
}
  \item{est.b.rater}{
Should the rater severities \eqn{b_r} be estimated?
}
  \item{est.a.item}{
Should the item slopes \eqn{a_i} be estimated?
}
  \item{est.a.rater}{
Should the rater slopes \eqn{a_r} be estimated?
}
  \item{max.b.increment}{
Maximum increment of item parameters during estimation
}
  \item{numdiff.parm}{
Numerical differentiation step width
}
  \item{maxdevchange}{
Maximum relative deviance change as a convergence criterion
}
  \item{globconv}{
Maximum parameter change
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{msteps}{
Maximum number of iterations during an M step
}
  \item{mstepconv}{
Convergence criterion in an M step
}
  \item{object}{
Object of class \code{rm.facets}
}
  \item{\dots}{
Further arguments to be passed
}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
This function models ratings \eqn{X_{pri}}
for person \eqn{p}, rater \eqn{r} and item \eqn{i}
and category \eqn{k}
\deqn{P( X_{pri} = k | \theta_p ) \propto
    exp( a_i a_r q_{ik} \theta_p - q_{ik} b_r -   \tau_{ik} ) \quad , 
	\quad \theta_p \sim N( 0 , \sigma^2 )}{%
    P( X_{pri} = k | \theta_p ) \propto
    exp( a_i a_r q_{ik} \theta_p - q_{ik} b_r - \tau_{ik} )    ,   
    		\theta_p ~ N( 0 , \sigma^2 )}
By default, the scores in the \eqn{Q} matrix are
\eqn{q_{ik}=k}. Item slopes \eqn{a_i} and rater slopes
\eqn{a_r} are standardized such that their product equals
one, i.e. \eqn{ \prod_i a_i = \prod_r a_r = 1}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria and number of parameters}
\item{item}{Data frame with item parameters}
\item{rater}{Data frame with rater parameters}
\item{person}{Data frame with person parameters: EAP and corresponding 
	standard errors}
\item{EAP.rel}{EAP reliability}
\item{sigma}{Standard deviation of the trait}
\item{tau.item}{Item parameters \eqn{\tau_{ik}}}
\item{se.tau.item}{Standard error of item parameters \eqn{\tau_{ik}}}
\item{a.item}{Item slopes \eqn{a_i}}
\item{se.a.item}{Standard error of item slopes \eqn{a_i}}
\item{b.rater}{Rater severity parameter \eqn{b_r}}
\item{se.b.rater}{Standard error of rater severity parameter \eqn{b_r}}
\item{a.rater}{Rater slope parameter \eqn{a_r}}
\item{se.a.rater}{Standard error of rater slope parameter \eqn{a_r}}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior distribution}
\item{probs}{Item probabilities at grid \code{theta.k}}
\item{n.ik}{Expected counts}
\item{maxK}{Maximum number of categories}
\item{procdata}{Processed data}
\item{iter}{Number of iterations}
\item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{
Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}. 
Chicago: MESA Press.
}


\author{
Alexander Robitzsch
}

\note{
If the trait standard deviation \code{sigma} strongly
differs from 1, then a user should investigate the sensitivity
of results using different theta integration points \code{theta.k}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{TAM} package for the estimation
of more complicated facet models.

See \code{\link{rm.hrm}} for estimating a hierarchical rater model.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Partial Credit Model and Generalized partial credit model 
#                   5 items - 1 rater
#############################################################################
data(data.ratings1)
dat <- data.ratings1

# select rater db01
dat <- dat[ paste(dat$rater) == "db01" , ]

# Model 1: Partial Credit Model
mod1 <- rm.facets( dat[ , paste0( "k",1:5) ] , pid=dat$idstud , maxiter=15)

# Model 2: Generalized Partial Credit Model
mod2 <- rm.facets( dat[ , paste0( "k",1:5) ] ,  pid=dat$idstud  , 
           est.a.item=TRUE , maxiter=15)

summary(mod1)
summary(mod2)

\dontrun{
#############################################################################
# EXAMPLE 2: Facets Model: 5 items, 7 raters
#############################################################################

data(data.ratings1)
dat <- data.ratings1

# Model 1: Partial Credit Model: no rater effects
mod1 <- rm.facets( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
             pid=dat$idstud  , est.b.rater=FALSE  , maxiter=15)

# Model 2: Partial Credit Model: intercept rater effects
mod2 <- rm.facets( dat[ , paste0( "k",1:5) ] , rater=dat$rater  , 
             pid=dat$idstud  , maxiter=15)

# Model 2a: compare results with TAM package
#   Results should be similar to Model 2
library(TAM)
mod2a <- TAM::tam.mml.mfr( resp= dat[ , paste0( "k",1:5) ] , 
             facets= dat[ , "rater" , drop=FALSE] ,
             pid= dat$pid , formulaA = ~ item*step + rater )

# Model 3: estimated rater slopes
mod3 <- rm.facets( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            est.a.rater=TRUE  , maxiter=15)

# Model 4: estimated item slopes
mod4 <- rm.facets( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
             pid=dat$idstud  ,  est.a.item=TRUE , maxiter=15)

# Model 5: estimated rater and item slopes
mod5 <- rm.facets( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
             pid=dat$idstud  , est.a.rater=TRUE , est.a.item=TRUE , maxiter=15)

summary(mod1)
summary(mod2)
summary(mod2a)
summary(mod3)
summary(mod4)
summary(mod5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rater model}
\keyword{Facets model}
\keyword{summary}
