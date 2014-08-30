\name{rm.hrm}
\alias{rm.hrm}
\alias{summary.rm.hrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierachical Rater Model Based on Signal Detection Theory (HRM-SDT)
}
\description{
This function estimates a version of the hierarchical rater
model (HRM) based on signal detection theory (HRM-SDT; DeCarlo, 2005; 
DeCarlo, Kim & Johnson, 2011).
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
rm.hrm(dat, pid, rater, Qmatrix = NULL, theta.k = seq(-9, 9, len = 30), 
    est.a.item = FALSE, est.c.rater = "n", est.d.rater = "n", 
    d.min = 0.5, d.max = 100, d.start = 3, max.increment = 1, numdiff.parm = 0.00001, 
    maxdevchange = 0.1, globconv = .001, maxiter = 1000, msteps = 4, mstepconv = 0.001)

\method{summary}{rm.hrm}(object,...)    
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
Rater identifier.
}
  \item{Qmatrix}{
An optional Q-matrix. If this matrix is not provided,
then by default the ordinary scoring of categories
(from 0 to the maximum score of \eqn{K}) is used.
}
  \item{theta.k}{
A grid of theta values for the ability distribution.
}
  \item{est.a.item}{
Should item parameters \eqn{a_i} be estimated?
}
  \item{est.c.rater}{
Type of estimation for item-rater parameters \eqn{c_{ir}}
in the signal detection model. Options are \code{'n'} (no estimation),
\code{'e'} (set all parameters equal to each other), 
\code{'i'} (item wise estmation), \code{'r'} (rater wise estimation)
and \code{'a'} (all parameters are estimated independently
from each other).
}
  \item{est.d.rater}{
Type of estimation of \eqn{d} parameters. Options are
the same as in \code{est.c.rater}.
}
  \item{d.min}{
Minimal \eqn{d} parameter to be estimated
}
  \item{d.max}{
Maximal \eqn{d} parameter to be estimated
}
  \item{d.start}{
Starting value of \eqn{d} parameters
}
  \item{max.increment}{
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
Object of class \code{rm.hrm}
}
  \item{\dots}{
Further arguments to be passed
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
The specification of the model follows DeCarlo et al. (2011).
The second level models the ideal rating (latent response) \eqn{\eta =0, ...,K}
of person \eqn{p} on item \eqn{i}
\deqn{ P( \eta_{pi} = \eta | \theta_p ) \propto
exp( a_{i} q_{ik} \theta_p - \tau_{ik} ) }
At the first level, the ratings \eqn{X_{pir}} for
person \eqn{p} on item \eqn{i} and rater \eqn{r}
are modelled as a signal detection model
\deqn{ P( X_{pir} \le k | \eta_{pi} ) =
G( c_{irk} - d_{ir} \eta_{pi} )}
where \eqn{G} is the logistic distribution function
and the categories are \eqn{k=1,\ldots , K+1}. 
The thresholds \eqn{c_{irk}} can be further restricted to
\eqn{c_{irk} = c_{k}} (\code{est.c.rater='e'}),
\eqn{c_{irk} = c_{ik}} (\code{est.c.rater='i'}) or
\eqn{c_{irk} = c_{ir}} (\code{est.c.rater='r'}). The same
holds for rater precision parameters \eqn{d_{ir}}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria and number of parameters}
\item{item}{Data frame with item parameters. The columns
\code{N} and \code{M} denote the number of oberved ratings and the
observed mean of all ratings, respectively. \cr
In addition to item parameters \eqn{\tau_{ik}} and \eqn{a_i}, the mean
for the latent response (\code{latM}) is computed as
\eqn{E( \eta_i ) = \sum_p P( \theta_p ) q_{ik} P( \eta_i = k | \theta_p ) } 
which provides an item parameter at the original metric of ratings. The latent standard
deviation (\code{latSD}) is computed in the same manner.
}
\item{rater}{Data frame with rater parameters.
Transformed \eqn{c} parameters
(\code{c_x.trans}) are computed as \eqn{c_{irk} / ( d_{ir} \cdot K )}.
}
\item{person}{Data frame with person parameters: EAP and corresponding standard 
errors}
\item{EAP.rel}{EAP reliability}
\item{sigma}{Standard deviation of the trait}
\item{tau.item}{Item parameters \eqn{\tau_{ik}}}
\item{se.tau.item}{Standard error of item parameters \eqn{\tau_{ik}}}
\item{a.item}{Item slopes \eqn{a_i}}
\item{se.a.item}{Standard error of item slopes \eqn{a_i}}
\item{c.rater}{Rater parameters \eqn{c_{irk}}}
\item{se.c.rater}{Standard error of rater severity parameter \eqn{c_{irk}}}
\item{d.rater}{Rater slope parameter \eqn{d_{ir}}}
\item{se.d.rater}{Standard error of rater slope parameter \eqn{d_{ir}}}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior distribution}
\item{probs}{Item probabilities at grid \code{theta.k}. Note that these
  probabilities are calculated on the pseudo items \eqn{i \times r}, 
  i.e. the interaction of item and rater.}
\item{prob.item}{Probabilities \eqn{P( \eta_i = \eta | \theta )} of 
	latent item responses evaluated at theta grid \eqn{\theta_p}.
}
\item{n.ik}{Expected counts}
\item{pi.k}{Estimated trait distribution \eqn{P(\theta_p)}.
}
\item{maxK}{Maximum number of categories}
\item{procdata}{Processed data}
\item{iter}{Number of iterations}
\item{\dots}{Further values}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\references{	
DeCarlo, L. T. (2005). A model of rater behavior in essay grading 
based on signal detection theory. 
\emph{Journal of Educational Measurement}, \bold{42}, 53-76.	

DeCarlo, L. T. (2010). \emph{Studies of a latent-class signal-detection model 
for constructed response scoring II: Incomplete and hierarchical designs}.
ETS Research Report ETS RR-10-08. Princeton NJ: ETS.

DeCarlo, T., Kim, Y., & Johnson, M. S. (2011).
A hierarchical rater model for constructed responses,
with a signal detection rater model.
\emph{Journal of Educational Measurement},
\bold{48}, 333-356.
}


\author{
Alexander Robitzsch
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The facets rater model can be estimated with \code{\link{rm.facets}}.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\examples{
#############################################################################
# EXAMPLE 1: Hierarchical rater model (HRM-SDT) data.ratings1
#############################################################################
data(data.ratings1)
dat <- data.ratings1

\dontrun{
# Model 1: Partial Credit Model: no rater effects
mod1 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="n" , est.d.rater="n" , 
            d.start=100 , maxiter=15)
summary(mod1)
            
# Model 2: Generalized Partial Credit Model: no rater effects
mod2 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="n" , est.d.rater="n" , 
            est.a.item =TRUE , d.start=100 , maxiter=15)
summary(mod2)
}
            
# Model 3: Equal effects in SDT
mod3 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="e" , est.d.rater="e" , maxiter=15)
summary(mod3)
  ##   Item Parameters 
  ##      item   N     M tau.Cat1 tau.Cat2 tau.Cat3 a  latM latSD
  ##   k1   k1 274 1.573   -2.164   -2.275    0.304 1 1.476 1.085
  ##   k2   k2 274 1.336   -1.883   -0.865    2.638 1 1.249 1.022
  ##   k3   k3 274 1.529   -4.414   -3.754   -0.274 1 1.481 0.872
  ##   k4   k4 274 1.372   -1.519   -1.366    3.544 1 1.234 0.974
  ##   k5   k5 274 1.412   -2.506   -2.213    2.986 1 1.300 0.894
  ##   -----------------------------------------------------------------
  ##   Rater Parameters 
  ##           item.rater  N     M     d   c_1   c_2   c_3 c_1.trans c_2.trans c_3.trans
  ##   k1-db01    k1-db01 41 0.854 3.593 0.645 5.475 9.723      0.06     0.508     0.902
  ##   k2-db01    k2-db01 41 0.951 3.593 0.645 5.475 9.723      0.06     0.508     0.902
  ##   k3-db01    k3-db01 41 1.122 3.593 0.645 5.475 9.723      0.06     0.508     0.902
  ##   k4-db01    k4-db01 41 1.024 3.593 0.645 5.475 9.723      0.06     0.508     0.902
  ##   k5-db01    k5-db01 41 1.146 3.593 0.645 5.475 9.723      0.06     0.508     0.902
  ##   k1-db02    k1-db02 37 1.351 3.593 0.645 5.475 9.723      0.06     0.508     0.902

\dontrun{
# Model 4: Rater effects in SDT
mod4 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="r" , est.d.rater="r" , maxiter=15)
summary(mod4)
}            
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rater model}
\keyword{Signal detection model}
\keyword{summary}
