\name{rm.hrm}
\alias{rm.hrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierachical Rater Model Based on Signal Detection Theory
}
\description{
This function estimates a version of the hierarchical rater
model based on signal detection theory (DeCarlo, Kim & Johnson, 2011).
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\usage{
rm.hrm(dat, pid, rater, Qmatrix = NULL, theta.k = seq(-9, 9, len = 30), 
    est.a.item = FALSE, est.c.rater = "n", est.d.rater = "n", 
    d.min = 0.5, d.max = 100, d.start = 3, max.increment = 1, numdiff.parm = 0.001, 
    maxdevchange = 0.1, globconv = .001, maxiter = 1000, msteps = 4, mstepconv = 0.001)
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
An optional \eqn{Q} matrix. If this matrix is not provided,
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
in the signal detection model. Options are \code{n} (no estimation),
\code{e} (set all parameters equal to each other), 
\code{i} (item wise estmation), \code{r} (rater wise estimation)
and \code{a} (all parameters are estimated independently
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
\deqn{ P( X_{pir} \leq k | \eta_{pi} ) =
G( c_{irk} - d_{ir} \eta_{pi} )}
where \eqn{G} is the logistic distribution function
and the categories are \eqn{k=1,\ldots , K+1}. 
The thresholds \eqn{c_{irk}} can be further restricted to
\eqn{c_{irk} = c_{k}} (\code{est.c.rater="e"}),
\eqn{c_{irk} = c_{ik}} (\code{est.c.rater="i"}) or
\eqn{c_{irk} = c_{ir}} (\code{est.c.rater="r"}). The same
holds for the rater precision parameters \eqn{d_{ir}}.
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
\item{person}{Data frame with person parameters: EAP and corresponding standard errors}
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
DeCarlo, T., Kim, Y. & Johnson, M. S. (2011).
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

######################################################
# EXAMPLE 1: HRM
######################################################
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
summary(mod2)}
            
# Model 3: Equal effects in SDT
mod3 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="e" , est.d.rater="e" , maxiter=15)
summary(mod3)

\dontrun{
# Model 4: Rater effects in SDT
mod4 <- rm.hrm( dat[ , paste0( "k",1:5) ] , rater=dat$rater , 
            pid=dat$idstud , est.c.rater="r" , est.d.rater="r" , maxiter=15)
summary(mod4)}            
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Rater model}
\keyword{Signal detection model}
