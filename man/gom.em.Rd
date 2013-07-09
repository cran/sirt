\name{gom.em}
\alias{gom.em}
\alias{summary.gom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discrete (Rasch) Grade of Membership Model
}
\description{
This function estimates the grade of membership model (Erosheva, Fienberg
& Joutard, 2007) by the EM algorithm assuming a discrete membership score
distribution.
}
\usage{
gom.em(dat, K=NULL, problevels=NULL, model="GOM", theta0.k=seq(-5, 5, len=15), 
    xsi0.k=exp(seq(-6, 3, len=15)), max.increment=0.3, numdiff.parm=0.001, 
    maxdevchange=10^(-5), globconv=0.001, maxiter=1000, msteps=4, mstepconv=0.001)

\method{summary}{gom}(object,...)    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with dichotomous responses
}
  \item{K}{
Number of classes (only applies for \code{model="GOM"})
}
  \item{problevels}{
Probability levels for membership functions (only applies for \code{model="GOM"})
}
  \item{model}{
The type of grade of membership model. The default \code{"GOM"} 
is the nonparametric grade of membership model. The probabilities and membership
functions specifications described in Details are called via \code{"GOMRasch"}.
}
  \item{theta0.k}{
Vector of \eqn{\tilde{\theta}_k} grid (applies only for \code{model="GOMRasch"})
}
  \item{xsi0.k}{
Vector of \eqn{\xi_p} grid (applies only for \code{model="GOMRasch"})
}
  \item{max.increment}{
Maximum increment
}
  \item{numdiff.parm}{
Numerical differentiation parameter
}
  \item{maxdevchange}{
Convergence criterion for change in relative deviance
}
  \item{globconv}{
Global convergence criterion for parameter change
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{msteps}{
Number of iterations within a M step
}
  \item{mstepconv}{
Convergence criterion within a M step
}
  \item{object}{
Object of class \code{gom}
}
  \item{\dots}{
Further arguments to be passed
}
}
\details{
The item response model of the grade of membership model 
(Erosheva, Fienberg & Joutard, 2007) with \eqn{K} classes
for dichotomous correct responses \eqn{X_{pi}}
of person \eqn{p} on item \eqn{i} is as follows (\code{model="GOM"})
\deqn{
    P(X_{pi}=1 | g_{p1}, \ldots , g_{pK} ) = \sum_k \lambda_{ik} g_{pk}
\quad , \quad \sum_k g_{pk} = 1 
                }{%
    P(X_{pi}=1 | g_{p1}, \ldots , g_{pK} ) = \sum_k \lambda_{ik} g_{pk}
    ,    \sum_k g_{pk} = 1 
            }
The membership function is assumed to be discretely represented (specified by
\code{problevels}).
           
The Rasch grade of membership model poses constraints on probabilities
\eqn{\lambda_{ik}} and membership functions \eqn{g_{pk}}. The membership 
function of person \eqn{p} is parametrized by a location parameter \eqn{\theta_p}
and a variability parameter \eqn{\xi_p}. Each class \eqn{k} is represented by
a location \eqn{\tilde{\theta}_k}. The membership function is defined as
\deqn{ g_{pk} \propto exp[ -  ( \theta_p - \tilde{\theta}_k )^2 / ( 2 \xi_p^2 )]}
The extremal class probabilities \eqn{\lambda_{ik}} follow the Rasch model
\deqn{ \lambda_{ik} = invlogit( \tilde{\theta}_k - b_i  ) }
            
            
}
\value{
A list with following entries:
\item{deviance}{Deviance}
\item{ic}{Information criteria}
\item{item}{Data frame with item parameters}
\item{person}{Data frame with person parameters}
\item{EAP.rel}{EAP reliability (only applies for \code{model="GOMRasch"})}
\item{MAP}{Maximum aposteriori estimate of the membership function 
	}
\item{classdesc}{Descriptives for class membership}	
\item{lambda}{Estimated response probabilities \eqn{\lambda_{ik}}}
\item{se.lambda}{Standard error for stimated response probabilities \eqn{\lambda_{ik}}}
\item{mu}{Mean of \eqn{(\theta_p , \xi_p)} 
	(only applies for \code{model="GOMRasch"})}
\item{Sigma}{Covariance matrix of \eqn{(\theta_p , \xi_p)} 
	(only applies for \code{model="GOMRasch"})}
\item{b}{Estimated item difficulties (only applies for \code{model="GOMRasch"})}
\item{se.b}{Standard error of estimated difficulties 
(only applies for \code{model="GOMRasch"})}
\item{f.yi.qk}{Individual likelihood}
\item{f.qk.yi}{Individual posterior}
\item{probs}{Array with response probabilities}
\item{n.ik}{Expected counts}
\item{iter}{Number of iterations}
\item{I}{Number of items}
\item{K}{Number of classes}
\item{TP}{Number of discrete integration points for \eqn{(g_{p1},...,g_{pK})}}
\item{theta.k}{Used grid of membership functions}
\item{\dots}{Further values}
}
\references{
Erosheva, E. A., Fienberg, S. E., & Joutard, C. (2007).
Describing disability through individual-level mixture models 
for multivariate binary data. \emph{Annals of Applied Statistics},
\bold{1}, 502-537.
}
\author{
Alexander Robitzsch
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For joint maximum likelihood estimation of the grade of membership model
see \code{\link{gom.jml}}.
}
\examples{
#########################################################
# EXAMPLE 1: PISA data mathematics
#########################################################
data(data.pisaMath)
dat <- data.pisaMath$data
dat <- dat[ , grep("M" , colnames(dat)) ]

#***
# Model 1: Discrete GOM with 3 classes and 5 probability levels
problevels <- seq( 0 , 1 , len=5 )
mod1 <- gom.em( dat , K=3 , problevels ,  model="GOM"  )            
summary(mod1)

\dontrun{
#***
# Model 2: Discrete GOM with 4 classes and 5 probability levels
problevels <- seq( 0 , 1 , len=5 )
mod2 <- gom.em( dat , K=4 , problevels ,  model="GOM"  )            
summary(mod2)}

#***
# Model 3: Rasch GOM
mod3 <- gom.em( dat , model="GOMRasch" , maxiter=20 )            
summary(mod3)

#***
# Model 4: 'Ordinary' Rasch model
mod4 <- rasch.mml2( dat )
summary(mod4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Grade of membership model}
\keyword{Rasch grade of membership model}
\keyword{summary}
