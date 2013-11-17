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
(Erosheva, Fienberg & Junker, 2002; 
Erosheva, Fienberg & Joutard, 2007) with \eqn{K} classes
for dichotomous correct responses \eqn{X_{pi}}
of person \eqn{p} on item \eqn{i} is as follows (\code{model="GOM"})
\deqn{
    P(X_{pi}=1 | g_{p1}, \ldots , g_{pK} ) = \sum_k \lambda_{ik} g_{pk}
\quad , \quad \sum_{k=1}^K g_{pk} = 1
\quad , \quad 0 \leq g_{pk} \leq 1 
                }

In most applications (e.g. Erosheva et al., 2007), the grade of
memebership function \eqn{\{g_{pk}\}} is assumed to follow a Dirichlet
distribution. In our \code{gom.em} implementation             
the membership function is assumed to be discretely represented 
by a grid \eqn{u=(u_1 , \ldots , u_L)} with entries between 0 and 1
(e.g. \code{seq(0,1,length=5)} with \eqn{L=5}). 
The values \eqn{g_{pk}} of the membership function can then
only take values in \eqn{\{ u_1 , \ldots , u_L \}} with the restriction
\eqn{\sum_k g_{pk} \sum_l \bold{1}(g_{pk} = u_l ) = 1}. 
The grid \eqn{u} is specified by using the argument \code{problevels}.
           
The Rasch grade of membership model (\code{model="GOMRasch"}) poses constraints 
on probabilities \eqn{\lambda_{ik}} and membership functions \eqn{g_{pk}}. 
The membership 
function of person \eqn{p} is parametrized by a location parameter \eqn{\theta_p}
and a variability parameter \eqn{\xi_p}. Each class \eqn{k} is represented by
a location parameter \eqn{\tilde{\theta}_k}. The membership function is defined as
\deqn{ g_{pk} \propto exp[ -  ( \theta_p - \tilde{\theta}_k )^2 / ( 2 \xi_p^2 )]}
The person parameter \eqn{\theta_p} indicates the usual 'ability', while
\eqn{\xi_p} describes the individual tendency to change between classes
\eqn{1,\ldots,K} and their corresponding locations
\eqn{\tilde{\theta}_1 , \ldots ,\tilde{\theta}_K}.
The extremal class probabilities \eqn{\lambda_{ik}} follow the Rasch model
\deqn{ \lambda_{ik} = invlogit( \tilde{\theta}_k - b_i  ) }
            
See Erosheva et al. (2002), Erosheva (2005) or Galyart (in press)
for a comparison of grade of memebership models with latent trait models 
and latent class models.
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
\item{se.lambda}{Standard error for stimated response probabilities 
\eqn{\lambda_{ik}}}
\item{mu}{Mean of the distribution of \eqn{(\theta_p , \xi_p)} 
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
Erosheva, E. A. (2005). Comparing latent structures of the grade of membership, 
Rasch, and latent class models. \emph{Psychometrika}, \bold{70}, 619-628.
	
Erosheva, E. A., Fienberg, S. E., & Junker, B. W. (2002). 
Alternative statistical models and representations for large sparse 
multi-dimensional contingency tables. 
\emph{Annales-Faculte Des Sciences Toulouse Mathematiques}, \bold{11},
485-505.

Erosheva, E. A., Fienberg, S. E., & Joutard, C. (2007).
Describing disability through individual-level mixture models 
for multivariate binary data. \emph{Annals of Applied Statistics},
\bold{1}, 502-537.

Galyardt, A. (in press). 
Interpreting mixed membership models: Implications of Erosheva's representation 
theorem. In E. M. Airoldi, D. Blei, E. A. Erosheva, & S. E. Fienberg:
\emph{Handbook of Mixed Membership Models}. Chapman & Hall.
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
#############################################################################
# EXAMPLE 1: PISA data mathematics
#############################################################################

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
