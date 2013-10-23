\name{invariance.alignment}
\alias{invariance.alignment}
\alias{summary.invariance.alignment}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment Procedure for Linking under Approximate Invariance
}
\description{
This function does alignment under approximate
invariance for \eqn{G} groups and \eqn{I} items
(Asparouhov & Muthen, 2013).
It is assumed that item loadings and intercepts are
previously estimated under the assumption of a factor
with zero mean and a variance of one. 
}
\usage{
invariance.alignment(lambda, nu, wgt, align.scale = c(1, 1), 
  align.pow = c(1, 1), eps = 1e-04, h = 0.001, max.increment = 0.2, 
  increment.factor = 1.001, maxiter = 3000, conv = 1e-04, 
  psi0.init = NULL, alpha0.init = NULL, progress = TRUE)

\method{summary}{invariance.alignment}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A \eqn{G \times I} matrix with item loadings
}
  \item{nu}{
A \eqn{G \times I}  matrix with item intercepts
}
  \item{wgt}{
A \eqn{G \times I} matrix for weighing groups 
for each item
}
  \item{align.scale}{
A vector of length two containing scale parameter
\eqn{a_\lambda} and \eqn{a_\nu} (see Details) 
}
  \item{align.pow}{
A vector of length two containing power
\eqn{p_\lambda} and \eqn{p_\nu} (see Details)
}
  \item{eps}{
A parameter in the optimization function
}
  \item{h}{
Numerical differentiation parameter 
}
  \item{max.increment}{
Maximum increment in each iteration
}
  \item{increment.factor}{
A numerical larger than one indicating the extent
of the decrease of \code{max.increment} in every
iteration.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{conv}{
Maximum parameter change of the optimization function
}
  \item{psi0.init}{
An optional vector of initial \eqn{\psi_0} parameters
}
  \item{alpha0.init}{
An optional vector of initial \eqn{\alpha_0} parameters
}
  \item{progress}{
An optional logical indicating whether computational
should be printed.
}
  \item{object}{
Object of class \code{invariance.alignment}
  }
  \item{\dots}{
Further optional arguments to be passed}
}

\details{
For \eqn{G} groups and \eqn{I} items, item loadings \eqn{\lambda_{ig0}}
and intercepts \eqn{\nu_{ig0}} are available and have been estimated
in a 1-dimensional factor analysis assuming a standardized factor.

The alignment procedure searches means \eqn{\alpha_{g0}}
and standard deviations \eqn{\psi_{g0}} using an alignment
optimization function \eqn{F}. This function is defined as
\deqn{F = \sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2} 
	f_\lambda( \lambda_{i g_1 ,1} - \lambda_{i g_2 ,1} )
	+ \sum_i \sum_{ g_1 < g_2} w_{i,g1} w_{i,g2} 
f_\nu( \nu_{i g_1 ,1} - \nu_{i g_2 ,1} ) }
where the aligned item parameters \eqn{\lambda_{i g ,1}}
and \eqn{\nu_{i g ,1}} are defined such that

\deqn{ \lambda_{i g ,1} = \lambda_{i g 0} / \psi_{g0} 
	\qquad \mbox{and} \qquad 
	\nu_{i g ,1} = \nu_{i g 0} -  \alpha_{g0} \lambda_{ig0} / \psi_{g0}
			}
and the optimization functions are defined as
\deqn{ f_\lambda (x) = [ ( x/ a_\lambda )^2 + \varepsilon ]^{p_\lambda}
	\qquad \mbox{and} \qquad
	f_\nu (x) = [ ( x/ a_\nu )^2 + \varepsilon ]^{p_\nu}
			}
using a small \eqn{ \varepsilon > 0} (e.g. .0001) to obtain
a differentiable optimization function. 

For identification reasons, the product \eqn{\Pi_g \psi_{g0}}
of all group standard deviations is set to one. The mean
\eqn{\alpha_{g0}} of the first group is set to zero.

Note that Asparouhov and Muthen (2013) use \eqn{a_\lambda=a_\nu=1}
(which can be modified in \code{align.scale})
and \eqn{p_\lambda = p_\nu = 1/4} (which can be modified in \code{align.pow}).
In case of \eqn{p_\lambda=1}, the penalty is approximately 
\eqn{f_\lambda(x) = x^2 }, in case of \eqn{p_\lambda=1/4}
it is approximately \eqn{f_\lambda(x) = \sqrt{|x|} }.

\eqn{R^2}-based effect sizes of approximate invariance have
been proposed by Asparouhov and Muthen (2013). These are
calculated separately for item loading and intercepts, resulting
in \eqn{R^2_\lambda} and \eqn{R^2_\nu} measures which are
included in the output \code{es.invariance}.
}

\value{
A list with following entries
\item{pars}{Aligned distribution parameters}
\item{itempars.aligned}{Aligned item parameters for all groups}
\item{es.invariance}{Effect sizes of approximate invariance}
\item{lambda.aligned}{Aligned \eqn{ \lambda_{i g ,1}} parameters}
\item{lambda.resid}{Residuals of \eqn{ \lambda_{i g ,1}} parameters}
\item{nu.aligned}{Aligned \eqn{ \nu_{i g ,1}} parameters}
\item{nu.resid}{Residuals of \eqn{ \nu_{i g ,1}} parameters}
\item{iter}{Number of iterations}
\item{miniter}{Iteration index with minimum optimizatiob value}
\item{fopt}{Minimum optimization value}
\item{fopt.history}{Vector of iteration history of optimization function}
\item{align.scale}{Used alignment scale parameters}
\item{align.pow}{Used alignment power parameters}
}


\references{
Asparouhov, T., & Muthen, B. (2013).
Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling}, \bold{xx}, xxx-xxx.
}

\author{
Alexander Robitzsch
}

\note{
In case of nonconvergence, decrease \code{max.increments}
or set \code{increment.factor} to a value larger than one
(maybe 1.05 or 1.1).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For IRT linking see also \code{\link{linking.haberman}}.

For modelling random item effects for loadings and intercepts
see \code{\link{mcmc.2pno.ml}}.
}

\examples{
#############################################################################
# EXAMPLE 1: Parameters cultural activities 
#############################################################################

data( data.activity.itempars )
lambda <- data.activity.itempars$lambda
nu <- data.activity.itempars$nu
Ng <-  data.activity.itempars$N
wgt <- matrix( sqrt(Ng) , length(Ng) , ncol(nu) )

#***
# Model 1: Alignment using a quadratic loss function
#   -> use the default of align.pow=c(1,1) and align.scale=c(1,1)
mod1 <- invariance.alignment( lambda , nu , wgt )
summary(mod1)

#****
# Model 2: Different powers for alignment
mod2 <- invariance.alignment( lambda , nu , wgt , 
        align.pow=c(.25,1/2) , align.scale=c(.95,.95) ,
        max.increment=.1)
summary(mod2)

# compare means from Models 1 and 2
plot( mod1$pars$alpha0 , mod2$pars$alpha0 , pch=16 ,
    xlab= "M (Model 1)" , ylab="M (Model 2)" , xlim=c(-.3,.3) , ylim=c(-.3,.3) )
lines( c(-1,1) , c(-1,1) , col="gray")    
round( cbind( mod1$pars$alpha0 , mod2$pars$alpha0 ) , 3 )
round( mod1$nu.resid , 3)
round( mod2$nu.resid ,3 )

#****
# Model 3: Low powers for alignment of scale and power
# Note that setting increment.factor larger than 1 seems necessary
mod3 <- invariance.alignment( lambda , nu , wgt , 
        align.pow=c(.25,.35) , align.scale=c(.55,.55) ,
        max.increment=.05 , increment.factor=1.1 , maxiter=500,
        psi0.init=mod1$psi0 , alpha0.init = mod1$alpha0 )
# plot iteration history
plot( 1:mod3$iter , mod3$fopt.history[ 1:mod3$iter ] , type="l" )
summary(mod3)

# compare mean and SD estimates of Models 1 and 3
plot( mod1$pars$alpha0 , mod3$pars$alpha0 , pch=16)
plot( mod1$pars$psi0 , mod3$pars$psi0 , pch=16)

# compare residuals for Models 1 and 3   
# lambda
plot( abs(as.vector(mod1$lambda.resid)) , abs(as.vector(mod3$lambda.resid)) ,
      pch=16 , xlab="Residuals lambda (Model 1)" ,
      ylab="Residuals lambda (Model 3)" , xlim=c(0,.1) , ylim=c(0,.1))
lines( c(-3,3),c(-3,3) , col="gray")    
# nu
plot( abs(as.vector(mod1$nu.resid)) , abs(as.vector(mod3$nu.resid)) ,
      pch=16 , xlab="Residuals nu (Model 1)" , ylab="Residuals nu (Model 3)" , 
      xlim=c(0,.4),ylim=c(0,.4))
lines( c(-3,3),c(-3,3) , col="gray")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linking}
\keyword{Approximate Invariance}
\keyword{Alignment}
\keyword{summary}
