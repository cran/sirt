\name{invariance.alignment}
\alias{invariance.alignment}
\alias{summary.invariance.alignment}
\alias{plot.invariance.alignment}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alignment Procedure for Linking under Approximate Invariance
}
\description{
This function does alignment under approximate
invariance for \eqn{G} groups and \eqn{I} items
(Asparouhov & Muthen, 2014).
It is assumed that item loadings and intercepts are
previously estimated under the assumption of a factor
with zero mean and a variance of one. 
}

%% \emph{Note that this function is in an early experimental stage.}


\usage{
invariance.alignment(lambda, nu, wgt, align.scale = c(1, 1), 
    align.pow = c(1, 1), eps = 1e-04, h = 0.001, max.increment = 0.2, 
    increment.factor = 1.001, maxiter = 3000, conv = 1e-04, 
    fac.oldpar=.85 , psi0.init = NULL, alpha0.init = NULL, progress = TRUE)

\method{summary}{invariance.alignment}(object,...)

# plot of convergence behavior of optimization functions
\method{plot}{invariance.alignment}(x,...)
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
\item{fac.oldpar}{Convergence acceleration parameter between 0 and 1. 
	This parameter defines the relative weight the previous parameter value
	for the calculation of the parameter update.
	The default is .85. But experiment with this value and study
	the obtained results.
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
  \item{x}{
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

Note that the standard deviations \eqn{\psi_g} are estimated
due to minimizing the sum of \eqn{f_\lambda} functions while
means \eqn{\alpha_g} are obtained by minimizing the \eqn{f_\nu}
part with fixed \eqn{\psi_g} parameters. Therefore, the original approach
of Asparouhov and Muthen (2014) is split into a two-step procedure.

Note that Asparouhov and Muthen (2014) use \eqn{a_\lambda=a_\nu=1}
(which can be modified in \code{align.scale})
and \eqn{p_\lambda = p_\nu = 1/4} (which can be modified in \code{align.pow}).
In case of \eqn{p_\lambda=1}, the penalty is approximately 
\eqn{f_\lambda(x) = x^2 }, in case of \eqn{p_\lambda=1/4}
it is approximately \eqn{f_\lambda(x) = \sqrt{|x|} }.

Effect sizes of approximate invariance based on \eqn{R^2} have
been proposed by Asparouhov and Muthen (2014). These are
calculated separately for item loading and intercepts, resulting
in \eqn{R^2_\lambda} and \eqn{R^2_\nu} measures which are
included in the output \code{es.invariance}. In addition, 
the average correlation of aligned item parameters among groups (\code{rbar})
is reported.

\emph{Metric invariance} means that all aligned item loadings \eqn{\lambda_{ig,1}}
are equal across groups and therefore \eqn{R^2_\lambda = 1}.
\emph{Scalar invariance} means that all aligned item loadings 
\eqn{\lambda_{ig,1}} and aligned item intercepts \eqn{\nu_{ig,1}} are
equal across groups and therefore \eqn{R^2_\lambda = 1} and \eqn{R^2_\nu = 1}
(see Vandenberg & Lance, 2000).
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
\item{Niter}{Number of iterations for \eqn{f_\lambda} and
\eqn{f_\nu} optimization functions}
\item{miniter}{Iteration index with minimum optimizatiob value}
\item{fopt}{Minimum optimization value}
\item{fopt.history}{Vector of iteration history of optimization function}
\item{align.scale}{Used alignment scale parameters}
\item{align.pow}{Used alignment power parameters}
}


\references{
Asparouhov, T., & Muthen, B. (2014).
Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling}, \bold{xx}, xxx-xxx.
\url{http://www.statmodel.com/Alignment.shtml}

Vandenberg, R. J., & Lance, C. E. (2000). A review and synthesis of the
measurement invariance literature: Suggestions, practices, and
recommendations for organizational research. \emph{Organizational Research
Methods}, \bold{3}, 4-70.
}

\author{
Alexander Robitzsch
}

\note{
In case of nonconvergence, decrease \code{max.increments}
or set \code{increment.factor} to a value larger than one
(maybe 1.01 or 1.05).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For IRT linking see also \code{\link{linking.haberman}}.

For modelling random item effects for loadings and intercepts
see \code{\link{mcmc.2pno.ml}}.
}

\examples{
#############################################################################
# EXAMPLE 1: Item parameters cultural activities 
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
##   Effect Sizes of Approximate Invariance 
##          loadings intercepts
##   R2       0.9266     0.9735
##   sqrtU2   0.2709     0.1628
##   rbar     0.9265     0.9738

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
        max.increment=.05 , increment.factor=1.01 , maxiter=500,
        psi0.init=mod1$psi0 , alpha0.init = mod1$alpha0 )
# plot iteration history
plot( mod3 )
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

\dontrun{
#############################################################################
# EXAMPLE 2: Comparison 4 groups | data.inv4gr 
#############################################################################

data(data.inv4gr)
dat <- data.inv4gr
library(semTools)

model1 <- "
    F =~ I01 + I02 + I03 + I04 + I05 + I06 + I07 + I08 + I09 + I10 + I11
    F ~~ 1*F
    "

res <- measurementInvariance(model1, std.lv =TRUE , data=dat , group="group")
##   Measurement invariance tests:
##   
##   Model 1: configural invariance:
##       chisq        df    pvalue       cfi     rmsea       bic 
##     162.084   176.000     0.766     1.000     0.000 95428.025 
##   
##   Model 2: weak invariance (equal loadings):
##       chisq        df    pvalue       cfi     rmsea       bic 
##     519.598   209.000     0.000     0.973     0.039 95511.835 
##   
##   [Model 1 versus model 2]
##     delta.chisq      delta.df delta.p.value     delta.cfi 
##         357.514        33.000         0.000         0.027 
##   
##   Model 3: strong invariance (equal loadings + intercepts):
##       chisq        df    pvalue       cfi     rmsea       bic 
##    2197.260   239.000     0.000     0.828     0.091 96940.676 
##   
##   [Model 1 versus model 3]
##     delta.chisq      delta.df delta.p.value     delta.cfi 
##        2035.176        63.000         0.000         0.172 
##   
##   [Model 2 versus model 3]
##     delta.chisq      delta.df delta.p.value     delta.cfi 
##        1677.662        30.000         0.000         0.144 
##   

# extract item parameters separate group analyses
ipars <- parameterEstimates(res$fit.configural)
# extract lambda's: groups are in rows, items in columns
lambda <- matrix( ipars[ ipars$op == "=~" , "est"] , nrow=4 ,  byrow=TRUE)
colnames(lambda) <- colnames(dat)[-1]
# extract nu's
nu <- matrix( ipars[ ipars$op == "~1"  & ipars$se != 0 , "est" ], nrow=4 ,  byrow=TRUE)
colnames(nu) <- colnames(dat)[-1]
# define wgt
wgt <- 1+0*nu

# Model 1: least squares optimization
mod1 <- invariance.alignment( lambda=lambda , nu=nu ,wgt=wgt)
summary(mod1)
##   Effect Sizes of Approximate Invariance 
##          loadings intercepts
##   R2       0.7007     0.8146
##   sqrtU2   0.5471     0.4305
##   rbar     0.6237     0.7822
##   -----------------------------------------------------------------
##   Group Means and Standard Deviations 
##     alpha0  psi0
##   1  0.000 0.965
##   2 -0.105 1.098
##   3 -0.081 1.011
##   4  0.171 0.935

# Model 2: sparse target function
mod2 <- invariance.alignment( lambda=lambda , nu=nu ,wgt=wgt , align.pow=c(1/4,1/4) , 
            increment.factor=1.01)
summary(mod2)
##   Effect Sizes of Approximate Invariance 
##          loadings intercepts
##   R2       0.6980     0.8063
##   sqrtU2   0.5495     0.4402
##   rbar     0.6237     0.7839
##   -----------------------------------------------------------------
##   Group Means and Standard Deviations 
##     alpha0  psi0
##   1  0.000 0.948
##   2 -0.042 1.106
##   3  0.003 0.990
##   4  0.089 0.964

# compare residuals
round( mod1$nu.resid , 2)
round( mod2$nu.resid , 2)
##   > round( mod1$nu.resid , 2)
##          I01   I02   I03   I04   I05   I06   I07   I08   I09   I10   I11
##   [1,]  0.16  0.08  0.00  0.08 -0.19 -0.08  0.08 -0.09 -0.07 -0.23  0.07
##   [2,]  0.09  0.14 -0.12 -0.03 -0.03 -0.04  0.12 -0.15 -0.17 -0.17  0.18
##   [3,] -0.17 -0.16  0.16  0.09  0.16  0.01 -0.16  0.08  0.25  0.12 -0.12
##   [4,] -0.07 -0.06 -0.04 -0.14  0.05  0.10 -0.04  0.15  0.00  0.28 -0.13
##   > round( mod2$nu.resid , 2)
##          I01   I02   I03   I04   I05   I06   I07   I08   I09   I10   I11
##   [1,]  0.17  0.09  0.01  0.09 -0.18 -0.07  0.09 -0.08 -0.06 -0.23  0.08
##   [2,]  0.06  0.13 -0.13 -0.05 -0.05 -0.06  0.09 -0.17 -0.19 -0.19  0.14
##   [3,] -0.22 -0.20  0.12  0.06  0.14 -0.02 -0.20  0.07  0.22  0.10 -0.16
##   [4,] -0.02 -0.02  0.01 -0.09  0.10  0.15  0.03  0.19  0.04  0.32 -0.06

summary( abs( as.vector( mod1$nu.resid )) )
summary( abs( as.vector( mod2$nu.resid )) )
##   > summary( abs( as.vector( mod1$nu.resid )) )
##        Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##   0.0002706 0.0732800 0.1080000 0.1116000 0.1596000 0.2775000 
##   > summary( abs( as.vector( mod2$nu.resid )) )
##       Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.006102 0.062720 0.090370 0.112300 0.171800 0.318900

#############################################################################
# EXAMPLE 3: European Social Survey data.ess2005 
#############################################################################

data(data.ess2005)
lambda <- data.ess2005$lambda
nu <- data.ess2005$nu
wgt <- 1+0*lambda

# Model 1: least squares optimization
mod1 <- invariance.alignment( lambda=lambda , nu=nu ,wgt=wgt)
summary(mod1)
##   Group Means and Standard Deviations 
##      alpha0  psi0
##   1   0.000 0.870
##   2  -0.566 1.034
##   3   0.121 0.970
##   4  -0.724 0.899
##   5   0.028 0.976
##   6   0.168 1.028

# Model 2: sparse target function
mod2 <- invariance.alignment( lambda=lambda , nu=nu ,wgt=wgt , 
            align.pow=c(1/4,1/4) , align.scale= c( .2 , .3) , 
            increment.factor=1.01 , max.increment = .05 )
plot(mod2)
summary(mod2)
## Group Means and Standard Deviations 
##    alpha0  psi0
## 1   0.000 0.894
## 2  -0.414 1.008
## 3   0.194 1.008
## 4  -0.584 0.930
## 5   0.037 1.010
## 6   0.234 1.028

plot( mod1$pars$alpha0 , mod2$pars$alpha0 , pch=16)
lines( c(-5,5) , c(-5,5) , col="gray" )

round( mod1$nu.resid[1:10,] , 2 )
round( mod2$nu.resid[1:10,] , 2 )
##   > round( mod1$nu.resid[1:10,] , 2 )
##         ipfrule ipmodst ipbhprp imptrad
##    [1,]    0.16   -0.26    0.00    0.01
##    [2,]   -0.19    0.24    0.10   -0.24
##    [3,]    0.22   -0.34    0.06   -0.02
##    [4,]    0.29   -0.04    0.11   -0.54
##    [5,]   -0.33    0.20    0.00    0.13
##    [6,]    0.06   -0.21    0.05    0.04
##   > round( mod2$nu.resid[1:10,] , 2 )
##         ipfrule ipmodst ipbhprp imptrad
##    [1,]    0.17   -0.25    0.01    0.02
##    [2,]   -0.27    0.14    0.00   -0.30
##    [3,]    0.18   -0.37    0.01   -0.05
##    [4,]    0.18   -0.13    0.00   -0.61
##    [5,]   -0.32    0.20    0.01    0.13
##    [6,]    0.01   -0.23    0.01    0.01
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linking}
\keyword{Approximate Invariance}
\keyword{Alignment}
\keyword{summary}
