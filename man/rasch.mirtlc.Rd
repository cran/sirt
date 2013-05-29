\name{rasch.mirtlc}
\alias{rasch.mirtlc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Latent Class 1PL and 2PL Model
}
\description{
This function estimates the multidimensional latent class
Rasch (1PL) and 2PL model (Bartolucci, 2007) for dichotomous data which
emerges from the original latent class model (Goodman, 1974).
}
\usage{
rasch.mirtlc(dat, Nclasses=NULL, modeltype="LC", dimensions=NULL ,
    group=NULL, weights=rep(1,nrow(dat)), theta.k=NULL, 
    distribution.trait= FALSE ,  range.b =c(-8,8), range.a=c(.2 , 6 ) ,  
    progress=TRUE, glob.conv=10^(-5), conv1=10^(-5), mmliter=1000, 
    mstep.maxit=3, seed=0, nstarts=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N} times \eqn{I} data frame
}
  \item{Nclasses}{
Number of latent classes. If the trait vector (or matrix)
\code{theta.k} is specified, then \code{Nclasses} is set
to the dimension of \code{theta.k}.
}
  \item{modeltype}{
Modeltype. \code{LC} is the latent class model of Goodman (1974).
\code{MLC1} is the multidimensional latent class Rasch model with
item discrimination parameter of 1. \code{MLC2} allows for the
estimation of item discriminations.
}
\item{dimensions}{
	Vector of dimension integers which allocate items to dimensions.
		}
  \item{group}{
A group identifier for multiple group estimation
}
  \item{weights}{
 Vector of sample weights
 }
  \item{theta.k}{
A grid of theta values can be specified if theta
should not be estimated. In the one-dimensional case, it
must be a vector, in the \eqn{D}-dimensional case it must
be a matrix of dimension \eqn{D}.
}
\item{distribution.trait}{A type of the assumed theta distribution can
be specified. One alternative is \code{normal}
for the normal distribution assumption.
The options \code{smooth2}, \code{smooth3} and
\code{smooth4} use the log-linear smoothing of
 Xu and von Davier (2008) to smooth the distribution
 up to two, three or four moments, respectively.
This function only works in unidimensional models.
    }
  \item{range.b}{
 Range of item difficulties which are allowed for estimation
 }
  \item{range.a}{
 Range of item slopes which are allowed for estimation
 }
  \item{progress}{
Display progress? Default is \code{TRUE}.
}
  \item{glob.conv}{
Global relative deviance convergence criterion
}
  \item{conv1}{
Item parameter convergence criterion
}
  \item{mmliter}{
Maximum number of iterations
}
  \item{mstep.maxit}{
Maximum number of iterations within an M step
}
  \item{seed}{
Set random seed for latent class estimation. A seed 
can be specified. If the seed is negative, then 
the function will generate a random seed.
}
\item{nstarts}{
	If a positive integer is provided, then a \code{nstarts}
	starts with different starting values are conducted.
}
}
\details{
The multidimensional latent class Rasch model (Bartolucci, 2007)
is an item response model which combines ideas from
latent class analysis and item response models with continuous variables.
With \code{modeltype="MLC2"} the following \eqn{D}-dimensional
item response model is estimated
\deqn{logit P(X_{pi} = 1 | \theta_p ) = a_i \theta_{pcd}- b_i}{%
logit P(X_{pi} = 1 | \theta_p ) = a_i \theta_{pcd}- b_i }
Besides the item thresholds \eqn{b_i} and item slopes \eqn{a_i}, 
for a prespecified number of latent classes \eqn{c=1,\ldots,C}{c=1,...,C}
a set of \eqn{C} \eqn{D}-dimensional \eqn{\{\theta_{cd} \}_{cd}}{{\theta_{cd} }_{cd}}
vectors are estimated.
These vectors represent the locations of latent classes. If the user
provides a grid of theta distribution \code{theta.k} as an argument in
\code{rasch.mirlc}, then the ability distribution is fixed.
}
\value{
A list with following entries
\item{pjk}{Item probabilities evaluated at discretized ability distribution}
\item{pi.k}{Estimated trait distribution}
\item{theta.k}{Discretized ability distribution}
\item{item}{Estimated item parameters}
\item{trait}{Estimated ability distribution (\code{theta.k} and \code{pi.k})}
\item{mean.trait}{Estimated mean of ability distribution}
\item{sd.trait}{Estimated standard deviation of ability distribution}
\item{skewness.trait}{Estimated skewnessof ability distribution}
\item{cor.trait}{Estimated correlation between abilities (only applies for 
	multidimensional models)}
\item{ic}{Information criteria}
\item{D}{Number of dimensions}
\item{G}{Number of groups}
\item{deviance}{Deviance}
\item{ll}{Log-likelihood}
\item{Nclasses}{Number of classes}
\item{modeltype}{Modeltype used}
\item{estep.res}{Result from E step: \code{f.qk.yi} is the individual posterior, 
	\code{f.yi.qk} is the individual likelihood}
\item{dat}{Original data frame}
\item{devL}{Vector of deviances if multiple random starts were conducted}
\item{seedL}{Vector of seed if multiple random starts were conducted}
\item{iter}{Number of iterations}
}
\references{
Bartolucci, F. (2007). A class of multidimensional
IRT models for testing unidimensionality and clustering
items. \emph{Psychometrika}, \bold{72}, 141-157.

Goodman, L. A. (1974). Exploratory latent structure 
analysis using both identifiable and unidentifiable 
models. \emph{Biometrika}, \bold{61}, 215-231.

Xu, X. & von Davier, M. (2008). \emph{Fitting the structured general 
diagnostic model to NAEP data}. ETS Research Report ETS RR-08-27. 
Princeton, ETS.
}
\author{
Alexander Robitzsch
}
\note{
For the estimation of latent class models, 
rerunning the model with starting values 
(different random seeds) is recommended. 

For fixed theta estimation in the multidimensional
case, there is a danger that too long vectors are
generated during estimation leading to memory overflow
in \R.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For an assessment of global model fit see \code{\link{modelfit.sirt}}.

The estimation of the multidimensional latent class item response
model for polytomous data can be conducted in the \pkg{MultiLCIRT}
package. Latent class analysis can be carried out with
\pkg{poLCA} and \pkg{randomLCA} packages.
}
\examples{
####################################
# Example 1: Reading data
####################################
data( data.read )
dat <- data.read

#******
# latent class models
# latent class model with 1 class
mod1 <- rasch.mirtlc( dat , Nclasses = 1 )
summary(mod1)

# latent class model with 2 classes
mod2 <- rasch.mirtlc( dat , Nclasses = 2 )
summary(mod2)

\dontrun{
# latent class model with 3 classes
mod3 <- rasch.mirtlc( dat , Nclasses = 3 , seed = - 30)  
summary(mod3)

# latent class model with 4 classes and
# 3 starts with different seeds
mod4 <- rasch.mirtlc( dat , Nclasses = 4 ,seed= -30 ,  nstarts=3 )   
# display different solutions
sort(mod1$devL)
summary(mod4)

# latent class multiple group model
# define group identifier
group <- rep( 1 , nrow(dat))
group[ 1:150 ] <- 2
mod5 <- rasch.mirtlc( dat , Nclasses = 3 , group = group )  
summary(mod5)

#*************
# Unidimensional IRT models with ordered trait

# 1PL model with 3 classes
mod11 <- rasch.mirtlc( dat , Nclasses = 3 , modeltype="MLC1" , mmliter=30)
summary(mod11)

# 1PL model with 11 classes
mod12 <- rasch.mirtlc( dat , Nclasses = 11 ,modeltype="MLC1", mmliter=30)
summary(mod12)

# 1PL model with 11 classes and fixed specified theta values
mod13 <- rasch.mirtlc( dat ,  modeltype="MLC1" , 
    theta.k = seq( -4 , 4 , len=11 ) , mmliter=100)
summary(mod13)

# 1PL model with fixed theta values and normal distribution
mod14 <- rasch.mirtlc( dat ,  modeltype="MLC1" , mmliter=30 , 
        theta.k = seq( -4 , 4 , len=11 ) , distribution.trait="normal")
summary(mod14)

# 1PL model with a smoothed trait distribution (3 moments)
mod15 <- rasch.mirtlc( dat ,  modeltype="MLC1" , mmliter=30 , 
        theta.k = seq( -4, 4 , len=11 ) ,  distribution.trait="smooth3")
summary(mod15)

# 2PL with 3 classes
mod16 <- rasch.mirtlc( dat , Nclasses=3 , modeltype="MLC2" , mmliter=30 )
summary(mod16)

# 2PL with fixed theta and smoothed distribution
mod17 <- rasch.mirtlc( dat, theta.k=seq(-4,4,len=12) , mmliter=30 ,
          modeltype="MLC2" , distribution.trait="smooth4"  )
summary(mod17)

# 1PL multiple group model with 8 classes
# define group identifier
group <- rep( 1 , nrow(dat))
group[ 1:150 ] <- 2
mod21 <- rasch.mirtlc( dat , Nclasses = 8 , modeltype="MLC1" , group=group )
summary(mod21)

#***************
# multidimensional latent class IRT models

# define vector of dimensions
dimensions <- rep( 1:3 , each = 4 )

# 3-dimensional model with 8 classes and seed 145
mod31 <- rasch.mirtlc( dat , Nclasses = 8 , mmliter=30 , 
        modeltype="MLC1" , seed = 145 , dimensions = dimensions )
summary(mod31)

# try the model above with different starting values
mod31s <- rasch.mirtlc( dat , Nclasses = 8 ,
        modeltype="MLC1" , seed = -30 , nstarts=30 , dimensions = dimensions )
summary(mod31s)     

# estimation with fixed theta vectors
# => 4^3 = 216 classes
theta.k <- seq(-4 , 4 , len=6 )
theta.k <- as.matrix( expand.grid( theta.k , theta.k , theta.k ) )
mod32 <- rasch.mirtlc( dat ,  dimensions=dimensions , 
        theta.k= theta.k , modeltype="MLC1"  )
summary(mod32)

# 3-dimensional 2PL model
mod33 <- rasch.mirtlc( dat ,  dimensions=dimensions , 
       theta.k= theta.k , modeltype="MLC2"  )
summary(mod33)

##############################################################
# Simulated Example 2: Skew trait distribution
##############################################################

N <- 1000   # number of persons
I <- 20     # number of items
theta <- sqrt( exp( rnorm( N ) ) )
theta <- theta - mean(theta )
# calculate skewness of theta distribution
mean( theta^3 ) / sd(theta)^3
# simulate item responses
dat <- sim.raschtype( theta , b=seq(-2,2,len=I ) )

# normal distribution
mod1 <- rasch.mirtlc( dat , theta.k=seq(-4,4,len=15)  , modeltype="MLC1",
            distribution.trait="normal" , mmliter=30)

# allow for skew distribution with smoothed distribution
mod2 <- rasch.mirtlc( dat , theta.k=seq(-4,4,len=15)  , modeltype="MLC1",
            distribution.trait="smooth3" , mmliter=30)

# nonparametric distribution
mod3 <- rasch.mirtlc( dat , theta.k=seq(-4,4,len=15)  , modeltype="MLC1",
            mmliter=30)

summary(mod1)
summary(mod2)            
summary(mod3)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Multidimensional Latent Class Rasch Model}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
