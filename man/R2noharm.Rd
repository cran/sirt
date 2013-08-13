\name{R2noharm}
\alias{R2noharm}
\alias{summary.R2noharm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of a NOHARM Analysis from within \R
}
\description{
This function enables the estimation of a NOHARM analysis
(Fraser & McDonald, 1988; McDonald, 1997) from within \R. 
NOHARM estimates a compensatory multidimensional factor analysis for dichotomous 
response data. Arguments of this function strictly follows
the rules of the NOHARM manual (see Fraser & McDonald, 2012).
}


\usage{
R2noharm(dat, model.type, dimensions = NULL, guesses = rep(0, ncol(dat)), 
	    noharm.path, F.pattern = NULL, F.init = NULL, P.pattern = NULL, 
	    P.init = NULL, digits.pm = 4, writename = NULL, display.fit = 5, 
	    dec = ".", display = TRUE)

\method{summary}{R2noharm}(object, logfile=NULL , \dots)	    
}
%- maybe also 'usage' for other objects documented here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{dat}{
An \eqn{N} times \eqn{I} data frame of item responses for \eqn{N} subjects
and \eqn{I} items
}
  \item{model.type}{
Can be \code{"EFA"} (exploratory factor analysis) or \code{"CFA"}
(confirmatory factor analysis).
}
  \item{dimensions}{
Number of dimensions in exploratory factor analysis
}
  \item{guesses}{
Fixed vector of guessing parameters of length \eqn{I}. 
The default is \code{rep(0,ncol(dat))}.
}
  \item{noharm.path}{
Local path where the NOHARM 4 command line 64-bit version is being located.
}
  \item{F.pattern}{
Pattern matrix for \eqn{F} (\eqn{I} times \eqn{D})
}
  \item{F.init}{
Initial matrix for \eqn{F} (\eqn{I} times \eqn{D})
}
  \item{P.pattern}{
Pattern matrix for \eqn{P} (\eqn{D} times \eqn{D})
}
  \item{P.init}{
Initial matrix for \eqn{P} (\eqn{D} times \eqn{D})
}
  \item{digits.pm}{
Number of digits after decimal separator which are used for
estimation
}
  \item{writename}{
Name for NOHARM input and output files
  }
  \item{display.fit}{
How many digits (after decimal separator) should be used
for printing results on the \R console?
}
  \item{dec}{
Decimal separator (\code{"."} or \code{","})
}
  \item{display}{
Display output?
}
  \item{object}{
Object of class \code{R2noharm}
}
\item{logfile}{File name if the summary should be sinked into a file}
\item{\dots}{Further arguments to be passed}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
NOHARM estimates a multidimensional compensatory
item response model with the probit link function \eqn{\Phi}. 
For item responses \eqn{X_{pi}} of person \eqn{p} on
item \eqn{i} the model equation is defined as
\deqn{P( X_{pi} = 1 | \bold{\theta}_p ) = c_i + ( 1 - c_i ) 
\Phi( f_{i0} + f_{i1} \theta_{p1} + ... +
f_{iD} \theta_{pD} ) }
where \eqn{F=(f_{id})} is a loading matrix and \eqn{P}
the covariance matrix of \eqn{\bold{\theta}_p}. The guessing
parameters \eqn{c_i} must be provided as fixed inputs.

For the definition of \eqn{F} and \eqn{P} matrices, please
consult the NOHARM manual.

This function needs the 64-bit command line version which can be downloaded
at 
\url{http://noharm.niagararesearch.ca/nh4cldl.html}.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\value{
A list with following entries
\item{tanaka}{Tanaka index}
\item{rmsr}{RMSR statistic}
\item{N.itempair}{Sample sizes of pairwise item observations}
\item{pm}{Product moment matrix}
\item{guesses}{Fixed guessing parameters}
\item{residuals}{Residual covariance matrix}
\item{final.constants}{Vector of final constants}
\item{thresholds}{Threshold parameters}
\item{uniquenesses}{Item uniquenesses}
\item{loadings.theta}{Matrix of loadings in theta parametrization
    (common factor parametrization)}
\item{factor.cor}{Covariance matrix of factors}
\item{difficulties}{Item difficulties (for unidimensional models)}
\item{discriminations}{Item discriminations (for unidimensional models)}
\item{loadings}{Loading matrix (latent trait parametrization)}
\item{model.type}{Used model type}
\item{Nobs}{Number of observations}
\item{Nitems}{Number of items}
\item{modtype}{Model type according to the NOHARM specification (see NOHARM manual)}
\item{F.init}{Initial loading matrix for \eqn{F}}
\item{F.pattern}{Pattern loading matrix for \eqn{F}}
\item{P.init}{Initial covariance matrix for \eqn{P}}
\item{P.pattern}{Pattern covariance matrix for \eqn{P}}
\item{dat}{Original data frame}
\item{systime}{System time}
\item{noharm.path}{Used NOHARM directory}
\item{digits.pm}{Number of digits in product moment matrix}
\item{dec}{Used decimal symbol}
\item{display.fit}{Number of digits for fit display}
\item{dimensions}{Number of dimensions}
\item{chisquare}{Statistic \eqn{\chi^2}}
\item{Nestpars}{Number of estimated parameters}
\item{df}{Degrees of freedom}
\item{chisquare_df}{Ratio \eqn{\chi^2 / df}}
\item{rmsea}{RMSEA statistic}
\item{p.chisquare}{Significance for \eqn{\chi^2} statistic}
}
\references{
Fraser, C., & McDonald, R. P. (1988). NOHARM: Least squares item factor analysis. 
\emph{Multivariate Behavioral Research}, \bold{23}, 267-269.	
	
Fraser, C., & McDonald, R. P. (2012). \emph{NOHARM 4 Manual}. 
\url{http://noharm.niagararesearch.ca/nh4man/nhman.html}

McDonald, R. P. (1997). Normal-ogive multidimensional model. 
In W. van der Linden & R. K. Hambleton (1997):
\emph{Handbook of modern item response theory} (pp. 257-269). 
New York: Springer.
}
\author{
Alexander Robitzsch
}
\note{
In general, \code{dec=","} corresponds to computers with European
systems and \code{dec="."} to computers with American systems.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For estimating standard errors see \code{\link{R2noharm.jackknife}}.
}
%% ~~objects to See Also as \code{\link{help}}, ~~~


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
\dontrun{
# load data
data( data.noharm )
dat <- data.noharm
I <- ncol(dat) # number of items

# locate noharm.path
noharm.path <- "c:/NOHARM"

#.............................................
# EXAMPLE 1
# 1-dimensional Rasch model (1-PL model)

# estimate one factor variance
P.pattern <- matrix( 1 , ncol=1 , nrow=1 )
P.init <- P.pattern
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0 , nrow=I , ncol=1 )
F.init <- 1 + 0*F.pattern       # 
# estimate model
mod <- R2noharm( dat = dat , model.type="CFA" ,  
		F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
        P.init = P.init , writename = "ex1__1dim_1pl" , 
		noharm.path = noharm.path , dec ="," )
# summary
summary(mod , logfile="ex1__1dim_1pl__SUMMARY")
# jackknife
jmod <- R2noharm.jackknife( mod , jackunits = 20 )
summary(jmod, logfile="ex1__1dim_1pl__JACKKNIFE")
# extract factor scores
emod <- R2noharm.EAP(mod)

#.............................................
# EXAMPLE 2
# 1-dimensional 2PL Model
# set trait variance equal to 1
P.pattern <- matrix( 0 , ncol=1 , nrow=1 )
P.init <- 1+0*P.pattern
# loading matrix
F.pattern <- matrix( 1 , nrow=I , ncol=1 )
F.init <- 1 + 0*F.pattern       

mod <- R2noharm( dat = dat , model.type="CFA" ,  
		F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
        P.init = P.init , writename = "ex2__1dim_2pl" , 
		noharm.path = noharm.path , dec = "," )

summary(mod)
jmod <- R2noharm.jackknife( mod , jackunits = 20 )
summary(jmod)

#.............................................
# EXAMPLE 3
# 1-dimensional 3PL Model with fixed guessing parameters

# set trait variance equal to 1
P.pattern <- matrix( 0 , ncol=1 , nrow=1 )
P.init <- 1+0*P.pattern
# loading matrix
F.pattern <- matrix( 1 , nrow=I , ncol=1 )
F.init <- 1 + 0*F.pattern       # 
# fix guessing parameters equal to .2 (for all items)
guesses <- rep( .1 , I )

mod <- R2noharm( dat = dat , model.type="CFA" ,  
		F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
        P.init = P.init , guesses = guesses , 
		writename = "ex3__1dim_3pl" , noharm.path = noharm.path , dec=","  )
summary(mod)
jmod <- R2noharm.jackknife( mod , jackunits = 20 )
summary(jmod)

#.............................................
# EXAMPLE 4
# 3-dimensional Rasch model

# estimate one factor variance
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- .8*P.pattern
diag(P.init) <- 1
# fix all entries in the loading matrix to 1
F.init <- F.pattern <- matrix( 0 , nrow=I , ncol=3 )
F.init[1:6,1] <- 1
F.init[7:12,2] <- 1
F.init[13:18,3] <- 1
F.init

mod <- R2noharm( dat = dat , model.type="CFA" ,  
		F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
        P.init = P.init , writename = "ex4__3dim_1pl" , 
		noharm.path = noharm.path  , dec ="," )
# write output from R console in a file
summary(mod , logfile="ex4__3dim_1pl__SUMMARY.Rout")

jmod <- R2noharm.jackknife( mod , jackunits = 20 )
summary(jmod)

# extract factor scores
emod <- R2noharm.EAP(mod)

#.............................................
# EXAMPLE 5
# 3-dimensional 2PL model

# estimate one factor variance
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- .8*P.pattern
diag(P.init) <- 0
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0 , nrow=I , ncol=3 )
F.pattern[1:6,1] <- 1
F.pattern[7:12,2] <- 1
F.pattern[13:18,3] <- 1
F.init <- F.pattern

mod <- R2noharm( dat = dat , model.type="CFA" ,  
		F.pattern = F.pattern , F.init = F.init , P.pattern = P.pattern ,
        P.init = P.init , writename = "ex5__3dim_2pl" , 
		noharm.path = noharm.path  , dec = "," )

summary(mod)
# use 50 jackknife units with 4 persons within a unit
jmod <- R2noharm.jackknife( mod , jackunits = seq( 1:50 , each = 4 ) )
summary(jmod)

#.............................................
# EXAMPLE 6
# Exploratory Factor Analysis with 3 factors

mod <- R2noharm( dat = dat , model.type="EFA" ,  dimensions = 3 , 
               writename = "ex6__3dim_efa" , noharm.path = noharm.path ,
			   dec = ",")
summary(mod)

jmod <- R2noharm.jackknife( mod , jackunits = 20 )
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NOHARM}
\keyword{summary}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
