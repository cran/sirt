\name{lsem.estimate}
\alias{lsem.estimate}
\alias{summary.lsem}
\alias{plot.lsem}
\alias{lsem.MGM.stepfunctions}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Structural Equation Models (LSEM)
}
\description{
Local structural equation models (LSEM) are structural equation models (SEM)
which are evaluated for each value of a pre-defined moderator variable
(Hildenbrandt, Wilhelm, & Robitzsch, 2009).
Like in nonparametric regression models, observations near a focal point - at
which the model is evaluated - obtain higher weights, far distant obervations
obtain lower weight. The LSEM can be specified by making use of \pkg{lavaan} syntax.
It is also possible to specify a discretized version of LSEM which 
means that values of the moderator are grouped into a small number of groups.
The LSEM can be tested by employing a permutation test, see
\code{\link{lsem.permutationTest}}.
The function \code{lsem.MGM.stepfunctions} outputs stepwise functions
for a multiple group model evaluated at a grid of focal points of the 
moderator, specified in \code{moderator.grid}.
}
\usage{
lsem.estimate(data, moderator, moderator.grid, lavmodel, type="LSEM", h = 1.1, 
    residualize=TRUE, fit_measures = c("rmsea", "cfi", "tli", "gfi", "srmr"), 
    eps = 1e-08, verbose = TRUE, ...)
    
\method{summary}{lsem}(object, file=NULL, digits=3, ...)

\method{plot}{lsem}(x , parindex=NULL , ask=TRUE , ci = TRUE , lintrend = TRUE , 
       parsummary = TRUE , ylim=NULL , xlab=NULL,  ylab=NULL , main=NULL , 
       digits=3, ...)
       
lsem.MGM.stepfunctions( object , moderator.grid )       
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame
}
  \item{moderator}{
Variable name of the moderator
}
  \item{moderator.grid}{
Focal points at which the LSEM should be evaluated. If \code{type="MGM"},
breaks are defined in this vector.
}
  \item{lavmodel}{
Specified SEM in \pkg{lavaan}. The function 
\code{\link[lavaan:sem]{sem}} (\pkg{lavaan}) is used.
}
\item{type}{Type of estimated model. The default is \code{type="LSEM"} which means
	that a local structural equation model is estimated.
	A multiple group model with a discretized moderator as the 
	grouping variable can be estimated with \code{type="MGM"}. In this
	case, the breaks must be defined in \code{moderator.grid}.
			}
  \item{h}{
Bandwidth factor
}
\item{residualize}{Logical indicating whether a residualization 
    should be applied. }
  \item{fit_measures}{
Vector with names of fit measures following the labels in \pkg{lavaan}
}
  \item{eps}{
Minimum number for weights
}
  \item{verbose}{
Optional logical printing information about computation progress.
}
\item{object}{Object of class \code{lsem}}
\item{file}{A file name in which the summary output will be written.}
\item{digits}{Number of digits.}
\item{x}{Object of class \code{lsem}.}
\item{parindex}{Vector of indices for parameters in plot function.}
\item{ask}{A logical which asks for changing the graphic for each parameter.}
\item{ci}{Logical indicating whether confidence intervals should be plotted.}
\item{lintrend}{Logical indicating whether a linear trend should be plotted.}
\item{parsummary}{Logical indicating whether a parameter summary
  should be displayed.}
\item{ylim}{Plot parameter \code{ylim}. Can be a list, see Examples.}
\item{xlab}{Plot parameter \code{xlab}. Can be a vector.}
\item{ylab}{Plot parameter \code{ylab}. Can be a vector.}
\item{main}{Plot parameter \code{main}. Can be a vector.}
  \item{\dots}{
Further arguments to be passed to \code{\link[lavaan:sem]{lavaan::sem}}.
}
}

%\details{
%EXPLAIN
%}

\value{
List with following entries
\item{parameters}{Data frame with all parameters estimated at focal points of 
	moderator}
\item{weights}{Data frame with weights at each focal point}
\item{bw}{Used bandwidth}
\item{h}{Used bandwidth factor}
\item{N}{Sample size}
\item{moderator.density}{Estimated frequencies and effective sample size for 
    moderator at focal points}
\item{moderator.stat}{Descriptive statistics for moderator}
\item{moderator}{Variable name of moderator}
\item{moderator.grid}{Used grid of focal points for moderator}
\item{moderator.grouped}{Data frame with informations about grouping of
moderator if \code{type="MGM"}.}
\item{residualized.intercepts}{Estmated intercept functions used for
residualization.}
\item{lavmodel}{Used lavaan model}
\item{data}{Used data frame, possibly residualized if \code{residualize=TRUE}}
}


\references{
Hildebrandt, A., Wilhelm, O., & Robitzsch, A. (2009). Complementary and 
competing factor analytic approaches for the investigation of measurement invariance. 
\emph{Review of Psychology}, \bold{16}, 87-102.
}


\author{
Alexander Robitzsch, Oliver Luedtke, Andrea Hildebrandt
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lsem.permutationTest}}
}

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.lsem01 | Age differentiation
#############################################################################	
	
data(data.lsem01)
dat <- data.lsem01

# specify lavaan model
lavmodel <- "
        F =~ v1+v2+v3+v4+v5
        F ~~ 1*F"
        
# define grid of moderator variable age
moderator.grid <- seq(4,23,1)

#********************************
#*** Model 1: estimate LSEM with bandwidth 2
mod1 <- lsem.estimate( dat , moderator="age" , moderator.grid=moderator.grid , 
               lavmodel=lavmodel , h=2 , std.lv=TRUE)
summary(mod1)
plot(mod1 , parindex=1:5)

# perform permutation test for Model 1
pmod1 <- lsem.permutationTest( mod1 , B=10 )    
          # only for illustrative purposes the number of permutations B is set 
          # to a low number of 10
summary(pmod1)
plot(pmod1, type="global")

#********************************
#*** Model 2: estimate multiple group model with 4 age groups

# define breaks for age groups
moderator.grid <- seq( 3.5 , 23.5 , len=5) # 4 groups 
# estimate model
mod2 <- lsem.estimate( dat , moderator="age" , moderator.grid=moderator.grid , 
           lavmodel=lavmodel , type="MGM" , std.lv=TRUE)
summary(mod2)

# output step functions
smod2 <- lsem.MGM.stepfunctions( object=mod2 , moderator.grid=seq(4,23,1) )
str(smod2)
}	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{LSEM}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
