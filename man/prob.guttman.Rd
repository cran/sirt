\name{prob.guttman}
\alias{prob.guttman}
\alias{summary.prob.guttman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probabilistic Guttman Model
}
\description{
This function estimates the probabilistic Guttman model which
is a special case of an ordered latent trait model (Hanson 2000,
Proctor, 1970).
}
\usage{
prob.guttman(dat, pid = NULL, guess.equal = FALSE,  slip.equal = FALSE, 
    itemlevel = NULL, conv1 = 0.001, glob.conv = 0.001, mmliter = 500)
    
\method{summary}{prob.guttman}(object,...)    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N} times \eqn{I} data frame of dichotomous item responses
}
  \item{pid}{
Optional vector of person identifiers
}
  \item{guess.equal}{
Should the same guessing parameters for all the items estimated?
}
  \item{slip.equal}{
Should the same slipping parameters for all the items estimated?
}
  \item{itemlevel}{
A vector of item levels of the Guttman scale for each item. If there
are \eqn{K} different item levels, then the Guttman scale possesses
\eqn{K} ordered trait values.
}
  \item{conv1}{
Convergence criterion for item parameters
}
  \item{glob.conv}{
Global convergence criterion for the deviance
}
  \item{mmliter}{
Maximum number of iterations
}
  \item{object}{
Object of class \code{prob.guttman}
}
  \item{\dots}{
Further arguments to be passed
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{prob.guttman}
\item{person}{Estimated person parameters}
\item{item}{Estimated item parameters}
\item{theta.k}{Ability levels}
\item{trait}{Estimated trait distribution}
\item{ic}{Information criteria}
\item{deviance}{Deviance}
\item{iter}{Number of iterations}
\item{itemdesign}{Specified allocation of items to trait levels}
}
\references{
Hanson, B. (2000). \emph{IRT parameter estimation using the EM algorithm}. 
Technical Report.

Proctor, C. H. (1970). A probabilistic formulation and statistical analysis 
for Guttman scaling. \emph{Psychometrika}, \bold{35}, 73-78.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
A summary can be obtained using \code{\link{summary.prob.guttman}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Dataset Reading
#############################################################################
data(data.read)
dat <- data.read

#***
# Model 1: estimate probabilistic Guttman model
mod1 <- prob.guttman( dat )
summary(mod1)

#***
# Model 2: Guttman model with equal guessing and slipping parameters
mod2 <- prob.guttman( dat , guess.equal=TRUE , slip.equal=TRUE)
summary(mod2)

#***
# Model 3: Guttman model with three a priori specified item levels
itemlevel <- rep(1,12)
itemlevel[ c(2,5,8,10,12) ] <- 2
itemlevel[ c(3,4,6) ] <- 3
mod3 <- prob.guttman( dat , itemlevel=itemlevel )
summary(mod3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Probabilistic Guttman model}
\keyword{summary}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
