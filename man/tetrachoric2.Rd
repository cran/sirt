\name{tetrachoric2}
\alias{tetrachoric2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tetrachoric Correlation Matrix
}
\description{
This function estimates a tetrachoric correlation matrix according to
the Tucker method (Method 2 of Froemel, 1971).
}
\usage{
tetrachoric2(dat, delta = 0.007, maxit = 1000000, cor.smooth=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A data frame of dichotomous response
}
  \item{delta}{
The step parameter. It is set by default to \eqn{2^{-7}} 
which is approximately .007.
}
  \item{maxit}{
Maximum number of iterations.
}
\item{cor.smooth}{
Should smoothing of the tetrachoric correlation matrix
be performed to ensure positive definiteness. Choosing
\code{cor.smooth=TRUE}, the function \code{cor.smooth}
from the \pkg{psych} package is used.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{tau}{Item thresholds}
\item{rho}{Tetrachoric correlation matrix}
}
\references{
Froemel, E. C.(1971). A comparison of computer routines for the 
calculation of the tetrachoric correlation coefficient.
\emph{Psychometrika}, \bold{36}, 165-174.                
}
\author{
Alexander Robitzsch

The code is adapted from an \R script of Cengiz Zopluoglu. 
See \url{http://www.tc.umn.edu/~zoplu001/Rscript.html}.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \code{tetrachoric} function in the \pkg{psych} package.
}
\examples{
data(data.read)	
t2 <- tetrachoric2( data.read )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Tetrachoric correlation}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
