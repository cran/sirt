\name{reliability.nonlinearSEM}
\alias{reliability.nonlinearSEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Reliability for Confirmatory Factor Analyses 
Based on Dichotomous Data
}
\description{
This function estimates a model based reliability using confirmatory factor
analysis (Green & Yang, 2009).
}
\usage{
reliability.nonlinearSEM(facloadings, thresh, cor.factors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{facloadings}{
Matrix of factor loadings
}
  \item{thresh}{
Vector of threshold
}
  \item{cor.factors}{
Optional matrix of covariances (correlations) between factors. The 
default is a diagonal matrix with variances of 1.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
Green, S.B., & Yang, Y. (2009). Reliability of summed item 
scores using structural equation modeling: An alternative to 
coefficient alpha. \emph{Psychometrika}, \bold{74}, 155-167. 
}
\author{
Alexander Robitzsch
}
\note{
This function needs the \pkg{mvtnorm} package.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
This function is used in \code{\link{greenyang.reliability}}.
}
%%\examples{
%% ####
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reliability}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
