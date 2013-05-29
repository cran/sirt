\name{pgenlogis}
\alias{pgenlogis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Probabilities in Generalized Logistic Item Response Models
}
\description{
Calculation of probabilities in the generalized logistic
item response model
}
\usage{
pgenlogis(x, alpha1 = 0, alpha2 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector
}
  \item{alpha1}{
Upper tail parameter \eqn{\alpha_1} in the generalized 
logistic item response model. The default is 0.
}
  \item{alpha2}{
Lower tail parameter \eqn{\alpha_2} parameter in the generalized 
logistic item response model. The default is 0.
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
Stukel, T. A. (1988). Generalized logistic models. 
\emph{Journal of the American Statistical Association},
\bold{83}, 426-431.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pgenlogis( x=c(-.3 , 0 , .25 , 1 ) , alpha1=0 , alpha2= .6 )
##   [1] 0.4185580 0.5000000 0.5621765 0.7310586
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Generalized logistic item response model}
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
