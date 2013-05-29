\name{dif.strata.variance}
\alias{dif.strata.variance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratified DIF Variance
}
\description{
Calculation of stratified DIF variance
}
\usage{
dif.strata.variance(dif, se.dif, itemcluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dif}{
Vector of uniform DIF effects
}
  \item{se.dif}{
Standard error of uniform DIF effects
}
  \item{itemcluster}{
Vector of item strata
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries:
\item{stratadif}{Summary statistics of DIF effects within item strata}
\item{weighted.DIFSD}{Weighted DIF standard deviation}
\item{unweigted.DIFSD}{DIF standard deviation}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{dif.logistic.regression}} for examples.
}
%\examples{
%%% see dif.logistic.regression
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Differential item functioning (DIF)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
