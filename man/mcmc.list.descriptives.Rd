\name{mcmc.list.descriptives}
\alias{mcmc.list.descriptives}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of Descriptive Statistics for a \code{mcmc.list} Object
}
\description{
Computation of descriptive statistics, Rhat convergence statistic
and MAP for a \code{mcmc.list} object. The Rhat statistic
is computed by splitting one Monte Carlo chain into three segments of equal 
length. The MAP is the mode estimate of the posterior distribution which is 
approximated by the mode of the kernel density estimate. 
}
\usage{
mcmc.list.descriptives( mcmcobj , quantiles=c(.025,.05,.1,.9,.95,.975) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcobj}{
Object of class \code{mcmc.list}
}
\item{quantiles}{Quantiles to be calculated for all parameters}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A data frame with descriptive statistics for all parameters in
the \code{mcmc.list} object.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
