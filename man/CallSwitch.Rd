\name{CallSwitch}
\alias{CallSwitch}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Switching between Calls of \R Function and \pkg{Rcpp} Function
}

\description{
Allows switching between calls of \R functions and \pkg{Rcpp} (or C)
functions within a package.
}


\usage{
CallSwitch(.NAME, ..., PACKAGE = "R")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.NAME}{
Name of the function
}
  \item{\dots}{
Function arguments
}
  \item{PACKAGE}{
Name of the package. Can be \code{"R"} for \R functions and
  (e.g.) )\code{"sirt"} for the \pkg{sirt} package.}
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

%\references{
%% ~put references to the literature/web site here ~
%}

%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[base:do.call]{base::do.call}},
\code{\link[base:.Call]{base::.Call}}
}

\examples{
\dontrun{
## The function is currently defined as	
function (.NAME, ..., PACKAGE = "R") 
{
    if (PACKAGE == "R") {
        args1 <- base::list(...)
        base::do.call(.NAME, args1)
    }
    else {
        base::.Call(.NAME, ..., PACKAGE = PACKAGE)
    }
  }
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
