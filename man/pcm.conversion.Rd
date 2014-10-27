\name{pcm.conversion}
\alias{pcm.conversion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conversion of the Parameterization of the Partial Credit Model
}
\description{
Converts a parameterization of the partial credit model (see Details).
}

\usage{
pcm.conversion(b)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{b}{
Matrix of item-category-wise intercepts \eqn{b_{ik}} (see Details).
}

}
\details{
Assume that the input matrix \code{b} containing parameters \eqn{b_{ik}}
is defined according to the following paramerization of the partial credit
model
\deqn{ P( X_{pi} = k | \theta_p ) \propto exp ( k \theta_p - b_{ik} ) }
if item \eqn{i} possesses \eqn{K_i} categories.
The transformed parameterization is defined as
\deqn{b_{ik} = k \delta + \sum_{k=1}^{K_i} \tau_k \quad
\mbox{with} \quad \sum_{k=1}^{K_i} \tau_k = 0 }
The function \code{pcm.conversions} has the \eqn{\delta} and \eqn{\tau}
parameters as values.
}



\value{
List with the following entries
  \item{delta}{Vector of \eqn{\delta} parameters}
  \item{tau}{Matrix of \eqn{\tau} parameters}
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


\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: Transformation PCM for data.mg
#############################################################################

library(CDM)
data(data.mg,package="CDM")
dat <- data.mg[ 1:1000 , paste0("I",1:11) ]

#*** Model 1: estimate partial credit model in parameterization "PCM"
mod1a <- TAM::tam.mml( dat , irtmodel="PCM")
# use parameterization "PCM2"
mod1b <- TAM::tam.mml( dat , irtmodel="PCM2")
summary(mod1a)
summary(mod1b)

# convert parameterization of Model 1a into parameterization of Model 1b
b <- mod2a$item[ , c("AXsi_.Cat1","AXsi_.Cat2","AXsi_.Cat3") ]
# compare results
pcm.conversion(b)
mod2b$xsi
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Partial credit model}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
