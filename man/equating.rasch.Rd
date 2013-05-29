\name{equating.rasch}
\alias{equating.rasch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Equating in the Generalized Logistic Rasch Model
}
\description{
This function does the linking in the generalized
logistic item response model. Only item difficulties (different \eqn{b}
item parameters) are allowed. Mean-mean linking and the methods 
of Haebara and Stocking-Lord are implemented (Kolen & Brennan, 2004).
}
\usage{
equating.rasch(x, y, theta = seq(-4, 4, len = 100), 
   alpha1 = 0, alpha2 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix with two columns: First column items, second column item
difficulties
}
  \item{y}{
Matrix with two columns: First columns item, second column item
difficulties
}
  \item{theta}{
Vector of theta values at which the linking functions
shall be evaluated
}
  \item{alpha1}{
Fixed \eqn{\alpha_1} parameter in the generalized item
response model
}
  \item{alpha2}{
Fixed \eqn{\alpha_2} parameter in the generalized item
response model
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{B.est}{Estimated linking constants according to the methods 
\code{Mean.Mean} (Mean-mean linking), \code{Haebara} (Haebara method)
and \code{Stocking.Lord} (Stocking-Lord method).}
\item{descriptives}{Descriptives of the linking. The linking error (\code{linkerror}) is 
calculated under the assumption of simple random sampling of items}
\item{anchor}{Original and transformed item parameters of anchor items}
\item{transf.par}{Original and transformed item parameters of all items}
}
\references{
Kolen, M. J. & Brennan, R. L. (2004). \emph{Test equating, Scaling, and Linking:
Methods and Practices}. New York: Springer.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For linking under more general item response models
see the \pkg{plink} package.
}
\examples{
##########################
# EXAMPLE 1: Linking item parameters of the PISA study

data(data.pisaPars)
pars <- data.pisaPars

# linking the tow studies with the Rasch model
mod <- equating.rasch(x=pars[,c("item","study1")] , 
			y=pars[,c("item","study2")])
##   Mean.Mean    Haebara Stocking.Lord
## 1   0.08828 0.08896269    0.09292838
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Equating}
\keyword{Linking}
