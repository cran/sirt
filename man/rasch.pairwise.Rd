%% File Name: rasch.pairwise.Rd
%% File Version: 0.356

\name{rasch.pairwise}
\alias{rasch.pairwise}
\alias{summary.rasch.pairwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Estimation Method of the Rasch Model
}
\description{
This function estimates the Rasch model with a minimum chi
square estimation method (cited in Fischer, 2007, p. 544)
which is a pairwise conditional likelihood estimation approach.
}
\usage{
rasch.pairwise(dat, weights=NULL, conv=1e-04, maxiter=3000, progress=TRUE,
        b.init=NULL, zerosum=FALSE, power=1, direct_optim=TRUE)

\method{summary}{rasch.pairwise}(object, digits=3, file=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses
}
\item{weights}{Optional vector of sampling weights}
  \item{conv}{
Convergence criterion
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{progress}{
Display iteration progress?
}
  \item{b.init}{
An optional vector of length \eqn{I} of item difficulties
}
\item{zerosum}{Optional logical indicating whether item difficulties
    should be centered in each iteration. The default is that
        no centering is conducted.}
\item{power}{Power used for computing pairwise response probabilities like
in row averaging approach}
\item{direct_optim}{Logical indicating whether least squares criterion
funcion should be minimized with \code{\link[stats:nlminb]{stats::nlminb}}
}
  \item{object}{
Object of class \code{rasch.pairwise}
}
\item{digits}{Number of digits after decimal for rounding}
\item{file}{Optional file name for summary output}
  \item{\dots}{
Further arguments to be passed
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{rasch.pairwise} with following entries
\item{b}{Item difficulties}
\item{eps}{Exponentiated item difficulties, i.e. \code{eps=exp(-b)}}
\item{iter}{Number of iterations}
\item{conv}{Convergence criterion}
\item{dat}{Original data frame}
\item{freq.ij}{Frequency table of all item pairs}
\item{item}{Summary table of item parameters}
}
\references{
  Fischer, G. H. (2007). Rasch models.
  In C. R. Rao and S. Sinharay (Eds.), \emph{Handbook of Statistics},
  Vol. 26 (pp. 515-585). Amsterdam: Elsevier.
}

%\author{
%Alexander Robitzsch
%}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{summary.rasch.pairwise}} for a summary.

A slightly different implementation of this conditional pairwise method
is implemented in \cr \code{\link{rasch.pairwise.itemcluster}}.

Pairwise marginal likelihood estimation (also labeled as pseudolikelihood
estimation) can be conducted with \code{\link{rasch.pml3}}.
}
\examples{
#############################################################################
# EXAMPLE 1: Reading data set | pairwise estimation Rasch model
#############################################################################

data(data.read)
dat <- data.read

#*** Model 1: no constraint on item difficulties
mod1 <- sirt::rasch.pairwise(dat)
summary(mod1)

#*** Model 2: sum constraint on item difficulties
mod2 <- sirt::rasch.pairwise(dat, zerosum=TRUE)
summary(mod2)

\dontrun{
#** obtain standard errors by bootstrap
mod2$item$b   # extract item difficulties

# Bootstrap of item difficulties
boot_pw <- function(data, indices ){
    dd <- data[ indices, ] # bootstrap of indices
    mod <- sirt::rasch.pairwise( dat=dd, zerosum=TRUE, progress=FALSE)
    return(mod$item$b)
}
set.seed(986)
library(boot)
bmod2 <- boot::boot(data=dat, statistic=boot_pw, R=999 )
print(bmod2)
summary(bmod2)
# quantiles for bootstrap sample (and confidence interval)
apply(bmod2$t, 2, stats::quantile, probs=c(.025, .5, .975) )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{Pairwise conditional maximum likelihood (PCML)}
%% \keyword{Minchi method}
%% \keyword{summary}

