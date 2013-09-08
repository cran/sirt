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
rasch.pairwise(dat, conv = 1e-04, maxiter = 3000, progress = TRUE, 
        b.init = NULL)

\method{summary}{rasch.pairwise}(object,...)        
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{N \times I} data frame of dichotomous item responses
}
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
  \item{object}{
Object of class \code{rasch.pairwise}
}
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
\author{
Alexander Robitzsch
}
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
data(data.read)
mod <- rasch.pairwise( data.read )
summary(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pairwise conditional maximum likelihood (PCML)}
\keyword{Minchi method}
\keyword{summary}

