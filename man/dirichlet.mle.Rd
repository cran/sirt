\name{dirichlet.mle}
\alias{dirichlet.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation of the Dirichlet Distribution
}
\description{
Maximum likelihood estimation of the parameters of the Dirichlet
distribution
}
\usage{
dirichlet.mle(x, eps = 10^(-16), convcrit = 1e-05 , maxit=1000, progress=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data frame with \eqn{N} observations and \eqn{K} variables
of a Dirichlet distribution
}
  \item{eps}{
Tolerance number which is added to prevent from logarithms of zero
}
  \item{convcrit}{
Convergence criterion
}
  \item{maxit}{
Maximum number of iterations
}
  \item{progress}{
Display iteration progress?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
  \item{alpha}{Vector of \eqn{\alpha} parameters}
  \item{alpha0}{The concentration parameter \eqn{\alpha_0 = \sum_k \alpha_k} }
  \item{xsi}{Vector of proportions \eqn{\xi_k = \alpha_k / \alpha_0} } 
}
\references{
Minka, T. P. (2012). \emph{Estimating a Dirichlet distribution}.
Technical Report.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For simulating Dirichlet vectors with matrixwise
\eqn{\bold{\alpha}} parameters see \code{\link{dirichlet.simul}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# (1) simulate data
set.seed(789)
N <- 200
probs <- c(.5 , .3 , .2 )
alpha0 <- .5
alpha <- alpha0*probs
alpha <- matrix( alpha , nrow=N , ncol=length(alpha) , byrow=TRUE  )
x <- dirichlet.simul( alpha )

# (2) estimate Dirichlet parameters
dirichlet.mle(x)
##   $alpha
##   [1] 0.24507708 0.14470944 0.09590745
##   $alpha0
##   [1] 0.485694
##   $xsi
##   [1] 0.5045916 0.2979437 0.1974648
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Dirichlet distribution}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
