\name{yen.q3}
\alias{yen.q3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the \eqn{Q_3} Statistic (Yen, 1984)
}
\description{
This function estimates the \eqn{Q_3} statistic according to Yen (1984). 
The statistic \eqn{Q_3} is calculated for every item pair \eqn{(i,j)} 
which is the correlation between item residuals after fitting the Rasch model.
}
\usage{
yen.q3(dat, theta, b, progress=TRUE)
}
\arguments{
  \item{dat}{
An \eqn{N} times \eqn{I} data frame of dichotomous item responses
}
  \item{theta}{
Vector of length \eqn{N} of person parameter estimates (e.g. obtained from 
	\code{\link{wle.rasch}})
}
  \item{b}{
Vector of length \eqn{I} (e.g. obtained from \code{\link{rasch.mml2}})
}
  \item{progress}{
Should iteration progress be displayed?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{q3.matrix}{An \eqn{I} times \eqn{I} matrix of \eqn{Q_3} statistics}
\item{q3.long}{Just the \code{q3.matrix} in long matrix format where every row
                  corresponds to an item pair}
\item{expected}{An \eqn{N} times \eqn{I} matrix of expected probabilities 
	by the Rasch model}
\item{residual}{An \eqn{N} times \eqn{I} matrix of residuals obtained after 
	fitting the Rasch model}
\item{Q3.stat}{Vector with descriptive statistics of \eqn{Q_3}}
%% ...
%%        res <- list( "q3.matrix" = q3.matr , "q3.long" = q3.long ,
%%		"expected" = expected , "residual" = residual )
}
\references{
Yen, W. M. (1984). Effects of local item dependence on the fit and equating 
performance of the three-parameter logistic model. 
\emph{Applied Psychological Measurement}, \bold{8}, 125-145.
}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For the estimation of average the \eqn{Q_3} statistic within testlets see 
\code{\link{testlet.yen.q3}}. 

For modelling testlet effects see \code{\link{rasch.testlet.glmer}}.

For handling local dependencies in IRT models see 
\code{\link{rasch.copula2}}, \code{\link{rasch.pml3}} or \cr
\code{\link{rasch.pairwise.itemcluster}}.
}
\examples{
########################
# use data.read
data(data.read)

# estimate the Rasch model
mod <- rasch.mml2( data.read)
mod$item

# estmate WLEs
mod.wle <- wle.rasch( dat = data.read , b = mod$item$b )

# Yen's Q3 statistic
mod.q3 <- yen.q3( dat = data.read , theta = mod.wle$theta , b = mod$item$b )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Testlets}
\keyword{Local dependence}
\keyword{Q3}

