\name{rasch.testlet.glmer}
\alias{rasch.testlet.glmer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the Rasch Testlet Model using Generalized
Linear Mixed Effects Models (\pkg{lme4})
}
%% or MCMCglmm: This is not yet implemented
\description{
This function is a wrapper to the \pkg{lme4} generalized linear mixed effects 
software in \R for estimating item response models (De Boeck et al., 2011). 
The estimation of unidimensional models and testlet models of
Rasch type (Wang & Wilson, 2005) is implemented.
}
%% \code{MCMCglmm} is planned to implement here!
%% [...]
\usage{
rasch.testlet.glmer( dat , testlet.matrix = NULL , 
    link = "logit" , verbose = TRUE , progress=TRUE  )
}
\arguments{
  \item{dat}{
An \eqn{N} times \eqn{I} data frame of dichotomous item responses
		}
  \item{testlet.matrix}{
A vector or a matrix. If it is a vector, the for all the items
in \code{dat} a testlet identifier should be used. If it is a matrix,
then it have two columns where the first columns contains testlet names
and the second column contains item names.
}
\item{link}{
Link function which is used in estimation of the generalized linear
mixed effects model (function \code{lmer} in \R package \pkg{lme4}).
The default is the logistic link \code{"logit"} which results in the
estimation of a Rasch model. Other options are the identity link 
(\code{"linear"}) which results in a linear probability model or the 
complementary log-log link (\code{"cloglog"}) or the log-log link function
(\code{"loglog"}) (see Goldstein 1980 for choosing a link function
in item response models).
}
\item{verbose}{
Print iteration history in \code{lmer}?
}
\item{progress}{
Print progress?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with the following entries
  \item{mod}{\code{lmer}-output of estimated model}
  \item{item}{item parameters}
  \item{sd.testlet}{standard deviations of trait and testlets}
  \item{dat1}{data matrix in long format}
  \item{link}{link function}
  \item{testlet.matrix}{used design matrix when a testlet model is estimated}
  \item{eap}{EAP estimates of person abilities. In testlet models
              only the EAP estimates of the primary dimension are provided.}
  \item{eap.rel}{Reliability of EAP estimates (Adams, 2005)} 
}
\references{
Adams, R. J. (2005). Reliability as a measurement design effect.
\emph{Studies in Educational Evaluation}, \bold{31}, 162-172.

De Boeck, P., Bakker, M., Zwitser, R., Nivard, M., Hofman, A., Tuerlinckx, F., 
& Partchev, I. (2011). The estimation of item response models with the lmer function 
from the lme4 package in R. \emph{Journal of Statistical Software}, \bold{39(12)}, 1-28.

Goldstein, H. (1980). Dimensionality, bias, independence and measurement scale problems 
in latent trait test score models.
\emph{British Journal of Mathematical and Statistical Psychology}, \bold{33}, 234-246.

Wang, W.-C., & Wilson, M. (2005). The Rasch testlet model.
\emph{Applied Psychological Measurement}, \bold{29}, 126-149.
}
\author{
Alexander Robitzsch
}
\note{
This function needs the \pkg{lme4} library.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For MCMC estimation of 3PNO testlet models see \code{\link{mcmc.3pno.testlet}}.	
	
For estimating marginal item response models which account for deviations
of local independence see \code{\link{rasch.copula2}}, \code{\link{rasch.pml3}}
or \code{\link{rasch.pairwise.itemcluster}}.
}
\examples{
\dontrun{
# load reading example data set
data(data.read)

# define testlet matrix
testlet.matrix <- data.frame( substring( colnames(data.read),1,1) , colnames(data.read))

# Rasch model
mod <- rasch.testlet.glmer(data.read )

# linear probability model
mod <- rasch.testlet.glmer(data.read , link = "linear")

# cloglog link function
mod <- rasch.testlet.glmer(data.read , link = "cloglog")

# Rasch testlet model
mod <- rasch.testlet.glmer(data.read , testlet.matrix = testlet.matrix )
##   [...]
##   Standard deviations
##       effect    SD
##   1    trait 1.060
##   2        A 0.682
##   3        B 0.000
##   4        C 1.570
##   5 Residual 1.814
##   
##   EAP reliability: 0.534 

# Rasch testlet models (with two testlets and 
#   some items without any testlet structure)
mod <- rasch.testlet.glmer(data.read , 
    testlet.matrix = testlet.matrix[ c(1:4 , 9:12) , ] )

# linear probability model with testlet structure
mod <- rasch.testlet.glmer(data.read , testlet.matrix = testlet.matrix , link="linear")

# loglog link function with testlet structure
mod <- rasch.testlet.glmer(data.read , testlet.matrix = testlet.matrix , link="loglog")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Testlet model}
\keyword{Testlets}
\keyword{Local dependence}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
