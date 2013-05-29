\name{dif.logistic.regression}
\alias{dif.logistic.regression}
\title{
Differential Item Functioning using Logistic Regression Analysis
}
\description{
This function estimates differential item
functioning using a logistic regression analysis
(Zumbo, 1999). 
}
\usage{
dif.logistic.regression(dat, group, score)
}
\arguments{
  \item{dat}{
Data frame with dichotomous item responses
}
  \item{group}{
Group identifier
}
  \item{score}{
Ability estimate, e.g. the WLE.
}
}
\value{
A data frame with following variables:
\item{item}{Item name}
\item{N}{Sample size per item}
\item{R}{Value of \code{group} variable for reference group}
\item{F}{Value of \code{group} variable for focal group}
\item{nR}{Sample size per item in reference group}
\item{nF}{Sample size per item in focal group}
\item{p}{Item \eqn{p} value}
\item{pR}{Item \eqn{p} value in reference group}
\item{pF}{Item \eqn{p} value in focal group}
\item{pdiff.adj}{Adjusted \eqn{p} value difference}
\item{uniformDIF}{Uniform DIF estimate}
\item{se.uniformDIF}{Standard error of uniform DIF}
\item{t.uniformDIF}{The \eqn{t} value for uniform DIF}
\item{sig.uniformDIF}{Significance label for uniform DIF}
\item{nonuniformDIF}{Nonuniform DIF estimate}
\item{se.nonuniformDIF}{Standard error of nonuniform DIF}
\item{t.nonuniformDIF}{The \eqn{t} value for nonuniform DIF}
\item{sig.nonuniformDIF}{Significance label for nonuniform DIF}
}
\references{
Zumbo, B. D. (1999). \emph{A Handbook on the Theory and Methods of Differential 
Item Functioning (DIF): Logistic Regression Modeling as a Unitary Framework 
for Binary and Likert-type (Ordinal) Item Scores}. 
Ottawa ON: Directorate of Human Resources Research and Evaluation, 
Department of National Defense. 
}
\author{
Alexander Robitzsch
}
\seealso{
For assessing DIF variance see \code{\link{dif.variance}} and
\code{\link{dif.strata.variance}}

See the \pkg{difR} package for a large collection of DIF detection
methods.
}
\examples{
#####################################
# EXAMPLE 1: Mathematics data

data( data.math )
dat <- data.math$data
items <- grep( "M" , colnames(dat))

# estimate item parameters and WLEs
mod <- rasch.mml2( dat[,items] )
wle <- wle.rasch( dat[,items] , b=mod$item$b )$theta

# assess DIF by logistic regression
mod1 <- dif.logistic.regression( dat=dat[,items] ,
         score=wle , group=dat$female )

# calculate DIF variance
dif1 <- dif.variance( dif=mod1$uniformDIF , se.dif = mod1$se.uniformDIF )
dif1$unweighted.DIFSD
## > dif1$unweighted.DIFSD
## [1] 0.1963958

# calculate stratified DIF variance
# stratification based on domains
dif2 <- dif.strata.variance( dif=mod1$uniformDIF , se.dif = mod1$se.uniformDIF ,
        itemcluster = data.math$item$domain )
## $unweighted.DIFSD
## [1] 0.1455916
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Differential item functioning (DIF)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
