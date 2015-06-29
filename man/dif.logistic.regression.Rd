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
dif.logistic.regression(dat, group, score,quant=1.645)
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
\item{quant}{Used quantile of the normal distribution 
for assessing statistical significance}
}

\details{
Items are classified into A (negligible DIF), B (moderate DIF) and 
C (large DIF) levels according to the
ETS classification system (Longford, Holland & Thayer, 1993, p. 175).
See also Monahan et al. (2007) for further DIF effect size
classifications.
}

\value{
A data frame with following variables:

\item{itemnr}{Numeric index of the item}
\item{sortDIFindex}{Rank of item with respect to the uniform DIF
	(from negative to positive values)
}
\item{item}{Item name}
\item{N}{Sample size per item}
\item{R}{Value of \code{group} variable for reference group}
\item{F}{Value of \code{group} variable for focal group}
\item{nR}{Sample size per item in reference group}
\item{nF}{Sample size per item in focal group}
\item{p}{Item \eqn{p} value}
\item{pR}{Item \eqn{p} value in reference group}
\item{pF}{Item \eqn{p} value in focal group}
\item{pdiff}{Item \eqn{p} value differences}
\item{pdiff.adj}{Adjusted \eqn{p} value difference}
\item{uniformDIF}{Uniform DIF estimate}
\item{se.uniformDIF}{Standard error of uniform DIF}
\item{t.uniformDIF}{The \eqn{t} value for uniform DIF}
\item{sig.uniformDIF}{Significance label for uniform DIF}
\item{DIF.ETS}{DIF classification according to the ETS classification
	system (see Details)
}
\item{uniform.EBDIF}{Empirical Bayes estimate of uniform DIF (Longford,
	Holland & Thayer, 1993) which takes degree of DIF standard error
        into account}
\item{DIF.SD}{Value of the DIF standard deviation}
\item{nonuniformDIF}{Nonuniform DIF estimate}
\item{se.nonuniformDIF}{Standard error of nonuniform DIF}
\item{t.nonuniformDIF}{The \eqn{t} value for nonuniform DIF}
\item{sig.nonuniformDIF}{Significance label for nonuniform DIF}
}


\references{
Longford, N. T., Holland, P. W., & Thayer, D. T. (1993).
Stability of the MH D-DIF statistics across populations.
In P. W. Holland & H. Wainer (Eds.). \emph{Differential
Item Functioning} (pp. 171-196). Hillsdale, NJ: Erlbaum.	

Monahan, P. O., McHorney, C. A., Stump, T. E., & Perkins, A. J. (2007). 
Odds ratio, delta, ETS classification, and standardization measures of 
DIF magnitude for binary logistic regression. 
\emph{Journal of Educational and Behavioral Statistics}, \bold{32}, 92-109.

Zumbo, B. D. (1999). \emph{A handbook on the theory and methods of differential 
item functioning (DIF): Logistic regression modeling as a unitary framework 
for binary and Likert-type (ordinal) item scores}. 
Ottawa ON: Directorate of Human Resources Research and Evaluation, 
Department of National Defense. 
}
\author{
Alexander Robitzsch
}
\seealso{
For assessing DIF variance see \code{\link{dif.variance}} and
\code{\link{dif.strata.variance}}

See also \code{\link{rasch.evm.pcm}} for assessing differential item
functioning in the partial credit model.

See the \pkg{difR} package for a large collection of DIF detection
methods.

For a download of the free \emph{DIF-Pack} software (SIBTEST, ...) see  
\url{http://psychometrictools.measuredprogress.org/home}.

}
\examples{
#############################################################################
# EXAMPLE 1: Mathematics data | Gender DIF
#############################################################################

data( data.math )
dat <- data.math$data
items <- grep( "M" , colnames(dat))

# estimate item parameters and WLEs
mod <- rasch.mml2( dat[,items] )
wle <- wle.rasch( dat[,items] , b=mod$item$b )$theta

# assess DIF by logistic regression
mod1 <- dif.logistic.regression( dat=dat[,items] , score=wle , group=dat$female)

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

\dontrun{
#****
# Likelihood ratio test and graphical model test in eRm package
miceadds::library_install("eRm")
# estimate Rasch model
res <- eRm::RM( dat[,items] )
summary(res)
# LR-test with respect to female
lrres <- eRm::LRtest(res, splitcr = dat$female)
summary(lrres)
# graphical model test
eRm::plotGOF(lrres) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Differential item functioning (DIF)}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
