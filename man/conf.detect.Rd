\name{conf.detect}
\alias{conf.detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confirmatory DETECT Analysis
}
\description{
This function computes the DETECT under a confirmatory specification
of item clusters (Stout, Habing, Douglas & Kim, 1996; Zhang & Stout,
1999a, 1999b).
}
\usage{
conf.detect(data, score, itemcluster, bwscale = 1.1, progress = TRUE, 
        thetagrid = seq(-3, 3, len = 200))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N} times \eqn{I} data frame of dichotomous responses.
Missing responses are allowed.
}
  \item{score}{
An ability estimate, e.g. the WLE
}
  \item{itemcluster}{
Item cluster for each item. The order of entries must correspond
to the columns in \code{data}.
}
  \item{bwscale}{
Bandwidth factor for calculation of conditional covariance 
(see \code{\link{ccov.np}})
}
  \item{progress}{
Display progress?
}
  \item{thetagrid}{
A vector which contains theta values where conditional 
covariances are evaluated.
}
}
\details{
The result of DETECT are the indices \code{DETECT}, \code{ASSI}
and \code{RATIO} (see Zhang 2007 for details) calculated 
for the options \code{unweighted} and \code{weighted}.
The option \code{unweighted} means that all conditional covariances of
item pairs are equally weighted, \code{weighted} means that
these covariances are weighted by the sample size of item pairs.
In case of multi matrix item designs, both types of indices can
differ. 

The classification scheme of these indices are as follows 
(Jang & Roussos, 2007; Zhang, 2007):

\tabular{lc}{
 Strong multidimensionality \tab DETECT > 1.00 \cr
Moderate multidimensionality \tab .40 < DETECT < 1.00 \cr
Weak multidimensionality \tab .20 < DETECT < .40 \cr
Essential unidimensionality \tab DETECT < .20 \cr
    }
\tabular{lcc}{
 Maximum value under simple structure \tab ASSI=1 \tab RATIO=1 \cr
Essential deviation from unidimensionality \tab ASSI < .25 \tab
RATIO > .36 \cr
Essential unidimensionality \tab  ASSI < .25 \tab
RATIO < .36 \cr
    }    
    
}
\value{
A list with following entries:
\item{detect}{Data frame with statistics DETECT, ASSI and RATIO}
\item{ccovtable}{Individual contributions to conditional covariance}
\item{ccov.matrix}{Evaluated conditional covariance}
}
\references{
Jang, E. E. & Roussos, L. (2007). An investigation into the dimensionality
of TOEFL using conditional covariance-based nonparametric approach.
\emph{Journal of Educational Measurement}, \bold{44}, 1-21.

Stout, W., Habing, B., Douglas, J. & Kim, H. R. (1996). 
Conditional covariance-based nonparametric multidimensionality assessment. 
\emph{Applied Psychological Measurement}, \bold{20}, 331-354.

Zhang, J. & Stout, W. (1999a). Conditional covariance structure 
of generalized compensatory multidimensional items, 
\emph{Psychometrika}, \bold{64}, 129-152.

Zhang, J. & Stout, W. (1999b). The theoretical DETECT index of 
dimensionality and its application to approximate simple structure, 
\emph{Psychometrika}, \bold{64}, 213-249. 

Zhang, J. (2007). Conditional covariance theory and DETECT for
polytomous items. \emph{Psychometrika}, \bold{72}, 69-91.
}
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
################################
# TIMSS mathematics data set
data(data.timss)

# extract data
dat <- data.timss$data
dat <- dat[ , substring( colnames(dat),1,1) == "M" ]

# extract iteminformations
iteminfo <- data.timss$item

# estimate Rasch model
mod1 <- rasch.mml2( dat )

# estimate WLEs
wle1 <- wle.rasch( dat , b = mod1$item$b )$theta

# DETECT Content Domains
detect1 <- conf.detect( data = dat , score = wle1 ,
                    itemcluster = iteminfo$Content.Domain )
##          unweighted weighted
##   DETECT      0.316    0.316
##   ASSI        0.273    0.273
##   RATIO       0.355    0.355

# DETECT Cognitive Domains
detect2 <- conf.detect( data = dat , score = wle1 ,
                    itemcluster = iteminfo$Cognitive.Domain )
##          unweighted weighted
##   DETECT      0.251    0.251
##   ASSI        0.227    0.227
##   RATIO       0.282    0.282

# DETECT Itemformat 
detect3 <- conf.detect( data = dat , score = wle1 ,
                    itemcluster = iteminfo$Format )
##          unweighted weighted
##   DETECT      0.056    0.056
##   ASSI        0.060    0.060
##   RATIO       0.062    0.062

# DETECT Block
detect4 <- conf.detect( data = dat , score = wle1 ,
                    itemcluster = iteminfo$Block )
##          unweighted weighted
##   DETECT      0.301    0.301
##   ASSI        0.193    0.193
##   RATIO       0.339    0.339

# Exploratory DETECT: Application of a cluster analysis
#   employing the Ward method
detect5 <- expl.detect( data = dat , score = wle1  , 
                nclusters = 10 , N.est = nrow(dat)  )

# Plot Cluster solution
pl <- plot( detect5$clusterfit , main = "Cluster solution" )
rect.hclust(detect5$clusterfit, k=4, border="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
