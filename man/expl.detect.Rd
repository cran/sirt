%% File Name: expl.detect.Rd
%% File Version: 1.274

\name{expl.detect}
\alias{expl.detect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exploratory DETECT Analysis
}
\description{
This function estimates the DETECT index (Stout, Habing, Douglas & Kim, 1996;
Zhang & Stout, 1999a, 1999b) in an exploratory way.
Conditional covariances of itempairs are transformed into a distance
matrix such that items are clustered by the hierarchical Ward
algorithm (Roussos, Stout & Marden, 1998). Note that the function will not
provide the same output as the original DETECT software.
}
\usage{
expl.detect(data, score, nclusters, N.est=NULL, seed=NULL, bwscale=1.1,
    smooth=TRUE, use_sum_score=FALSE, hclust_method="ward.D", estsample=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An \eqn{N \times I} data frame of dichotomous or polytomous responses.
Missing responses are allowed.
}
  \item{score}{
An ability estimate, e.g. the WLE, sum score or mean score
}
  \item{nclusters}{
Maximum number of clusters used in the exploratory analysis
}
  \item{N.est}{
Number of students in a (possible) validation of the DETECT index.
\code{N.est} students are drawn at random from \code{data}.
}
  \item{seed}{
Random seed
}
  \item{bwscale}{
Bandwidth scale factor
}
\item{smooth}{Logical indicating whether smoothing should be
applied for conditional covariance estimation}
\item{use_sum_score}{Logical indicating whether sum score should be used.
With this option, the bias corrected conditional covariance of Zhang and
Stout (1999) is used.}
\item{hclust_method}{Clustering method used as the argument
\code{method} in \code{\link[stats:hclust]{stats::hclust}}.
}
\item{estsample}{Optional vector of subject indices that defines the
estimation sample}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with following entries
\item{detect.unweighted}{Unweighted DETECT statistics}
\item{detect.weighted}{Weighted DETECT statistics. Weighting is done
proportionally to sample sizes of item pairs.}
\item{clusterfit}{Fit of the cluster method}
\item{itemcluster}{Cluster allocations}
use_sum_score
}
\references{
Roussos, L. A., Stout, W. F., & Marden, J. I. (1998). Using new proximity
measures with hierarchical cluster analysis to detect multidimensionality.
\emph{Journal of Educational Measurement,
35}, 1-30.

Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996).
Conditional covariance-based nonparametric multidimensionality assessment.
\emph{Applied Psychological Measurement, 20}, 331-354.

Zhang, J., & Stout, W. (1999a). Conditional covariance structure
of generalized compensatory multidimensional items,
\emph{Psychometrika, 64}, 129-152.

Zhang, J., & Stout, W. (1999b). The theoretical DETECT index of
dimensionality and its application to approximate simple structure,
\emph{Psychometrika, 64}, 213-249.
}


%\author{
%Alexander Robitzsch
%}


\seealso{
For examples see \code{\link{conf.detect}}.
}
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% \keyword{DETECT}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
