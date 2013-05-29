\name{sirt-package}
\alias{sirt-package}
\alias{sirt}
\docType{package}
\title{
Supplementary Item Response Theory Models
}
\description{
Supplementary item response theory models to complement existing functions in \R. 
}
\details{
\tabular{ll}{
Package: \tab \pkg{sirt} \cr
Type: \tab Package\cr
Version: \tab 0.31 \cr
Publication Year: \tab 2013\cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{https://sites.google.com/site/alexanderrobitzsch/software} \cr
} This package enables the estimation of following models:  
\itemize{   

%% M-dim generalized item response model
\item Multidimensional marginal maximum likelihood estimation (MML)
      of generalized logistic Rasch type models using the 
      generalized logistic link function (Stukel, 1988) can be conducted 
      with \code{\link{rasch.mml2}} and the argument \code{itemtype="raschtype"}.
      This model also allows the estimation of the 4PL item
      response model (Loken & Rulison, 2010).      
      Multiple group estimation, latent regression models and
      plausible value imputation are supported.

%% M-dim noncompensatory and compensatory IRT model
\item Multidimensional noncompensatory and compensatory item response models
for dichotomous item responses (Reckase, 2009) can be estimated with the 
\code{\link{smirt}} function and the options \code{irtmodel="noncomp"} 
and \code{irtmodel="comp"}.
      
%% 1-dim Ramsay type model
\item The unidimensional quotient model (Ramsay, 1989) can be estimated
using \code{\link{rasch.mml2}} with \code{itemtype="ramsay.qm"}.

%% 1-dim nonparametric IRT models
\item Unidimensional nonparametric item response models can be estimated 
 employing MML estimation (Rossi, Wang & Ramsay, 2002) by making use of 
 \code{\link{rasch.mml2}} with \code{itemtype="npirt"}. 
 Kernel smoothing for item response function estimation (Ramsay, 1991) 
 is implemented in \code{\link{np.dich}}.

%% 1-dim Copula model
\item The unidimensional IRT copula model (Braeken, 2011) can be applied 
    for handling local dependencies, see \code{\link{rasch.copula2}}.

%% 1-dim JML
\item Unidimensional joint maximum likelihood estimation (JML) of the Rasch model 
is possible with the \code{\link{rasch.jml}} function. Bias correction methods
    for item parameters are included in \code{\link{jackknife.rasch.jml}}
    and \code{\link{rasch.jml.biascorr}}.

%% M-dim LC Rasch model
\item The multidimensional latent class Rasch and 2PL model (Bartolucci, 2007)
which employs a discrete trait distribution can be estimated
with \code{\link{rasch.mirtlc}}.

%% Rater Models
\item The unidimensional 2PL rater facets model (Lincare, 1994) can be estimated 
with \code{\link{rm.facets}}. A hierarchical rater model based on
signal detection theory (DeCarlo, Kim & Johnson, 2011) can be conducted
with \code{\link{rm.hrm}}.

%% Grade of membership model
\item The discrete grade of membership model (Erosheva, Fienberg & Joutard, 2007)
and the Rasch grade of membership model can be estimated by \code{\link{gom.em}}.

%% 1-dim PCML
\item Unidimensional pairwise conditional likelihood estimation 
    (PCML; Zwinderman, 1995) is implemented in \code{\link{rasch.pairwise}} or 
    \code{\link{rasch.pairwise.itemcluster}}.

%% 1-dim PMML
\item Unidimensional pairwise marginal likelihood estimation 
   (PMML; Renard, Molenberghs & Geys, 2004) 
    can be conducted using \code{\link{rasch.pml3}}. In this function
    local dependence can be handled by imposing residual error structure
    or omitting item pairs within a dependent item cluster from the
    estimation.

%% 1-dim Rasch model variational approximation
\item The Rasch model can be estimated by variational approximation
(Rijmen & Vomlel, 2008) using \code{\link{rasch.va}}.
    
%% 1-dim Guttman model
\item The unidimensional probabilistic Guttman model (Proctor, 1970) can be
    specified with \code{\link{prob.guttman}}.

%% jackknife WLE
\item A jackknife method for the estimation of standard errors of the 
    weighted likelihood trait estimate (Warm, 1989) is available in
    \code{\link{wle.rasch.jackknife}}.

%% reliability
\item Model based reliability for dichotomous data can be calculated by the method 
of Green and Yang (2009) with \code{\link{greenyang.reliability}} and the
marginal true score method of Dimitrov (2003) using the function 
\code{\link{marginal.truescore.reliability}}.

%% DETECT
\item Essential unidimensionality can be assessed by the DETECT
    index (Stout, Habing, Douglas & Kim, 1996), see the function 
    \code{\link{conf.detect}}.  

%% Person Fit
\item Some person fit statistics in the Rasch model (Meijer & Sijtsma, 2001) 
    are included in \code{\link{personfit.stat}}.

%% LSDM
\item An alternative to the linear logistic test model (LLTM), the
    so called least squares distance model for cognitive diagnosis 
    (LSDM; Dimitrov, 2007), can be estimated with the function 
    \code{\link{lsdm}}.
    }
    
}
\author{
Alexander Robitzsch \cr 
Federal Institute for Education Research, Innovation and Development of the 
Austrian School System (BIFIE Salzburg), Austria
\cr \cr
Maintainer: Alexander Robitzsch <a.robitzsch@bifie.at>
}


\references{
Bartolucci, F. (2007). A class of multidimensional
IRT models for testing unidimensionality and clustering
items. \emph{Psychometrika}, \bold{72}, 141-157.

Braeken, J. (2011). A boundary mixture approach to violations of conditional 
independence. \emph{Psychometrika}, \bold{76}, 57-76.

DeCarlo, T., Kim, Y. & Johnson, M. S. (2011).
A hierarchical rater model for constructed responses, with a signal detection 
rater model. \emph{Journal of Educational Measurement}, \bold{48}, 333-356.

Dimitrov, D. (2003). Marginal true-score measures and reliability
for binary items as a function of their IRT parameters.
\emph{Applied Psychological Measurement}, \bold{27}, 440-458.

Dimitrov, D. M. (2007). Least squares distance method of cognitive validation 
and analysis for binary items using their item response theory parameters. 
\emph{Applied Psychological Measurement}, \bold{31}, 367-387.

Erosheva, E. A., Fienberg, S. E. & Joutard, C. (2007).
Describing disability through individual-level mixture models 
for multivariate binary data. \emph{Annals of Applied Statistics},
\bold{1}, 502-537.

Green, S.B., & Yang, Y. (2009). Reliability of summed item 
scores using structural equation modeling: An alternative to 
coefficient alpha. \emph{Psychometrika}, \bold{74}, 155-167. 

Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}. 
Chicago: MESA Press.

Loken, E. & Rulison, K. L. (2010). Estimation of a four-parameter 
item response theory model. \emph{British Journal of Mathematical 
and Statistical Psychology}, \bold{63}, 509-525.

Meijer, R. R., & Sijtsma, K. (2001). Methodology
review: Evaluating person fit. \emph{Applied Psychological
Measurement}, \bold{25}, 107-135.

Proctor, C. H. (1970). A probabilistic formulation and statistical analysis 
for Guttman scaling. \emph{Psychometrika}, \bold{35}, 73-78.

Ramsay, J. O. (1989). A comparison of three simple test theory models. 
\emph{Psychometrika}, \bold{54}, 487-499.

Ramsay, J. O. (1991). Kernel smoothing approaches to
nonparametric item characteristic curve estimation.
\emph{Psychometrika}, \bold{56}, 611-630.

Reckase, M.D. (2009). \emph{Multidimensional item response theory}.
New York: Springer.

Rijmen, F. & Vomlel, J. (2008). Assessing the performance of
variational methods for mixed logistic regression models.
\emph{Journal of Statistical Computation and Simulation},
\bold{78}, 765-779.

Renard, D., Molenberghs, G. & Geys, H. (2004). A pairwise likelihood
approach to estimation in multilevel probit models. \emph{Computational Statistics
& Data Analysis}, \bold{44}, 649-667.

Rossi, N., Wang, X. & Ramsay, J. O. (2002). Nonparametric item response function 
estimates with the EM algorithm. 
\emph{Journal of Educational and Behavioral Statistics}, 
\bold{27}, 291-317. 

Stout, W., Habing, B., Douglas, J. & Kim, H. R. (1996). 
Conditional covariance-based nonparametric multidimensionality assessment. 
\emph{Applied Psychological Measurement}, \bold{20}, 331-354.

Stukel, T. A. (1988). Generalized logistic models. 
\emph{Journal of the American Statistical Association},
\bold{83}, 426-431.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
\emph{Psychometrika}, \bold{54}, 427-450.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
\emph{Applied Psychological Measurement}, \bold{19}, 369-375.

}
\keyword{ package }


\seealso{
For estimating multidimensional models for polytomous item resonses 
see the \pkg{mirt} and \pkg{TAM} (\url{https://sites.google.com/site/irttam1/})
packages. For dichotomous data,the free NOHARM software 
(\url{http://noharm.niagararesearch.ca/nh4cldl.html})
is a reliable alternative, see the \pkg{R2NOHARM} package 
(\url{https://sites.google.com/site/alexanderrobitzsch/software}) for running
NOHARM from within \R.

For conditional maximum likelihood estimation see the \pkg{eRm}  package.

For pairwise estimation likelihood methods (also known as composite
likelihood methods) see \pkg{pln} or \pkg{lavaan}.

The estimation of cognitive diagnostic models is possible
using the \pkg{CDM} package.

For the multidimensional latent class IRT model see the \pkg{MultiLCIRT}
package which also allows for polytomous item responses.

Latent class analysis can be carried out with
\pkg{poLCA} and \pkg{randomLCA} packages.
}
%%\examples{
%~~ simple examples of the most important functions ~~
%%}
\section{\R Function Versions }{
\code{ccov.np__0.06.Rd}, 
\code{class.accuracy.rasch__0.06.Rd}, 
\code{conf.detect__1.06.Rd}, \cr
\code{data.long__0.03.Rd}, 
\code{data.math__0.04.Rd}, 
\code{data.pisaMath__0.04.Rd}, 
\code{data.pisaPars__0.02.Rd}, \cr
\code{data.pisaRead__0.03.Rd}, 
\code{data.ratings1__0.02.Rd}, 
\code{data.read__0.05.Rd}, \cr
\code{data.timss__0.03.Rd}, 
\code{detect.index__0.07.Rd}, 
\code{dif.logistic.regression__0.07.Rd}, \cr
\code{dif.strata.variance__0.04.Rd}, 
\code{dif.variance__0.05.Rd}, 
\code{dirichlet.mle__0.06.Rd}, \cr
\code{dirichlet.simul__0.02.Rd}, 
\code{eigenvalues.manymatrices__0.04.Rd}, \cr
\code{equating.rasch.jackknife__0.07.Rd}, 
\code{equating.rasch__1.07.Rd}, 
\code{expl.detect__1.04.Rd}, \cr
\code{gom.em__1.09.Rd}, 
\code{gom.jml__0.06.Rd}, 
\code{greenyang.reliability__1.08.Rd}, \cr
\code{jackknife.rasch.jml__1.07.Rd}, 
\code{latent.regression.em.raschtype__1.11.Rd}, 
\code{lsdm__1.03.Rd}, \cr
\code{marginal.truescore.reliability__0.04.Rd}, 
\code{modelfit.sirt__0.14.Rd}, 
\code{np.dich__0.08.Rd}, \cr
\code{pbivnorm2__0.03.Rd}, 
\code{person.parameter.rasch.copula__1.03.Rd}, 
\code{personfit.stat__0.03.Rd}, \cr
\code{pgenlogis__0.05.Rd}, 
\code{plausible.value.imputation.raschtype__1.05.Rd}, \cr
\code{plot.np.dich__0.08.Rd}, 
\code{prmse.subscores.scales__0.08.Rd}, 
\code{prob.guttman__1.04.Rd}, \cr
\code{qmc.nodes__0.03.Rd}, 
\code{rasch.copula__1.20.Rd}, 
\code{rasch.jml.biascorr__0.07.Rd}, \cr
\code{rasch.jml__1.10.Rd}, 
\code{rasch.mirtlc__2.26.Rd}, 
\code{rasch.mml__3.08.Rd}, \cr
\code{rasch.pairwise.itemcluster__0.08.Rd}, 
\code{rasch.pairwise__0.09.Rd}, 
\code{rasch.pml__2.23.Rd}, \cr
\code{rasch.prox__1.04.Rd}, 
\code{rasch.testlet.glmer__0.14.Rd}, 
\code{rasch.va__0.02.Rd}, \cr
\code{reliability.nonlinearSEM__0.06.Rd}, 
\code{rm.facets__0.10.Rd}, 
\code{rm.hrm__0.08.Rd}, \cr
\code{sim.qm.ramsay__0.09.Rd}, 
\code{sim.rasch.dep__0.13.Rd}, 
\code{sim.raschtype__0.05.Rd}, \cr
\code{sirt-package__1.30.Rd}, 
\code{smirt__1.14.Rd}, 
\code{stratified.cronbach.alpha__0.08.Rd}, \cr
\code{summary.gom__0.06.Rd}, 
\code{summary.lsdm__0.07.Rd}, 
\code{summary.prob.guttman__0.07.Rd}, \cr
\code{summary.rasch.copula__0.08.Rd}, 
\code{summary.rasch.jml__0.08.Rd}, \cr
\code{summary.rasch.mirtlc__0.06.Rd}, 
\code{summary.rasch.mml__0.07.Rd}, \cr
\code{summary.rasch.pairwise__0.06.Rd}, 
\code{summary.rasch.pml__0.07.Rd}, \cr
\code{summary.rm.facets__0.06.Rd}, 
\code{summary.rm.hrm__0.06.Rd}, 
\code{summary.smirt__0.07.Rd}, \cr
\code{testlet.yen.q3__0.08.Rd}, 
\code{tetrachoric2__1.03.Rd}, 
\code{wle.rasch.jackknife__1.06.Rd}, \cr
\code{wle.rasch__1.04.Rd}, 
\code{yen.q3__0.11.Rd}, 
}
\section{\emph{Rcpp} Function Versions }{
\code{rm.arraymult__0.03.cpp}, 
\code{rm.probraterfct__0.02.cpp}, 
\code{smirt_calcpost__0.02.cpp}, \cr
\code{smirt_calcprob_comp__0.02.cpp}, 
\code{smirt_calcprob_noncomp__0.02.cpp}, 
}
\section{\emph{Rd} Documentation Versions }{
\code{ccov.np__0.06.Rd}, 
\code{class.accuracy.rasch__0.06.Rd}, 
\code{conf.detect__1.06.Rd}, \cr
\code{data.long__0.03.Rd}, 
\code{data.math__0.04.Rd}, 
\code{data.pisaMath__0.04.Rd}, 
\code{data.pisaPars__0.02.Rd}, \cr
\code{data.pisaRead__0.03.Rd}, 
\code{data.ratings1__0.02.Rd}, 
\code{data.read__0.05.Rd}, \cr
\code{data.timss__0.03.Rd}, 
\code{detect.index__0.07.Rd}, 
\code{dif.logistic.regression__0.07.Rd}, \cr
\code{dif.strata.variance__0.04.Rd}, 
\code{dif.variance__0.05.Rd}, 
\code{dirichlet.mle__0.06.Rd}, \cr
\code{dirichlet.simul__0.02.Rd}, 
\code{eigenvalues.manymatrices__0.04.Rd}, \cr
\code{equating.rasch.jackknife__0.07.Rd}, 
\code{equating.rasch__1.07.Rd}, 
\code{expl.detect__1.04.Rd}, \cr
\code{gom.em__1.09.Rd}, 
\code{gom.jml__0.06.Rd}, 
\code{greenyang.reliability__1.08.Rd}, \cr
\code{jackknife.rasch.jml__1.07.Rd}, 
\code{latent.regression.em.raschtype__1.11.Rd}, 
\code{lsdm__1.03.Rd}, \cr
\code{marginal.truescore.reliability__0.04.Rd}, 
\code{modelfit.sirt__0.14.Rd}, 
\code{np.dich__0.08.Rd}, \cr
\code{pbivnorm2__0.03.Rd}, 
\code{person.parameter.rasch.copula__1.03.Rd}, 
\code{personfit.stat__0.03.Rd}, \cr
\code{pgenlogis__0.05.Rd}, 
\code{plausible.value.imputation.raschtype__1.05.Rd}, \cr
\code{plot.np.dich__0.08.Rd}, 
\code{prmse.subscores.scales__0.08.Rd}, 
\code{prob.guttman__1.04.Rd}, \cr
\code{qmc.nodes__0.03.Rd}, 
\code{rasch.copula__1.20.Rd}, 
\code{rasch.jml.biascorr__0.07.Rd}, \cr
\code{rasch.jml__1.10.Rd}, 
\code{rasch.mirtlc__2.26.Rd}, 
\code{rasch.mml__3.08.Rd}, \cr
\code{rasch.pairwise.itemcluster__0.08.Rd}, 
\code{rasch.pairwise__0.09.Rd}, 
\code{rasch.pml__2.23.Rd}, \cr
\code{rasch.prox__1.04.Rd}, 
\code{rasch.testlet.glmer__0.14.Rd}, 
\code{rasch.va__0.02.Rd}, \cr
\code{reliability.nonlinearSEM__0.06.Rd}, 
\code{rm.facets__0.10.Rd}, 
\code{rm.hrm__0.08.Rd}, \cr
\code{sim.qm.ramsay__0.09.Rd}, 
\code{sim.rasch.dep__0.13.Rd}, 
\code{sim.raschtype__0.05.Rd}, \cr
\code{sirt-package__1.31.Rd}, 
\code{smirt__1.14.Rd}, 
\code{stratified.cronbach.alpha__0.08.Rd}, \cr
\code{summary.gom__0.06.Rd}, 
\code{summary.lsdm__0.07.Rd}, 
\code{summary.prob.guttman__0.07.Rd}, \cr
\code{summary.rasch.copula__0.08.Rd}, 
\code{summary.rasch.jml__0.08.Rd}, \cr
\code{summary.rasch.mirtlc__0.06.Rd}, 
\code{summary.rasch.mml__0.07.Rd}, \cr
\code{summary.rasch.pairwise__0.06.Rd}, 
\code{summary.rasch.pml__0.07.Rd}, \cr
\code{summary.rm.facets__0.06.Rd}, 
\code{summary.rm.hrm__0.06.Rd}, 
\code{summary.smirt__0.07.Rd}, \cr
\code{testlet.yen.q3__0.08.Rd}, 
\code{tetrachoric2__1.03.Rd}, 
\code{wle.rasch.jackknife__1.06.Rd}, \cr
\code{wle.rasch__1.04.Rd}, 
\code{yen.q3__0.11.Rd}, 
}
