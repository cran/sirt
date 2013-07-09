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
Version: \tab 0.33 \cr
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
    for item parameters are included in \code{\link{rasch.jml.jackknife1}}
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

%% MCMC estimation multilevel IRT models    
\item Some hierarchical IRT models and random item models for dichotomous
and normally distributed data (van den Noortgate, de Boeck & Meulders, 2003;
Fox & Verhagen, 2010) can be estimated with \code{\link{mcmc.2pno.ml}}.

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

%% MCMC estimation of some models    
\item Some item response models in \pkg{sirt} can be estimated via
Markov Chain Monte Carlo (MCMC) methods. In \code{\link{mcmc.2pno}}
the two-parameter normal ogive model can be estimated. A hierarchical 
version of this model (Janssen, Tuerlinckx, Meulders & de Boeck, 2000)
is implemented in \code{\link{mcmc.2pnoh}}. The 3PNO testlet model
(Wainer, Bradlow & Wang, 2007; Glas, 2012) can be estimated with 
\code{\link{mcmc.3pno.testlet}}.
Some hierarchical IRT models and random item models 
(van den Noortgate, de Boeck & Meulders, 2003) can be estimated
with \code{\link{mcmc.2pno.ml}}.

%% NOHARM
\item For dichotomous response data,the free NOHARM software 
(McDonald, 1997) estimates the multidimensional compensatory 3PL model and the function 
\code{\link{R2noharm}} runs NOHARM from within \R. Note that NOHARM must be
downloaded from \url{http://noharm.niagararesearch.ca/nh4cldl.html}
at first.

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
Maintainer: Alexander Robitzsch <a.robitzsch@bifie.at> \cr \cr
URL: \url{https://sites.google.com/site/alexanderrobitzsch/}
}


\references{
Bartolucci, F. (2007). A class of multidimensional
IRT models for testing unidimensionality and clustering
items. \emph{Psychometrika}, \bold{72}, 141-157.

Braeken, J. (2011). A boundary mixture approach to violations of conditional 
independence. \emph{Psychometrika}, \bold{76}, 57-76.

DeCarlo, T., Kim, Y., & Johnson, M. S. (2011).
A hierarchical rater model for constructed responses, with a signal detection 
rater model. \emph{Journal of Educational Measurement}, \bold{48}, 333-356.

Dimitrov, D. (2003). Marginal true-score measures and reliability
for binary items as a function of their IRT parameters.
\emph{Applied Psychological Measurement}, \bold{27}, 440-458.

Dimitrov, D. M. (2007). Least squares distance method of cognitive validation 
and analysis for binary items using their item response theory parameters. 
\emph{Applied Psychological Measurement}, \bold{31}, 367-387.

Erosheva, E. A., Fienberg, S. E., & Joutard, C. (2007).
Describing disability through individual-level mixture models 
for multivariate binary data. \emph{Annals of Applied Statistics},
\bold{1}, 502-537.

Fox, J.-P., & Verhagen, A.-J. (2010). Random item effects modeling for 
cross-national survey data. 
In E. Davidov & P. Schmidt, and J. Billiet (Eds.), 
\emph{Cross-cultural Analysis: Methods and Applications} 
(pp. 467-488), London: Routeledge Academic.

Glas, C. A. W. (2012). \emph{Estimating and testing the extended testlet model.}
LSAC Research Report Series, RR 12-03.

Green, S.B., & Yang, Y. (2009). Reliability of summed item 
scores using structural equation modeling: An alternative to 
coefficient alpha. \emph{Psychometrika}, \bold{74}, 155-167. 

Janssen, R., Tuerlinckx, F., Meulders, M., & de Boeck, P. (2000).
A hierarchical IRT model for criterion-referenced measurement.
\emph{Journal of Educational and Behavioral Statistics}, \bold{25}, 285-306.

Linacre, J. M. (1994). \emph{Many-Facet Rasch Measurement}. 
Chicago: MESA Press.

Loken, E. & Rulison, K. L. (2010). Estimation of a four-parameter 
item response theory model. \emph{British Journal of Mathematical 
and Statistical Psychology}, \bold{63}, 509-525.

McDonald, R. P. (1997). Normal-ogive multidimensional model. 
In W. van der Linden & R. K. Hambleton (1997):
\emph{Handbook of modern item response theory} (pp. 257-269). 
New York: Springer.

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

Rijmen, F., & Vomlel, J. (2008). Assessing the performance of
variational methods for mixed logistic regression models.
\emph{Journal of Statistical Computation and Simulation},
\bold{78}, 765-779.

Renard, D., Molenberghs, G., & Geys, H. (2004). A pairwise likelihood
approach to estimation in multilevel probit models. \emph{Computational Statistics
& Data Analysis}, \bold{44}, 649-667.

Rossi, N., Wang, X. & Ramsay, J. O. (2002). Nonparametric item response function 
estimates with the EM algorithm. 
\emph{Journal of Educational and Behavioral Statistics}, 
\bold{27}, 291-317. 

Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996). 
Conditional covariance-based nonparametric multidimensionality assessment. 
\emph{Applied Psychological Measurement}, \bold{20}, 331-354.

Stukel, T. A. (1988). Generalized logistic models. 
\emph{Journal of the American Statistical Association},
\bold{83}, 426-431.

van den Noortgate, W., De Boeck, P., & Meulders, M. (2003). 
Cross-classification multilevel logistic models in psychometrics. 
\emph{Journal of Educational and Behavioral Statistics}, \bold{28}, 369-386.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
\emph{Psychometrika}, \bold{54}, 427-450.

Wainer, H., Bradlow, E. T., & Wang, X. (2007). 
\emph{Testlet response theory and its applications}. 
Cambridge: Cambridge University Press.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
\emph{Applied Psychological Measurement}, \bold{19}, 369-375.
}



\seealso{
For estimating multidimensional models for polytomous item resonses 
see the \pkg{mirt} and \pkg{TAM} (\url{https://sites.google.com/site/irttam1/})
packages. 

For conditional maximum likelihood estimation see the \pkg{eRm} package.

For pairwise estimation likelihood methods (also known as composite
likelihood methods) see \pkg{pln} or \pkg{lavaan}.

The estimation of cognitive diagnostic models is possible
using the \pkg{CDM} package.

For the multidimensional latent class IRT model see the \pkg{MultiLCIRT}
package which also allows for polytomous item responses.

Latent class analysis can be carried out with \pkg{covLCA}, \pkg{poLCA} or 
\pkg{randomLCA} packages.

Markov Chain Monte Carlo estimation for item response models can also
be found in the \pkg{MCMCpack} package (see the \code{MCMCirt} functions
therein).
}
%%\examples{
%~~ simple examples of the most important functions ~~
%%}

\keyword{package}
\section{\R Function Versions }{
\code{class.accuracy.rasch__0.03.R}, 
\code{data.prep__1.03.R}, 
\code{detect__0.02.R}, \cr
\code{dif.logisticregression__0.02.R}, 
\code{dif.variance__0.04.R}, 
\code{dirichlet__0.04.R}, \cr
\code{eigenvalues.manymatrices__0.03.R}, 
\code{equating.rasch__0.02.R}, 
\code{gom.em.alg__5.01.R}, \cr
\code{gom.em__5.01.R}, 
\code{gom.jml__0.08.R}, 
\code{gom.jml_alg__0.04.R}, 
\code{greenyang.reliability__0.03.R}, \cr
\code{latent.regression__1.06.R}, 
\code{lsdm__0.04.R}, 
\code{marginal.truescore.reliability__0.01.R}, \cr
\code{matrix_functions__0.01.R}, 
\code{mcmc.2pno.ml__3.06.R}, 
\code{mcmc.2pno.ml_alg__3.12.R}, \cr
\code{mcmc.2pno.ml_output__1.02.R}, 
\code{mcmc.2pno__1.13.R}, 
\code{mcmc.2pno_alg__1.07.R}, \cr
\code{mcmc.2pnoh__1.01.R}, 
\code{mcmc.2pnoh_alg__0.05.R}, 
\code{mcmc.3pno.testlet__2.02.R}, \cr
\code{mcmc.3pno.testlet_alg__2.05.R}, 
\code{mcmc.3pno.testlet_output__1.01.R}, 
\code{mcmc.aux__0.02.R}, \cr
\code{mcmc.list.descriptives__0.02.R}, 
\code{mcmclist2coda__0.01.R}, 
\code{modelfit.cor__2.05.R}, \cr
\code{normal2.cw__0.04.R}, 
\code{np.dich__0.08.R}, 
\code{nr.numdiff__0.01.R}, 
\code{pbivnorm2__1.03.R}, \cr
\code{personfit.stat__0.01.R}, 
\code{personfit__0.02.R}, 
\code{plausible.values.raschtype__2.03.R}, \cr
\code{plot.mcmc.sirt__0.06.R}, 
\code{prmse.subscores__0.02.R}, 
\code{prob.guttman__1.02.R}, \cr
\code{qmc.nodes__0.01.R}, 
\code{R2noharm-utility__1.01.R}, 
\code{R2noharm.EAP__0.03.R}, \cr
\code{R2noharm.jackknife__1.01.R}, 
\code{R2noharm__2.01.R}, 
\code{rasch.copula__0.992.R}, \cr
\code{rasch.copula2__4.01.R}, 
\code{rasch.jml.biascorr__0.03.R}, 
\code{rasch.jml__3.12.R}, \cr
\code{rasch.mirtlc__8.01.R}, 
\code{rasch.mirtlc_aux__8.03.R}, 
\code{rasch.mml.npirt__1.02.R}, \cr
\code{rasch.mml.ramsay__2.03.R}, 
\code{rasch.mml.raschtype__2.16.R}, 
\code{rasch.mml__2.01.R}, \cr
\code{rasch.mml2__5.03.R}, 
\code{rasch.pairwise__0.11.R}, 
\code{rasch.pml__2.11.R}, 
\code{rasch.pml_aux__1.03.R}, \cr
\code{rasch.pml2__4.07.R}, 
\code{rasch.pml2_aux__3.17.R}, 
\code{rasch.pml3__5.03.R}, 
\code{rasch.pml3_aux__5.01.R}, \cr
\code{rasch.prox__1.03.R}, 
\code{rasch.testlet__0.09.R}, 
\code{rasch.va__0.01.R}, 
\code{rm.facets__2.02.R}, \cr
\code{rm.facets_alg__2.03.R}, 
\code{rm.hrm__7.04.R}, 
\code{rm.hrm_alg__7.06.R}, 
\code{rm_proc__0.01.R}, \cr
\code{sim.rasch.dep__0.03.R}, 
\code{smirt__5.30.R}, 
\code{smirt_alg_comp__1.02.R}, \cr
\code{smirt_alg_noncomp__2.05.R}, 
\code{smirt_preproc__0.04.R}, 
\code{stratified.cronbach.alpha__0.02.R}, \cr
\code{summary.gom.em__0.04.R}, 
\code{summary.mcmc.sirt__1.02.R}, \cr
\code{summary.R2noharm.jackknife__1.01.R}, 
\code{summary.R2noharm__0.01.R}, \cr
\code{summary.rasch.copula__1.01.R}, 
\code{summary.rasch.mirtlc__7.03.R}, \cr
\code{summary.rasch.pml__0.02.R}, 
\code{summary.rm.facets__0.02.R}, 
\code{summary.rm.hrm__0.02.R}, \cr
\code{summary.smirt__0.04.R}, 
\code{tetrachoric2__0.06.R}, 
\code{wle.rasch__1.03.R}, 
\code{yen.q3__0.04.R}, \cr
\code{zzz__0.07.R}, 
}
\section{\emph{Rcpp} Function Versions }{
\code{mml2_calcpost_v1__0.01.cpp}, 
\code{mml2_raschtype_counts__0.01.cpp}, 
\code{rm.arraymult__0.03.cpp}, \cr
\code{rm.probraterfct__0.02.cpp}, 
\code{rm_calcpost__0.01.cpp}, 
\code{smirt_calcpost__0.02.cpp}, \cr
\code{smirt_calcprob_comp__0.02.cpp}, 
\code{smirt_calcprob_noncomp__0.02.cpp}, 
}
\section{\emph{Rd} Documentation Versions }{
\code{ccov.np__0.07.Rd}, 
\code{class.accuracy.rasch__0.06.Rd}, 
\code{conf.detect__1.06.Rd}, \cr
\code{data.long__0.03.Rd}, 
\code{data.math__0.04.Rd}, 
\code{data.ml1__0.02.Rd}, 
\code{data.ml2__0.01.Rd}, \cr
\code{data.noharm__0.03.Rd}, 
\code{data.pisaMath__0.04.Rd}, 
\code{data.pisaPars__0.02.Rd}, \cr
\code{data.pisaRead__0.03.Rd}, 
\code{data.ratings1__0.02.Rd}, 
\code{data.read__0.05.Rd}, \cr
\code{data.timss__0.03.Rd}, 
\code{detect.index__0.07.Rd}, 
\code{dif.logistic.regression__0.07.Rd}, \cr
\code{dif.strata.variance__0.05.Rd}, 
\code{dif.variance__0.06.Rd}, 
\code{dirichlet.mle__0.06.Rd}, \cr
\code{dirichlet.simul__0.03.Rd}, 
\code{eigenvalues.manymatrices__0.04.Rd}, \cr
\code{equating.rasch.jackknife__0.08.Rd}, 
\code{equating.rasch__1.08.Rd}, 
\code{expl.detect__1.05.Rd}, \cr
\code{gom.em__1.11.Rd}, 
\code{gom.jml__0.07.Rd}, 
\code{greenyang.reliability__1.09.Rd}, \cr
\code{latent.regression.em.raschtype__1.13.Rd}, 
\code{lsdm__1.10.Rd}, \cr
\code{marginal.truescore.reliability__0.07.Rd}, 
\code{mcmc.2pno.ml__0.10.Rd}, 
\code{mcmc.2pno__1.16.Rd}, \cr
\code{mcmc.2pnoh__0.10.Rd}, 
\code{mcmc.3pno.testlet__0.09.Rd}, 
\code{mcmc.list.descriptives__0.01.Rd}, \cr
\code{mcmclist2coda__0.03.Rd}, 
\code{modelfit.sirt__0.17.Rd}, 
\code{np.dich__0.08.Rd}, 
\code{pbivnorm2__0.04.Rd}, \cr
\code{person.parameter.rasch.copula__1.04.Rd}, 
\code{personfit.stat__0.04.Rd}, 
\code{pgenlogis__0.05.Rd}, \cr
\code{plausible.value.imputation.raschtype__1.09.Rd}, 
\code{plot.mcmc.sirt__0.04.Rd}, \cr
\code{plot.np.dich__0.09.Rd}, 
\code{prmse.subscores.scales__0.09.Rd}, 
\code{prob.guttman__1.05.Rd}, \cr
\code{qmc.nodes__0.04.Rd}, 
\code{R2noharm.EAP__0.05.Rd}, 
\code{R2noharm.jackknife__1.05.Rd}, \cr
\code{R2noharm__1.14.Rd}, 
\code{rasch.copula__1.22.Rd}, 
\code{rasch.jml.biascorr__0.09.Rd}, \cr
\code{rasch.jml.jackknife1__2.02.Rd}, 
\code{rasch.jml__1.17.Rd}, 
\code{rasch.mirtlc__2.30.Rd}, \cr
\code{rasch.mml__3.18.Rd}, 
\code{rasch.pairwise.itemcluster__0.09.Rd}, 
\code{rasch.pairwise__0.10.Rd}, \cr
\code{rasch.pml__2.27.Rd}, 
\code{rasch.prox__1.05.Rd}, 
\code{rasch.testlet.glmer__0.16.Rd}, \cr
\code{rasch.va__0.03.Rd}, 
\code{reliability.nonlinearSEM__0.06.Rd}, 
\code{rm.facets__0.12.Rd}, \cr
\code{rm.hrm__0.13.Rd}, 
\code{sim.qm.ramsay__0.11.Rd}, 
\code{sim.rasch.dep__0.13.Rd}, \cr
\code{sim.raschtype__0.06.Rd}, 
\code{sirt-package__1.43.Rd}, 
\code{smirt__1.25.Rd}, \cr
\code{stratified.cronbach.alpha__0.08.Rd}, 
\code{summary.mcmc.sirt__0.05.Rd}, \cr
\code{testlet.yen.q3__0.08.Rd}, 
\code{tetrachoric2__1.03.Rd}, 
\code{wle.rasch.jackknife__1.08.Rd}, \cr
\code{wle.rasch__1.04.Rd}, 
\code{yen.q3__0.11.Rd}, 
}
