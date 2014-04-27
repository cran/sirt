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
Version: \tab 0.44 \cr
Publication Year: \tab 2014\cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{https://sites.google.com/site/alexanderrobitzsch/software} \cr
} 


This package enables the estimation of following models: 

\itemize{   

%% M-dim generalized item response model
\item Multidimensional marginal maximum likelihood estimation (MML)
      of generalized logistic Rasch type models using the 
      generalized logistic link function (Stukel, 1988) can be conducted 
      with \code{\link{rasch.mml2}} and the argument \code{itemtype="raschtype"}.
      This model also allows the estimation of the 4PL item
      response model (Loken & Rulison, 2010).      
      Multiple group estimation, latent regression models and
      plausible value imputation are supported. In addition, pseudo-likelihood
      estimation for fractional item response data can be conducted.

%% M-dim noncompensatory and compensatory IRT model
\item Multidimensional noncompensatory, compensatory and partially compensatory
item response models for dichotomous item responses (Reckase, 2009) can be estimated 
with the \code{\link{smirt}} function and the options \code{irtmodel="noncomp"} 
, \code{irtmodel="comp"} and \code{irtmodel="partcomp"}.
      
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
\item The multidimensional IRT copula model (Braeken, 2011) can be applied 
    for handling local dependencies, see \code{\link{rasch.copula3}}.

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
with \code{\link{rm.hrm}}. A simple latent class model for two exchangeable
raters is implemented in \code{\link{lc.2raters}}.

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
    estimation. \cr    
    The function \code{\link{rasch.evm.pcm}} estimates the mutiple group
    partial credit model based on the pairwise eigenvector approach
    which avoids iterative estimation.

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
\item For dichotomous response data, the free NOHARM software 
(McDonald, 1997) estimates the multidimensional compensatory 3PL model and the function 
\code{\link{R2noharm}} runs NOHARM from within \R. Note that NOHARM must be
downloaded from \url{http://noharm.niagararesearch.ca/nh4cldl.html}
at first. A pure \R implementation of the NOHARM model with some extensions
can be found in \code{\link{noharm.sirt}}.

%% Nonparametric item response theory / ISOP model
\item The measurement theoretic founded nonparametric item
response models of Scheiblechner (1995, 1999) -- the ISOP
and the ADISOP model -- can be estimated with 
\code{\link{isop.dich}} or \code{\link{isop.poly}}.
Item scoring within this theory can be conducted with
\code{\link{isop.scoring}}.

%% Functional unidimensional item response model
\item The functional unidimensional item response model
(Ip et al., 2013) can be estimated with \code{\link{f1d.irt}}.

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

%% linking / alignment
\item Item parameters from several studies can be linked using the Haberman
method (Haberman, 2009) in \code{\link{linking.haberman}}. See also
\code{\link{equating.rasch}} and \code{\link{linking.robust}}. 
The alignment procedure (Asparouhov & Muthen, 2013)
\code{\link{invariance.alignment}} is originally for comfirmatory factor 
analysis and aims at obtaining approximate invariance.
    
    
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
Asparouhov, T., & Muthen, B. (2014).
Multiple-group factor analysis alignment.
\emph{Structural Equation Modeling}, \bold{xx}, xxx-xxx.	

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
In E. Davidov, P. Schmidt, & J. Billiet (Eds.), 
\emph{Cross-cultural Analysis: Methods and Applications} 
(pp. 467-488), London: Routledge Academic.

Glas, C. A. W. (2012). \emph{Estimating and testing the extended testlet model.}
LSAC Research Report Series, RR 12-03.

Green, S.B., & Yang, Y. (2009). Reliability of summed item 
scores using structural equation modeling: An alternative to 
coefficient alpha. \emph{Psychometrika}, \bold{74}, 155-167. 

Haberman, S. J. (2009). \emph{Linking parameter estimates derived
from an item respone model through separate calibrations}. 
ETS Research Report ETS RR-09-40. Princeton, ETS.

Ip, E. H., Molenberghs, G., Chen, S. H., Goegebeur, Y., & 
De Boeck, P. (2013). Functionally unidimensional item 
response models for multivariate binary data. 
\emph{Multivariate Behavioral Research}, \bold{48}, 534-562.

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

Rusch, T., Mair, P., & Hatzinger, R. (2013). 
\emph{Psychometrics with \R: A Review of CRAN Packages for Item Response Theory.}
\url{http://epub.wu.ac.at/4010/1/resrepIRThandbook.pdf}.

Scheiblechner, H. (1995). Isotonic ordinal
probabilistic models (ISOP). \emph{Psychometrika},
\bold{60}, 281-304.

Scheiblechner, H. (1999). Additive conjoint isotonic
probabilistic models (ADISOP). \emph{Psychometrika},
\bold{64}, 295-316.

Stout, W., Habing, B., Douglas, J., & Kim, H. R. (1996). 
Conditional covariance-based nonparametric multidimensionality assessment. 
\emph{Applied Psychological Measurement}, \bold{20}, 331-354.

Stukel, T. A. (1988). Generalized logistic models. 
\emph{Journal of the American Statistical Association},
\bold{83}, 426-431.

Uenlue, A., & Yanagida, T. (2011). \R you ready for \R?: 
The CRAN psychometrics task view.
\emph{British Journal of Mathematical and Statistical Psychology},
\bold{64}, 182-186.

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
see the \pkg{mirt} and \pkg{TAM} packages. 

For conditional maximum likelihood estimation see the \pkg{eRm} package.

For pairwise estimation likelihood methods (also known as composite
likelihood methods) see \pkg{pln} or \pkg{lavaan}.

The estimation of cognitive diagnostic models is possible
using the \pkg{CDM} package.

For the multidimensional latent class IRT model see the \pkg{MultiLCIRT}
package which also allows the estimation IRT models with polytomous item responses.

Latent class analysis can be carried out with \pkg{covLCA}, \pkg{poLCA} or 
\pkg{randomLCA} packages.

Markov Chain Monte Carlo estimation for item response models can also
be found in the \pkg{MCMCpack} package (see the \code{MCMCirt} functions
therein).

See Rusch, Mair and Hatzinger (2013) and Uenlue and Yanagida (2011)
for reviews of psychometrics packages in \R.
}


\examples{
##  
##   |-----------------------------------------------------------------|
##   | sirt 0.40-4 (2013-11-26)                                        |
##   | Supplementary Item Response Theory                              |
##   | Maintainer: Alexander Robitzsch <a.robitzsch at bifie.at >      |
##   | https://sites.google.com/site/alexanderrobitzsch/software       |
##   |-----------------------------------------------------------------|
##     
##                       _/              _/   
##              _/_/_/      _/  _/_/  _/_/_/_/
##           _/_/      _/  _/_/        _/     
##              _/_/  _/  _/          _/      
##         _/_/_/    _/  _/            _/_/   
##   
}

\keyword{package}
\section{\R Function Versions }{
\emph{attach.environment.sirt__0.01.R}, 
\emph{automatic.recode__1.04.R}, \cr
\emph{class.accuracy.rasch__0.03.R}, 
\emph{data.prep__1.03.R}, 
\emph{data.wide2long__0.12.R}, \cr
\emph{detect__0.04.R}, 
\emph{dif.logisticregression__1.01.R}, 
\emph{dif.variance__0.05.R}, \cr
\emph{dirichlet__0.04.R}, 
\emph{eigenvalues.manymatrices__0.03.R}, 
\emph{eigenvalues.sirt__0.05.R}, \cr
\emph{equating.rasch__0.02.R}, 
\emph{f1d.irt__1.03.R}, 
\emph{fit.adisop__2.09.R}, \cr
\emph{fit.gradedresponse__1.05.R}, 
\emph{fit.gradedresponse_alg__1.06.R}, 
\emph{fit.isop__2.05.R}, \cr
\emph{fit.logistic__2.03.R}, 
\emph{fit.logistic_alg__0.04.R}, 
\emph{fuzcluster__0.07.R}, \cr
\emph{fuzcluster_alg__0.10.R}, 
\emph{fuzdiscr__0.02.R}, 
\emph{fuzirt__0.01.R}, 
\emph{fuzlc__0.01.R}, \cr
\emph{gom.em.alg__5.04.R}, 
\emph{gom.em__5.06.R}, 
\emph{gom.jml__0.08.R}, 
\emph{gom.jml_alg__0.04.R}, \cr
\emph{greenyang.reliability__1.03.R}, 
\emph{invariance.alignment.aux__1.03.R}, \cr
\emph{invariance.alignment__2.22.R}, 
\emph{isop.dich__3.08.R}, 
\emph{isop.poly__2.06.R}, \cr
\emph{isop.scoring__1.03.R}, 
\emph{isop.test__0.05.R}, 
\emph{isop_alg__0.03.R}, \cr
\emph{latent.regression.em.normal__2.02.R}, 
\emph{latent.regression.em.raschtype__2.46.R}, \cr
\emph{lc.2raters.aux__0.02.R}, 
\emph{lc.2raters__0.11.R}, 
\emph{linking.haberman.aux__0.02.R}, \cr
\emph{linking.haberman__2.11.R}, 
\emph{linking.robust__1.06.R}, 
\emph{lllsa__0.01.R}, 
\emph{lsdm__1.07.R}, \cr
\emph{lsdm_aux__0.01.R}, 
\emph{marginal.truescore.reliability__0.01.R}, 
\emph{matrix_functions__0.02.R}, \cr
\emph{matrixfunctions_sirt__0.06.R}, 
\emph{mcmc.2pno.ml__3.08.R}, 
\emph{mcmc.2pno.ml_alg__3.12.R}, \cr
\emph{mcmc.2pno.ml_output__1.02.R}, 
\emph{mcmc.2pno__1.14.R}, 
\emph{mcmc.2pno_alg__1.08.R}, \cr
\emph{mcmc.2pnoh__1.02.R}, 
\emph{mcmc.2pnoh_alg__0.05.R}, 
\emph{mcmc.3pno.testlet__4.05.R}, \cr
\emph{mcmc.3pno.testlet_alg__2.11.R}, 
\emph{mcmc.3pno.testlet_output__1.06.R}, 
\emph{mcmc.aux__0.02.R}, \cr
\emph{mcmc.list.descriptives__0.04.R}, 
\emph{mcmclist2coda__0.01.R}, 
\emph{md.pattern.sirt__0.03.R}, \cr
\emph{mle.pcm.group__0.04.R}, 
\emph{modelfit.cor__2.20.R}, 
\emph{monoreg.rowwise__0.03.R}, \cr
\emph{noharm.sirt.est.aux__4.04.R}, 
\emph{noharm.sirt.preprocess__0.10.R}, 
\emph{noharm.sirt__0.35.R}, \cr
\emph{normal2.cw__0.04.R}, 
\emph{np.dich__0.09.R}, 
\emph{nr.numdiff__0.01.R}, 
\emph{pbivnorm2__1.06.R}, \cr
\emph{personfit.stat__0.01.R}, 
\emph{personfit__0.04.R}, 
\emph{pgenlogis__1.01.R}, \cr
\emph{plausible.values.raschtype__2.08.R}, 
\emph{plot.invariance.alignment__0.01.R}, \cr
\emph{plot.isop__1.05.R}, 
\emph{plot.mcmc.sirt__0.07.R}, 
\emph{plot.rasch.mml__0.06.R}, 
\emph{polychoric2__0.04.R}, \cr
\emph{prmse.subscores__0.02.R}, 
\emph{prob.guttman__1.02.R}, 
\emph{qmc.nodes__0.02.R}, \cr
\emph{R2noharm-utility__1.04.R}, 
\emph{R2noharm.EAP__0.15.R}, 
\emph{R2noharm.jackknife__1.03.R}, \cr
\emph{R2noharm__2.14.R}, 
\emph{rasch.copula__0.993.R}, 
\emph{rasch.copula2__6.12.R}, \cr
\emph{rasch.copula2_aux__1.05.R}, 
\emph{rasch.copula3.covariance__0.04.R}, 
\emph{rasch.copula3__6.33.R}, \cr
\emph{rasch.copula3_aux__6.09.R}, 
\emph{rasch.evm.pcm__1.09.R}, 
\emph{rasch.evm.pcm_aux__0.02.R}, \cr
\emph{rasch.jml.biascorr__0.03.R}, 
\emph{rasch.jml__3.12.R}, 
\emph{rasch.mirtlc__91.23.R}, \cr
\emph{rasch.mirtlc_aux__91.12.R}, 
\emph{rasch.mml.npirt__2.03.R}, 
\emph{rasch.mml.ramsay__2.03.R}, \cr
\emph{rasch.mml.raschtype__2.37.R}, 
\emph{rasch.mml__2.02.R}, 
\emph{rasch.mml2.missing1__0.05.R}, \cr
\emph{rasch.mml2__6.817.R}, 
\emph{rasch.pairwise.itemcluster__0.02.R}, 
\emph{rasch.pairwise__0.14.R}, \cr
\emph{rasch.pml__2.14.R}, 
\emph{rasch.pml_aux__1.03.R}, 
\emph{rasch.pml2__4.10.R}, 
\emph{rasch.pml2_aux__3.17.R}, \cr
\emph{rasch.pml3__5.07.R}, 
\emph{rasch.pml3_aux__5.01.R}, 
\emph{rasch.prox__1.03.R}, 
\emph{rasch.va__0.01.R}, \cr
\emph{reliability.nonlinear.sem__1.04.R}, 
\emph{rm.facets__3.13.R}, 
\emph{rm.facets_alg__3.12.R}, \cr
\emph{rm.hrm__7.14.R}, 
\emph{rm.hrm_alg__7.06.R}, 
\emph{rm_proc__0.01.R}, 
\emph{sim.rasch.dep__0.06.R}, \cr
\emph{smirt__7.13.R}, 
\emph{smirt_alg_comp__1.05.R}, 
\emph{smirt_alg_noncomp__2.27.R}, \cr
\emph{smirt_alg_partcomp__0.05.R}, 
\emph{smirt_postproc__0.01.R}, 
\emph{smirt_preproc__1.04.R}, \cr
\emph{stratified.cronbach.alpha__0.02.R}, 
\emph{summary.fuzcluster__0.03.R}, \cr
\emph{summary.gom.em__0.05.R}, 
\emph{summary.invariance.alignment__0.06.R}, 
\emph{summary.isop__0.04.R}, \cr
\emph{summary.latent.regression__0.01.R}, 
\emph{summary.mcmc.sirt__1.02.R}, \cr
\emph{summary.noharm.sirt__1.05.R}, 
\emph{summary.R2noharm.jackknife__1.01.R}, \cr
\emph{summary.R2noharm__0.04.R}, 
\emph{summary.rasch.copula__2.03.R}, \cr
\emph{summary.rasch.evm.pcm__0.04.R}, 
\emph{summary.rasch.mirtlc__7.03.R}, \cr
\emph{summary.rasch.mml2__1.02.R}, 
\emph{summary.rasch.pml__0.05.R}, 
\emph{summary.rm.facets__0.02.R}, \cr
\emph{summary.rm.hrm__0.02.R}, 
\emph{summary.smirt__0.06.R}, 
\emph{testlet.marginalized__0.04.R}, \cr
\emph{tetrachoric2__1.12.R}, 
\emph{truescore.irt__0.08.R}, 
\emph{unidim.csn__0.05.R}, 
\emph{wle.rasch__1.04.R}, \cr
\emph{yen.q3__1.01.R}, 
\emph{zzz__1.05.R}, 
}
\section{\emph{Rcpp} Function Versions }{
\emph{eigenvaluessirt__3.06.cpp}, 
\emph{evm_comp_matrix_poly__1.25.cpp}, 
\emph{evm_eigaux_fcts__4.07.h}, \cr
\emph{evm_eigenvals2__0.01.h}, 
\emph{first_eigenvalue_sirt__2.18.h}, 
\emph{gooijer_isop__4.02.cpp}, \cr
\emph{gooijercsntableaux__1.08.h}, 
\emph{matrixfunctions_sirt__1.07.cpp}, \cr
\emph{mle_pcm_group_c__1.02.cpp}, 
\emph{noharm_sirt_auxfunctions__2.09.cpp}, \cr
\emph{pbivnorm_rcpp_aux__0.29.h}, 
\emph{polychoric2_tetrachoric2_rcpp_aux__2.02.cpp}, \cr
\emph{probs_multcat_items_counts_csirt__2.03.cpp}, 
\emph{rm_smirt_mml2_code__4.05.cpp}, 
}
\section{\emph{Rd} Documentation Versions }{
\emph{automatic.recode__0.04.Rd}, 
\emph{ccov.np__0.08.Rd}, 
\emph{class.accuracy.rasch__0.09.Rd}, \cr
\emph{conf.detect__1.18.Rd}, 
\emph{data.activity.itempars__0.04.Rd}, 
\emph{data.big5__0.20.Rd}, \cr
\emph{data.ess2005__0.03.Rd}, 
\emph{data.inv4gr__0.03.Rd}, 
\emph{data.liking.science__0.03.Rd}, \cr
\emph{data.long__0.11.Rd}, 
\emph{data.math__0.07.Rd}, 
\emph{data.mcdonald__0.10.Rd}, 
\emph{data.mixed1__0.04.Rd}, \cr
\emph{data.ml__1.06.Rd}, 
\emph{data.noharm__2.06.Rd}, 
\emph{data.pars1.rasch__0.07.Rd}, \cr
\emph{data.pirlsmissing__0.03.Rd}, 
\emph{data.pisaMath__0.07.Rd}, 
\emph{data.pisaPars__0.05.Rd}, \cr
\emph{data.pisaRead__0.05.Rd}, 
\emph{data.ratings1__0.09.Rd}, 
\emph{data.raw1__0.03.Rd}, 
\emph{data.read__1.22.Rd}, \cr
\emph{data.reck__0.14.Rd}, 
\emph{data.si__0.15.Rd}, 
\emph{data.timss__0.05.Rd}, \cr
\emph{data.timss07.G8.RUS__0.02.Rd}, 
\emph{data.wide2long__0.09.Rd}, 
\emph{detect.index__0.09.Rd}, \cr
\emph{dif.logistic.regression__0.14.Rd}, 
\emph{dif.strata.variance__0.05.Rd}, \cr
\emph{dif.variance__0.06.Rd}, 
\emph{dirichlet.mle__0.06.Rd}, 
\emph{dirichlet.simul__0.03.Rd}, \cr
\emph{eigenvalues.manymatrices__0.06.Rd}, 
\emph{eigenvalues.sirt__0.03.Rd}, \cr
\emph{equating.rasch.jackknife__0.10.Rd}, 
\emph{equating.rasch__1.17.Rd}, 
\emph{expl.detect__1.06.Rd}, \cr
\emph{f1d.irt__1.12.Rd}, 
\emph{fit.isop__1.08.Rd}, 
\emph{fuzcluster__0.08.Rd}, 
\emph{fuzdiscr__0.08.Rd}, \cr
\emph{gom.em__1.37.Rd}, 
\emph{gom.jml__0.10.Rd}, 
\emph{greenyang.reliability__1.14.Rd}, \cr
\emph{invariance.alignment__1.12.Rd}, 
\emph{isop.scoring__1.11.Rd}, 
\emph{isop.test__0.08.Rd}, \cr
\emph{isop__3.13.Rd}, 
\emph{latent.regression.em.raschtype__1.25.Rd}, 
\emph{lc.2raters__0.07.Rd}, \cr
\emph{linking.haberman__0.25.Rd}, 
\emph{linking.robust__0.11.Rd}, 
\emph{lsdm__2.02.Rd}, \cr
\emph{marginal.truescore.reliability__0.11.Rd}, 
\emph{matrixfunctions.sirt__1.08.Rd}, \cr
\emph{mcmc.2pno.ml__0.24.Rd}, 
\emph{mcmc.2pno__1.18.Rd}, 
\emph{mcmc.2pnoh__0.13.Rd}, \cr
\emph{mcmc.3pno.testlet__1.11.Rd}, 
\emph{mcmc.list.descriptives__0.06.Rd}, 
\emph{mcmclist2coda__0.06.Rd}, \cr
\emph{md.pattern.sirt__0.08.Rd}, 
\emph{mle.pcm.group__0.11.Rd}, 
\emph{modelfit.sirt__0.34.Rd}, \cr
\emph{monoreg.rowwise__0.05.Rd}, 
\emph{noharm.sirt__0.12.Rd}, 
\emph{np.dich__0.10.Rd}, 
\emph{pbivnorm2__0.07.Rd}, \cr
\emph{person.parameter.rasch.copula__1.05.Rd}, 
\emph{personfit.stat__0.07.Rd}, 
\emph{pgenlogis__0.12.Rd}, \cr
\emph{plausible.value.imputation.raschtype__1.15.Rd}, 
\emph{plot.mcmc.sirt__0.04.Rd}, \cr
\emph{plot.np.dich__0.10.Rd}, 
\emph{polychoric2__0.03.Rd}, 
\emph{prmse.subscores.scales__0.12.Rd}, \cr
\emph{prob.guttman__1.08.Rd}, 
\emph{qmc.nodes__0.04.Rd}, 
\emph{R2noharm.EAP__0.08.Rd}, \cr
\emph{R2noharm.jackknife__1.07.Rd}, 
\emph{R2noharm__2.16.Rd}, 
\emph{rasch.copula__1.35.Rd}, \cr
\emph{rasch.evm.pcm__0.17.Rd}, 
\emph{rasch.jml.biascorr__0.11.Rd}, 
\emph{rasch.jml.jackknife1__2.07.Rd}, \cr
\emph{rasch.jml__1.19.Rd}, 
\emph{rasch.mirtlc__2.56.Rd}, 
\emph{rasch.mml__3.60.Rd}, \cr
\emph{rasch.pairwise.itemcluster__0.19.Rd}, 
\emph{rasch.pairwise__0.13.Rd}, 
\emph{rasch.pml__2.43.Rd}, \cr
\emph{rasch.prox__1.07.Rd}, 
\emph{rasch.va__0.04.Rd}, 
\emph{reliability.nonlinearSEM__0.10.Rd}, \cr
\emph{rm.facets__0.18.Rd}, 
\emph{rm.hrm__0.27.Rd}, 
\emph{sim.qm.ramsay__0.15.Rd}, 
\emph{sim.rasch.dep__0.16.Rd}, \cr
\emph{sim.raschtype__0.09.Rd}, 
\emph{sirt-package__1.76.Rd}, 
\emph{smirt__2.09.Rd}, \cr
\emph{stratified.cronbach.alpha__0.13.Rd}, 
\emph{summary.mcmc.sirt__0.05.Rd}, \cr
\emph{testlet.marginalized__0.07.Rd}, 
\emph{testlet.yen.q3__1.01.Rd}, 
\emph{tetrachoric2__1.22.Rd}, \cr
\emph{truescore.irt__0.09.Rd}, 
\emph{unidim.test.csn__1.05.Rd}, 
\emph{wle.rasch.jackknife__1.11.Rd}, \cr
\emph{wle.rasch__1.06.Rd}, 
\emph{yen.q3__1.04.Rd}, 
}
