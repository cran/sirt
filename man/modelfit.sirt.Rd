%% File Name: modelfit.sirt.Rd
%% File Version: 0.624

\name{modelfit.sirt}
\alias{modelfit.sirt}
\alias{modelfit.cor.poly}
\alias{IRT.modelfit.sirt}

\title{
Assessing Model Fit and Local Dependence by Comparing Observed and Expected
Item Pair Correlations
}

\description{
This function computes several measures of absolute model fit and local
dependence indices for dichotomous item responses which are
based on comparing observed and expected frequencies of item pairs
(Chen, de la Torre & Zhang, 2013; see \code{\link[CDM:modelfit.cor]{modelfit.cor}}
for more details).
}
\usage{
modelfit.sirt(object)

modelfit.cor.poly( data, probs, theta.k, f.qk.yi)

\method{IRT.modelfit}{sirt}(object, mod, ...)
}
\arguments{
  \item{object}{
An object generated by \code{\link{rasch.mml2}},
\code{\link{rasch.mirtlc}}, \code{\link{rasch.pml3}} (\code{\link{rasch.pml2}}),
\code{\link{smirt}}, \code{\link{R2noharm}}, \code{\link{noharm.sirt}},
\code{\link{gom.em}}, \code{\link[TAM:tam.mml]{TAM::tam.mml}},
\code{\link[TAM:tam.mml.2pl]{TAM::tam.mml.2pl}},
\code{\link[TAM:tam.fa]{TAM::tam.fa}},
  \code{\link[mirt:mirt]{mirt::mirt}}
}
\item{data}{Dataset with polytomous item responses}
\item{probs}{Item response probabilities at grid \code{theta.k}}
\item{theta.k}{Grid of theta vector}
\item{f.qk.yi}{Individual posterior}
\item{mod}{Model name}
\item{\ldots}{Further arguments to be passed}
}

\value{
A list with following entries:
  \item{modelfit}{Model fit statistics:

   \code{MADcor}: mean of absolute deviations in observed and expected correlations
        (DiBello et al., 2007)

   \code{SRMSR}: standardized mean square root of squared residuals
        (Maydeu-Olivares, 2013; Maydeu-Olivares & Joe, 2014)

    \code{MX2}: Mean of \eqn{\chi^2} statistics of all item pairs
            (Chen & Thissen, 1997)

    \code{MADRESIDCOV}: Mean of absolute deviations of residual
            covariances (McDonald & Mok, 1995)

    \code{MADQ3}: Mean of absolute values of \eqn{Q_3} statistic (Yen, 1984)

    \code{MADaQ3}: Mean of absolute values of centered \eqn{Q_3} statistic
            }

  \item{itempairs}{Fit of every item pair}
}
\references{
Chen, W., & Thissen, D. (1997). Local dependence indexes for item pairs
 using item response theory. \emph{Journal of Educational and Behavioral Statistics,
 22}, 265-289.

  DiBello, L. V., Roussos, L. A., & Stout, W. F. (2007) Review of
  cognitively diagnostic assessment and a summary of psychometric models.
  In C. R. Rao and S. Sinharay (Eds.), \emph{Handbook of Statistics},
  Vol. 26 (pp. 979--1030). Amsterdam: Elsevier.

Maydeu-Olivares, A. (2013). Goodness-of-fit assessment of item response
theory models (with discussion).
\emph{Measurement: Interdisciplinary Research and Perspectives,
11}, 71-137.

Maydeu-Olivares, A., & Joe, H. (2014). Assessing approximate fit in categorical
data analysis. \emph{Multivariate Behavioral Research, 49}, 305-328.

  McDonald, R. P., & Mok, M. M.-C. (1995). Goodness of fit in item response models.
  \emph{Multivariate Behavioral Research, 30}, 23-40.

Yen, W. M. (1984). Effects of local item dependence on the fit and equating
performance of the three-parameter logistic model.
\emph{Applied Psychological Measurement, 8}, 125-145.
}


%\author{
%Alexander Robitzsch
%}




\note{
The function \code{modelfit.cor.poly} is just a wrapper to
\code{\link[TAM:tam.modelfit]{TAM::tam.modelfit}} in the \pkg{TAM} package.
}

\seealso{
Supported classes: \code{\link{rasch.mml2}},
\code{\link{rasch.mirtlc}}, \code{\link{rasch.pml3}} (\code{\link{rasch.pml2}}),
\code{\link{smirt}}, \code{\link{R2noharm}}, \code{\link{noharm.sirt}},
\code{\link{gom.em}},
\code{\link[TAM:tam.mml]{TAM::tam.mml}},
\code{\link[TAM:tam.mml.2pl]{TAM::tam.mml.2pl}},
\code{\link[TAM:tam.fa]{TAM::tam.fa}},
  \code{\link[mirt:mirt]{mirt::mirt}}

For more details on fit statistics of this function
see \code{\link[CDM:modelfit.cor]{CDM::modelfit.cor}}.
 }

\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: Reading data
#############################################################################
data(data.read)
dat <- data.read
I <- ncol(dat)

#*** Model 1: Rasch model
mod1 <- sirt::rasch.mml2(dat)
fmod1 <- sirt::modelfit.sirt( mod1 )
summary(fmod1)

#*** Model 1b: Rasch model in TAM package
library(TAM)
mod1b <- TAM::tam.mml(dat)
fmod1b <- sirt::modelfit.sirt( mod1b )
summary(fmod1b)

#*** Model 2: Rasch model with smoothed distribution
mod2 <- sirt::rasch.mml2( dat, distribution.trait="smooth3" )
fmod2 <- sirt::modelfit.sirt( mod2 )
summary(fmod2 )

#*** Model 3: 2PL model
mod3 <- sirt::rasch.mml2( dat, distribution.trait="normal", est.a=1:I )
fmod3 <- sirt::modelfit.sirt( mod3 )
summary(fmod3 )

#*** Model 3: 2PL model in TAM package
mod3b <- TAM::tam.mml.2pl( dat )
fmod3b <- sirt::modelfit.sirt(mod3b)
summary(fmod3b)
# model fit in TAM package
tmod3b <- TAM::tam.modelfit(mod3b)
summary(tmod3b)
# model fit in mirt package
library(mirt)
mmod3b <- sirt::tam2mirt(mod3b)   # convert to mirt object
mirt::M2(mmod3b$mirt)         # global fit statistic
mirt::residuals( mmod3b$mirt, type="LD")  # local dependence statistics

#*** Model 4: 3PL model with equal guessing parameter
mod4 <- TAM::rasch.mml2( dat, distribution.trait="smooth3", est.a=1:I, est.c=rep(1,I) )
fmod4 <- sirt::modelfit.sirt( mod4 )
summary(fmod4 )

#*** Model 5: Latent class model with 2 classes
mod5 <- sirt::rasch.mirtlc( dat, Nclasses=2 )
fmod5 <- sirt::modelfit.sirt( mod5 )
summary(fmod5 )

#*** Model 6: Rasch latent class model with 3 classes
mod6 <- sirt::rasch.mirtlc( dat, Nclasses=3, modeltype="MLC1", mmliter=100)
fmod6 <- sirt::modelfit.sirt( mod6 )
summary(fmod6 )

#*** Model 7: PML estimation
mod7 <- sirt::rasch.pml3( dat )
fmod7 <- sirt::modelfit.sirt( mod7 )
summary(fmod7 )

#*** Model 8: PML estimation
#      Modelling error correlations:
#          joint residual correlations for each item cluster
error.corr <- diag(1,ncol(dat))
itemcluster <- rep( 1:4,each=3 )
for ( ii in 1:3){
    ind.ii <- which( itemcluster==ii )
    error.corr[ ind.ii, ind.ii ] <- ii
        }
mod8 <- sirt::rasch.pml3( dat, error.corr=error.corr )
fmod8 <- sirt::modelfit.sirt( mod8 )
summary(fmod8 )

#*** Model 9: 1PL in smirt
Qmatrix <- matrix( 1, nrow=I, ncol=1 )
mod9 <- sirt::smirt( dat, Qmatrix=Qmatrix )
fmod9 <- sirt::modelfit.sirt( mod9 )
summary(fmod9 )

#*** Model 10: 3-dimensional Rasch model in NOHARM
noharm.path <- "c:/NOHARM"
Q <- matrix( 0, nrow=12, ncol=3 )
Q[ cbind(1:12, rep(1:3,each=4) ) ] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("A","B","C")
# covariance matrix
P.pattern <- matrix( 1, ncol=3, nrow=3 )
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
# loading matrix
F.pattern <- 0*Q
F.init <- Q
# estimate model
mod10 <- sirt::R2noharm( dat=dat, model.type="CFA", F.pattern=F.pattern,
            F.init=F.init, P.pattern=P.pattern, P.init=P.init,
            writename="ex4e", noharm.path=noharm.path, dec="," )
fmod10 <- sirt::modelfit.sirt( mod10 )
summary(fmod10)

#*** Model 11: Rasch model in mirt package
library(mirt)
mod11 <- mirt::mirt(dat, 1, itemtype="Rasch",verbose=TRUE)
fmod11 <- sirt::modelfit.sirt( mod11 )
summary(fmod11)
# model fit in mirt package
mirt::M2(mod11)
mirt::residuals(mod11)
}
}
%% \keyword{Model fit}

