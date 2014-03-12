\name{data.read}
\alias{data.read}
\docType{data}
\title{
Dataset Reading
}
\description{
This dataset contains \eqn{N=328} students and \eqn{I=12} items measuring reading 
competence. All 12 items are arranged into 3 testlets (items with common 
text stimulus) labeled as 
A, B and C. The allocation of items to testlets is indicated by their 
variable names.
}
\usage{data(data.read)}
\format{
  A data frame with 328 persons on the following 12 variables.
  Rows correspond to persons and columns to items. The following items are 
  included in \code{data.read}:
  
  Testlet A: \code{A1}, \code{A2}, \code{A3}, \code{A4} 
  
  Testlet B: \code{B1}, \code{B2}, \code{B3}, \code{B4} 
  
  Testlet C: \code{C1}, \code{C2}, \code{C3}, \code{C4} 
  
}



\section{Example Index}{
	
\code{\link{class.accuracy.rasch}} (Example 1),
\code{\link{fit.isop}} (Example 1),
\code{\link{greenyang.reliability}} (Example 1),
\code{\link{isop.dich}} (Example 1),
\code{\link{isop.scoring}} (Example 1),
\code{\link{marginal.truescore.reliability}} (Example 3),
\code{\link{mcmc.2pno}} (Example 1),
\code{\link{mcmc.3pno.testlet}} (Example 1),
\code{\link{mcmclist2coda}} (Example 1),
\code{\link{mle.pcm.group}} (Example 2),
\code{\link{modelfit.sirt}} (Example 1),
\code{\link{np.dich}} (Example 1),
\code{\link{person.parameter.rasch.copula}} (Example 1),
\code{\link{personfit.stat}} (Example 1),
\code{\link{prmse.subscores.scales}} (Example 1),
\code{\link{prob.guttman}} (Example 1),
\code{\link{rasch.copula2}} (Example 1),
\code{\link{rasch.jml.biascorr}} (Example 1),
\code{\link{rasch.mirtlc}} (Example 1),
\code{\link{rasch.mml2}} (Example 1),
\code{\link{rasch.pairwise}} (Example 1),
\code{\link{rasch.pml2}} (Example 1),
\code{\link{rasch.prox}} (Example 1),
\code{\link{smirt}} (Example 2),
\code{\link{stratified.cronbach.alpha}} (Example 1),
\code{\link{testlet.yen.q3}} (Example 1),
\code{\link{tetrachoric2}} (Example 1),
\code{\link{truescore.irt}} (Example 2),
\code{\link{unidim.test.csn}} (Example 1),
\code{\link{wle.rasch}} (Example 1),
\code{\link{wle.rasch.jackknife}} (Example 1),
\code{\link{yen.q3}} (Example 1),
}

\examples{
\dontrun{
data(data.read)
dat <- data.read
I <- ncol(dat)

library(eRm); library(ltm); library(TAM); library(mRm)
library(CDM); library(mirt)

#****************
# Model 1: Rasch model
#****************

# M1a: rasch.mml2 (in sirt)
mod1a <- rasch.mml2(dat)
summary(mod1a)

# M1b: smirt (in sirt)
Qmatrix <- matrix(1,nrow=I , ncol=1)
mod1b <- smirt(dat,Qmatrix=Qmatrix)
summary(mod1b)

# M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
mod1c <- gdm(dat,theta.k=theta.k,irtmodel="1PL", skillspace="normal")
summary(mod1c)

# M1d: tam.mml (in TAM)
mod1d <- tam.mml( resp=dat )
summary(mod1d)

# M1e: RM (in eRm)
mod1e <- RM( dat )
summary(mod1e)

# M1f: mrm (in mRm)
mod1f <- mrm( dat , cl=1)
mod1f$beta  # item parameters

# M1g: mirt (in mirt)
mod1g <- mirt( dat , model=1 , itemtype="1PL" )
summary(mod1g)
coef(mod1g)

# M1h: ltm (in ltm)
mod1h <- ltm( dat ~ z1 , control=list(verbose=TRUE ) )
summary(mod1h)
coef(mod1h)

#****************
# Model 2: 2PL model
#****************

# M2a: rasch.mml2 (in sirt)
mod2a <- rasch.mml2(dat , est.a=1:I)
summary(mod2a)

# M2b: smirt (in sirt)
mod2b <- smirt(dat,Qmatrix=Qmatrix,est.a="2PL")
summary(mod2b)

# M2c: gdm (in CDM)
mod2c <- gdm(dat,theta.k=theta.k,irtmodel="2PL", skillspace="normal")
summary(mod2c)

# M2d: tam.mml (in TAM)
mod2d <- tam.mml.2pl( resp=dat )
summary(mod2d)

# M2e: mirt (in mirt)
mod2e <- mirt( dat , model=1 , itemtype="2PL" )
summary(mod2e)
coef(mod2e)

# M2f: ltm (in ltm)
mod2f <- ltm( dat ~ z1 , control=list(verbose=TRUE ) )
summary(mod2f)
coef(mod2f)
plot(mod2f)

# M2g: R2noharm (in NOHARM, running from within R using sirt package)
# define noharm.path where 'NoharmCL.exe' is located
noharm.path <- "c:/NOHARM"
# covariance matrix
P.pattern <- matrix( 1 , ncol=1 , nrow=1 )
P.init <- P.pattern
P.init[1,1] <- 0
# loading matrix
F.pattern <- matrix(1,I,1)
F.init <- F.pattern
# estimate model
mod2g <- R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
    F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
    writename = "ex2g" , noharm.path = noharm.path  , dec ="," )
summary(mod2g)

#****************
# Model 3: 3PL model
#****************

# M3a: rasch.mml2 (in sirt)
mod3a <- rasch.mml2(dat , est.a=1:I, est.c=1:I)
summary(mod3a)

# M3b: smirt (in sirt)
mod3b <- smirt(dat,Qmatrix=Qmatrix,est.a="2PL" , est.c=1:I)
summary(mod3b)

# M3c: mirt (in mirt)
mod3c <- mirt( dat , model=1 , itemtype="3PL" )
summary(mod3c)
coef(mod3c)

# M3d: ltm (in ltm)
mod3d <- tpm( dat , control=list(verbose=TRUE ) , max.guessing=.3)
summary(mod3d)
coef(mod3d) # => numerical instabilities

#****************
# Model 4: 3-dimensional Rasch model
#****************

# define Q-matrix
Q <- matrix( 0 , nrow=12 , ncol=3 )
Q[ cbind(1:12 , rep(1:3,each=4) ) ] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("A","B" , "C")

# define nodes
theta.k <- seq(-6,6,len=13 )

# M4a: smirt (in sirt)
mod4a <- smirt(dat,Qmatrix=Q,irtmodel="comp" , theta.k=theta.k , maxiter=30)
summary(mod4a)

# M4b: rasch.mml2 (in sirt)
mod4b <- rasch.mml2(dat,Q=Q,theta.k=theta.k , mmliter=30)
summary(mod4b)

# M4c: gdm (in CDM)
mod4c <- gdm( dat , irtmodel="1PL" , theta.k=theta.k , skillspace="normal" , 
            Qmatrix=Q , maxiter=30 , centered.latent=TRUE )
summary(mod4c)

# M4d: tam.mml (in TAM)
mod4d <- tam.mml( resp=dat , Q=Q , control=list(nodes=theta.k , maxiter=30) )
summary(mod4d)

# M4e: R2noharm (in NOHARM, running from within R using sirt package)
noharm.path <- "c:/NOHARM"
# covariance matrix
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
# loading matrix
F.pattern <- 0*Q
F.init <- Q
# estimate model
mod4e <- R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
    F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
    writename = "ex4e" , noharm.path = noharm.path  , dec ="," )
summary(mod4e)

# M4f: mirt (in mirt)
cmodel <- mirt.model(
    'F1 = 1-4
     F2 = 5-8 
     F3 = 9-12
     # equal item slopes correspond to the Rasch model
     CONSTRAIN = (1-4, a1), (5-8, a2) , (9-12,a3)
     COV = F1*F2, F1*F3 , F2*F3 ' )
mod4f <- mirt(dat, cmodel , verbose=TRUE)
summary(mod4f)

#****************
# Model 5: 3-dimensional 2PL model
#****************

# M5a: smirt (in sirt)
mod5a <- smirt(dat,Qmatrix=Q,irtmodel="comp" , est.a="2PL" , theta.k=theta.k , 
           maxiter=30)
summary(mod5a)

# M5b: rasch.mml2 (in sirt)
mod5b <- rasch.mml2(dat,Q=Q,theta.k=theta.k ,est.a=1:12, mmliter=30)
summary(mod5b)

# M5c: gdm (in CDM)
mod5c <- gdm( dat , irtmodel="2PL" , theta.k=theta.k , skillspace="loglinear" , 
            Qmatrix=Q , maxiter=30 , centered.latent=TRUE ,
            standardized.latent=TRUE)
summary(mod5c)

# M5d: tam.mml (in TAM)
mod5d <- tam.mml.2pl( resp=dat , Q=Q , control=list(nodes=theta.k , maxiter=30) )
summary(mod5d)

# M5e: R2noharm (in NOHARM, running from within R using sirt package)
noharm.path <- "c:/NOHARM"
# covariance matrix
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
diag(P.pattern) <- 0
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
# loading matrix
F.pattern <- Q
F.init <- Q
# estimate model
mod5e <- R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
    F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
    writename = "ex5e" , noharm.path = noharm.path  , dec ="," )
summary(mod5e)

# M5f: mirt (in mirt)
cmodel <- mirt.model('
   F1 = 1-4
   F2 = 5-8
   F3 = 9-12
   COV = F1*F2, F1*F3 , F2*F3 '  )
mod5f <- mirt(dat, cmodel , verbose=TRUE)
summary(mod5f)
}
}
\keyword{datasets}
