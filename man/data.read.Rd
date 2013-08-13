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
\examples{
\dontrun{
data(data.read)
dat <- data.read
I <- ncol(dat)

library(eRm); library(ltm); library(TAM); library(mRm)
library(CDM); library(mirt)

#***
# Model 1: Rasch model

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

#***
# Model 2: 2PL model

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

#***
# Model 3: 3PL model

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
}
}
\keyword{datasets}
