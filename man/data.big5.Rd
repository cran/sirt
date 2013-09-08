\name{data.big5}
\alias{data.big5}
\docType{data}
\title{
Dataset Big 5 from \pkg{qgraph} Package
}
\description{
This is a Big 5 dataset from the \pkg{qgraph} package (Dolen, Oorts,
Stoel, Wicherts, 2009). It contains 500 subjects on 240 items.
}
\usage{data(data.big5)}

\format{
The format is: \cr \cr
\code{ num [1:500, 1:240] 1 0 0 0 0 1 1 2 0 1 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:240] "N1" "E2" "O3" "A4" ...} \cr
}


\details{
In the dataset, there exist 48 items for each dimension. The Big 5
dimensions are Neuroticism (\code{N}), Extraversion (\code{E}), 
Openness (\code{O}), Agreeableness (\code{A}) and 
Conscientiousness (\code{C}). Note that the original data was recoded into
three categories 0,1 and 2.
}


\source{
See \code{big5} in \pkg{qgraph} package.
}


\references{
Dolan, C. V., Oort, F. J., Stoel, R. D., & Wicherts, J. M. (2009). 
Testing measurement invariance in the target rotates multigroup exploratory 
factor model. \emph{Structural Equation Modeling}, \bold{16}, 295-314.
}



\examples{
\dontrun{
data(data.big5)	
# extract first 10 openness items
items <- which( substring( colnames(data.big5) , 1 , 1 ) == "O"  )[1:10]
dat <- data.big5[ , items ]
I <- ncol(dat)
summary(dat)

# load needed packages
library(TAM); library(eRm); library(CDM); library(mirt);
library(ltm); library(mokken)

#****************
# Model 1: Partial Credit Model
#****************

# M1a: rm.facets (in sirt)
m1a <- rm.facets( dat )
summary(m1a)

# M1b: tam.mml (in TAM)
m1b <- tam.mml( resp=dat )
summary(m1b)

# M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
m1c <- gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="normal")
summary(m1c)
# compare results with loglinear skillspace
m1c2 <- gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="loglinear")
summary(m1c2)

# M1d: PCM (in eRm)
m1d <- PCM( dat )
summary(m1d)

# M1e: gpcm (in ltm)
m1e <- gpcm( dat , constraint = "1PL" , control=list(verbose=TRUE))
summary(m1e)

# M1f: mirt (in mirt)
m1f <- mirt( dat , model=1 , itemtype="1PL" , verbose=TRUE)
summary(m1f)
coef(m1f)

#****************
# Model 2: Generalized partial credit model
#****************

# M2a: rm.facets (in sirt)
m2a <- rm.facets( dat , est.a.item=TRUE)
summary(m2a)
# Note that in rm.facets the mean of item discriminations is fixed to 1

# M2b: tam.mml.2pl (in TAM)
m2b <- tam.mml.2pl( resp=dat , irtmodel="GPCM")
summary(m2b)

# M2c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
m2c <- gdm( dat , irtmodel="2PL" ,theta.k=theta.k ,
            skillspace="normal" , standardized.latent=TRUE)
summary(m2c)

# M2d: gpcm (in ltm)
m2d <- gpcm( dat , control=list(verbose=TRUE))
summary(m2d)

# M2e: mirt (in mirt)
m2e <- mirt( dat , model=1 ,  verbose=TRUE)
summary(m2e)
coef(m2e)

#****************
# Model 3: Nonparametric Item Response Model
#****************

# M3a: ISOP and ADISOP model - isop.poly (in sirt)
m3a <- isop.poly( dat )
summary(m3a)
plot(m3a)

# M3b: Mokken scale analysis (in mokken) 
# Scalability coefficients
coefH(dat)
# Assumption of monotonicity
monotonicity.list <- check.monotonicity(dat)
summary(monotonicity.list)
plot(monotonicity.list)
# Assumption of non-intersecting ISRFs using method restscore
restscore.list <- check.restscore(dat)
summary(restscore.list)
plot(restscore.list) }
}


\keyword{datasets}
