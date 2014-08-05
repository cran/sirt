\name{data.big5}
\alias{data.big5}
\alias{data.big5.qgraph}


\docType{data}
\title{
Dataset Big 5 from \pkg{qgraph} Package
}

\description{
This is a Big 5 dataset from the \pkg{qgraph} package (Dolen, Oorts,
Stoel, Wicherts, 2009). It contains 500 subjects on 240 items.
}

\usage{
data(data.big5)
data(data.big5.qgraph)
}

\format{
\itemize{
\item The format of \code{data.big5} is: \cr 
\code{ num [1:500, 1:240] 1 0 0 0 0 1 1 2 0 1 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:240] "N1" "E2" "O3" "A4" ...} \cr

\item The format of \code{data.big5.qgraph} is: \cr  

\code{ num [1:500, 1:240] 2 3 4 4 5 2 2 1 4 2 ...} \cr
\code{ - attr(*, "dimnames")=List of 2} \cr
\code{  ..$ : NULL} \cr
\code{  ..$ : chr [1:240] "N1" "E2" "O3" "A4" ...} \cr
}
}


\details{
In these datasets, there exist 48 items for each dimension. The Big 5
dimensions are Neuroticism (\code{N}), Extraversion (\code{E}), 
Openness (\code{O}), Agreeableness (\code{A}) and 
Conscientiousness (\code{C}). Note that the \code{data.big5} differs from
\code{data.big5.qgraph} in a way that original items were recoded into
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





\section{Example Index}{

%%%%%%% data.big5
\emph{Dataset} \code{data.big5} 	
	
\code{\link{conf.detect}} (Example 2)

%%%%%%% data.big5.qgraph
\emph{Dataset} \code{data.big5.qgraph} 	

--
%%\code{\link{mcmc.2pno.ml}} (Example 2)
}

\examples{
\dontrun{
data(data.big5)	
# extract first 10 openness items
items <- which( substring( colnames(data.big5) , 1 , 1 ) == "O"  )[1:10]
dat <- data.big5[ , items ]
I <- ncol(dat)
summary(dat)

# list of needed packages for the following examples
packages <- scan(what="character")
     sirt   TAM   eRm   CDM   mirt  ltm   mokken  psychotools  psychomix
     psych
     
# install packages     
if (FALSE){  # default is FALSE
   install.packages(packages)
   	  }
# load packages
for (pack in packages){ 
    library(pack, character.only=TRUE) 
        }

# descriptive statistics
psych::describe(dat)
        
#****************
# Model 1: Partial Credit Model
#****************

#-- M1a: rm.facets (in sirt)
m1a <- sirt::rm.facets( dat )
summary(m1a)

#-- M1b: tam.mml (in TAM)
m1b <- TAM::tam.mml( resp=dat )
summary(m1b)

#-- M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
m1c <- CDM::gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="normal")
summary(m1c)
# compare results with loglinear skillspace
m1c2 <- CDM::gdm( dat , irtmodel="1PL" ,theta.k=theta.k , skillspace="loglinear")
summary(m1c2)

#-- M1d: PCM (in eRm)
m1d <- eRm::PCM( dat )
summary(m1d)

#-- M1e: gpcm (in ltm)
m1e <- ltm::gpcm( dat , constraint = "1PL" , control=list(verbose=TRUE))
summary(m1e)

#-- M1f: mirt (in mirt)
m1f <- mirt::mirt( dat , model=1 , itemtype="1PL" , verbose=TRUE)
summary(m1f)
coef(m1f)

#-- M1g: PCModel.fit (in psychotools)
mod1g <- psychotools::PCModel.fit(dat)
summary(mod1g)
plot(mod1g)

#****************
# Model 2: Generalized Partial Credit Model
#****************

#-- M2a: rm.facets (in sirt)
m2a <- sirt::rm.facets( dat , est.a.item=TRUE)
summary(m2a)
# Note that in rm.facets the mean of item discriminations is fixed to 1

#-- M2b: tam.mml.2pl (in TAM)
m2b <- TAM::tam.mml.2pl( resp=dat , irtmodel="GPCM")
summary(m2b)

#-- M2c: gdm (in CDM)
m2c <- CDM::gdm( dat , irtmodel="2PL" ,theta.k=seq(-6,6,len=21) ,
                   skillspace="normal" , standardized.latent=TRUE)
summary(m2c)

#-- M2d: gpcm (in ltm)
m2d <- ltm::gpcm( dat , control=list(verbose=TRUE))
summary(m2d)

#-- M2e: mirt (in mirt)
m2e <- mirt::mirt( dat , model=1 ,  verbose=TRUE)
summary(m2e)
coef(m2e)

#****************
# Model 3: Nonparametric Item Response Model
#****************

#-- M3a: ISOP and ADISOP model - isop.poly (in sirt)
m3a <- sirt::isop.poly( dat )
summary(m3a)
plot(m3a)

#-- M3b: Mokken scale analysis (in mokken) 
# Scalability coefficients
mokken::coefH(dat)
# Assumption of monotonicity
monotonicity.list <- mokken::check.monotonicity(dat)
summary(monotonicity.list)
plot(monotonicity.list)
# Assumption of non-intersecting ISRFs using method restscore
restscore.list <- mokken::check.restscore(dat)
summary(restscore.list)
plot(restscore.list) 
}
}


\keyword{datasets}
