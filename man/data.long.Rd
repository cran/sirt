\name{data.long}
\alias{data.long}
\docType{data}
\title{
Longitudinal Dataset
}
\description{
This dataset contains 200 observations on
12 items. 6 items (\code{I1T1}, ... ,\code{I6T1})
were administered at measurement occasion T1
and 6 items at T2 (\code{I3T2}, ... , \code{I8T2}). There were 4 anchor items
which were presented at both time points.
The first column in the dataset contains the student identifier.
}
\usage{
data(data.long)
}

\format{
The format of the dataset is

\code{'data.frame':   200 obs. of  13 variables:} \cr
\code{ $ idstud: int  1001 1002 1003 1004 1005 1006 1007 1008 1009 1010 ...} \cr
\code{ $ I1T1  : int  1 1 1 1 1 1 1 0 1 1 ...} \cr
\code{ $ I2T1  : int  0 0 1 1 1 1 0 1 1 1 ...} \cr
\code{ $ I3T1  : int  1 0 1 1 0 1 0 0 0 0 ...} \cr
\code{ $ I4T1  : int  1 0 0 1 0 0 0 0 1 1 ...} \cr
\code{ $ I5T1  : int  1 0 0 1 0 0 0 0 1 0 ...} \cr
\code{ $ I6T1  : int  1 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ $ I3T2  : int  1 1 0 0 1 1 1 1 0 1 ...} \cr
\code{ $ I4T2  : int  1 1 0 0 1 1 0 0 0 1 ...} \cr
\code{ $ I5T2  : int  1 0 1 1 1 1 1 0 1 1 ...} \cr
\code{ $ I6T2  : int  1 1 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ I7T2  : int  1 0 0 0 0 0 0 0 0 1 ...} \cr
\code{ $ I8T2  : int  0 0 0 0 1 0 0 0 0 0 ...} \cr
}


\section{Example Index}{
\code{\link{linking.haberman}} (Example 2),
\code{\link{rasch.mml2}} (Example 2)
}


\examples{
\dontrun{
data(data.long)
dat <- data.long
dat <- dat[,-1]
I <- ncol(dat)

#****************
# Model 1: 2-dimensional Rasch model
#****************
# define Q-matrix
Q <- matrix(0,I,2)
Q[1:6,1] <- 1
Q[7:12,2] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("T1","T2")

# vector with same items
itemnr <- as.numeric( substring( colnames(dat) ,2,2) )
# fix mean at T2 to zero
mu.fixed <- cbind( 2,0 )

#--- M1a: rasch.mml2 (in sirt)
mod1a <- sirt::rasch.mml2(dat, Q=Q , est.b=itemnr , mu.fixed=mu.fixed)
summary(mod1a)

#--- M1b: smirt (in sirt)
mod1b <- sirt::smirt(dat, Qmatrix=Q , irtmodel="comp" , est.b= itemnr , 
                  mu.fixed=mu.fixed )

#--- M1c: tam.mml (in TAM)

# assume equal item difficulty of I3T1 and I3T2, I4T1 and I4T2, ...
# create draft design matrix and modify it
A <- TAM::designMatrices(resp=dat)$A
dimnames(A)[[1]] <- colnames(dat)
##   > str(A)
##    num [1:12, 1:2, 1:12] 0 0 0 0 0 0 0 0 0 0 ...
##    - attr(*, "dimnames")=List of 3
##     ..$ : chr [1:12] "Item01" "Item02" "Item03" "Item04" ...
##     ..$ : chr [1:2] "Category0" "Category1"
##     ..$ : chr [1:12] "I1T1" "I2T1" "I3T1" "I4T1" ...
A1 <- A[  ,  , c(1:6 , 11:12 ) ]
A1[7,2,3] <- -1     # difficulty(I3T1) = difficulty(I3T2)
A1[8,2,4] <- -1     # I4T1 = I4T2
A1[9,2,5] <- A1[10,2,6] <- -1
dimnames(A1)[[3]] <- substring( dimnames(A1)[[3]],1,2)
##   > A1[,2,]
##        I1 I2 I3 I4 I5 I6 I7 I8
##   I1T1 -1  0  0  0  0  0  0  0
##   I2T1  0 -1  0  0  0  0  0  0
##   I3T1  0  0 -1  0  0  0  0  0
##   I4T1  0  0  0 -1  0  0  0  0
##   I5T1  0  0  0  0 -1  0  0  0
##   I6T1  0  0  0  0  0 -1  0  0
##   I3T2  0  0 -1  0  0  0  0  0
##   I4T2  0  0  0 -1  0  0  0  0
##   I5T2  0  0  0  0 -1  0  0  0
##   I6T2  0  0  0  0  0 -1  0  0
##   I7T2  0  0  0  0  0  0 -1  0
##   I8T2  0  0  0  0  0  0  0 -1

# estimate model
# set intercept of second dimension (T2) to zero
beta.fixed <- cbind( 1 , 2 , 0 )
mod1c <- TAM::tam.mml( resp=dat , Q=Q , A=A1 , beta.fixed=beta.fixed)
summary(mod1c)

#***
# Model 2: 2-dimensional 2PL model

# set variance at T2 to 1
variance.fixed <- cbind( 2,2,1)

# M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat, Q=Q , est.b=itemnr , est.a=itemnr , mu.fixed=mu.fixed,
             variance.fixed=variance.fixed , mmliter=100)
summary(mod2a)
}
}
\keyword{datasets}
