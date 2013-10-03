\name{data.reck}

\alias{data.reck}
\alias{data.reck21}
\docType{data}
\title{
Datasets from Reckase' Book \emph{Multidimensional Item Response Theory}
}
\description{
Some simulated datasets from Reckase (2009).
}
\usage{
data(data.reck21)
}

\format{
\itemize{
\item The format of the \code{data.reck21} (Table 2.1, p. 45) is:

\code{List of 2} \cr
\code{ $ data: num [1:2500, 1:50] 0 0 0 1 1 0 0 0 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:50] "I0001" "I0002" "I0003" "I0004" ...} \cr
\code{ $ pars:'data.frame':} \cr
\code{  ..$ a: num [1:50] 1.83 1.38 1.47 1.53 0.88 0.82 1.02 1.19 1.15 0.18 ...} \cr
\code{  ..$ b: num [1:50] 0.91 0.81 0.06 -0.8 0.24 0.99 1.23 -0.47 2.78 -3.85 ...} \cr
\code{  ..$ c: num [1:50] 0 0 0 0.25 0.21 0.29 0.26 0.19 0 0.21 ...} \cr

}
  }
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Simulated datasets
}
\references{
Reckase, M. (2009). \emph{Multidimensional item response theory}. 
New York: Springer.
}

\examples{
\dontrun{	
#############################################################################
# EXAMPLE 1: data.reck21 dataset, Table 2.1, p. 45	
#############################################################################
data(data.reck21)

dat <- data.reck21$dat      # extract dataset

# items with zero guessing parameters
guess0 <- c( 1 , 2 , 3 , 9 ,11 ,27 ,30 ,35 ,45 ,49 ,50 )
I <- ncol(dat)

#***
# Model 1: 3PL estimation using rasch.mml2
est.c <- est.a <- 1:I
est.c[ guess0 ] <- 0
mod1 <- rasch.mml2( dat , est.a=est.a , est.c=est.c , mmliter= 300 )
summary(mod1)

#***
# Model 2: 3PL estimation using smirt
Q <- matrix(1,I,1)
mod2 <- smirt( dat , Qmatrix=Q , est.a= "2PL"  , est.c=est.c , increment.factor=1.01)
summary(mod2)

#***
# Model 3: estimation in mirt package
library(mirt)
itemtype <- rep("3PL" , I )
itemtype[ guess0 ] <- "2PL"
mod3 <- mirt(dat, 1, itemtype = itemtype , verbose=TRUE)
summary(mod3)

c3 <- unlist( coef(mod3) )[ 1:(4*I) ]
c3 <- matrix( c3 , I , 4 , byrow=TRUE )
# compare estimates of rasch.mml2, smirt and true parameters
round( cbind( mod1$item$c , mod2$item$c ,c3[,3]  ,data.reck21$pars$c ) , 2 )
round( cbind( mod1$item$a , mod2$item$a.Dim1 ,c3[,1], data.reck21$pars$a ) , 2 )
round( cbind( mod1$item$b , mod2$item$b.Dim1 / mod2$item$a.Dim1 , - c3[,2] / c3[,1] , 
            data.reck21$pars$b ) , 2 )
}
}
\keyword{datasets}
