\name{data.reck}

\alias{data.reck}
\alias{data.reck21}
\alias{data.reck61DAT1}
\alias{data.reck61DAT2}

\docType{data}


\title{
Datasets from Reckase' Book \emph{Multidimensional Item Response Theory}
}
\description{
Some simulated datasets from Reckase (2009).
}
\usage{
data(data.reck21)
data(data.reck61DAT1)
data(data.reck61DAT2)
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

\item The format of the datasets \code{data.reck61DAT1} and
\code{data.reck61DAT2} (Table 6.1, p. 153) is

\code{List of 4} \cr
\code{ $ data : num [1:2500, 1:30] 1 0 0 1 1 0 0 1 1 0 ...} \cr
\code{  ..- attr(*, "dimnames")=List of 2} \cr
\code{  .. ..$ : NULL} \cr
\code{  .. ..$ : chr [1:30] "A01" "A02" "A03" "A04" ...} \cr
\code{ $ pars :'data.frame':} \cr
\code{  ..$ a1: num [1:30] 0.747 0.46 0.861 1.014 0.552 ...} \cr
\code{  ..$ a2: num [1:30] 0.025 0.0097 0.0067 0.008 0.0204 0.0064 0.0861 ...} \cr
\code{  ..$ a3: num [1:30] 0.1428 0.0692 0.404 0.047 0.1482 ...} \cr
\code{  ..$ d : num [1:30] 0.183 -0.192 -0.466 -0.434 -0.443 ...} \cr
\code{ $ mu   : num [1:3] -0.4 -0.7 0.1} \cr
\code{ $ sigma: num [1:3, 1:3] 1.21 0.297 1.232 0.297 0.81 ...} \cr



The dataset \code{data.reck61DAT2} has correlated dimensions while
\code{data.reck61DAT1} has uncorrelated dimensions.

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
            
#############################################################################
# EXAMPLE 2: data.reck61 dataset, Table 6.1, p. 153	
#############################################################################            

data(data.reck61DAT1)
dat <- data.reck61DAT1$data

#***
# Model 1: Exploratory factor analysis

#-- Model 1a: tam.fa in TAM
library(TAM)
mod1a <- tam.fa( dat , irtmodel="efa" , nfactors=3 )
# varimax rotation
varimax(mod1a$B.stand)

# Model 1b: EFA in NOHARM (Promax rotation)
mod1b <- R2noharm( dat = dat , model.type="EFA" ,  dimensions = 3 , 
           writename = "reck61__3dim_efa", noharm.path = "c:/NOHARM" ,dec = ",")
summary(mod1b)

#***
# Model 2: Confirmatory factor analysis

#-- Model 2a: tam.fa in TAM
dims <- c( rep(1,10) , rep(3,10) , rep(2,10)  )
Qmatrix <- matrix( 0 , nrow=30 , ncol=3 )
Qmatrix[ cbind( 1:30 , dims) ] <- 1
mod2a <- tam.mml.2pl( dat ,Q=Qmatrix , 
    control=list( snodes=1000, QMC=TRUE , maxiter=200) )
summary(mod2a)

#-- Model 2b: smirt in sirt
mod2b <- smirt( dat ,Qmatrix =Qmatrix , est.a="2PL" , maxiter=20 , qmcnodes=1000 )
summary(mod2b)

#-- Model 2c: rasch.mml2 in sirt
mod2c <- rasch.mml2( dat ,Qmatrix =Qmatrix , est.a= 1:30 , 
           mmliter =200 , theta.k = seq(-5,5,len=11) )
summary(mod2c)

#-- Model 2d: mirt in mirt
cmodel <- mirt.model(
    'F1 = 1-10
     F2 = 21-30
     F3 = 11-20
     COV = F1*F2, F1*F3 , F2*F3 ' )
mod2d <- mirt(dat, cmodel , verbose=TRUE)
summary(mod2d)
coef(mod2d)

#-- Model 2e: CFA in NOHARM
# specify covariance pattern
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- .4*P.pattern
diag(P.pattern) <- 0
diag(P.init) <- 1
# fix all entries in the loading matrix to 1
F.pattern <- matrix( 0 , nrow=30 , ncol=3 )
F.pattern[1:10,1] <- 1
F.pattern[21:30,2] <- 1
F.pattern[11:20,3] <- 1
F.init <- F.pattern
# estimate model
mod2e <- R2noharm( dat = dat , model.type="CFA" , P.pattern=P.pattern, 
            P.init=P.init , F.pattern=F.pattern, F.init=F.init ,
            writename = "reck61__3dim_cfa", noharm.path = "c:/NOHARM" ,dec = ",")
summary(mod2e)
}
}
\keyword{datasets}
