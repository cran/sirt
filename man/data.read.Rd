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
\code{\link{noharm.sirt}} (Example 1),
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

# list of needed packages for the following examples
packages <- scan(what="character")
     eRm  ltm  TAM mRm  CDM  mirt psychotools  IsingFit  igraph  qgraph  pcalg 
     poLCA  randomLCA psychomix
     
# install packages     
if (FALSE){  # default is FALSE
   install.packages(packages)
   	  }
# load packages
for (pack in packages){ 
    library(pack, character.only=TRUE) 
        }

#*****************************************************
# Model 1: Rasch model
#*****************************************************

#-- M1a: rasch.mml2 (in sirt)
mod1a <- sirt::rasch.mml2(dat)
summary(mod1a)

#-- M1b: smirt (in sirt)
Qmatrix <- matrix(1,nrow=I , ncol=1)
mod1b <- sirt::smirt(dat,Qmatrix=Qmatrix)
summary(mod1b)

#-- M1c: gdm (in CDM)
theta.k <- seq(-6,6,len=21)
mod1c <- CDM::gdm(dat,theta.k=theta.k,irtmodel="1PL", skillspace="normal")
summary(mod1c)

#-- M1d: tam.mml (in TAM)
mod1d <- TAM::tam.mml( resp=dat )
summary(mod1d)

#-- M1e: RM (in eRm) 
mod1e <- eRm::RM( dat )
  # eRm uses Conditional Maximum Likelihood (CML) as the estimation method.
summary(mod1e)
eRm::plotPImap(mod1e)

#-- M1f: mrm (in mRm)
mod1f <- mRm::mrm( dat , cl=1)   # CML estimation
mod1f$beta  # item parameters

#-- M1g: mirt (in mirt)
mod1g <- mirt::mirt( dat , model=1 , itemtype="Rasch" )
summary(mod1g)
coef(mod1g)

#-- M1h: rasch (in ltm)
mod1h <- ltm::rasch( dat  , control=list(verbose=TRUE ) )
summary(mod1h)
coef(mod1h)

#-- M1i: RaschModel.fit (in psychotools)
mod1i <- psychotools::RaschModel.fit(dat)  # CML estimation
summary(mod1i)
plot(mod1i)

#-- M1j: noharm.sirt (in sirt)
Fpatt <- matrix( 0 , I , 1 )
Fval <- 1 + 0*Fpatt
Ppatt <- Pval <- matrix(1,1,1)
mod1j <- sirt::noharm.sirt( dat=dat , Ppatt=Ppatt,Fpatt=Fpatt , Fval=Fval , Pval=Pval )
summary(mod1j)
  #   Normal-ogive model, multiply item discriminations with constant D=1.7. 
  #   The same holds for other examples with noharm.sirt and R2noharm.
plot(mod1j) 

#-- M1k: rasch.pml3 (in sirt)
mod1k <- sirt::rasch.pml3( dat=dat)
  #         pairwise marginal maximum likelihood estimation
summary(mod1k)

#*****************************************************
# Model 2: 2PL model
#*****************************************************

#-- M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat , est.a=1:I)
summary(mod2a)

#-- M2b: smirt (in sirt)
mod2b <- sirt::smirt(dat,Qmatrix=Qmatrix,est.a="2PL")
summary(mod2b)

#-- M2c: gdm (in CDM)
mod2c <- CDM::gdm(dat,theta.k=theta.k,irtmodel="2PL", skillspace="normal")
summary(mod2c)

#-- M2d: tam.mml (in TAM)
mod2d <- TAM::tam.mml.2pl( resp=dat )
summary(mod2d)

#-- M2e: mirt (in mirt)
mod2e <- mirt::mirt( dat , model=1 , itemtype="2PL" )
summary(mod2e)
coef(mod2e)

#-- M2f: ltm (in ltm)
mod2f <- ltm::ltm( dat ~ z1 , control=list(verbose=TRUE ) )
summary(mod2f)
coef(mod2f)
plot(mod2f)

#-- M2g: R2noharm (in NOHARM, running from within R using sirt package)
  # define noharm.path where 'NoharmCL.exe' is located
noharm.path <- "c:/NOHARM"
  # covariance matrix
P.pattern <- matrix( 1 , ncol=1 , nrow=1 )
P.init <- P.pattern
P.init[1,1] <- 1
  # loading matrix
F.pattern <- matrix(1,I,1)
F.init <- F.pattern
  # estimate model
mod2g <- sirt::R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
             F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
             writename = "ex2g" , noharm.path = noharm.path  , dec ="," )
summary(mod2g)

#-- M2h: noharm.sirt (in sirt)
mod2h <- sirt::noharm.sirt( dat=dat , Ppatt=P.pattern,Fpatt=F.pattern , 
              Fval=F.init , Pval=P.init )
summary(mod2h)
plot(mod2h)

#-- M2i: rasch.pml2 (in sirt)
mod2i <- sirt::rasch.pml2(dat, est.a=1:I)
summary(mod2i)

#*****************************************************
# Model 3: 3PL model (note that results are very unstable!)
#*****************************************************

#-- M3a: rasch.mml2 (in sirt)
mod3a <- sirt::rasch.mml2(dat , est.a=1:I, est.c=1:I)
summary(mod3a)

#-- M3b: smirt (in sirt)
mod3b <- sirt::smirt(dat,Qmatrix=Qmatrix,est.a="2PL" , est.c=1:I)
summary(mod3b)

#-- M3c: mirt (in mirt)
mod3c <- mirt::mirt( dat , model=1 , itemtype="3PL" )
summary(mod3c)
coef(mod3c)

#-- M3d: ltm (in ltm)
mod3d <- ltm::tpm( dat , control=list(verbose=TRUE ) , max.guessing=.3)
summary(mod3d)
coef(mod3d) # => numerical instabilities

#*****************************************************
# Model 4: 3-dimensional Rasch model
#*****************************************************

# define Q-matrix
Q <- matrix( 0 , nrow=12 , ncol=3 )
Q[ cbind(1:12 , rep(1:3,each=4) ) ] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("A","B","C")

# define nodes
theta.k <- seq(-6,6,len=13 )

#-- M4a: smirt (in sirt)
mod4a <- sirt::smirt(dat,Qmatrix=Q,irtmodel="comp" , theta.k=theta.k , maxiter=30)
summary(mod4a)

#-- M4b: rasch.mml2 (in sirt)
mod4b <- sirt::rasch.mml2(dat,Q=Q,theta.k=theta.k , mmliter=30)
summary(mod4b)

#-- M4c: gdm (in CDM)
mod4c <- CDM::gdm( dat , irtmodel="1PL" , theta.k=theta.k , skillspace="normal" , 
            Qmatrix=Q , maxiter=30 , centered.latent=TRUE )
summary(mod4c)

#-- M4d: tam.mml (in TAM)
mod4d <- TAM::tam.mml( resp=dat , Q=Q , control=list(nodes=theta.k , maxiter=30) )
summary(mod4d)

#-- M4e: R2noharm (in NOHARM, running from within R using sirt package)
noharm.path <- "c:/NOHARM"
  # covariance matrix
P.pattern <- matrix( 1 , ncol=3 , nrow=3 )
P.init <- 0.8+0*P.pattern
diag(P.init) <- 1
  # loading matrix
F.pattern <- 0*Q
F.init <- Q
  # estimate model
mod4e <- sirt::R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
    F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
    writename = "ex4e" , noharm.path = noharm.path  , dec ="," )
summary(mod4e)

#-- M4f: mirt (in mirt)
cmodel <- mirt::mirt.model("
     F1 = 1-4
     F2 = 5-8 
     F3 = 9-12
     # equal item slopes correspond to the Rasch model
     CONSTRAIN = (1-4, a1), (5-8, a2) , (9-12,a3)
     COV = F1*F2, F1*F3 , F2*F3 
     " )
mod4f <- mirt::mirt(dat, cmodel , verbose=TRUE)
summary(mod4f)

#*****************************************************
# Model 5: 3-dimensional 2PL model
#*****************************************************

#-- M5a: smirt (in sirt)
mod5a <- sirt::smirt(dat,Qmatrix=Q,irtmodel="comp" , est.a="2PL" , theta.k=theta.k , 
                 maxiter=30)
summary(mod5a)

#-- M5b: rasch.mml2 (in sirt)
mod5b <- sirt::rasch.mml2(dat,Q=Q,theta.k=theta.k ,est.a=1:12, mmliter=30)
summary(mod5b)

#-- M5c: gdm (in CDM)
mod5c <- CDM::gdm( dat , irtmodel="2PL" , theta.k=theta.k , skillspace="loglinear" , 
            Qmatrix=Q , maxiter=30 , centered.latent=TRUE ,
            standardized.latent=TRUE)
summary(mod5c)

#-- M5d: tam.mml (in TAM)
mod5d <- TAM::tam.mml.2pl( resp=dat , Q=Q , control=list(nodes=theta.k , maxiter=30) )
summary(mod5d)

#-- M5e: R2noharm (in NOHARM, running from within R using sirt package)
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
mod5e <- sirt::R2noharm( dat = dat , model.type="CFA" , F.pattern = F.pattern , 
    F.init = F.init , P.pattern = P.pattern , P.init = P.init , 
    writename = "ex5e" , noharm.path = noharm.path  , dec ="," )
summary(mod5e)

#-- M5f: mirt (in mirt)
cmodel <- mirt::mirt.model("
   F1 = 1-4
   F2 = 5-8
   F3 = 9-12
   COV = F1*F2, F1*F3 , F2*F3 
   "  )
mod5f <- mirt::mirt(dat, cmodel , verbose=TRUE)
summary(mod5f)

#*****************************************************
# Model 6: Network models (Graphical models)
#*****************************************************

#-- M6a: Ising model using the IsingFit package (undirected graph)
#        - fit Ising model using the "OR rule" (AND=FALSE)
mod6a <- IsingFit::IsingFit(x=dat, family="binomial" , AND=FALSE)
summary(mod6a)
##           Network Density:                 0.29 
##    Gamma:                  0.25 
##    Rule used:              Or-rule 
# plot results
qgraph(mod6a$weiadj,fade = FALSE)

#**-- graph estimation using pcalg package

# some packages from Bioconductor must be downloaded at first (if not yet done)
if (FALSE){  # set 'if (TRUE)' if packages should be downloaded 
     source("http://bioconductor.org/biocLite.R")
     biocLite("RBGL")
     biocLite("Rgraphviz")
	}

#-- M6b: graph estimation based on Pearson correlations
V <- colnames(dat)
n <- nrow(dat)
mod6b <- pcalg::pc(suffStat = list(C = cor(dat), n = n ),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha=0.05, labels = V, verbose = TRUE)
plot(mod6b)
# plot in qgraph package
qgraph::qgraph(mod6b , label.color= rep( c( "red" , "blue","darkgreen" ) , each=4 ) ,
         edge.color="black")
summary(mod6b)

#-- M6c: graph estimation based on tetrachoric correlations
mod6c <- pcalg::pc(suffStat = list(C = tetrachoric2(dat)$rho, n = n ),
             indepTest = gaussCItest, alpha=0.05, labels = V, verbose = TRUE)
plot(mod6c)
summary(mod6c)

#-- M6d: Statistical implicative analysis (in sirt)
mod6d <- sirt::sia.sirt(dat , significance=.85 )
  # plot results with igraph and qgraph package
plot( mod6d$igraph.obj  , vertex.shape="rectangle" , vertex.size=30 )
qgraph::qgraph( mod6d$adj.matrix )

#*****************************************************
# Model 7: Latent class analysis with 3 classes
#*****************************************************

#-- M7a: randomLCA (in randomLCA)
  #        - use two trials of starting values
mod7a <- randomLCA::randomLCA(dat,nclass=3, notrials=2, verbose=TRUE)
summary(mod7a)
plot(mod7a,type="l" , xlab="Item")

#-- M7b: rasch.mirtlc (in sirt)
mod7b <- sirt::rasch.mirtlc( dat , Nclasses = 3 ,seed= -30 ,  nstarts=2 )   
summary(mod7b)
matplot( t(mod7b$pjk) , type="l" , xlab="Item" )

#-- M7c: poLCA (in poLCA)
  #   define formula for outcomes
f7c <- paste0( "cbind(" , paste0(colnames(dat),collapse=",") , ") ~ 1 " )
dat1 <- as.data.frame( dat + 1 ) # poLCA needs integer values from 1,2,..
mod7c <- poLCA::poLCA( as.formula(f7c),dat1,nclass=3 , verbose=TRUE)
plot(mod7c)

#-- M7d: gom.em (in sirt) 
  #    - the latent class model is a special grade of membership model
mod7d <- sirt::gom.em( dat , K=3 , problevels=c(0,1) ,  model="GOM"  )            
summary(mod7d)

#*****************************************************
# Model 8: Mixed Rasch model with two classes
#*****************************************************

#-- M8a: raschmix (in psychomix)
mod8a <- psychomix::raschmix(data= as.matrix(dat) , k = 2, scores = "saturated")
summary(mod8a)

#*****************************************************
# Model 9: Mixed 2PL model with two classes
#*****************************************************

#--- TO BE INCLUDED.

#*****************************************************
# Model 10: Rasch testlet model
#*****************************************************

#-- M10a: tam.fa (in TAM))
dims <- substring( colnames(dat),1,1 )  # define dimensions
mod10a <- TAM::tam.fa( resp=dat , irtmodel="bifactor1" , dims=dims ,
                control=list(maxiter=60) )
summary(mod10a)

#-- M10b: bfactor (in mirt)
cmodel <- mirt::mirt.model("
        G = 1-12 
        A = 1-4
        B = 5-8
        C = 9-12
        CONSTRAIN = (1-12,a1), (1-4, a2), (5-8, a3) , (9-12,a4)    
      ")
mod10b <- mirt::mirt(dat, model=cmodel , verbose=TRUE)
summary(mod10b)
coef(mod10b)
mod10b@logLik   # equivalent is slot( mod10b , "logLik")

#-- M10c: smirt (in sirt)
# define Q-matrix
Qmatrix <- matrix(0,12,4)
Qmatrix[,1] <- 1
Qmatrix[ cbind( 1:12 , match( dims , unique(dims)) +1 ) ]  <- 1
  # uncorrelated factors
variance.fixed <- cbind( c(1,1,1,2,2,3) , c(2,3,4,3,4,4) , 0 )
  # estimate model
mod10c <- sirt::smirt( dat , Qmatrix=Qmatrix , irtmodel="comp" ,
              variance.fixed=variance.fixed , qmcnodes=1000 , maxiter=60)
summary(mod10c)

#*****************************************************
# Model 11: Bifactor model
#*****************************************************

#-- M11a: tam.fa (in TAM)
dims <- substring( colnames(dat),1,1 )  # define dimensions
mod11a <- TAM::tam.fa( resp=dat , irtmodel="bifactor2" , dims=dims ,
                 control=list(maxiter=60) )
summary(mod11a)

#-- M11b: bfactor (in mirt)
dims1 <- match( dims , unique(dims) )
mod11b <- mirt::bfactor(dat, model=dims1 , verbose=TRUE)
summary(mod11b)
coef(mod11b)
mod11b@logLik

#-- M11c: smirt (in sirt)
  # define Q-matrix
Qmatrix <- matrix(0,12,4)
Qmatrix[,1] <- 1
Qmatrix[ cbind( 1:12 , match( dims , unique(dims)) +1 ) ]  <- 1
  # uncorrelated factors
variance.fixed <- cbind( c(1,1,1,2,2,3) , c(2,3,4,3,4,4) , 0 )
  # estimate model
mod11c <- sirt::smirt( dat , Qmatrix=Qmatrix , irtmodel="comp" , est.a="2PL" ,  
                variance.fixed=variance.fixed , qmcnodes=1000 , maxiter=60)
summary(mod11c)
}
}
\keyword{datasets}
