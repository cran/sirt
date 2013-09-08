\name{linking.haberman}
\alias{linking.haberman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linking in the 2PL/Generalized Partial Credit Model
}
\description{
This function does the linking of serval studies which are calibrated
using the 2PL or the generalized item response model according to
Haberman (2009). This method is a generalization of log-mean-mean
linking from tow to several studies.
}
\usage{
linking.haberman(itempars, conv = 1e-05, maxiter = 1000, progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{itempars}{
A data frame with four or five column. The first four columns contain
in the order: study name, item name, \eqn{a} parameter, \eqn{b} parameter.
The fifth column is an optional weight for every item and every study.
}
\item{conv}{
Convergence criterion.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{progress}{
An optional logical indicating whether computational progress
should be displayed.
}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



\details{
For \eqn{t=1,\ldots,T} studies, item difficulties \eqn{b_{it}} and
item slopes \eqn{a_{it}} are available. For dichotomous responses, these
parameters are defined by the 2PL response equation
\deqn{ log P(X_{pi}=1| \theta_p ) = a_i ( \theta_p - b_i ) }
while for polytomous responses the generalized partial credit model holds
\deqn{ logit \frac{P(X_{pi}=k| \theta_p )}{P(X_{pi}=k-1| \theta_p )}
= a_i ( \theta_p - b_i + d_{ik} ) }

The parameters \eqn{ \{ a_{it} , b_{it} \}} of all items and studies are
linearly transformed using equations \eqn{a_{it} \approx A_t \cdot a_i} and
\eqn{b_{it} / A_t \approx B_t + b_i}. For identification reasons
we define \eqn{A_1=1} and \eqn{B_1} = 0.

The optimization function (which is a least squares criterion;
see Haberman 2009) seeks the transformation parameters \eqn{A_t} and 
\eqn{B_t} with an alternating least squares 
method. Note that every item \eqn{i} and every study \eqn{t} can
be weighted (specified in the fifth column of \code{itempars}).
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\value{
A list with following entries
  \item{transf.itempars}{Data frame with linear transformation functions
  for item parameters}
  \item{transf.personpars}{Data frame with linear transformation functions
  for person parameters}  
  \item{joint.itempars}{Estimated joint item parameters \eqn{a_i} and \eqn{b_i}
  	}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\references{
Haberman, S. J. (2009). \emph{Linking parameter estimates derived
from an item respone model through separate calibrations}. 
ETS Research Report ETS RR-09-40. Princeton, ETS.
}

\author{
Alexander Robitzsch
}


%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \pkg{plink} package for a wide diversity of linking methods.

Mean-mean linking, Stocking-Lord and Haebara linking in the generalized
logistic item response model can be conducted with 
\code{\link{equating.rasch}}.
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
#############################################################################
# EXAMPLE 1: Item parameters data.pars1.rasch and data.pars1.2pl
#############################################################################

# Model 1: Linking three studies calibrated by the Rasch model
data(data.pars1.rasch)
mod1 <- linking.haberman( itempars=data.pars1.rasch )

# Model 1b: Linking these studies but weigh these studies by
#     proportion weights 3 : 0.5 : 1  (see below).
#     All weights are the same for each item but they could also
#     be item specific.
itempars <- data.pars1.rasch
itempars$wgt <- 1
itempars[ itempars$study == "study1","wgt"] <- 3
itempars[ itempars$study == "study2","wgt"] <- .5
mod1b <- linking.haberman( itempars=itempars )

# Model 2: Linking three studies calibrated by the 2PL model
data(data.pars1.2pl)
mod2 <- linking.haberman( itempars=data.pars1.2pl )

\dontrun{
#############################################################################
# EXAMPLE 2: Linking longitudinal data
#############################################################################
data(data.long)

#******
# Model 1: Scaling with the 1PL model

# scaling at T1
dat1 <- data.long[ , grep("T1" , colnames(data.long) ) ]
resT1 <- rasch.mml2( dat1 )
itempartable1 <- data.frame( "study"="T1" , resT1$item[ , c("item" , "a" , "b" ) ] )
# scaling at T2
dat2 <- data.long[ , grep("T2" , colnames(data.long) ) ]
resT2 <- rasch.mml2( dat2 )
summary(resT2)
itempartable2 <- data.frame( "study"="T2" , resT2$item[ , c("item" , "a" , "b" ) ] )
itempartable <- rbind( itempartable1 , itempartable2 )
itempartable[,2] <- substring( itempartable[,2] , 1, 2 )
# estimate linking parameters
mod1 <- linking.haberman( itempars= itempartable )

#******
# Model 2: Scaling with the 2PL model

# scaling at T1
dat1 <- data.long[ , grep("T1" , colnames(data.long) ) ]
resT1 <- rasch.mml2( dat1 , est.a=1:6)
itempartable1 <- data.frame( "study"="T1" , resT1$item[ , c("item" , "a" , "b" ) ] )

# scaling at T2
dat2 <- data.long[ , grep("T2" , colnames(data.long) ) ]
resT2 <- rasch.mml2( dat2 , est.a=1:6)
summary(resT2)
itempartable2 <- data.frame( "study"="T2" , resT2$item[ , c("item" , "a" , "b" ) ] )
itempartable <- rbind( itempartable1 , itempartable2 )
itempartable[,2] <- substring( itempartable[,2] , 1, 2 )
# estimate linking parameters
mod2 <- linking.haberman( itempars= itempartable )
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linking}
\keyword{Equating}
