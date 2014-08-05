\name{sim.qm.ramsay}
\alias{sim.qm.ramsay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate from Ramsay's Quotient Model
}
\description{
This function simulates dichotomous item response data
according to Ramsay's quotient model (Ramsay, 1989).
}
\usage{
sim.qm.ramsay(theta, b, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
Vector of of length \eqn{N} person parameters (must be positive!)
}
  \item{b}{
Vector of length \eqn{I} of item difficulties (must be positive)
}
  \item{K}{
Vector of length \eqn{I} of guessing parameters (must be positive)
}
}
\details{
Ramsay's quotient model (Ramsay, 1989) is defined by the equation 
\deqn{P(X_{pi} = 1 | \theta_p ) = \frac{ \exp { ( \theta_p / b_i ) } }
	{ K_i + \exp { ( \theta_p / b_i ) } }}
}
\value{
An \eqn{N \times I}  data frame with dichotomous item responses.
}
\references{
Ramsay, J. O. (1989). A comparison of three simple test theory models. 
\emph{Psychometrika}, \bold{54}, 487-499.

van der Maas, H. J. L., Molenaar, D., Maris, G., Kievit, R. A., & 
Borsboom, D. (2011).
Cognitive psychology meets psychometric theory: On the relation between 
process models for decision making and latent variable models for 
individual differences.
\emph{Psychological Review}, \bold{318}, 339-356.
}
\author{
Alexander Robitzsch
}
%\section{Specific Objective Warning}{
%    If the same arguments like in the Rasch model with respect to the scale level 
%    would apply for Ramsay's quotient model, this model could provide
%    ratio scaled person abilities. This is prevalent in the model equation
%    because person abilities are assumed to be nonnegative, i.e. \eqn{\theta_p >0}.
%    }
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{rasch.mml2}} for estimating Ramsay's quotient model.

See \code{\link{sim.raschtype}} for simulating response data from
the generalized logistic item response model.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#***
set.seed(657)
# simulate data according to the Ramsay model
N <- 1000       # persons
I <- 11         # items
theta <- exp( rnorm( N ) )  # person ability
b <- exp( seq(-2,2,len=I))  # item difficulty
K <- rep( 3 , I )           # K parameter (=> guessing)

# apply simulation function
dat <- sim.qm.ramsay( theta , b , K )

#***
# analysis
mmliter <- 50       # maximum number of iterations
I <- ncol(dat)
fixed.K <- rep( 3 , I ) 

# Ramsay QM with fixed K parameter (K=3 in fixed.K specification)
mod1 <- rasch.mml2( dat , mmliter = mmliter , irtmodel = "ramsay.qm", 
           numdiff.parm = .0001 ,  fixed.K = fixed.K )                  
summary(mod1)

# Ramsay QM with joint estimated K parameters
mod2 <- rasch.mml2( dat , mmliter = mmliter , irtmodel = "ramsay.qm" , 
               numdiff.parm = .0001 ,  est.K = rep(1,I)  )                  
summary(mod2)

\dontrun{
# Ramsay QM with itemwise estimated K parameters
mod3 <- rasch.mml2( dat , mmliter = mmliter , irtmodel = "ramsay.qm" , 
             numdiff.parm = .0001 ,   est.K = 1:I  )                  
summary(mod3)

# Rasch model
mod4 <- rasch.mml2( dat )
summary(mod4)

# generalized logistic model
mod5 <- rasch.mml2( dat , est.alpha = TRUE )
summary(mod5) 

# 2PL model
mod6 <- rasch.mml2( dat , est.a = rep(1,I) )
summary(mod6) 

# Difficulty + Guessing (b+c) Model
mod7 <- rasch.mml2( dat , est.c = rep(1,I) )
summary(mod7) 

# estimate separate guessing (c) parameters
mod8 <- rasch.mml2( dat , est.c = 1:I  )
summary(mod8) 
} 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simulating IRT models}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
