\name{data.wide2long}
\alias{data.wide2long}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converting a Data Frame from Wide Format in a Long Format
}
\description{
Converts a data frame in wide format into long format.
}
\usage{
data.wide2long(dat, id = NULL, X = NULL, Q = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses and a person identifier if
\code{id != NULL}.
}
  \item{id}{
An optional string with the variable name of the person identifier.
}
  \item{X}{Data frame with person covariates for inclusion in the
  	  data frame of long format
}
  \item{Q}{Data frame with item predictors. Item labels must be included
  	  as a column named by \code{"item"}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Data frame in long format
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alexander Robitzsch
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#############################################################################
# EXAMPLE 1: data.pisaRead
#############################################################################

data(data.pisaRead) 
dat <- data.pisaRead$data 
Q <- data.pisaRead$item   # item predictors

# define items
items <- colnames(dat)[ substring( colnames(dat) , 1 , 1 ) == "R" ]
dat1 <- dat[ , c( "idstud" , items ) ]
# matrix with person predictors
X <- dat[ , c("idschool" , "hisei" , "female" , "migra") ]

# create dataset in long format
dat.long <- data.wide2long( dat=dat1 , id="idstud" , X=X , Q=Q )

#***
# Model 1: Rasch model
mod1 <- lmer( resp ~ 0 + ( 1 | idstud ) + as.factor(item) , data = dat.long ,
         family="binomial" , verbose=TRUE)
summary(mod1)

#***
# Model 2: Rasch model and inclusion of person predictors
mod2 <- lmer( resp ~ 0 + ( 1 | idstud ) + as.factor(item) + 
        female + hisei + migra  , data = dat.long ,
         family="binomial" , verbose=TRUE)
summary(mod2)

#***
# Model 3: LLTM
mod3 <- lmer(resp ~ (1|idstud) + as.factor(ItemFormat) + as.factor(TextType), 
            data = dat.long , family="binomial" , verbose=TRUE)
summary(mod3)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{R utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
