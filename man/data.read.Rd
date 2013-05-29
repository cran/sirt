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
data(data.read)
colMeans(data.read)

# estimation of the Rasch model with rasch.mml2
mod <- rasch.mml2( data.read )
summary(mod)
}
\keyword{datasets}
