\name{data.ratings}
\alias{data.ratings}
\alias{data.ratings1}
\docType{data}
\title{
Rating Dataset 1
}
\description{
This data set contains ratings for 135 students from
7 raters on 5 items.
}
\usage{
data(data.ratings1)
}



\format{
	
\itemize{	
\item Dataset \code{data.ratings1}: \cr

Data frame with 274 observations containing 5 criteria (\code{k1}, ..., \code{k5}),
135 students and 7 raters.
  
\code{'data.frame':   274 obs. of  7 variables:} \cr
\code{ $ idstud: int  100020106 100020106 100070101 100070101 100100109  ...} \cr
\code{ $ rater : Factor w/ 16 levels "db01","db02",..: 3 15 5 10 2 1 5 4 1 5 ...} \cr
\code{ $ k1    : int  1 1 0 1 2 0 1 3 0 0 ...} \cr
\code{ $ k2    : int  1 1 1 1 1 0 0 3 0 0 ...} \cr
\code{ $ k3    : int  1 1 1 1 2 0 0 3 1 0 ...} \cr
\code{ $ k4    : int  1 1 1 2 1 0 0 2 0 1 ...} \cr
\code{ $ k5    : int  2 2 1 2 0 1 0 3 1 0 ...} \cr

Data from a 2009 Austrian survey of national educational
standards for 8th graders in German language writing.
Variables \code{k1} to \code{k5} denote several rating
criteria of writing compentency. \cr




}
}

%\source{
%xxx
%}





\section{Example Index}{

%%% data.ratings1
\emph{Dataset} \code{data.ratings1} 

\code{\link{rm.facets}} (Example 1, Example 2),
\code{\link{rm.hrm}} (Example 1),

%%% data.si02
%\emph{Dataset} \code{data.si02} 

%\code{\link{rasch.mirtlc}} (Example 3)

}

%\examples{
%# include an example of rating data here	
%}	

\keyword{datasets}
