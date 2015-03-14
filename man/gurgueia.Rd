\name{gurgueia}
\alias{gurgueia}
\docType{data}
\title{
Daily river quota of Gurgueia River.
}
\description{
These data are the monthly maximum river quota of Gurgueia River, in Brazil, from 1975 to 2012.
}
\usage{data(gurgueia)}
\format{
  A data frame with 415 observations on the following 2 variables.
  \describe{
    \item{\code{date}}{month/year}
    \item{\code{maximum}}{a numeric vector with monthly maximum}
  }
}
\examples{
data(gurgueia)
hist(gurgueia[,2])
}
\keyword{datasets}
