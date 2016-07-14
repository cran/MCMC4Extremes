\name{ibovespa}
\alias{ibovespa}
\docType{data}
\title{
Daily returns of ibovespa
}
\description{
These data are the daily returns of ibovespa from 2000 to 2009.
}
\usage{data(ibovespa)}
\format{
  A data frame with 2369 observations on the following 4 variables.
  \describe{
    \item{\code{month}}{a numeric vector with month}
    \item{\code{day}}{a numeric vector with day}
    \item{\code{year}}{a numeric vector with year}
    \item{\code{returns}}{a numeric vector with returns}
  }
}
\examples{
data(ibovespa)
hist(ibovespa[,4], main=NULL)
}
\keyword{datasets}
