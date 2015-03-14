\name{retlevel.gev}
\alias{retlevel.gev}
\title{
Expected Return Level in each Block Times of GEV
}
\description{
Calculates the expected return level in each t block times.
}
\usage{
retlevel.gev(vector, t)
}
\arguments{
  \item{vector}{
  a list object returned by posterior.gev
}
  \item{t}{
  specifies the particular return level to be estimated
}
}
\value{
The program returns the posterior mean of expected return level in each t
periods of time.
}
\examples{
# Calculate the return level of each 10 years in River nidd data# Predictive for River nidd data
data(nidd.annual)
out=posterior.gev(nidd.annual,1,1000)
r10=retlevel.gev(out,10)
# Expected return level for each t=20 maxima for ibovespa
data(ibovespa)
postibv=posterior.gev(ibovespa[,4],15,300)
retlevel.gev(postibv,20)
}
