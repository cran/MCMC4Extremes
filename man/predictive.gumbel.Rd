\name{predictive.gumbel}
\alias{predictive.gumbel}
\title{ Predictive Density of Gumbel
}
\description{
Predictive density resulting of posterior distribution of Gumbel
parameters.
}
\usage{
predictive.gumbel(vector, data)
}
\arguments{
  \item{vector}{
a list object returned by posterior.gumbel
}
  \item{data}{
data vector
}
}
\value{
The program draws the predictive distribution of datased fitted using the posterior distribution of the Gumbel parameters.
}
\examples{
# Vector of maxima return for each 15 days for ibovespa data
data(ibovespa)
postibv=posterior.gumbel(ibovespa[,4],15,1000)
datai=gev(ibovespa[,4],15)$data
predictive.gumbel(postibv,datai)
}