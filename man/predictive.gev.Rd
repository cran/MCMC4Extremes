\name{predictive.gev}
\alias{predictive.gev}
\title{
Predictive Density of GEV
}
\description{
Predictive density resulting of posterior distribution of GEV parameters.
}
\usage{
predictive.gev(vector, data)
}
\arguments{
  \item{vector}{
a list object returned by posterior.gev
}
  \item{data}{
data vector
}
}
\value{
The program draws the predictive distribution of dataset fitted using
the posterior distribution of the GEV parameteres.
}
\examples{
# Fitting the predictive distribution on simulated data
x=rgev(250,xi=0.1,mu=10,sigma=5)
ajuste=posterior.gev(x,1,500)
predictive.gev(ajuste,x)
# Predictive for montly maxima of River nidd
data(nidd.annual)
out=posterior.gev(nidd.annual,1,300)
predictive.gev(out,nidd.annual)
}
