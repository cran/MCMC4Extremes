\name{predictive.gpd}
\alias{predictive.gpd}
\title{
Predictive Density of GPD
}
\description{
Predictive density resulting of posterior distribution of GPD parameters.
}
\usage{
predictive.gpd(vector, data, threshold)
}
\arguments{
  \item{vector}{
a list object returned by posterior.gpd
}
  \item{data}{
data vector
}
  \item{threshold}{
a threshold value
}
}
\value{
The program draws the predictive distribution of dataset fitted using
the posterior distribution of the GPD parameteres.
}
\examples{
# Fitting the predictive distribution on simulated data
x=rgpd(500,xi=0.3,mu=9,beta=2)
ajuste=posterior.gpd(x,9,800)
predictive.gpd(ajuste,x,9)
# Danish data
data(danish)
out=posterior.gpd(danish,10,500)
predictive.gpd(out,danish,10)
}
