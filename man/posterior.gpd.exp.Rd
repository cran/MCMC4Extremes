\name{posterior.gpd.exp}
\alias{posterior.gpd.exp}
\title{
Posterior Distribution with Parameters of GPD, when xi=0
}
\description{
MCMC runs of posterior distribution of data with parameters of
Generalized Pareto Distribution (GPD), with parameters sigma, when xi=0.
}
\usage{
posterior.gpd.exp(data, threshold, int)
}
\arguments{
  \item{data}{
data vector
}
  \item{threshold}{
a threshold value. When threshold=0, we have the particular case of exponential.
}
  \item{int}{
number of iteractions selected in MCMC. The program selects
1 in each 10 iteraction, then thin=10. The first thin*int/3 iteractions
is used as burn-in. After that, is runned thin*int iteraction, in which
1 of thin is selected for the final MCMC chain, resulting the number of int
iteractions.
}
}
\value{
The program has as result a vector with dimensions int, which represents points of posterior points of sigma parameter of GPD.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgpd(1000,xi=0.0001,mu=9,beta=2) # in this case beta is the scale parameter sigma
# Obtaning 1000 points of posterior distribution
ajuste=posterior.gpd.exp(x,9,1000)
# Histogram of posterior distribution of the parameters,with 95\% credibility intervals
hist(ajuste,freq=FALSE,main="mu")
abline(v=quantile(ajuste,0.025),lwd=2)
abline(v=quantile(ajuste,0.975),lwd=2)
# Absolute returns of Ibovespa
data(ibovespa)
absibv=abs(ibovespa[,4])
postibv=posterior.gpd.exp(absibv,0,1000)
}
