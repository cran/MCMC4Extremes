\name{posterior.normal}
\alias{posterior.normal}
\title{
Posterior Distribution with Normal Density
}
\description{
MCMC runs of posterior distribution of data with Normal(mu,1/tau) density, where tau is the inverse of variance.
}
\usage{
posterior.normal(data, int)
}
\arguments{
  \item{data}{
data vector
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
The program has as result a matrix with dimensions int x 2, where the first
represents points of posterior distribution of mu parameter of Normal,
and the second column is posterior distribution of tau parameter. The
non-informative prior distribution of these parameters are Normal(0,10000000)
for the parameter mu and Gamma(0.001,0.001) for the parameter tau.
During the MCMC runs, screen shows the proportion
of iteractions made.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rnorm(300,2,sqrt(10))
# Obtaning 1000 points of posterior distribution
ajuste=posterior.normal(x,1000)
# Histogram of posterior distribution of the parameters,with 95% credibility intervals
par(mfrow=c(1,2))
hist(ajuste[,1],freq=FALSE,main="mu")
abline(v=quantile(ajuste[,1],0.025),lwd=2)
abline(v=quantile(ajuste[,1],0.975),lwd=2)
hist(ajuste[,2],freq=FALSE,main="tau")
abline(v=quantile(ajuste[,2],0.025),lwd=2)
abline(v=quantile(ajuste[,2],0.975),lwd=2)
# Posterior distribution of variance sigma^2
hist(1/ajuste[,2],freq=FALSE,main="sigma^2")
# Posterior distribution of river Nile dataset
data(Nile)
postnile=posterior.normal(Nile,1000)
}
