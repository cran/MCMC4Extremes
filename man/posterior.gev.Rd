\name{posterior.gev}
\alias{posterior.gev}
\title{
Posterior Distribution with Parameters of GEV
}
\description{
MCMC runs of posterior distribution of data with parameters of
Generalized Extreme Value (GEV) density, with parameters mu,sigma and xi.
}
\usage{
posterior.gev(data, block, int)
}
\arguments{
  \item{data}{
data vector
}
  \item{block}{the block size. A numeric value is interpreted as the number of data values in
each successive block. All the data is used, so the last block may not contain
block observations.
}
  \item{int}{
Number of iteractions selected in MCMC. The program selects
1 in each 10 iteraction, then thin=10. The first thin*int/3 iteractions
is used as burn-in. After that, is runned thin*int iteraction, in which
1 of thin is selected for the final MCMC chain, resulting the number of int
iteractions.
}
}
\value{
The program has as result a matrix with dimensions int x 3, where the first
represents points of posterior points of mu parameter of GEV,
the second column is posterior points of sigma parameter, and the third
column is the posterior of xi . The
non-informative prior distribution of these parameters are Normal(0,1000)
for the parameter mu, Gamma(0.001,0.001) for the parameter sigma and N(0,100)
for parameter xi. During the MCMC runs, screen shows the proportion
of iteractions made.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgev(300,xi=0.4,mu=10,sigma=5)
# Obtaning 600 points of posterior distribution
ajuste=posterior.gev(x,1,600)
# Histogram of posterior distribution of the parameters,with 95% credibility intervals
par(mfrow=c(1,3))
hist(ajuste[,1],freq=FALSE,main="mu")
abline(v=quantile(ajuste[,1],0.025),lwd=2)
abline(v=quantile(ajuste[,1],0.975),lwd=2)
hist(ajuste[,2],freq=FALSE,main="sigma")
abline(v=quantile(ajuste[,2],0.025),lwd=2)
abline(v=quantile(ajuste[,2],0.975),lwd=2)
hist(ajuste[,3],freq=FALSE,main="xi")
abline(v=quantile(ajuste[,3],0.025),lwd=2)
abline(v=quantile(ajuste[,3],0.975),lwd=2)
# Return level of river nidd data
data(nidd.annual)
out=posterior.gev(nidd.annual,1,300)
# Vector of maxima return for each 15 days for ibovespa data
data(ibovespa)
postibv=posterior.gev(ibovespa[,4],15,300)
plot.ts(postibv)
}
