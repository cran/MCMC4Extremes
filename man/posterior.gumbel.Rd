\name{posterior.gumbel}
\alias{posterior.gumbel}
\title{
Posterior Distribution with GEV, where xi=0
}
\description{
MCMC runs of posterior distribution of data with parameters of Generalized Extreme Value (GEV) density, in the particular case where xi=0 with parameters mu,sigma.
}
\usage{
posterior.gumbel(data, block, int)
}
\arguments{
  \item{data}{
data vector
}
  \item{block}{
the block size. A numeric value is interpreted as the number of data values in
each successive block. All the data is used, so the last block may not contain
block observations.
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
represents points of posterior points of mu parameter of GEV, and
the second column is posterior points of sigma parameter. The
non-informative prior distribution of these parameters are Normal(0,1000)
for the parameter mu and Gamma(0.001,0.001) for the parameter sigma. During the MCMC runs, screen shows the proportion of iteractions made.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgev(200,xi=0.0001,mu=10,sigma=5)
# Obtaning 600 points of posterior distribution
ajuste=posterior.gumbel(x,1,600)
# Histogram of posterior distribution of the parameters,with 95\% credibility intervals
par(mfrow=c(1,2))
hist(ajuste[,1],freq=FALSE,main="mu")
abline(v=quantile(ajuste[,1],0.025),lwd=2)
abline(v=quantile(ajuste[,1],0.975),lwd=2)
hist(ajuste[,2],freq=FALSE,main="sigma")
abline(v=quantile(ajuste[,2],0.025),lwd=2)
abline(v=quantile(ajuste[,2],0.975),lwd=2)
# Maxima of each month in river nidd data
data(nidd.annual)
out=posterior.gumbel(nidd.annual,1,500)
# Predictive distribution for 15 day maxima ibovespa returns
data(ibovespa)
postibv=posterior.gumbel(ibovespa[,4],15,500)
datai=gev(ibovespa[,4],15)$data
predictive.gumbel(postibv,datai)
}