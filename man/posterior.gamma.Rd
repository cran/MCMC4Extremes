\name{posterior.gamma}
\alias{posterior.gamma}

\title{
Posterior Distribution with Gamma Density
}
\description{
MCMC runs of posterior distribution of data with Gamma(alpha,beta) density.
}
\usage{
posterior.gamma(data, int)
}
\arguments{
  \item{data}{
data vector
}
  \item{int}{
number of iteractions selected in MCMC. The program selects 1 in each 10 iteractions, then thin=10. The first thin*int/3 iteractions
is used as burn-in. After that, is runned thin*int iteraction, in which 1 of thin is selected for the final MCMC chain, resulting the number of int
iteractions.
}
}

\value{
The program has as result a matrix with dimensions int x 2, where the first
represents points of posterior distribution of alpha parameter of Gamma,
and the second column is posterior distribution of beta parameter. The
non-informative prior distribution of these parameters are both
Gamma(0.0001,0.0001). During the MCMC runs, screen shows the proportion
of iteractions made.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgamma(300,3,5)
# Obtaning 1000 points of posterior distribution
ajuste=posterior.gamma(x,1000)
# Histogram of posterior distribution of the parameters,with 95% credibility intervals
par(mfrow=c(1,2))
hist(ajuste[,1],freq=FALSE,main="alpha")
abline(v=quantile(ajuste[,1],0.025),lwd=2)
abline(v=quantile(ajuste[,1],0.975),lwd=2)
hist(ajuste[,2],freq=FALSE,main="beta")
# Vector of maxima return for each 10 days for ibovespa data
data(ibovespa)
ibmax=gev(ibovespa[,4],10)$data
ibovpost=posterior.gamma(ibmax,1000)
}