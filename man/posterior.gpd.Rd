\name{posterior.gpd}
\alias{posterior.gpd}
\title{
Posterior Distribution with Parameters of GPD
}
\description{
MCMC runs of posterior distribution of data with parameters of Generalized Pareto Distribution (GPD), with parameters sigma and xi.
}
\usage{
posterior.gpd(data, threshold, int)
}
\arguments{
  \item{data}{
data vector
}
  \item{threshold}{
a threshold value
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
The program has as result a matrix with dimensions int x 2, where the first represents points of posterior points of sigma parameter of GPD and the second column is posterior points of xi parameter. The joint priordistribution for these parameters is the Jeffreys prior Given as Castellanos and Cabras (2007).
}
\references{
Castellanos, M. A. and Cabras, S. (2007). A default Bayesian procedure for the
generalized Pareto distribution, Journal of Statistical Planning and Inference, 137, 473-483.
}
\examples{
# Obtaining posterior distribution of a vector of simulated points 
x=rgpd(300,xi=0.3,mu=9,beta=2) # in this case beta is the scale parameter sigma
# Obtaning 1000 points of posterior distribution
ajuste=posterior.gpd(x,9,1000)
# Histogram of posterior distribution of the parameters,with 95\% credibility intervals
par(mfrow=c(1,2))
hist(ajuste[,1],freq=FALSE,main="mu")
abline(v=quantile(ajuste[,1],0.025),lwd=2)
abline(v=quantile(ajuste[,1],0.975),lwd=2)
hist(ajuste[,2],freq=FALSE,main="sigma")
abline(v=quantile(ajuste[,2],0.025),lwd=2)
abline(v=quantile(ajuste[,2],0.975),lwd=2)
# Danish data for evir package, modelling losses over 10
data(danish)
out=posterior.gpd(danish,10,300)
}
