\name{gpdp}
\alias{gpdp}
\title{
Posterior Distribution with Parameters of GPD
}
\description{
MCMC runs of posterior distribution of data with parameters of Generalized Pareto Distribution
(GPD), with parameters \code{sigma} and \code{xi} .
}
\usage{
gpdp(data, threshold, int=1000)
}
\arguments{
  \item{data}{
data vector
}
  \item{threshold}{
a threshold value
}
  \item{int}{
number of iteractions selected in MCMC. The program selects 1 in each 10
iteraction, then \code{thin=10}. The first \code{thin*int/3} iteractions is used as burn-in. After
that, is runned \code{thin*int} iteraction, in which 1 of thin is selected for the final
MCMC chain, resulting the number of int iteractions.
}
}
\value{
An object of class \code{gpdp} that gives a list containing the points of posterior distributions of \code{sigma} and \code{xi} of the gpd distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\references{
Castellanos, M. A. and Cabras, S. (2007). A default Bayesian procedure for the generalized Pareto
distribution, \emph{Journal of Statistical Planning and Inference}, 137, 473-483.
}
\note{
The joint priordistribution for these parameters is the Jeffreys prior Given as Castellanos and Cabras
(2007).
}
\seealso{
\code{\link{plot.gpdp}}, \code{\link{summary.gpdp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points 
x=rgpd(300,xi=0.1,mu=9,beta=2) # in this case beta is the scale parameter sigma

# Obtaning 1000 points of posterior distribution
ajuste=gpdp(x,9, 200)

# Histogram of posterior distribution of the parameters,with 95% credibility intervals
# Danish data for evir package, modelling losses over 10
## Not run data(danish)
## Not run out=gpdp(danish,10,300)
}
