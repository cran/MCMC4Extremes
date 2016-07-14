\name{gammap}
\alias{gammap}
\title{
Posterior Distribution with Gamma Density
}
\description{
MCMC runs of posterior distribution of data with \code{Gamma(alpha,beta)} density.
}
\usage{
gammap(data, int=1000)
}
\arguments{
  \item{data}{
data vector
}
  \item{int}{
number of iteractions selected in MCMC. The program selects 1 in each 10
iteractions, then \code{thin=10}. The first \code{thin*int/3} iteractions is used as burn-in. After
that, is runned \code{thin*int} iteraction, in which 1 of thin is selected for the final
MCMC chain, resulting the number of int iteractions.
}
}
\value{
An object of class \code{gammap} that gives a list containing the points of posterior distributions of \code{alpha} and \code{beta} of the gamma distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\note{
The non-informative prior distribution of these parameters are both
\code{Gamma(0.0001,0.0001)}. During the MCMC runs, screen shows the proportion of iteractions made
}
\examples{
# Vector of maxima return for each 10 days for ibovespa data
data(ibovespa)
ibmax=gev(ibovespa[,4],10)$data
# obtaining 500 points of posterior distribution 
ibovpost=gammap(ibmax,300)
}
