\name{gevp}
\alias{gevp}
\title{
Posterior Distribution with Parameters of GEV
}
\description{
MCMC runs of posterior distribution of data with parameters of Generalized Extreme Value (GEV)
density, with parameters \code{mu}, \code{sigma} and \code{xi}.
}
\usage{
gevp(data, block, int=1000)
}
\arguments{                                                                               X
  \item{data}{
data vector
}
  \item{block}{
the block size. A numeric value is interpreted as the number of data values in
each successive block. All the data is used, so the last block may not contain
block observations
}
  \item{int}{
Number of iteractions selected in MCMC. The program selects 1 in each 10
iteraction, then \code{thin=10}. The first \code{thin*int/3} iteractions is used as burn-in. After
that, is runned \code{thin*int} iteraction, in which 1 of thin is selected for the final
MCMC chain, resulting the number of int iteractions.
}
}
\value{
An object of class \code{gevp} that gives a list containing the points of posterior distributions of \code{mu}, \code{sigma} and \code{xi} of the gev distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\note{
 The non-informative prior distribution of these parameters
are \code{Normal(0,1000)} for the parameter \code{mu}, \code{Gamma(0.001,0.001)} for the parameter \code{sigma}
and \code{Normal(0,100)} for parameter \code{xi}. During the MCMC runs, screen shows the proportion of iteractions
made.
}
\seealso{
\code{\link{plot.gevp}}, \code{\link{summary.gevp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgev(300,xi=0.1,mu=10,sigma=5)

# Obtaning 600 points of posterior distribution
ajuste=gevp(x,1,200)

# Obtaining 300 points of ponterior distribution of river nidd data
## Not run: data(nidd.annual)
## Not run: out=gevp(nidd.annual,1,300)

# Vector of maxima return for each 15 days for ibovespa data
## Not run: data(ibovespa)
## Not run: postibv=gevp(ibovespa[,4],15,300)
## Not run: plot.ts(postibv$posterior)
}
