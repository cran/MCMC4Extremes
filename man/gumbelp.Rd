\name{gumbelp}
\alias{gumbelp}
\title{
Posterior Distribution with GEV, where xi=0
}
\description{
MCMC runs of posterior distribution of data with parameters of Generalized Extreme Value (GEV)
density, in the particular case where \code{xi=0} with parameters \code{mu}, \code{sigma}.
}
\usage{
gumbelp(data, block, int=1000)
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
number of iteractions selected in MCMC. The program selects 1 in each 10
iteraction, then \code{thin=10}. The first \code{thin*int/3} iteractions is used as burn-in. After
that, is runned \code{thin*int} iteraction, in which 1 of thin is selected for the final
MCMC chain, resulting the number of int iteractions
}
}
\value{
An object of class \code{gumbelp} that gives a list containing the points of posterior distributions of \code{mu} and \code{sigma} of the gev distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\note{
The non-informative prior distribution of these parameters are \code{Normal(0,1000)} for the
parameter \code{mu} and \code{Gamma(0.001,0.001)} for the parameter \code{sigma}. During the MCMC runs, screen
shows the proportion of iteractions made.
}
\seealso{
\code{\link{plot.gumbelp}}, \code{\link{summary.gumbelp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rgev(200,xi=0.0001,mu=10,sigma=5)
# Obtaning 600 points of posterior distribution
ajuste=gumbelp(x,1,600)

# Maxima of each month in river nidd data
## Not run: data(nidd.annual)
## Not run: out=gumbelp(nidd.annual,1,500)

# Predictive distribution for 15 day maxima ibovespa returns
## Not run: data(ibovespa)
## Not run: postibv=gumbelp(ibovespa[,4],15,500)
}
