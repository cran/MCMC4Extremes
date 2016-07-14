\name{normalp}
\alias{normalp}
\title{
Posterior Distribution with Normal Density
}
\description{
MCMC runs of posterior distribution of data with \code{Normal(mu,1/tau)} density, where \code{tau} is the inverse
of variance.
}
\usage{
normalp(data, int=1000)
}
\arguments{
  \item{data}{
data vector
}
  \item{int}{
number of iteractions selected in MCMC. The program selects 1 in each 10
iteraction, then \code{thin=10}. The first \code{thin*int/3} iteractions is used as burn-in. After
that, is runned \code{thin*int} iteraction, in which 1 of thin is selected for the final
MCMC chain, resulting the number of int iteractions
}
}
\value{
An object of class \code{gumbelp} that gives a list containing the points of posterior distributions of \code{mu} and \code{tau} of the normal distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\note{
The non-informative prior distribution of these parameters are \code{Normal(0,10000000)}
for the parameter mu and \code{Gamma(0.001,0.001)} for the parameter \code{tau} . During the MCMC runs,
screen shows the proportion of iteractions made.
}
\seealso{
\code{\link{plot.normalp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
x=rnorm(300,2,sqrt(10))

# Obtaning 1000 points of posterior distribution
ajuste=normalp(x, 200)

# Posterior distribution of river Nile dataset
## Not run: data(Nile)
## Not run: postnile=normalp(Nile,1000)
}
