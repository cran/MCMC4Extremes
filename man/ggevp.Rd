\name{ggevp}
\alias{ggevp}
\title{
Posterior Distribution with Parameters of Dual Gamma Generalized Extreme Value Distribution
}
\description{
MCMC runs of posterior distribution of data with parameters of Dual Gamma Generalized Extreme Value Distribution 
density, with parameters \code{mu}, \code{sigma} and \code{xi}.
}
\usage{
ggevp(data, block, int=1000, delta)
}
\arguments{
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
  \item{delta}{
additional shape parameter of GGEV extension
}
}
\value{
An object of class \code{ggevp} that gives a list containing the points of posterior distributions of \code{mu}, \code{sigma} and \code{xi} of the dual gamma generalized extreme value distribution, the data, mean posterior, median posterior and the credibility interval of the parameters.
}
\references{
Nascimento, F. F.; Bourguigon, M. ; Leao, J. S. (2015). Extended generalized extreme value distribution with applications in environmental data. \emph{HACET J MATH STAT}.
}
\seealso{
\code{\link{plot.ggevp}}, \code{\link{summary.ggevp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
w=rggev(300,0.1,10,5,0.5)

# Obtaning 500 points of posterior distribution with delta=0.5
ajust=ggevp(w,1,200,0.5)
}
