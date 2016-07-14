\name{summary.gumbelp}
\alias{summary.gumbelp}
\title{
Summarizing Posterior Distribution with Parameters of Gumbel
}
\description{
summary method for class "gumbelp"
}
\usage{
\method{summary}{gumbelp}(object, \dots)
}
\arguments{
  \item{object}{
an object of class \code{"gumbelp"}, usually, a result of a call to \code{\link{gumbelp}}.
}
\item{\dots}{	
further arguments passed to or from other methods.
}
}
\value{
The function \code{summary.gumbelp} computes and returns a list of summary statistics of the posterior distribution given in \code{object}.
\item{postmean}{mean posterior}
\item{postmedian}{median posterior}
\item{postCI}{credibility interval}
\item{fitm}{fit measures for standard Gumbel model}
}
\seealso{
\code{\link{gumbelp}}
}
\examples{
# Example with simulated datapoints
x=rgev(300,0.01,10,5)
fit=gumbelp(x,1,300)
fitgum=summary(fit)

# Compare if the fit measures of gumbel is better than measures using GEV
## Not run: fit2=gevp(x,1,300)
## Not run: fitgev=summary(fit2)
# the best model is that with lowest fit measures
}
