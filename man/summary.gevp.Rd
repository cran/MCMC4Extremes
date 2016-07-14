\name{summary.gevp}
\alias{summary.gevp}
\title{
Summarizing Posterior Distribution with Parameters of GEV
}
\description{
summary method for class "gevp"
}
\usage{
\method{summary}{gevp}(object, \dots)
}
\arguments{
  \item{object}{
an object of class \code{"gevp"}, usually, a result of a call to \code{\link{gevp}}.
}
\item{\dots}{	
further arguments passed to or from other methods.
}
}
\value{
The function \code{summary.gevp} computes and returns a list of summary statistics of the posterior distribution given in \code{object}.
\item{postmean}{mean posterior}
\item{postmedian}{median posterior}
\item{postCI}{credibility interval}
\item{fitm}{fit measures for standard GEV model}
}

\seealso{
\code{\link{gevp}}
}
\examples{
# Return level of river nidd data
data(nidd.annual)
out=gevp(nidd.annual,1,300)
a=summary(out)
a
}

