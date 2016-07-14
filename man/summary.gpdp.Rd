\name{summary.gpdp}
\alias{summary.gpdp}
\title{
Summarizing Posterior Distribution with Parameters of GPD
}
\description{
summary method for class "gpdp"
}
\usage{
\method{summary}{gpdp}(object, \dots)
}
\arguments{
  \item{object}{
an object of class \code{"gpdp"}, usually, a result of a call to \code{\link{gpdp}}.
}
\item{\dots}{	
further arguments passed to or from other methods.
}
}
\value{
The function \code{summary.ggevp} computes and returns a list of summary statistics of the posterior distribution given in \code{object}.
\item{postmean}{mean posterior}
\item{postmedian}{median posterior}
\item{postCI}{credibility interval}
\item{fitm}{fit measures for standard GPD model}
}
\seealso{
\code{\link{gpdp}}
}
\examples{
data(danish)
out=gpdp(danish,10,300)
a=summary(out)
a
}

