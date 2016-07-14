\name{plot.gevp}
\alias{plot.gevp}
\title{
Plot Fitted GEV Model
}
\description{
  The plot method plot.gevp provides three differents plots: a histogram of the gev parameters, a plot of predictive density resulting of posterior distribution of GEV parameters, and a return level plot of GEV distribution.
}
\usage{
\method{plot}{gevp}(x, type = c("histogram", "predictive", "retlevel"), t=2, k=100, \dots)
}
\arguments{
  \item{x}{
a \code{gevp} object
}
  \item{type}{
which chosen plot
}
  \item{t}{
start return level
}
  \item{k}{
end return level
}
  \item{\dots}{
other graphics parameters
}
}
\seealso{
\code{\link{gevp}}
}
\examples{
# Return level of river nidd data
data(nidd.annual)
out=gevp(nidd.annual,1,300)
## Not run: plot(out,"histogram")
plot(out,"predictive")
## Not run: plot(out,"retlevel", 10, 50)
}
