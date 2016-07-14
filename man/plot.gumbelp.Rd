\name{plot.gumbelp}
\alias{plot.gumbelp}
\title{
Plot Fitted Gumbel Model
}
\description{
The plot method plot.gumbelp provides three differents plots: a histogram of the gumbel parameters, a plot of predictive density resulting of posterior distribution of gumbel parameters, and a return level plot of gumbel distribution.
}
\usage{
\method{plot}{gumbelp}(x, type = c("histogram", "predictive", "retlevel"), t=2, k=100, \dots)
}
\arguments{
  \item{x}{
a \code{gumbelp} object
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
\code{\link{gumbelp}}
}
\examples{
data(nidd.annual)
out=gumbelp(nidd.annual,1,500)
## Not run: plot(out,"histogram")
## Not run: plot(out,"predictive")
plot(out,"retlevel", 10)
}
