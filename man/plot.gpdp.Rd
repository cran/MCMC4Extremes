\name{plot.gpdp}
\alias{plot.gpdp}
\title{
Plot Fitted GPD Model
}
\description{
The plot method plot.gpdp provides three differents plots: a histogram of the GPD parameters, a plot of predictive density resulting of posterior distribution of GPD parameters, and a return level plot of GPD distribution.
}
\usage{
\method{plot}{gpdp}(x, type = c("histogram", "predictive", "retlevel"), t=2, k=100, \dots)
}
\arguments{
  \item{x}{
a \code{gpdp} object
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
\code{\link{gpdp}}
}
\examples{
data(danish)
out=gpdp(danish,10,300)
## Not run: plot(out,"histogram")
## Not run: plot(out,"predictive")
plot(out,"retlevel", 10, 50)
}