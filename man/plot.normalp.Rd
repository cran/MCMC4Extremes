\name{plot.normalp}
\alias{plot.normalp}
\title{
Plot Fitted Normal Model
}
\description{
The plot method plot.normalp provides three differents plots: a histogram of the normal parameters, a plot of predictive density resulting of posterior distribution of normal parameters, and a return level plot of normal distribution.
}
\usage{
\method{plot}{normalp}(x, type = c("histogram"), \dots)
}
\arguments{
  \item{x}{
a \code{normalp} object
}
  \item{type}{
which chosen plot
}
  \item{\dots}{
other graphics parameters
}
}
\seealso{
\code{\link{normalp}}
}
\examples{
data(Nile)
p=normalp(Nile,600)
plot(p,"histogram")  
}
