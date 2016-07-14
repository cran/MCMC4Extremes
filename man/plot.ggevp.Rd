\name{plot.ggevp}
\alias{plot.ggevp}
\title{
Plot Fitted for the Dual Gamma Generalized Extreme Value Distribution (GGEV) Model
}
\description{
The plot method plot.ggevp provides three differents plots: a histogram of the GGEV parameters, a plot of predictive density resulting of posterior distribution of GGEV parameters, and a return level plot of GGEV distribution.
}
\usage{
\method{plot}{ggevp}(x, type = c("histogram", "predictive", "retlevel"), t=2, k = 100, \dots)
}
\arguments{
  \item{x}{
a \code{ggevp} object
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
\references{
Nascimento, F. F.; Bourguigon, M. ; Leao, J. S. (2015). Extended generalized extreme value distribution with applications in environmental data. \emph{HACET J MATH STAT}.
}
\seealso{
\code{\link{ggevp}}
}
\examples{
# Obtaining posterior distribution of a vector of simulated points
w=rggev(300,0.4,10,5,0.5)

# Obtaning 300 points of posterior distribution with delta=0.5
fit=ggevp(w,1,200,0.5)
## Not run: plot(fit,"histogram")
plot(fit,"predictive")
## Not run: plot(fit,"retlevel", 10, 50)
}