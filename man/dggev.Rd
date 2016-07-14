\name{dggev}
\alias{dggev}
\alias{pggev}
\alias{qggev}
\alias{rggev}
\title{
Dual Gamma Generalized Extreme Value Distribution
}
\description{
Cumulative probability, quantiles, density and random generation from the dual gamma  generalized extreme value distribution. 
}
\usage{
pggev(q, xi, mu, sigma, delta)
qggev(p, xi, mu, sigma, delta)
dggev(x, xi, mu, sigma, delta)
rggev(n, xi, mu, sigma, delta)
}

\arguments{
\item{q}{
vector of quantiles
}
\item{p}{
vector of probabilities
}
  \item{x}{
vector of values at which to evaluate density
}
\item{n}{
sample size
}
  \item{xi}{
shape parameter
}
  \item{mu}{
location parameter
}
  \item{sigma}{
scale parameter
}
  \item{delta}{
additional shape parameter of GGEV extension
}
}
\value{
Probability (\code{pggev}), quantile (\code{qggev}), density (\code{dggev}) or random sample (\code{rggev}) for the GGEV distribution.
}
\references{
Nascimento, F. F.; Bourguigon, M. ; Leao, J. S. (2015). Extended generalized extreme value distribution with applications in environmental data. \emph{HACET J MATH STAT}.
}
\seealso{
\code{\link{ggevp}}
}
