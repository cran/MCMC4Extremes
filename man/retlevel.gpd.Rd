\name{retlevel.gpd}
\alias{retlevel.gpd}
\title{
Expected Return Level in each Block Times of GPD
}
\description{
Calculates the expected return level in each t times.
}
\usage{
retlevel.gpd(vector, t, data, threshold)
}
\arguments{
  \item{vector}{
a list object returned by posterior.gpd
}
  \item{t}{
  specifies the particular return level to be estimated
}
  \item{data}{
  data vector
}
  \item{threshold}{
  a threshold value
}
}
\value{
The program returns the posterior mean of expected return level in each t
periods of time.
}
\examples{
# Calculate the return level of each 12 periods of time in danish data
data(danish)
out=posterior.gpd(danish,10,1000)
r12=retlevel.gpd(out,12,danish,10)
}
