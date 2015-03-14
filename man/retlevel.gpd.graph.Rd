\name{retlevel.gpd.graph}
\alias{retlevel.gpd.graph}
\title{
Return level plot of GPD
}
\description{
Return level plot of GPD distribution.
}
\usage{
retlevel.gpd.graph(vector, data, threshold)
}
\arguments{
  \item{vector}{
  a list object returned by posterior.gpd
}
  \item{data}{
  data vector
}
  \item{threshold}{
a threshold value
}
}
\value{
The program returns return level plot from 1 until 100 periods in time,
with 95\% credibility intervals.
}
\examples{
# Return level plot for danish data
data(danish)
out=posterior.gpd(danish,10,1000)
retlevel.gpd.graph(out,danish,10)
}
