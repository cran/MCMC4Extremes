\name{retlevel.gev.graph}
\alias{retlevel.gev.graph}
\title{
Return level plot of GEV
}
\description{
Return level plot of GEV distribution.
}
\usage{
retlevel.gev.graph(vector)
}
\arguments{
  \item{vector}{
  a list object returned by posterior.gev
}
}
\value{
The program returns return level plot from 1 until 100 periods in time,
with 95\% credibility intervals
}
\examples{
# Return level plot for River nidd data
data(nidd.annual)
out=posterior.gev(nidd.annual,1,1000)
retlevel.gev.graph(out)
# Return level plot for ibovespa 15 day maxima
data(ibovespa)
postibv=posterior.gev(ibovespa[,4],15,300)
retlevel.gev.graph(postibv)
}
