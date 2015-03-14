\name{retlevel.gumbel.graph}
\alias{retlevel.gumbel.graph}
\title{
Return level plot of Gumbel
}
\description{
Return level plot for gumbel distribution.
}
\usage{
retlevel.gumbel.graph(vector)
}
\arguments{
  \item{vector}{
  a list object returned by posterior.gumbel
}
}
\value{
The program returns the return level plot from 1 until 100 periods of time, with 95\% credibility intervals.
}
\examples{
# Return level plot for 15 day ibovespa maxima returns
data(ibovespa)
postibv=posterior.gumbel(ibovespa[,4],15,1000)
retlevel.gumbel.graph(postibv)
}
