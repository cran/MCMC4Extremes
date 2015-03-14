\name{retlevel.gumbel}
\alias{retlevel.gumbel}
\title{
Expected Return Level in each Block Times of Gumbel
}
\description{
Calculates the expected return level in each t block times.
}
\usage{
retlevel.gumbel(vector, t)
}
\arguments{
  \item{vector}{
  a list object returned by posterior.gumbel
}
  \item{t}{
  specifies the particular return level to be estimated
}
}
\value{
The program returns the posterior mena of expected return level in each t periods of time.
}
\examples{
# Expected return level for 15 day maxima ibovespa return
data(ibovespa)
postibv=posterior.gumbel(ibovespa[,4],15,1000)
retlevel.gumbel(postibv,20)
}
