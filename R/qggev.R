qggev <-
function(p,xi,mu,sigma,delta)
      {q=mu+(sigma/xi)*(qgamma(1-p,delta,1)^(-xi)-1)
       return(q)  }
