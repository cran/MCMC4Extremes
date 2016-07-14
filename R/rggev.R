rggev <-
function(n,xi,mu,sigma,delta)
     {u=runif(n)
      x=mu+(sigma/xi)*( (qgamma(1-u,delta,1))^(-xi) - 1 ) 
      return(x)}
