dggev <-
function(x,xi,mu,sigma,delta)
      {x=(1/(sigma*gamma(delta)))*(1+xi*(x-mu)/sigma)^(-1-delta/xi)*exp(-(1+xi*(x-mu)/sigma)^(-1/xi))
      return(x)}
