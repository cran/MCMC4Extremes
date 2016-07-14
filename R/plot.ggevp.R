plot.ggevp <-
function(x, type = c("histogram","predictive", "retlevel"), t=2, k=100, ...){

    if(type=="histogram"){
      par(mfrow=c(1,3))
      hist(x$posterior[,1],freq=FALSE,main=expression(mu), xlab=NULL)
      abline(v=quantile(x$posterior[,1],0.025),lwd=2)
      abline(v=quantile(x$posterior[,1],0.975),lwd=2)
      hist(x$posterior[,2],freq=FALSE,main=expression(sigma), xlab=NULL)
      abline(v=quantile(x$posterior[,2],0.025),lwd=2)
      abline(v=quantile(x$posterior[,2],0.975),lwd=2)
      hist(x$posterior[,3],freq=FALSE,main=expression(xi), xlab=NULL)
      abline(v=quantile(x$posterior[,3],0.025),lwd=2)
      abline(v=quantile(x$posterior[,3],0.975),lwd=2)
      par(mfrow=c(1,1))
    }
    if(type=="predictive"){
       mu=x$posterior[,1]
       sigma=x$posterior[,2]
       xi=x$posterior[,3]
       linf=min(x$data)
       lsup=max(x$data)
       dat1=seq(linf,lsup,(lsup-linf)/300)
       n=length(dat1)
       int=length(x$posterior[,1])
       res=array(0,c(n))
       for (i in 1:n)
          {for (j in 1:int)
              {res[i]=res[i]+(1/int)*(1/(sigma[j]*gamma(x$delta)))*(1+xi[j]*(dat1[i]-mu[j])/sigma[j])^(-1-x$delta/xi[j])*exp(-(1+xi[j]*(dat1[i]-mu[j])/sigma[j])^(-1/xi[j]))}}
      hist(x$data,freq=F,nclass=30, xlab="data", ylim=c(min(res),max(res)))
      lines(dat1,res)
      out<-list(pred=res)
      return(out)
    }
    
    if(type=="retlevel"){
      mu=x$posterior[,1]
      sigma=x$posterior[,2]
      xi=x$posterior[,3]
      qu=1-1/t
      sampl=mu+(sigma/xi)*(qgamma(1-qu,x$delta,1)^(-xi)-1)
      res=quantile(sampl,0.5)
  
      ta=seq(1,k,1)
      n=length(ta)
      li=array(0,c(n))
      ls=array(0,c(n))
      pred=array(0,c(n))
      mu=x$posterior[,1]
      sigma=x$posterior[,2]
      xi=x$posterior[,3]
      for (s in 1:n)
            {qu=1-1/s
             sampl=mu+(sigma/xi)*(qgamma(1-qu,x$delta,1)^(-xi)-1)
             li[s]=quantile(sampl,0.025)
             ls[s]=quantile(sampl,0.975)
             pred[s]=quantile(sampl,0.5)
            }
      plot(ta,pred,type="l", xlab="t", ylim=c(li[2],max(ls)))
      lines(ta,li,lty=2)
      lines(ta,ls,lty=2)
      
      out<-list(retmedian=res, retpred=pred)
      return(out)
    }

}
