plot.gumbelp <-
function(x, type = c("histogram","predictive", "retlevel"), t=2, k=100, ...){

    if(type=="histogram"){
      par(mfrow=c(1,2))
      hist(x$posterior[,1],freq=FALSE,main=expression(mu), xlab=NULL)
      abline(v=quantile(x$posterior[,1],0.025),lwd=2)
      abline(v=quantile(x$posterior[,1],0.975),lwd=2)
      hist(x$posterior[,2],freq=FALSE,main=expression(sigma), xlab=NULL)
      abline(v=quantile(x$posterior[,2],0.025),lwd=2)
      abline(v=quantile(x$posterior[,2],0.975),lwd=2)
      par(mfrow=c(1,1))
    }
    if(type=="predictive"){
       linf=max(min(x$data)-1,0)
       lsup=11*max(x$data)/10
       dat1=seq(linf,lsup,(lsup-linf)/300)
       n=length(dat1)
       int=length(x$posterior[,1])
       mu=x$posterior[,1]
       sigma=x$posterior[,2]
       res=array(0,c(n))
       for (i in 1:n)
          {for (j in 1:int)
              {
               res[i]=res[i]+(1/int)*(1/sigma[j])*exp(-(dat1[i]-mu[j])/sigma[j])*exp(-exp(-(dat1[i]-mu[j])/sigma[j]))}}
      
      hist(x$data,freq=F,ylim=c(min(res),max(res)),xlab="data",ylab="density", ...)
      lines(dat1,res)
      out<-list(pred=res)
      return(out)
    }
    
    if(type=="retlevel"){
      sampl=x$posterior[,1]-x$posterior[,2]*log(-log(1-1/t))
      res=quantile(sampl,0.5)
  
      t=seq(1,k,1)
      n=length(t)
      li=array(0,c(n))
      ls=array(0,c(n))
      pred=array(0,c(n))
      for (s in 1:n)
            {sampl=x$posterior[,1]-x$posterior[,2]*log(-log(1-1/s))
             li[s]=quantile(sampl,0.025)   
             ls[s]=quantile(sampl,0.975)
             pred[s]=quantile(sampl,0.5)
            }
      plot(t,pred,type="l",ylim=c(li[2],max(ls)),ylab="returns", ...)
      lines(t,li,lty=2)
      lines(t,ls,lty=2)
      out<-list(retmedian=res, retpred=pred)
  
      return(out)
    }

}
