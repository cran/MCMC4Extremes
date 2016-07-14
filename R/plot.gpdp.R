plot.gpdp <-
function(x, type = c("histogram","predictive", "retlevel"), t=2, k=100, ...){

    if(type=="histogram"){
      par(mfrow=c(1,2))
      hist(x$posterior[,1],freq=FALSE,main=expression(sigma),xlab=NULL)
      abline(v=quantile(x$posterior[,1],0.025),lwd=2)
      abline(v=quantile(x$posterior[,1],0.975),lwd=2)
      hist(x$posterior[,2],freq=FALSE,main=expression(xi),xlab=NULL)
      abline(v=quantile(x$posterior[,2],0.025),lwd=2)
      abline(v=quantile(x$posterior[,2],0.975),lwd=2)
      par(mfrow=c(1,1))
    }
    if(type=="predictive"){
       data=x$data[x$data>x$threshold]
       linf=min(data)
       lsup=max(data)
       dat1=seq(linf,lsup,(lsup-linf)/300)
       n=length(dat1)
       int=length(x$posterior[,1])
       res=array(0,c(n))
       for (i in 1:n)
          {for (j in 1:int)
              {res[i]=res[i]+(1/int)*dgpd(dat1[i],x$posterior[j,2],x$threshold,x$posterior[j,1])}}
      hist(data,freq=F,ylim=c(min(res),max(res)),breaks=seq(x$threshold,max(data),(max(data)-x$threshold)/20),main=NULL,xlab="data",ylab="density", ...)
      lines(dat1,res)
      
      out<-list(pred=res)
      return(out)
    }
    
    if(type=="retlevel"){
      p=1-1/t
      Nu=length(x$data[x$data>x$threshold])
      N=length(x$data)
      sampl=qgpd(1-(1/t)*N/Nu,x$posterior[,2],x$threshold,x$posterior[,1])
      res=quantile(sampl,0.5)
      res
  
      Nu=length(x$data[x$data>x$threshold])
      N=length(x$data)
      plimiar=1-Nu/N
      ta=seq(trunc(1/(1-plimiar)+1),k,1)
      n=length(ta)
      li=array(0,c(n))
      ls=array(0,c(n))
      pred=array(0,c(n))
      for (s in 1:n)
            {sampl=qgpd(1-(1/(s+min(t)))*N/Nu,x$posterior[,2],x$threshold,x$posterior[,1])
             li[s]=quantile(sampl,0.025)   
             ls[s]=quantile(sampl,0.975)
             pred[s]=quantile(sampl,0.5)
            }
      plot(ta,pred,type="l", xlab="t", ylim=c(max(min(li),0),max(ls)),ylab="returns", ...)
      lines(ta,li,lty=2)
      lines(ta,ls,lty=2)
  
      out<-list(retmedian=res, retpred=pred)
  
      return(out)
    }

}
