gpdp <-
function(data,threshold,int=1000)
  {thin=10;burnin=int*thin/2
   data=data[data>threshold];n=length(data)
   lpost=function(sigma,xi)
      {logpost=-n*log(sigma)-(1+1/xi)*sum(log(1+xi*(data-threshold)/sigma))
       logpost=logpost-log(sigma)-log(1+xi)-0.5*log(1+2*xi)
       logpost}
      sigmamc=array(0,c(burnin+int,1));ximc=array(0,c(burnin+int,1))
      sigmamc[1]=sd(data);ximc[1]=0.1
      Vsigma=sqrt(sigmamc[1]);Vxi=0.2
      while ((min(1+ximc[1]*(data-threshold)/sigmamc[1])<0)  | (ximc[1]< -0.5))
              {ximc[1]=rnorm(1,ximc[1],Vxi)
              sigmamc[1]=rgamma(1,sigmamc[1]^2/Vsigma,sigmamc[1]/Vsigma)}
      for (i in 2:burnin)
         {xiest=rnorm(1,ximc[i-1],Vxi)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
          while ((min(1+xiest*(data-threshold)/sigmaest)<0)  | (xiest< -0.5)  )
              {xiest=rnorm(1,ximc[i-1],Vxi)
              sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)}
      alpha=exp(lpost(sigmaest,xiest)-lpost(sigmamc[i-1],ximc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0}
          u=runif(1)
          if (u<alpha)
            {sigmamc[i]=sigmaest
             ximc[i]=xiest}
          else
             {sigmamc[i]=sigmamc[i-1]
             ximc[i]=ximc[i-1]}
            if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  sigmamcb=array(0,c(int));ximcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {xiest=rnorm(1,ximc[i-1],Vxi)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
          while ((min(1+xiest*(data-threshold)/sigmaest)<0)  | (xiest< -0.5)  )
              {xiest=rnorm(1,ximc[i-1],Vxi)
              sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)}
      alpha=exp(lpost(sigmaest,xiest)-lpost(sigmamc[i-1],ximc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0} 
          u=runif(1)
          if (u<alpha)
            {sigmamc[i]=sigmaest
             ximc[i]=xiest}
          else
             {sigmamc[i]=sigmamc[i-1]
             ximc[i]=ximc[i-1]}
          if ((i%%thin)==0)
              {sigmamcb[j]=sigmamc[i]
               ximcb[j]=ximc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
          estim=cbind(sigmamcb,ximcb)
          
    ests<-c(mean(estim[,1]), mean(estim[,2]))
    ests1<-c(median(estim[,1]), median(estim[,2]))
    ests2=array(0,c(2,2))
    ests2[1,]<-c(quantile(estim[,1],0.025), quantile(estim[,2],0.025))
    ests2[2,]<-c(quantile(estim[,1],0.975), quantile(estim[,2],0.975))
   
    out<-list(posterior=estim, data=data, postmean = ests, postmedian = ests1, postCI=ests2, threshold=threshold)
    names(out$postmean) <- c("sigma", "xi")
    names(out$postmedian) <- c("sigma", "xi")
    dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("sigma", "xi"))
      
    class(out) <- "gpdp"
    out                 
}
