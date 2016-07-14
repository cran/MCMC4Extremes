gumbelp <-
function(data,block,int=1000)
  {thin=10;burnin=int*thin/2  
   fit=gev(data,block)
   data=fit$data;n=length(data)
   lpost=function(mu,sigma)
      {logpost=-n*log(sigma)-sum(((data-mu)/sigma))-sum(exp(-(data-mu)/sigma))
       logpost=logpost+(0.001-1)*log(sigma)-0.001*sigma-mu^2/2000
       logpost}
      mumc=array(0,c(burnin+int,1));sigmamc=array(0,c(burnin+int,1))
      mumc[1]=fit$par.ests[3];sigmamc[1]=fit$par.ests[2]
      Vu=(sigmamc[1]/10)
      Vsigma=(sigmamc[1]/25)^2
      for (i in 2:burnin)
         {muest=rnorm(1,mumc[i-1],Vu)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
      alpha=exp(lpost(muest,sigmaest)-lpost(mumc[i-1],sigmamc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0}  
          u=runif(1)
          if (u<alpha)
            {mumc[i]=muest
             sigmamc[i]=sigmaest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]}
            if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  mumcb=array(0,c(int));sigmamcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {muest=rnorm(1,mumc[i-1],Vu)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
      alpha=exp(lpost(muest,sigmaest)-lpost(mumc[i-1],sigmamc[i-1]))
      alpha=alpha*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
      alpha=alpha/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
       if (is.nan(alpha)){alpha=0}
          u=runif(1)
          if (u<alpha)
            {mumc[i]=muest
             sigmamc[i]=sigmaest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]}
          if ((i%%thin)==0)
              {mumcb[j]=mumc[i]
               sigmamcb[j]=sigmamc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
          cbind(mumcb,sigmamcb)

          estim=cbind(mumcb,sigmamcb)
      
          ests<-c(mean(estim[,1]), mean(estim[,2]))
          ests1<-c(median(estim[,1]), median(estim[,2]))
          ests2=array(0,c(2,2))
          ests2[1,]<-c(quantile(estim[,1],0.025), quantile(estim[,2],0.025))
          ests2[2,]<-c(quantile(estim[,1],0.975), quantile(estim[,2],0.975))
         
          out<-list(posterior=estim, data=data, postmean = ests, postmedian = ests1, postCI=ests2 , block=block)
          names(out$postmean) <- c("mu", "sigma")
          names(out$postmedian) <- c("mu", "sigma")
          dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("mu", "sigma"))
            
          class(out) <- "gumbelp"
          out           

}
