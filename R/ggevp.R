ggevp <-
function(data,block,int=1000,delta)
  {thin=10
   burnin=thin*int/2
   fit=gev(data,block)
   dat=fit$data;n=length(dat)
   lpost1=function(mu,sigma,xi)
      {
       logpost=-n*log(sigma)-sum((1+xi*(dat-mu)/sigma)^(-1/xi))
       logpost=logpost-(1+delta/xi)*sum(log((1+xi*(dat-mu)/sigma)))
       logpost=logpost+(0.001-1)*log(sigma)-0.001*sigma-mu^2/2000-xi^2/2
       logpost}
      mumc=array(0,c(burnin+int,1));sigmamc=array(0,c(burnin+int,1))
      ximc=array(0,c(burnin+int,1))
      mumc[1]=fit$par.ests[3];sigmamc[1]=fit$par.ests[2]
      ximc[1]=fit$par.ests[1];
      Vu=0.1;Vsigma=0.5;Vxi=0.1;
      while (min(1+ximc[1]*(dat-mumc[1])/sigmamc[1])<0)
              {mumc[1]=rnorm(1,mumc[1],Vu)
              ximc[1]=rnorm(1,ximc[1],Vxi)
              sigmamc[1]=rgamma(1,sigmamc[1]^2/Vsigma,sigmamc[1]/Vsigma)}
      for (i in 2:burnin)
         {muest=rnorm(1,mumc[i-1],Vu)
          xiest=rnorm(1,ximc[i-1],Vxi)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
          while (min(1+xiest*(dat-muest)/sigmaest)<0)
              {muest=rnorm(1,mumc[i-1],Vu)
               xiest=rnorm(1,ximc[i-1],Vxi)
               sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)}
  alphaa=exp(lpost1(muest,sigmaest,xiest)-lpost1(mumc[i-1],sigmamc[i-1],ximc[i-1]))
  alphaa=alphaa*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
  alphaa=alphaa/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
   if (is.nan(alphaa)){alphaa=0}
          u=runif(1)
          if (u<alphaa)
            {mumc[i]=muest
             sigmamc[i]=sigmaest
             ximc[i]=xiest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]
             ximc[i]=ximc[i-1]}
            if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  mumcb=array(0,c(int));sigmamcb=array(0,c(int));ximcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {muest=rnorm(1,mumc[i-1],Vu)
          xiest=rnorm(1,ximc[i-1],Vxi)
          sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)
          while (min(1+xiest*(dat-muest)/sigmaest)<0)
              {muest=rnorm(1,mumc[i-1],Vu)
               xiest=rnorm(1,ximc[i-1],Vxi)
               sigmaest=rgamma(1,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma)}
  alphaa=exp(lpost1(muest,sigmaest,xiest)-lpost1(mumc[i-1],sigmamc[i-1],ximc[i-1]))
  alphaa=alphaa*dgamma(sigmamc[i-1],sigmaest^2/Vsigma,sigmaest/Vsigma)
  alphaa=alphaa/(dgamma(sigmaest,sigmamc[i-1]^2/Vsigma,sigmamc[i-1]/Vsigma))
   if (is.nan(alphaa)){alphaa=0}
          u=runif(1)
          if (u<alphaa)
            {mumc[i]=muest
             sigmamc[i]=sigmaest
             ximc[i]=xiest}
          else
             {mumc[i]=mumc[i-1]
             sigmamc[i]=sigmamc[i-1]
             ximc[i]=ximc[i-1]}
          if ((i%%thin)==0)
              {mumcb[j]=mumc[i]
               sigmamcb[j]=sigmamc[i]
               ximcb[j]=ximc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int)) }
    estim=cbind(mumcb,sigmamcb,ximcb)
  
    ests<-c(mean(estim[,1]), mean(estim[,2]), mean(estim[,3]))
    ests1<-c(median(estim[,1]), median(estim[,2]), mean(estim[,3]))
    ests2=array(0,c(2,3))
    ests2[1,]<-c(quantile(estim[,1],0.025), quantile(estim[,2],0.025), quantile(estim[,3],0.025))
    ests2[2,]<-c(quantile(estim[,1],0.975), quantile(estim[,2],0.975), quantile(estim[,3],0.975))
   
    out<-list(posterior=estim, data=dat, postmean = ests, postmedian = ests1, postCI=ests2, block=block, delta=delta)
    names(out$postmean) <- c("mu","sigma", "xi")
    names(out$postmedian) <- c("mu","sigma", "xi")
    dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("mu","sigma", "xi"))
      
    class(out) <- "ggevp"
    out      

}

