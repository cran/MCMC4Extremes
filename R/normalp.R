normalp <-
function(data,int=1000)
{n=length(data);thin=10;burnin=int*thin/2
mumc=array(0,c(burnin+int,1));taumc=array(0,c(burnin+int,1))
      mumc[1]=mean(data);taumc[1]=1/var(data);meandata=mean(data)
      for (i in 2:burnin)
         {tau1=n*taumc[i-1]+1/10000000
          mumc[i]=rnorm(1,(1/tau1)*(n*meandata*taumc[i-1]),1/sqrt(tau1))
          taumc[i]=rgamma(1,n/2+0.001,0.001+0.5*sum((data-mumc[i])^2))
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
  mumcb=array(0,c(int));taumcb=array(0,c(int));j=1
      for (i in (burnin+1):(burnin+thin*int))
         {tau1=n*taumc[i-1]+1/1000000
          mumc[i]=rnorm(1,(1/tau1)*(n*meandata*taumc[i-1]),1/sqrt(tau1))
          taumc[i]=rgamma(1,n/2+0.001,0.001+0.5*sum((data-mumc[i])^2))
          if ((i%%thin)==0)
              {mumcb[j]=mumc[i]
               taumcb[j]=taumc[i]
               j=j+1}
          if ((i%%100)==0)
             print(i/(burnin+thin*int))}
          estim=cbind(mumcb,taumcb)
          
          ests<-c(mean(estim[,1]), mean(estim[,2]))
          ests1<-c(median(estim[,1]), median(estim[,2]))
          ests2=array(0,c(2,2))
          ests2[1,]<-c(quantile(estim[,1],0.025), quantile(estim[,2],0.025))
          ests2[2,]<-c(quantile(estim[,1],0.975), quantile(estim[,2],0.975))
         
          out<-list(posterior=estim, data=data, postmean = ests, postmedian = ests1, postCI=ests2)
          names(out$postmean) <- c("mu", "tau")
          names(out$postmedian) <- c("mu", "tau")
          dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("mu", "tau"))
            
          class(out) <- "normalp"
          out          
          
}
