summary.gumbelp <-
function(object, ...){
     mu=object$posterior[,1]
     sigma=object$posterior[,2]
     I=length(mu)
     x=object$data
     n=length(object$data)

     med=0
     for (i in 1:I)
         {med=med-n*log(sigma[i])-sum((x-mu[i])/sigma[i])-sum(exp(-(x-mu[i])/sigma[i]))
          }
     med=-2*med/I
     mediaDIC=-n*log(mean(sigma))-sum((x-mean(mu))/mean(sigma))-sum(exp(-(x-mean(mu))/mean(sigma)))
     mediaDIC=-2*mediaDIC
     BIC=med+4*log(n)
     pD=med-mediaDIC
     DIC=med+pD
     AIC=med+4    
     fitm=cbind(AIC,BIC,pD,DIC) 
     
     ests<-c(mean(object$posterior[,1]), mean(object$posterior[,2]))
     ests1<-c(median(object$posterior[,1]), median(object$posterior[,2]))
     ests2=array(0,c(2,2))
     ests2[1,]<-c(quantile(object$posterior[,1],0.025), quantile(object$posterior[,2],0.025))
     ests2[2,]<-c(quantile(object$posterior[,1],0.975), quantile(object$posterior[,2],0.975))
     
     out<-list(postmean = ests, postmedian = ests1, postCI=ests2, block=object$block, fitm=fitm)
     
     names(out$postmean) <- c("mu", "sigma")
     names(out$postmedian) <- c("mu", "sigma")
     dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("mu", "sigma"))
     
     class(out)<-"summary.gumbelp"
     return(out)
    }
