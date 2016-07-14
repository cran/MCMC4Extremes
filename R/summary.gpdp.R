summary.gpdp <-
function(object, ...){
     data=object$data[object$data>object$threshold];n=length(data)
     n=length(data)
     sigmamc=object$posterior[,1]
     ximc=object$posterior[,2]
     inte=length(sigmamc)
     sigmabar=mean(sigmamc)
     xibar=mean(ximc)
     dbar=-2*(-n*log(sigmabar)-((1+xibar)/xibar)*sum(log(1+xibar*(data-object$threshold)/sigmabar)))
     bard=0
     for (i in 1:inte)
        {
        bard=bard-2*(-n*log(sigmamc[i])-((1+ximc[i])/ximc[i])*sum(log(1+ximc[i]*(data-object$threshold)/sigmamc[i])))       
        }
     bard=bard/inte
     pD=bard-dbar
     DIC=bard+pD
     AIC=bard+4 
     BIC=bard+4*log(n)
     fitm=cbind(AIC,BIC,pD,DIC)
     
     ests<-c(mean(object$posterior[,1]), mean(object$posterior[,2]))
     ests1<-c(median(object$posterior[,1]), median(object$posterior[,2]))
     ests2=array(0,c(2,2))
     ests2[1,]<-c(quantile(object$posterior[,1],0.025), quantile(object$posterior[,2],0.025))
     ests2[2,]<-c(quantile(object$posterior[,1],0.975), quantile(object$posterior[,2],0.975))
     
     out<-list(postmean = ests, postmedian = ests1, postCI=ests2, threshold=object$threshold, fitm=fitm)
     
     names(out$postmean) <- c("sigma", "xi")
     names(out$postmedian) <- c("sigma", "xi")
     dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("sigma", "xi"))
     
     class(out)<-"summary.gevp"
     return(out)
    }
