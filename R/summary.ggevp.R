summary.ggevp <-
function(object, ...){
    mu=object$posterior[,1]
    sigma=object$posterior[,2]
    xi=object$posterior[,3]
    I=length(mu)
    x=object$data
    n=length(object$data)

    med=0
    for (i in 1:I)
        {med=med-n*log(sigma[i]*gamma(object$delta))-sum((1+xi[i]*(x-mu[i])/sigma[i])^(-1/xi[i]))-(1+object$delta/xi[i])*sum(log((1+xi[i]*(x-mu[i])/sigma[i])))
         }
    med=-2*med/I
    mediaDIC=-n*log(mean(sigma)*gamma(object$delta))-sum((1+mean(xi)*(x-mean(mu))/mean(sigma))^(-1/mean(xi)))-(1+object$delta/mean(xi))*sum(log((1+mean(xi)*(x-mean(mu))/mean(sigma))))
    mediaDIC=-2*mediaDIC
    BIC=med+6*log(n)
    pD=med-mediaDIC
    DIC=med+pD
    AIC=med+6        
    
    fitm=cbind(AIC,BIC,pD,DIC)
         
     ests<-c(mean(object$posterior[,1]), mean(object$posterior[,2]), mean(object$posterior[,3]))
     ests1<-c(median(object$posterior[,1]), median(object$posterior[,2]), median(object$posterior[,3]))
     ests2=array(0,c(2,3))
     ests2[1,]<-c(quantile(object$posterior[,1],0.025), quantile(object$posterior[,2],0.025), quantile(object$posterior[,3],0.025))
     ests2[2,]<-c(quantile(object$posterior[,1],0.975), quantile(object$posterior[,2],0.975), quantile(object$posterior[,3],0.975))
     
     out<-list(postmean = ests, postmedian = ests1, postCI=ests2, block=object$block, delta=object$delta, fitm=fitm)
     
     names(out$postmean) <- c("mu", "sigma", "xi")
     names(out$postmedian) <- c("mu", "sigma", "xi")
     dimnames(out$postCI) <- list(c("lower bound", "upper bound"),c("mu", "sigma", "xi"))
     
     class(out)<-"summary.ggevp"
     return(out)
    }
