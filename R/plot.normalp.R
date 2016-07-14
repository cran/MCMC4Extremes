plot.normalp <-
function(x, type = c("histogram"), ...){

    if(type=="histogram"){
      par(mfrow=c(1,2))
      hist(x$posterior[,1],freq=FALSE,main=expression(mu), xlab=NULL, ...)
      abline(v=quantile(x$posterior[,1],0.025),lwd=2)
      abline(v=quantile(x$posterior[,1],0.975),lwd=2)
      hist(x$posterior[,2],freq=FALSE,main=expression(tau), xlab=NULL, ...)
      abline(v=quantile(x$posterior[,2],0.025),lwd=2)
      abline(v=quantile(x$posterior[,2],0.975),lwd=2)
      par(mfrow=c(1,1))
    }
}
