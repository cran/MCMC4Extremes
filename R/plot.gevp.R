plot.gevp <-
function(x, type = c("histogram","predictive", "retlevel"), t=2, k=100, ...){

    if(type=="histogram"){
      par(mfrow=c(1,3))
      hist(x$posterior[,1],freq=FALSE,main=expression(mu), xlab=NULL)
      abline(v=quantile(x$posterior[,1],0.025),lwd=2)
      abline(v=quantile(x$posterior[,1],0.975),lwd=2)
      hist(x$posterior[,2],freq=FALSE,main=expression(sigma), xlab=NULL)
      abline(v=quantile(x$posterior[,2],0.025),lwd=2)
      abline(v=quantile(x$posterior[,2],0.975),lwd=2)
      hist(x$posterior[,3],freq=FALSE,main=expression(xi), xlab=NULL)
      abline(v=quantile(x$posterior[,3],0.025),lwd=2)
      abline(v=quantile(x$posterior[,3],0.975),lwd=2)
      par(mfrow=c(1,1))
    }
    if(type=="predictive"){
      dat=x$data
      linf = max(min(dat) - 1, 0)
      lsup = 11 * max(dat)/10
      dat1 = seq(linf, lsup, (lsup - linf)/70)
      n = length(dat1)
      int = length(x$posterior[, 1])
      res = array(0, c(n))
      for (i in 1:n) {
          for (j in 1:int) {
              if ((x$posterior[j, 3] > 0) && (dat1[i] > (x$posterior[j, 1] - 
                  x$posterior[j, 2]/x$posterior[j, 3]))) 
                  res[i] = res[i] + (1/int) * dgev(dat1[i], x$posterior[j, 
                    3], x$posterior[j, 1], x$posterior[j, 2])
              if ((x$posterior[j, 3] < 0) && (dat1[i] < (x$posterior[j, 1] - 
                  x$posterior[j, 2]/x$posterior[j, 3]))) 
                  res[i] = res[i] + (1/int) * dgev(dat1[i], x$posterior[j, 
                    3], x$posterior[j, 1], x$posterior[j, 2])
          }
      }
      hist(dat, freq = F, ylim = c(min(res), max(res)), main=NULL, xlab = "data", ylab = "density", ...)
      lines(dat1, res)
      out<-list(pred=res)
      return(out)
    }
    
    if(type=="retlevel"){
      sampl = qgev(1 - 1/t, x$posterior[, 3], x$posterior[, 1], x$posterior[, 2])
      res = quantile(sampl, 0.5)
  
      ta=seq(1,k,1)
      n = length(ta)
      li = array(0, c(n))
      ls = array(0, c(n))
      pred = array(0, c(n))
      for (s in 1:n) {
          sampl = qgev(1 - 1/s, x$posterior[, 3], x$posterior[, 1], x$posterior[, 
              2])
          li[s] = quantile(sampl, 0.025)
          ls[s] = quantile(sampl, 0.975)
          pred[s] = quantile(sampl, 0.5)
      }
      plot(ta, pred, type = "l", xlab="t", ylim = c(li[2], max(ls)), ylab = "returns", ...)
      lines(ta, li, lty = 2)
      lines(ta, ls, lty = 2)
      out<-list(retmedian=res, retpred=pred)
  
      return(out)
    }

}
