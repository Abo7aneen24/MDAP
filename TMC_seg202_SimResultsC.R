#  TMC_seg202_SimResultsC.R
#
#  Produces summarizing figures C for simulation data 
#  Is called by TMC_seg202_SimResults.R

#  19.04.2021 202.050 changed to plot of cdf
#  18.04.2021 Additional output files 202.040, 202.050
#  17.04.2021 sum.reps renamed to scen.sum
#  15.12.2020 Output file renamed to 202.030
#  11.12.2020 Start
# =============================================================================

if (print.log.message) { cat("%%%   TMC_seg202_SimResultsC  Start\n") }

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.030) & (!is.na(lambda.c.gen)))
{   
  dev.set(fig202.030)
  bringToTop(fig202.030)

  #  Identical scale on both axes is needed
  xplotmin <- min(c(FPR1, FPR2), na.rm=TRUE)
  xplotmax <- max(c(FPR1, FPR2), na.rm=TRUE)
  yplotmin <- xplotmin
  yplotmax <- xplotmax

  for (imeth in 1:meth.list.202.n)
  { 
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    meth.col <- meth.dsg[which(meth.dsg[ ,1]==meth) ,2]
    meth.pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1]==meth) ,3])
    if (imeth == 1)
    {
      plot(FPR1[ok], FPR2[ok], type="p", col=meth.col, pch=meth.pch,
           xlim=c(xplotmin, xplotmax), ylim=c(yplotmin, yplotmax),
           xlab="FPR1 (%)", ylab="FPR2 (%)", 
           main="Absolute False Positive Rates (%) per replicate",
           sub=subtitle.202, cex.sub=0.8)
      points(mean(FPR1[ok], na.rm=TRUE), mean(FPR2[ok], na.rm=TRUE), 
            col="black", pch=16, cex=1.5) 
    } else
    {
      lines(FPR1[ok], FPR2[ok], type="p", col=meth.col, pch=meth.pch)
      points(mean(FPR1[ok], na.rm=TRUE), mean(FPR2[ok], na.rm=TRUE), 
             col=meth.col, pch=16, cex=1.5) 
    }
  }
  abline(v=RL1.p, col=gencol)
  abline(h=1-RL2.p, col=gencol)
  abline(c(RL1.p+1-RL2.p, -1), col=gencol, lty=2)

  legend("topright",meth.dsg[meth.list.202.num,1],
         col=meth.dsg[meth.list.202.num,2],
         pch=as.numeric(meth.dsg[meth.list.202.num,3]), lwd=3, cex=0.8)

  savePlot(file=fig202.030.file,type=figtype)
}

# ----------------------------------------------------------------------------
# Density of FPR sum

if ((nrep1 > 1) & !is.na(fig202.040) & (!is.na(lambda.c.gen)))
{   
  dev.set(fig202.040)
  bringToTop(fig202.040)

  xplotmin <- min(FPR.sum, na.rm=TRUE)
  xplotmax <- max(FPR.sum, na.rm=TRUE)

  FPR.kde <- matrix(NA, nrow=512, ncol=meth.list.202.n)
  for (imeth in 1:meth.list.202.n)
  {    
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    kde      <- density(FPR.sum[ok], from=xplotmin, to=xplotmax, 
                                 n=512)
    FPR.kde[ , imeth] <- kde$y
  }

  yplotmin <- 0
  yplotmax <- max(FPR.kde, na.rm=TRUE)

  for (imeth in 1:meth.list.202.n)
  { 
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    meth.col <- meth.dsg[which(meth.dsg[ ,1]==meth) ,2]
    meth.pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1]==meth) ,3])
    if (imeth == 1)
    {
      plot(kde$x, FPR.kde[ ,imeth], type="l", col=meth.col,
           xlim=c(xplotmin, xplotmax), ylim=c(yplotmin, yplotmax),
           xlab="FPR.sum (%)", ylab="KDE", 
           main="PDF of False Positive Rates(%)",
           sub=subtitle.202, cex.sub=0.8)
      abline(v=mean(FPR.sum[ok]), col=meth.col)
      abline(v=median(FPR.sum[ok]), col=meth.col, lty=2)
    } else
    {
      lines(kde$x, FPR.kde[ ,imeth], type="l", col=meth.col)
      abline(v=mean(FPR.sum[ok]), col=meth.col)
      abline(v=median(FPR.sum[ok]), col=meth.col, lty=2)
    }
  }
  abline(h=0)
  abline(v=RL1.p + 1 - RL2.p, col="magenta")
  legend("topright",meth.dsg[meth.list.202.num,1],
         col=meth.dsg[meth.list.202.num,2],
         lwd=3,cex=0.8)

  savePlot(file=fig202.040.file,type=figtype)
}

# ----------------------------------------------------------------------------
#  delta FPR. Show cdf, not pdf

if ((nrep1 > 1) & !is.na(fig202.050) & (!is.na(lambda.c.gen)))
{   
  dev.set(fig202.050)
  bringToTop(fig202.050)

  qua.list <- seq(0.05, 1, by=0.05)
  delta.FPR.qua <- matrix(NA, nrow=meth.list.202.n, ncol=length(qua.list))
  
  xplotmin <- 0
  xplotmax <- max(delta.FPR, na.rm=TRUE)

  delta.FPR.x   <- matrix(NA, nrow=nrep+1, ncol=meth.list.202.n)
  delta.FPR.cdf <- matrix(NA, nrow=nrep+1, ncol=meth.list.202.n)
  for (imeth in 1:meth.list.202.n)
  {    
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    #kde      <- density(delta.FPR[ok], from=xplotmin, to=xplotmax, 
    #                             n=512)
    cdf.x     <- sort(delta.FPR[ok])
    cdf.x     <- cdf.x[!is.na(cdf.x)]
    cdf.x.n   <- length(cdf.x)
    cdf.y     <- (1:cdf.x.n)/cdf.x.n
    cdf.x     <- c(0, cdf.x)
    cdf.y     <- c(0, cdf.y)
    cdf.x.n   <- cdf.x.n + 1  
    delta.FPR.x[1:cdf.x.n, imeth]  <- cdf.x
    delta.FPR.cdf[1:cdf.x.n, imeth] <- cdf.y
    delta.FPR.qua[imeth, ] <- Quantile(delta.FPR[ok], probs=qua.list)
  }

  #  Weird action to bring 'meth' into the data frame. Using meth as row names
  #  produces a csv file with the first col name missing. 
  delta.FPR.qua <- data.frame(meth=meth.list.202, delta.FPR.qua, 
                              stringsAsFactors=FALSE)
  colnames(delta.FPR.qua) <- c("meth", paste(formatC(100*qua.list, width=3, 
                                             digits=0,
                                             format="f"), "%", sep=""))
  #  Write assessment matrix to file
  write.table(delta.FPR.qua, file=outname.scen.asm.csv, row.names=FALSE,
              col.names=TRUE, sep=";", dec=".") 

  yplotmin <- 0
  # yplotmax <- max(delta.FPR.kde, na.rm=TRUE)
  yplotmax <- 1

  for (imeth in 1:meth.list.202.n)
  { 
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    imeth.n  <- sum(!is.na(delta.FPR.x[ , imeth]))

    meth.col <- meth.dsg[which(meth.dsg[ ,1]==meth) ,2]
    meth.pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1]==meth) ,3])
    if (imeth == 1)
    {
      plot(delta.FPR.x[1:imeth.n, imeth], delta.FPR.cdf[1:imeth.n, imeth], 
           type="s", col=meth.col, lwd=2,
           xlim=c(xplotmin, xplotmax), ylim=c(yplotmin, yplotmax),
           xlab="delta.FPR (%)", ylab="cdf(delta.FPR)", 
           main="CDF of absolute sum of FPR errors (%)",
           sub=subtitle.202, cex.sub=0.8)
      abline(v=mean(delta.FPR[ok]), col=meth.col)
      # abline(v=median(delta.FPR[ok]), col=meth.col, lty=2)
    } else
    {
      lines(delta.FPR.x[1:imeth.n, imeth], delta.FPR.cdf[1:imeth.n, imeth], 
            type="s", col=meth.col, lwd=2)
      abline(v=mean(delta.FPR[ok]), col=meth.col)
      # abline(v=median(delta.FPR[ok]), col=meth.col, lty=2)
    }
  }
  #dimnames(delta.FPR.qua) <- list(meth.list.202, )

  axis(2, at=seq(0.05, 0.95, by=0.05), labels=FALSE)

  legend("bottomright",meth.dsg[meth.list.202.num,1],
         col=meth.dsg[meth.list.202.num,2],
         lwd=3,cex=0.8)

  abline(h=c(0, 1))
  abline(h=seq(0.1, 0.9, by=0.1), lty=2)
  abline(v=0)
  savePlot(file=fig202.050.file,type=figtype)
}
  
if (print.log.message) { cat("%%%   TMC_seg202_SimResultsC  End\n") }
