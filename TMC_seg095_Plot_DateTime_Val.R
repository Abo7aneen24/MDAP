#  TMC_seg095_Plot_DateTime_Val.R

#  (c) wwosniok@math.uni-bremen.de

#  Plot of values (mean+CI, median+mad)  against weekday and time of day

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg095_Plot_DateTime_Val.R  Start\n") }
# ===========================================================================

# ===========================================================================
#  Namen der Ausgabedatei in TMC_seg080_NamesPerFile.R festgelegt
# ===========================================================================

#  Graph windows

source("TMC_seg094_OpenWindows.R")

# ----------------------------------------------------------------------------

#  Execution only if date and time are given in the data

if (
     !is.na(spalte.dz) &&
     (!all(is.na(dataset[ ,"ana.wday"])) &  !all(is.na(dataset[ ,"ana.hour"])) ) 
   )
{ 
  #  Arrange data by weekday and time of day
  #  Find day of week

  #ana.day <-   31
  #ana.mon <-    5
  #ana.yea <- 2019
  #ana.date <- format(as.Date(paste(ana.yea,"-", ana.mon,"-", ana.day, sep="")),
  #                    "%d %b %Y" )
  #ana.wday <- format(as.Date(paste(ana.yea,"-", ana.mon,"-", ana.day, sep="")),
  #                   "%a" )   #  weekday, character

  xplot        <- dataset[ ,"ana.wday"] + dataset[ ,"ana.hour"]/24   
                                                          #  + ana.minu/1440
  xplot.list   <- sort(unique(xplot))
  xplot.list.n <- length(xplot.list)

  #  Go through xplot list and calculate statistics
  yplot.mea  <- rep(NA, times=xplot.list.n)
  yplot.sd   <- rep(NA, times=xplot.list.n)
  yplot.med  <- rep(NA, times=xplot.list.n)
  yplot.mad  <- rep(NA, times=xplot.list.n)
  yplot.n    <- rep(NA, times=xplot.list.n)
  yplot1.lo  <- rep(NA, times=xplot.list.n)
  yplot1.hi  <- rep(NA, times=xplot.list.n)
  yplot2.lo  <- rep(NA, times=xplot.list.n)
  yplot2.hi  <- rep(NA, times=xplot.list.n)

  #  Calculate mean, median ... per weekday and hour
  for (iplot in 1:xplot.list.n)
  {
    ok               <- (xplot == xplot.list[iplot])
    yplot.n[iplot]   <- sum(ok)
    yplot.mea[iplot] <- mean(dataset[ok, spalte.w], na.rm=TRUE)
    yplot.sd[iplot]  <-   sd(dataset[ok, spalte.w], na.rm=TRUE)   
    yplot.med[iplot] <- median(dataset[ok, spalte.w], na.rm=TRUE)
    yplot.mad[iplot] <- median(abs(dataset[ok, spalte.w]-yplot.med[iplot]),
                               na.rm=TRUE)
  }
  yplot.sd[is.na(yplot.sd)] <- 0   # if only 1 value

  #  If only 1 observation: set sd and MAD to 0
  n.eq.1     <- yplot.n == 1
  if (sum(n.eq.1) > 0)
  {
    yplot.sd[n.eq.1]  <- 0
    yplot.mad[n.eq.1] <- 0
  }
  yplot1.lo  <- yplot.mea - 1.96*yplot.sd/sqrt(yplot.n)
  yplot1.hi  <- yplot.mea + 1.96*yplot.sd/sqrt(yplot.n)
  yplot2.lo  <- yplot.med - 1.96*yplot.mad/sqrt(yplot.n)
  yplot2.hi  <- yplot.med + 1.96*yplot.mad/sqrt(yplot.n)

  #data.frame(xplot.list, yplot.n, yplot1.lo, yplot.mea, yplot1.hi)
  
  if (print.daytime.contours)
  {
    cat("\n\nMeasurand       ",xlabel,
        "\nMean   values with SD per weekday and time of day",
        "\nMedian values with median absolute deviation per weekday and time of day",
        "\nDateTime = 0: Monday, 00:00",
        "\n-----------------------------------------------------",
        "\n")
    print(data.frame(DateTime=xplot.list, Count=yplot.n, 
                    Mean=yplot.mea,
                   SD=yplot.sd, 
                   Median=yplot.med, 
                   MAD=yplot.mad))
    cat("\n-------------------------------------------------------\n")
  }

  # .......................................................................
  #  Plot mean with 95% CI per day and hour
  #  Colors: purple, powderblue

  xplotmin <- 0
  xplotmax <- 7
  yplotmin <- max(0.98*min(yplot1.lo), 0)
  yplotmax <- 1.02*max(yplot1.hi, na.rm=TRUE)

  if (plot.fig095.010)
  { 
    dev.set(fig095.010)
    par(tcl=par.tcl)
    par(las=par.las)

    par(xaxt="n")
    plot(xplot.list, yplot.mea, 
       type="s", lwd=1, lty=1, col="purple",
       xlim=c(xplotmin,xplotmax), ylim=c(yplotmin,yplotmax),
       xlab="Time in week", ylab=xlabel, 
       main="Mean + CI  per weekday and time of day",
       sub=paste(date()," / ",outname," / n =",x.n),cex.sub=0.7)

    par(xaxt="s")

    axis(1,at=seq(0,7,by=0.25), labels=FALSE, col="gray50" )
    axis(1,at=0:7, labels=c("Mo 00:00", "Tu 00:00", "We 00:00", 
       "Th 00:00", "Fr 00:00", "Sa 00:00", "Su 00:00", " "), cex.axis=0.8 )
    lines(c(tail(xplot.list,1), xplotmax),  
        c(tail(yplot.mea,1), tail(yplot.mea,1)),lwd=1, lty=1, col="purple")

    xplot1.poly <- rep(c(xplot.list,xplotmax), each=2) 
    xplot1.poly <- xplot1.poly[-1]
    xplot1.poly <- xplot1.poly[-length(xplot1.poly)]

    yplot1.lo.poly <- rep(yplot1.lo, each=2)
    yplot1.hi.poly <- rev(rep(yplot1.hi, each=2))

    polygon(c(xplot1.poly, rev(xplot1.poly)),
          c(yplot1.lo.poly, yplot1.hi.poly),
          border="powderblue",col="powderblue")

    lines(xplot.list, yplot.mea, type="s", lwd=1, lty=1, col="purple")
    lines(c(tail(xplot.list,1), xplotmax),  
        c(tail(yplot.mea,1), tail(yplot.mea,1)),lwd=1, lty=1, col="purple")
    # abline(v=seq(0,7,by=0.25),lty=2, col="gray50")
    abline(v=seq(0,7, by=1), col="gray50")

    savePlot(file=fig095.010.file,type=figtype)

  }  #  if (plot ...

  # .......................................................................
  #  Plot median with 2*mad per day and hour
  #  Original colors: darkgreen, darkolivegreen1
  #  Colors now: coral4, darkgoldenrod1

  yplotmin <- max(0.98*min(yplot2.lo), 0)
  yplotmax <- 1.02*max(yplot2.hi, na.rm=TRUE)

  if (plot.fig095.020)
  {   
    dev.set(fig095.020)
    par(tcl=par.tcl)
    par(las=par.las)

    par(xaxt="n")
    plot(xplot.list, yplot.med, type="s", lwd=1, lty=1, col="coral4",
       xlim=c(xplotmin,xplotmax), ylim=c(yplotmin,yplotmax),
       xlab="Time in week", ylab=xlabel, 
       main="Median +/- 2*MAD/sqrt(n)  per weekday and time of day",
       sub=paste(date()," / ",outname," / n =",x.n),cex.sub=0.7)

    par(xaxt="s")
    axis(1,at=seq(0,7,by=0.25), labels=FALSE, col="gray50" )
    axis(1,at=0:7, labels=c("Mo 00:00", "Tu 00:00", "We 00:00", 
       "Th 00:00", "Fr 00:00", "Sa 00:00", "Su 00:00", " "), cex.axis=0.8 )
    lines(c(tail(xplot.list,1), xplotmax),  
        c(tail(yplot.med,1), tail(yplot.med,1)),lwd=1, lty=1, col="coral4")
  
    yplot2.lo.poly <- rep(yplot2.lo, each=2)
    yplot2.hi.poly <- rev(rep(yplot2.hi, each=2))

    polygon(c(xplot1.poly, rev(xplot1.poly)),
          c(yplot2.lo.poly, yplot2.hi.poly),
          border="darkgoldenrod1",col="darkgoldenrod1")
    lines(xplot.list, yplot.med, type="s", lwd=1, lty=1, col="coral4")
    lines(c(tail(xplot.list,1), xplotmax),  
        c(tail(yplot.med,1), tail(yplot.med,1)),lwd=1, lty=1, col="coral4")

    # abline(v=seq(0,7,by=0.25),lty=2, col="gray50")
    abline(v=seq(0,7, by=1), col="gray50")

    savePlot(file=fig095.020.file,type=figtype)
  
  }  #  if (plot ...
 
  # .......................................................................
  #  Plot mean +/- 1.96*sd,   per hour, collapsed from Mo to Fr only
  #  Colors: purple, powderblue

  xplot3        <- dataset[ ,"ana.hour"]
  mo.fr         <- dataset[ ,"ana.wday"] < 5
  xplot3        <- xplot3[mo.fr]
  xplot3.list   <- sort(unique(xplot3))
  xplot3.list.n <- length(xplot3.list)         # 
  yplot3        <- dataset[mo.fr, spalte.w]    # 640

  if (plot.fig095.030)
  { 
    #  Go through xplot list and calculate statistics
    yplot3.mea <- rep(NA, times=xplot3.list.n)
    yplot3.sd  <- rep(NA, times=xplot3.list.n)
    yplot3.med <- rep(NA, times=xplot3.list.n)
    yplot3.mad <- rep(NA, times=xplot3.list.n)
    yplot3.n   <- rep(NA, times=xplot3.list.n)
    yplot3.lo  <- rep(NA, times=xplot3.list.n)
    yplot3.hi  <- rep(NA, times=xplot3.list.n)
    yplot4.lo  <- rep(NA, times=xplot3.list.n)
    yplot4.hi  <- rep(NA, times=xplot3.list.n)

    for (iplot in 1:xplot3.list.n)
    {
      ok <- xplot3 == xplot3.list[iplot]
      yplot3.mea[iplot] <- mean(yplot3[ok], na.rm=TRUE)
      yplot3.sd[iplot]  <- sd(yplot3[ok], na.rm=TRUE)
      if (is.na(yplot3.sd[iplot])) { yplot3.sd[iplot] <- 0 }

      yplot3.med[iplot] <- median(yplot3[ok], na.rm=TRUE)
      yplot3.mad[iplot] <- median(abs(yplot3[ok]-yplot3.med[iplot]),
                               na.rm=TRUE)
      yplot3.n[iplot]   <- sum(ok)
      yplot3.lo[iplot]  <- yplot3.mea[iplot] - 1.96*yplot3.sd[iplot]/sqrt(yplot3.n[iplot])
      yplot3.hi[iplot]  <- yplot3.mea[iplot] + 1.96*yplot3.sd[iplot]/sqrt(yplot3.n[iplot])
      yplot4.lo[iplot]  <- yplot3.med[iplot] - 2*yplot3.mad[iplot]/sqrt(yplot3.n[iplot])
      yplot4.hi[iplot]  <- yplot3.med[iplot] + 2*yplot3.mad[iplot]/sqrt(yplot3.n[iplot])
    }

    # Write table to tab/2_file/...
    tab095.3.file <- paste(path.tab.file,outname,"-Val_DayTime.txt", sep="")
     
    sink(file=tab095.3.file, append=FALSE, split=print.daytime.contours)
 
    cat("\n\nMeasurand      ", xlabel,
        "\nData file        ", outname,
        "\nAnalysis date    ", format(date()),
        "\nTable stored in  ",tab095.3.file, 
        "\nMean   values with SD and 95% CI, per time of day",
        "\nMedian values with median absolute deviation (MAD) and",
        "\n       median +/- 2*MAD,  per time of day",
        "\nData from Monday to Friday",
        "\n--------------------------------------------------------------",
        "\n")
    print(data.frame(Time=xplot3.list, 
                     Count=yplot3.n,  
                     Mean=formatC(yplot3.mea, format="f", digits=2),
                     SD=formatC(yplot3.sd, format="f", digits=3), 
                     Mean.lo=formatC(yplot3.lo, format="f", digits=2), 
                     Mean.hi=formatC(yplot3.hi, format="f", digits=2), 
                     Median=formatC(yplot3.med, format="f", digits=2), 
                     MAD=formatC(yplot3.mad, format="f", digits=3),
                     Median.lo=formatC(yplot4.lo, format="f", digits=2), 
                     Median.hi=formatC(yplot4.hi, format="f", digits=2) ))
    cat("\n----------------------------------------------------------------\n")

    sink() 

    xplotmin <- 0
    xplotmax <- 24
    yplotmin <- max(0.98*min(yplot3.lo), 0)
    yplotmax <- 1.02*max(yplot3.hi, na.rm=TRUE)

    dev.set(fig095.030)
    par(tcl=par.tcl)
    par(las=par.las)

    par(xaxt="n")
    plot(xplot3.list, yplot3.mea, type="s", lwd=1, lty=1, col="purple",
       xlim=c(xplotmin,xplotmax), ylim=c(yplotmin,yplotmax),
       xlab="Time of day", ylab=xlabel, 
       main="Mean +/- 1.96*SD/sqrt(n)  per time of day (Mo-Fr only)",
       sub=paste(date()," / ",outname," / n =",x.n),cex.sub=0.7)

    par(xaxt="s")
    axis(1,at=seq(0,24,by=1), labels=FALSE, col="gray50" )
    axis(1,at=seq(0,24,by=6), labels=c("00:00", "06:00", "12:00", 
       "18:00", "24:00"), cex.axis=0.8 )
    lines(c(tail(xplot3.list,1), tail(xplot3.list.n,1)),  
        c(tail(yplot3.mea,1), tail(yplot3.mea,1)),lwd=1, lty=1, col="purple")

    xplot3.poly <- rep(c(xplot3.list,xplot3.list.n), each=2)
    xplot3.poly <- xplot3.poly[-1]
    xplot3.poly <- xplot3.poly[-length(xplot3.poly)]

    yplot3.lo.poly <- rep(yplot3.lo, each=2)
    yplot3.hi.poly <- rev(rep(yplot3.hi, each=2))
    length(yplot3.lo.poly)
    length(yplot3.hi.poly)

    polygon(c(xplot3.poly, rev(xplot3.poly)),
          c(yplot3.lo.poly, yplot3.hi.poly),
          border="powderblue",col="powderblue")

    lines(xplot3.list, yplot3.mea, type="s", lwd=1, lty=1, col="purple")
    lines(c(tail(xplot3.list,1), tail(xplot3.list.n,1)),  
        c(tail(yplot3.mea,1), tail(yplot3.mea,1)),lwd=1, lty=1, col="purple")
    # abline(v=seq(0,7,by=0.25),lty=2, col="gray50")
    abline(v=seq(0,24, by=6), col="gray50")

    savePlot(file=fig095.030.file,type=figtype) 
  }  #  if (plot ...

  # .......................................................................
  #  Plot median with 2*mad  per hour, collapsed from Mo to Fr only
  #  Original colors: darkgreen, darkolivegreen1
  #  Nw: coral4, darkgoldenrod1

  xplotmin <- 0
  xplotmax <- 24
  yplotmin <- max(0.98*min(yplot4.lo), 0)
  yplotmax <- 1.02*max(yplot4.hi, na.rm=TRUE)

  if (plot.fig095.040)
  { 
    dev.set(fig095.040)
    par(tcl=par.tcl)
    par(las=par.las)
    par(xaxt="n")
    plot(xplot3.list, yplot3.med, type="s", lwd=1, lty=1, col="coral4",
       xlim=c(xplotmin,xplotmax), ylim=c(yplotmin,yplotmax),
       xlab="Time of day (Mo-Fr)", ylab=xlabel, 
       main="Median +/- 2*MAD/sqrt(n), per time of day, Mo-Fr only",
       sub=paste(date()," / ",outname," / n =",x.n),cex.sub=0.7)

    par(xaxt="s")
    axis(1,at=seq(0,24,by=1), labels=FALSE, col="gray50" )
    axis(1,at=seq(0,24,by=6), labels=c("00:00", "06:00", "12:00", 
       "18:00", "24:00"), cex.axis=0.8 )
    lines(c(tail(xplot3.list,1), tail(xplot3.list.n,1)),  
        c(tail(yplot3.med,1), tail(yplot3.med,1)),lwd=1, lty=1, col="coral4")

    xplot3.poly <- rep(c(xplot3.list,xplot3.list.n), each=2) 
    xplot3.poly <- xplot3.poly[-1]
    xplot3.poly <- xplot3.poly[-length(xplot3.poly)]

    yplot4.lo.poly <- rep(yplot4.lo, each=2)
    yplot4.hi.poly <- rev(rep(yplot4.hi, each=2))
  
    polygon(c(xplot3.poly, rev(xplot3.poly)),
          c(yplot4.lo.poly, yplot4.hi.poly),
          border="darkgoldenrod1",col="darkgoldenrod1")

    lines(xplot3.list, yplot3.med, type="s", lwd=1, lty=1, col="coral4")
    lines(c(tail(xplot3.list,1), tail(xplot3.list.n,1)),  
        c(tail(yplot3.med,1), tail(yplot3.med,1)),lwd=1, lty=1, col="coral4")
    # abline(v=seq(0,7,by=0.25),lty=2, col="gray50")
    abline(v=seq(0,24, by=6), col="gray50")

    savePlot(file=fig095.040.file,type=figtype)
  }  # if (plot...
}  

# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg095_Plot_DateTime_Val.R  End\n") }
# ===========================================================================


