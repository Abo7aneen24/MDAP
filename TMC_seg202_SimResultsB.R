#  TMC_seg202_SimResultsB.R
#
#  Produces summarizing figure B for simulation data 
#  Is called by TMC_seg202_SimResultsA.R

#  17.04.2021 sum.reps renamed to scen.sum
#  15.12.2020 Output file renamed to 202.020
#  11.12.2020 Start
# =============================================================================

if (print.log.message) { cat("%%%   TMC_seg202_SimResultsB  Start\n") }

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.020) & ("npa" %in% meth.list.202))
{   
  dev.set(fig202.020)
  bringToTop(fig202.020)

  i.ifc  <- which(meth.list.202 == "npa")
  ok.ifc <- stab1[ ,"method"]==meth.list.202[i.ifc]
  RL1.ifc <- as.numeric(stab1[ok.ifc ,"RL1"])
  RL2.ifc <- as.numeric(stab1[ok.ifc ,"RL2"])

  RL1.diff.pct <- as.numeric(stab1[ ,"RL1"])
  RL2.diff.pct <- as.numeric(stab1[ ,"RL2"])

  for (imeth in 1:meth.list.202.n)
  { 
    ok <- stab1[ ,"method"]==meth.list.202[imeth]
    RL1.diff.pct[ok] <- 100*(RL1.diff.pct[ok] - RL1.ifc)/RL1.ifc
    RL2.diff.pct[ok] <- 100*(RL2.diff.pct[ok] - RL2.ifc)/RL2.ifc
  }
  
  #cat("[seg201] Kontrolle der Differenzberechnung\n")
  #print(data.frame(stab1[ ,c("method","RL1","RL2")],
  #               RL1.diff.pct,RL2.diff.pct)) 

  xplotmin <- min(RL1.diff.pct,na.rm=TRUE)
  xplotmax <- max(RL1.diff.pct,na.rm=TRUE)
  yplotmin <- min(RL2.diff.pct,na.rm=TRUE)
  yplotmax <- max(RL2.diff.pct,na.rm=TRUE)

  for (imeth in 1:meth.list.202.n)
  {
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    meth.col <- meth.dsg[which(meth.dsg[ ,1]==meth) ,2]
    meth.pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1]==meth) ,3])

    RL1.diff.pct.mea  <- mean(RL1.diff.pct[ok],na.rm=TRUE)
    RL1.diff.pct.sd   <- sd(RL1.diff.pct[ok],na.rm=TRUE)
    RL1.diff.pct.025  <- Quantile(RL1.diff.pct[ok],probs=0.025)
    RL1.diff.pct.med  <- Quantile(RL1.diff.pct[ok],probs=0.500)
    RL1.diff.pct.975  <- Quantile(RL1.diff.pct[ok],probs=0.975)
    RL2.diff.pct.mea  <- mean(RL2.diff.pct[ok],na.rm=TRUE)
    RL2.diff.pct.sd   <- sd(RL2.diff.pct[ok],na.rm=TRUE)
    RL12.diff.pct.cov <- cov(RL1.diff.pct[ok], RL2.diff.pct[ok])
  
    if (imeth == 1)
    {
      plot(RL1.diff.pct[ok],RL2.diff.pct[ok],type="p",
           col=meth.col,pch=meth.pch,
           xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax),
           xlab="100*(Estimate RL1 - RL1(IFCC))/RL1(IFCC)",
           ylab="100*(Estimate RL2 - RL2(IFCC))/RL2(IFCC)",
           sub=subtitle.202,
           cex.sub=0.8)

        points(RL1.diff.pct.mea,RL2.diff.pct.mea,
               type="p", col="black", pch=19, cex=1.5)
      #  Write parameter information into the plot
      ok.info <- Infoblock(RunId, yplotmin, yplotmax, 
                     outname.file, x.n1, subsample.n, oh.sym, dev.sym, 
                     round.unit, smooth.hist1, smooth.hist2, 
                     n.per.bin.min, bins.n.min, bins.n.max,  
                     lambda.min, lambda.max,   
                     l.fact, p.fact, r.fact, s.fact, w.fact, 
                     x.tr.prop.min, x.tr.prop.max, p.fit.min)
    }  else
    {
      points(RL1.diff.pct[ok],RL2.diff.pct[ok],type="p",
             col=meth.col,pch=meth.pch)
      points(RL1.diff.pct.mea,RL2.diff.pct.mea,
             type="p",col=meth.col,pch=19,cex=1.5)
    }

    #  Confidence ellipse for (RL1.diff.pct, RL2.diff.pct) if meth != ifc
    if (imeth != i.ifc)
    { cov.elli <- PlotConfElli(c(RL1.diff.pct.mea,RL2.diff.pct.mea),
                    matrix(c(RL1.diff.pct.sd^2,-RL12.diff.pct.cov,
                             -RL12.diff.pct.cov,RL2.diff.pct.sd^2),
                           byrow=TRUE,ncol=2),
                           0.05,NA)
      lines(cov.elli,type="l",col=meth.col)
    }
  }

  legend("topright",meth.dsg[meth.list.202.num,1],
          col=meth.dsg[meth.list.202.num,2],
          pch=as.numeric(meth.dsg[meth.list.202.num,3]), lwd=3, cex=0.8)

  abline(v=0,col=gencol)
  abline(h=0,col=gencol)

  savePlot(file=fig202.020.file,type=figtype)
}

if (print.log.message) { cat("%%%   TMC_seg202_SimResultsB  End\n") }
  

  