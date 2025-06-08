#  TMC_seg202_SimResultsA.R
#
#  Produces summarizing figure A for simulation data
#  Is called by TMC_seg202_SimResults.R

#  17.04.2021 sum.reps renamed to scen.sum
#  05.04.2021 Error in processing RL1.clip.min corrected
#  15.12.2020 Output file renamed to 202.010
#  11.12.2020 Start
# =============================================================================

if (print.log.message) { cat("%%%   TMC_seg202_SimResultsA  Start\n") }

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.010))
{
  # ---------------------------------------------------------------------------
  #  Show values of RL1 and RL2

  dev.set(fig202.010)
  bringToTop(fig202.010)
  
  if (is.na(RL1.clip.min))
  { #  No user specification (if present, must be present for all axes)
    xplotmin <- min(as.numeric(stab1[ ,"RL1"]),na.rm=TRUE)
    xplotmax <- max(as.numeric(stab1[ ,"RL1"]),na.rm=TRUE)
    yplotmin <- min(as.numeric(stab1[ ,"RL2"]),na.rm=TRUE)
    yplotmax <- max(as.numeric(stab1[ ,"RL2"]),na.rm=TRUE)
  }  else
  {
    xplotmin <- RL1.clip.min     
    xplotmax <- RL1.clip.max     
    yplotmin <- RL2.clip.min     
    yplotmax <- RL2.clip.max     
  }

  for (imeth in 1:meth.list.202.n)
  {
    meth     <- meth.list.202[imeth]
    ok       <- stab1[ ,"method"]==meth
    meth.col <- meth.dsg[which(meth.dsg[ ,1]==meth) ,2]
    meth.pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1]==meth) ,3])
        
    # Extract sub-dataset for meth for easier referencing
    stab1.sub <- data.frame(
                           sol.type=stab1[ok,"rc"],
                           errcode=stab1[ok,"errcode"],
                           iter=as.numeric(stab1[ok,"iter"]),
                           RL1=as.numeric(stab1[ok,"RL1"]),
                           RL2=as.numeric(stab1[ok,"RL2"]), 
                           prev.l=as.numeric(stab1[ok,"prev.l"]), 
                           prev.c=as.numeric(stab1[ok,"prev.c"]), 
                           prev.r=as.numeric(stab1[ok,"prev.r"]), 
                           lambda=as.numeric(stab1[ok,"lambda"]), 
                           mue=as.numeric(stab1[ok,"mue"]), 
                           sigma=as.numeric(stab1[ok,"sigma"]),
                           n.per.bin=as.numeric(stab1[ok,"n.per.bin.min.eff"]),  
                           prop=as.numeric(stab1[ok,"prop"]),  
                           p.fit=as.numeric(stab1[ok,"p.fit"]) ) 
      cat("\n  [202A] imeth = ", imeth, meth,  "  stab1.sub \n")
      print(stab1.sub)

      if (imeth == 1)
      {
        #  Establish the plot
        #  No axes if  RL1.clip, RL2.clip are specified
        if (!is.na(RL1.clip.min))
        { 
          par(xaxt="n") 
          par(yaxt="n")
        }
        plot(stab1.sub[ ,"RL1"],stab1.sub[ ,"RL2"],type="p",
           col=meth.col,pch=meth.pch,
           xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax),
           xlab="Estimated RL1",ylab="Estimated RL2",
           sub=subtitle.202,
           cex.sub=0.8)

        if (!is.na(RL1.clip.min))
        { # more axis details according to user specification
          par(xaxt="s") 
          par(yaxt="s")

          axis(1, at=seq(RL1.clip.min, RL1.clip.max, by=RL1.clip.by1) )  
          axis(1, at=seq(RL1.clip.min, RL1.clip.max, by=RL1.clip.by2), 
               labels=FALSE )  
          axis(2, at=seq(RL2.clip.min, RL2.clip.max, by=RL2.clip.by1) )  
          axis(2, at=seq(RL2.clip.min, RL2.clip.max, by=RL2.clip.by2), 
               labels=FALSE )            
        }

        #  Write parameter information into the plot
        ok.info <- Infoblock(RunId, yplotmin, yplotmax, 
                     outname.file, x.n1, subsample.n, oh.sym, dev.sym, 
                     round.unit, smooth.hist1, smooth.hist2, 
                     n.per.bin.min, bins.n.min, bins.n.max,
                     lambda.min, lambda.max,   
                     l.fact, p.fact, r.fact, s.fact, w.fact, 
                     x.tr.prop.min, x.tr.prop.max, p.fit.min)

      } else    # not first method in the list
      {
        points(stab1.sub[ ,"RL1"],stab1.sub[ ,"RL2"],type="p",
             col=meth.col,pch=meth.pch)
      }
 
      #  Calculate statistics

      #  For each method
      scen.sum[meth,"RL1.mw"]   <- mean(stab1.sub[ ,"RL1"],na.rm=TRUE) 
      scen.sum[meth,"RL1.med"]  <- median(stab1.sub[ ,"RL1"],na.rm=TRUE) 
      scen.sum[meth,"RL1.sd"]   <-   sd(stab1.sub[ ,"RL1"],na.rm=TRUE) 
      scen.sum[meth,"RL1.min"]  <-  min(stab1.sub[ ,"RL1"],na.rm=TRUE) 
      scen.sum[meth,"RL1.max"]  <-  max(stab1.sub[ ,"RL1"],na.rm=TRUE) 
      scen.sum[meth,"RL2.mw"]   <- mean(stab1.sub[ ,"RL2"],na.rm=TRUE) 
      scen.sum[meth,"RL2.med"]  <- median(stab1.sub[ ,"RL2"],na.rm=TRUE) 
      scen.sum[meth,"RL2.sd"]   <-   sd(stab1.sub[ ,"RL2"],na.rm=TRUE)    
      scen.sum[meth,"RL2.min"]  <-  min(stab1.sub[ ,"RL2"],na.rm=TRUE) 
      scen.sum[meth,"RL2.max"]  <-  max(stab1.sub[ ,"RL2"],na.rm=TRUE) 
      scen.sum[meth,"cov.RL1RL2"] <- cov(stab1.sub[ ,"RL1"],
                                       stab1.sub[ ,"RL2"])
      #  Method tmc & tmu only
      if (meth == "tmc" | meth == "tmu")
      { 
        scen.sum[meth,"n.per.bin.min"] <- n.per.bin.min 
  
        scen.sum[meth,"prop.mw"]  <- mean(stab1.sub[ ,"prop"],na.rm=TRUE) 
        scen.sum[meth,"prop.sd"]  <- sd(stab1.sub[ ,"prop"],na.rm=TRUE) 
        scen.sum[meth,"prop.min"] <- min(stab1.sub[ ,"prop"],na.rm=TRUE) 
        scen.sum[meth,"prop.max"] <- max(stab1.sub[ ,"prop"],na.rm=TRUE) 

        scen.sum[meth,"prev.l.mw"]  <- mean(stab1.sub[ ,"prev.l"],na.rm=TRUE) 
        scen.sum[meth,"prev.l.sd"]  <-   sd(stab1.sub[ ,"prev.l"],na.rm=TRUE) 
        scen.sum[meth,"prev.l.min"] <-  min(stab1.sub[ ,"prev.l"],na.rm=TRUE) 
        scen.sum[meth,"prev.l.max"] <-  max(stab1.sub[ ,"prev.l"],na.rm=TRUE) 

        scen.sum[meth,"prev.c.mw"]  <- mean(stab1.sub[ ,"prev.c"],na.rm=TRUE) 
        scen.sum[meth,"prev.c.sd"]  <-   sd(stab1.sub[ ,"prev.c"],na.rm=TRUE) 
        scen.sum[meth,"prev.c.min"] <-  min(stab1.sub[ ,"prev.c"],na.rm=TRUE) 
        scen.sum[meth,"prev.c.max"] <-  max(stab1.sub[ ,"prev.c"],na.rm=TRUE) 

        scen.sum[meth,"prev.r.mw"]  <- mean(stab1.sub[ ,"prev.r"],na.rm=TRUE) 
        scen.sum[meth,"prev.r.sd"]  <-   sd(stab1.sub[ ,"prev.r"],na.rm=TRUE) 
        scen.sum[meth,"prev.r.min"] <-  min(stab1.sub[ ,"prev.r"],na.rm=TRUE) 
        scen.sum[meth,"prev.r.max"] <-  max(stab1.sub[ ,"prev.r"],na.rm=TRUE) 

        scen.sum[meth,"p.fit.mw"]  <-  mean(stab1.sub[ ,"p.fit"],na.rm=TRUE) 
        scen.sum[meth,"p.fit.sd"]  <-   sd(stab1.sub[ ,"p.fit"],na.rm=TRUE) 
        scen.sum[meth,"p.fit.min"] <-  min(stab1.sub[ ,"p.fit"],na.rm=TRUE) 
        scen.sum[meth,"p.fit.max"] <-  max(stab1.sub[ ,"p.fit"],na.rm=TRUE)
      }

      #  Methods tmc & tml 
      if (meth == "tmc" | meth== "tml")
      { 
        scen.sum[meth,"lambda.mw"]  <- mean(stab1.sub[ ,"lambda"],na.rm=TRUE) 
        scen.sum[meth,"lambda.sd"]  <-   sd(stab1.sub[ ,"lambda"],na.rm=TRUE) 
        scen.sum[meth,"lambda.min"] <-  min(stab1.sub[ ,"lambda"],na.rm=TRUE) 
        scen.sum[meth,"lambda.max"] <-  max(stab1.sub[ ,"lambda"],na.rm=TRUE) 

        scen.sum[meth,"mue.mw"]  <- mean(stab1.sub[ ,"mue"],na.rm=TRUE) 
        scen.sum[meth,"mue.sd"]  <-   sd(stab1.sub[ ,"mue"],na.rm=TRUE) 
        scen.sum[meth,"mue.min"] <-  min(stab1.sub[ ,"mue"],na.rm=TRUE) 
        scen.sum[meth,"mue.max"] <-  max(stab1.sub[ ,"mue"],na.rm=TRUE) 

        scen.sum[meth,"sig.mw"]  <- mean(stab1.sub[ ,"sigma"],na.rm=TRUE) 
        scen.sum[meth,"sig.sd"]  <-   sd(stab1.sub[ ,"sigma"],na.rm=TRUE) 
        scen.sum[meth,"sig.min"] <-  min(stab1.sub[ ,"sigma"],na.rm=TRUE) 
        scen.sum[meth,"sig.max"] <-  max(stab1.sub[ ,"sigma"],na.rm=TRUE) 
      }

      if (meth=="npa")
      { 
        points(scen.sum[meth,"RL1.mw"],scen.sum[meth,"RL2.mw"],
               type="p",col="black",pch=19,cex=1.5)
      }  else
      { 
        points(scen.sum[meth,"RL1.mw"],scen.sum[meth,"RL2.mw"],
               type="p",col=meth.col,pch=19,cex=1.5)
      }
      abline(v=scen.sum[meth,"RL1.med"], col=meth.col, lty=2)
      abline(h=scen.sum[meth,"RL2.med"], col=meth.col, lty=2)              
   
      #  Confidence ellipse for (RL1,RL2)
   
      cov.elli <- PlotConfElli(c(scen.sum[meth,"RL1.mw"],scen.sum[meth,"RL2.mw"]),
                  matrix(c(scen.sum[meth,"RL1.sd"]^2,-scen.sum[meth,"cov.RL1RL2"],
                           -scen.sum[meth,"cov.RL1RL2"],scen.sum[meth,"RL2.sd"]^2),
                         byrow=TRUE,ncol=2),
                          0.05,NA)
      lines(cov.elli,type="l",col=meth.col)

      #  Mark extreme results of tmc
      #  Mark invalid results of tmc
      if (meth=="tmc")
      { 
        tmc.RL1.min.idx <- which.min(stab1.sub[ ,"RL1"])
        tmc.RL1.max.idx <- which.max(stab1.sub[ ,"RL1"])
        tmc.RL2.min.idx <- which.min(stab1.sub[ ,"RL2"])
        tmc.RL2.max.idx <- which.max(stab1.sub[ ,"RL2"])
        text(stab1.sub[tmc.RL1.min.idx,"RL1"],stab1.sub[tmc.RL1.min.idx,"RL2"],
           labels=tmc.RL1.min.idx,col=meth.col,cex=0.8,pos=4)  
        text(stab1.sub[tmc.RL1.max.idx,"RL1"],stab1.sub[tmc.RL1.max.idx,"RL2"],
           labels=tmc.RL1.max.idx,col=meth.col,cex=0.8,pos=2)  
        text(stab1.sub[tmc.RL2.min.idx,"RL1"],stab1.sub[tmc.RL2.min.idx,"RL2"],
           labels=tmc.RL2.min.idx,col=meth.col,cex=0.8,pos=3)  
        text(stab1.sub[tmc.RL2.max.idx,"RL1"],stab1.sub[tmc.RL2.max.idx,"RL2"],
           labels=tmc.RL2.max.idx,col=meth.col,cex=0.8,pos=1)  
        
        # invalid <- (stab1.sub[ ,"sol.type"] != 2)
        invalid <- (stab1.sub[ ,"p.fit"] < p.fit.min)
        if (sum(invalid) > 0)
        {
          points(stab1.sub[invalid, "RL1"],stab1.sub[invalid, "RL2"],
                 type="p", col="cyan", pch=meth.pch)
        }
      }
    }

    #  Construct legend. Needs all methods addressed in meth.dsg. 
    meth.list.202.num <- rep(NA, times=length(meth.list.202))
    for (i in 1:meth.list.202.n)
    { 
      meth.list.202.num[i] <- which(meth.dsg[ ,1] == meth.list.202[i])
    }
    legend("topright",meth.dsg[meth.list.202.num,1],col=meth.dsg[meth.list.202.num,2],
          pch=as.numeric(meth.dsg[meth.list.202.num,3]),lwd=3, cex=0.8)

    abline(v=xc.RL1.gen,col=gencol)
    abline(h=xc.RL2.gen,col=gencol)

    savePlot(file=fig202.010.file,type=figtype)
  }  #  if (plot...

if (print.log.message) { cat("%%%   TMC_seg202_SimResultsA  End\n") }
