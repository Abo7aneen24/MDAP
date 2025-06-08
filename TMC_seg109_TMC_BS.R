#  TMC_seg109_TMC_BS.R

  if (go.on & (nbs > 0))
  {  
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 910\n") }

    #  Generate plot windows, if requested
    if (plot.fig100.115)          # Criteria for prop.tmc selection when 
                                  # bootstrapping
    { 
      if (!exists("fig100.115"))
      { #  Open graph new window
        win.graph()
        fig100.115 <- dev.cur()
      }
      fig100.115.file <- paste(path.fig.stra,outname.stra,"-F100.115.",
                             figtype,sep="")
      if (file.exists(fig100.115.file)) { file.remove(fig100.115.file) }
    }

    if (plot.fig100.125)          # Estimated RLs vs prop.tmc when bootstrapping
    { if (!exists("fig100.125"))
      { #  Open graph new window
        win.graph()
        fig100.125 <- dev.cur()
      }
      fig100.125.file <- paste(path.fig.stra,outname.stra,"-F100.125.",
                               figtype,sep="")
      if (file.exists(fig100.125.file)) { file.remove(fig100.125.file) }
    }

    if (plot.fig100.130)          #  Convergence of RL estimates when bootstrapping 
    { if (!exists("fig100.130"))
      { #  Open graph new window
        win.graph()
        fig100.130 <- dev.cur()
      }
      fig100.130.file <- paste(path.fig.stra,outname.stra,"-F100.130.",
                               figtype,sep="")
      if (file.exists(fig100.130.file)) { file.remove(fig100.130.file) }
    }
     
    if (plot.fig100.140)          #  Bootstrapped histograms and densities

    { if (!exists("fig100.140"))
      { #  Open graph new window
        win.graph()
        fig100.140 <- dev.cur()
      }
      fig100.140.file <- paste(path.fig.stra,outname.stra,"-F100.140.",
                               figtype,sep="")
      if (file.exists(fig100.140.file)) { file.remove(fig100.140.file) }
    }
     
    if (plot.fig100.145)          #  Densities of bootstrapped results
    { if (!exists("fig100.145"))
      { #  Open graph new window
        win.graph()
        fig100.145 <- dev.cur()
      }
      fig100.145.file <- paste(path.fig.stra,outname.stra,"-F100.145.",
                               figtype,sep="")
      if (file.exists(fig100.145.file)) { file.remove(fig100.145.file) }
    }
     
    #  Berechnung von Toleranzintervallen per Bootstrap ist erwünscht
    #  Tabelle für Einzelergebnisse anlegen
    tab.bs.names     <- c("lambda","mue","sigma","RL1","RL2","prev",
                          "prop", "x.tr.lo", "x.tr.hi", "p.fit")
    tab.bs           <- matrix(NA,nrow=nbs,ncol=length(tab.bs.names))
    colnames(tab.bs) <- tab.bs.names

    theta.ini <- c(lambda.tmc, mue.tmc, sigma.tmc)

    if (plot.fig100.140)
    { dev.set(fig100.140)
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(x), col="black", lwd=3,
           xlim=c(0,40),
           main="Data and BS densities",
           sub=subtitle, cex.sub=0.7)
    }

    for (ibs in 1:nbs)
    {
        t1 <- Sys.time() 

        #  Bootstrap sample herstellen
  
        xbs <- sample(x, x.n, replace=TRUE)   #  empirical BS

        #xbs <- rxref(x.kde.pdf$x,x.cdf,x.n) # nichtparametrischer bootstrap,
        #                                # beruhend auf der Dichteschätzung der 
        #                                # Daten

        xbs <- sort(xbs)
        #save(xbs, file="SmoothHistTest.RData")
        #save(x.breaks, file="SmoothHistBreaks.RData")

        #cat("   ibs:",ibs, "mean(x) ", mean(x), " mean(xbs) ", mean(xbs), "\n")
        cat("   ibs:",ibs, "mean(x) ", mean(x), " mean(xbs) ", mean(xbs), "\x0D")
        #cat("   ibs:",ibs,  "\x0D")
         
        #  Generate histogram and density estimate for the xbs sample
        #  Parameter 'freq' not used here
        xbs.hist <- hist(xbs,breaks=x.breaks, right=FALSE,
                         plot=FALSE,
                         warn.unused=FALSE)
        xbs.kde  <- density(xbs,n=kernel.n,adjust=x.adjust,from=x.from,to=x.to)

        #  Smooth the histogram
        xbs.hist.smo <- SmoothHist(xbs.hist, fig=FALSE)

        save(xbs.hist, file="xbs.hist.RData")
        save(xbs.hist.smo, file="xbs.hist.smo.RData")

        if (plot.fig100.140)
        { dev.set(fig100.140)       
          lines(xbs.kde, col="green3", lwd=1)
          par(las=par.las)    
          par(tcl=par.tcl)    
          plot(xbs.hist, add=TRUE, border="cyan")
        } 

        #  Find mode from kde
        xbs.kde.mode.idx <- which.max(xbs.kde$y)

        #  If there is more than 1 maximum: message
        if (length(xbs.kde.mode.idx) > 1)
        {
          cat(
               "\n\n +++ [Analysis] BS density estimate has > 1 maximum - first is taken\n\n" )
          xbs.kde.mode.idx <- xbs.kde.mode.idx[1] 
        } 
        xbs.kde.mode     <- xbs.kde$x[xbs.kde.mode.idx]

        xbs.breaks.mode.idx <- which((x.breaks[1:(x.breaks.n-1)] <= xbs.kde.mode) &
                                     (xbs.kde.mode < x.breaks[2:x.breaks.n] ) ) 

        #  Plot history: show used time
        #t2 <- Sys.time() 
        #delta21 <- difftime(t2,t1)
        ##  Verbrauchte Zeit darstellen
        #if (ibs==1)
        #{ xplotmin11 <- 1
        #  xplotmax11 <- nbs
        #  yplotmin11 <- 0      
        #  yplotmax11 <- 10 * as.numeric(delta21)

        #  plot(ibs,delta21,type="h",col="blue",
        #         xlim=c(xplotmin11,xplotmax11),ylim=c(yplotmin11,yplotmax11),
        #         main="Bootstrap progress",
        #         xlab="# of BS sample",ylab="Time used")
        #}  else
        #{  
        #  lines(ibs,delta21,type="h",col="blue")
        #}

        #  Free choice of prop.tmc produces extreme fluctuation of the 
        #  truncation intervel and thus in the estimates.   
        TMC.bs <- tmc.master(df.est, df.con, ErrMsg, figtype, 
                       outname.stra, p.fact, path.fig.stra,
                       pathol.position, 
                       plot.fig100.115, plot.fig100.125,
                       fig100.115, fig100.125,
                       fig100.115.file, fig100.125.file,   
                       prop.tmc.seq0, prop.tmc.seq.n,
                       RL1.p, RL2.p, theta.ini, w.fact,
                       xbs.breaks.mode.idx, xbs.hist.smo, xbs.kde,  
                       x.n, x.n.min, x.tr.bins.min,
                       print.log.message, fastnull, fastnull.chi)

        tab.bs[ibs,"lambda"] <- TMC.bs[["est"]]["lambda"]
        tab.bs[ibs,"mue"]    <- TMC.bs[["est"]]["mue"]
        tab.bs[ibs,"sigma"]  <- TMC.bs[["est"]]["sigma"]
        tab.bs[ibs,"RL1"]    <- TMC.bs[["est"]]["x.RL1"]
        tab.bs[ibs,"RL2"]    <- TMC.bs[["est"]]["x.RL2"]
        tab.bs[ibs,"prev"]   <- TMC.bs[["fit"]]["prev.tmc"]

        tab.bs[ibs,"prop"]      <- TMC.bs[["est"]]["prop.tmc"]
        tab.bs[ibs,"x.tr.lo"]   <- TMC.bs[["est"]]["x.tr.lo"]
        tab.bs[ibs,"x.tr.hi"]   <- TMC.bs[["est"]]["x.tr.hi"]
        tab.bs[ibs,"p.fit"]     <- TMC.bs[["est"]]["p.fit"]

        # show estimated density 
        lines(x.kde.pdf$x, (1-tab.bs[ibs,"prev"])*
                       pdf.PN(x.kde.pdf$x, tab.bs[ibs,"lambda"],
                                       tab.bs[ibs,"mue"],
                                       tab.bs[ibs,"sigma"]), 
              col=boscol)

        #t3 <- Sys.time()
        #delta32 <- difftime(t3,t2)
        #print(delta32)
        #lines(c(ibs,ibs),c(delta21,delta21+delta32),col="green3")
        #cat("   ibs:",ibs,  difftime(t3,t1),"\x0D")

        # Alternative to elapsed time: show convergence of RLs
        if (plot.fig100.130)
        { dev.set(fig100.130)
          if (ibs==1)
          { xplotmin11 <- 1
            xplotmax11 <- nbs
            yplotmin11 <- 0      
            yplotmax11 <- 2 * tab.bs[ibs,"RL2"]

            par(las=par.las)    
            par(tcl=par.tcl)    
            plot(ibs,tab.bs[ibs,"RL1"],type="p",col=boscol, pch=3, 
               xlim=c(xplotmin11,xplotmax11),ylim=c(yplotmin11,yplotmax11),
               main="Bootstrap progress",
               xlab="BS sample #",ylab="Mean RIL",
               sub=subtitle, cex.sub=0.7)
            points(ibs,tab.bs[ibs,"RL2"],type="p",col=boscol, pch=3)  
          }  else
          {  
            points(ibs,tab.bs[ibs, "RL1"],type="p",col="gray", pch=20)
            points(ibs,tab.bs[ibs, "RL2"],type="p",col="gray", pch=20)
            points(ibs,median(tab.bs[1:ibs, "RL1"]),type="p",col=boscol, pch=3)
            points(ibs,median(tab.bs[1:ibs, "RL2"]),type="p",col=boscol, pch=3)
          }
        }

    }  # if (plot ...
    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 920\n") }

    if ((nbs > 0) & plot.fig100.130)
    { #  Add tmc estimate
      dev.set(fig100.130)
      abline(h=x.RL1.tmc, col=tmccol)
      abline(h=x.RL2.tmc, col=tmccol)
      abline(h=TMC[["ci"]][3], col=tmccol, lty=3)
      abline(h=TMC[["ci"]][4], col=tmccol, lty=3)

      abline(h=median(tab.bs[1:nbs, "RL1"]), col=boscol)
      abline(h=median(tab.bs[1:nbs, "RL2"]), col=boscol)
      abline(h=Quantile(tab.bs[1:nbs, "RL2"], probs=0.025), col=boscol, lty=3)
      abline(h=Quantile(tab.bs[1:nbs, "RL2"], probs=0.975), col=boscol, lty=3)
    }
    savePlot(file=fig100.130.file, type=figtype)

    #cat("Single nonparametric bootstrap results for tolerance interval\n") 
    #print(tab.bs[order(tab.bs[ ,"RL1"]), ])

    if (plot.fig100.145)
    { dev.set(fig100.145)
      par(mfrow=c(2,2))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(x),main="Estimated density of x")
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"RL1"]),main="BS density of RL1",
         sub=paste("nbs =",nbs))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"RL2"]),main="BS density of RL2",
         sub=paste("nbs =",nbs))
      par(las=par.las)    
      par(tcl=par.tcl)    
      plot(density(tab.bs[ ,"prev"]),main="BS density of prev",
         sub=paste("nbs =",nbs))
    }  # if (plot...
    
    x.RL.est["tmc.tmc","Prevtilo"] <- Quantile(tab.bs[ ,"prev"],probs=0.025) 
    x.RL.est["tmc.tmc","Prevtihi"] <- Quantile(tab.bs[ ,"prev"],probs=0.975) 
    x.RL.est["tmc.tmc","x.RL1bsmea"] <- mean(tab.bs[ ,"RL1"]) 
    x.RL.est["tmc.tmc","x.RL1tilo"] <- Quantile(tab.bs[ ,"RL1"],probs=0.025) 
    x.RL.est["tmc.tmc","x.RL1tihi"] <- Quantile(tab.bs[ ,"RL1"],probs=0.975) 
    x.RL.est["tmc.tmc","x.RL2bsmea"] <- mean(tab.bs[ ,"RL2"]) 
    x.RL.est["tmc.tmc","x.RL2tilo"] <- Quantile(tab.bs[ ,"RL2"],probs=0.025) 
    x.RL.est["tmc.tmc","x.RL2tihi"] <- Quantile(tab.bs[ ,"RL2"],probs=0.975) 
 
    #cat("[seg100] Bootstrap-Ergebnis für tmc\n") 
    #print(x.RL.est["tmc.tmc", ])
  
    #  Nutzung der CIs in SplineVerlauf()

    if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 930\n") }
  }     #   if (go.on & nbs > 0) 

