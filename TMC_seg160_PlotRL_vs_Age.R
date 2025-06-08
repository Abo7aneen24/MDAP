#  TMC_seg160_PlotRL_vs_Age.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Plot RL vs age class, smoothed.
#  Only 1 method
#  Data stratified by M / F  
#  Age position is mean within age class
#  Colors and plot characters mostly taken from style file, with some
#  exceptions at the beginning of this segment
#  Plot scales are fixed in TMC_seg50_master.R

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  09.04.2021 oh.sym, dev.sym added
#  02.12.2020 x.RL1, ... changed to RL1, ..., call for infoblock updated
#  16.07.2020 Common par parameters for plotting
#  27.02.2020 Infoblock replaced by function
#  24.01.2017 No details into plot if no spline was generated
#  08.01.2020 sexcolFMFM, sexpchFMFM specified in seg030_Style
#             Var names modified for TMC4
#  13.05.     Legend, gridlines controlled by switch
#  01.05.     Check for no points in input data
#  18.03.2019 Start from WWW17_seg48  (directory ProgWW17_RH)
# ===========================================================================

if (print.log.message) { cat("%%%   TMC_seg160_PlotRL_vs_Age.R  Start\n") }

# ===========================================================================

source("TMC_seg159_OpenWindows.R")

if (PlotId)
{ 
  #  gtab auf die gefragte Methode und subset.type reduzieren 
  gtab.M1 <- gtab[ (gtab[ ,"method"]==Verlauf.Meth1) &
                   (gtab[ ,"subset.type"]==3)   , ]

  #  Weiter, wenn gtab.M1 nicht leer
  gtab.M1.n <- nrow(gtab.M1)

  if (gtab.M1.n > 1)
  { #  Is a matrix with > 1 rows
     
    #  RL.label == NA: leave construction to RL.age.spline()
    if (RL.age.meth == "tmc" & !is.na(RL.label.tmc) )
    { RL.label <- RL.label.tmc }   #  user-defined label

    #  Find existing sex categories
    sexcode.list0 <- sort(unique(gtab.M1[ ,"Sex"]))
    sexcode.list0.n <- length(sexcode.list0)

    #  Code mit "+" suchen
    sexcode.list  <- "dummy"

    for (i in 1:sexcode.list0.n)
    { mit.plus <- FALSE
      s.code <- sexcode.list0[i] 
      for (j in 1:nchar(s.code))
      {  
        mit.plus <- mit.plus |(substr(s.code,j,j) == "+") 
      }
      if (mit.plus) { sexcode.list <- c(sexcode.list,s.code)}
    }
    sexcode.list   <- sexcode.list[-1]      #  get rid of dummy entry
    sexcode.list.n <- length(sexcode.list)

    #  There should be at most one code containing a '+'
    #  Error message, if not
    if (sexcode.list.n > 1)
    { cat("\n ++++++ [seg160]  More than 1 sex code contains a '+'",
          "\n                  Only the first will be processed",
          "\n")
      print(sexcode.list )
      sexcode.list <- sexcode.list[1]
    }

    group.shift <- 0       # no groups to plot
 
    #  Process the sex code. Plot scale was determined in seg150. 
    isex <- 0
    for (sexcode in sexcode.list)
    { isex <- isex + 1
      ok <- gtab.M1[ ,"Sex"] == sexcode
      gtab.M1.subset <- gtab.M1[ok, ]

      ok <-  !is.na(as.numeric(gtab.M1.subset[ ,"x.RL1"]))

      #  Plot only if at least 2 support points (then plot points)
      #  If >= 4 points, plot splines also
      if (sum(ok) > 1)
      {
        gtab.M1.subset <- gtab.M1.subset[ok, ] 
 
        gtab.M1.spl <- RL.age.spline(PlotId,
                      gtab.M1.subset,
                      sexcode,Verlauf.Meth1,group.shift*isex,
                      age.class[1, "lo"], age.class[age.class.n, "hi"],
                      yplotmin150.010,yplotmax150.010, 
                      yplotby1150.010,yplotby2150.010,   
                      age.label, RL.label, RL.age.meth,
                      age.clip.min, age.clip.max, age.clip.by1, age.clip.by2,
                      par.tcl, par.las,
                      datafile,sexcolFM[isex],sexpchFM[isex],
                      arrow.length=0.04,newplot=TRUE)


        #  Permisible difference is calculated here from the smoothed 
        #  RL-age relation. Not from estimated RL per age class, as in Tab2c

        if (plot.perm.diff.band & gtab.M1.spl$spline.ok)
        {
          #  Permissible difference band
          pd1 <- PermissibleDiff(gtab.M1.spl$RL.smo["RL1.hut"],
                                 gtab.M1.spl$RL.smo["RL2.hut"],
                                 gtab.M1.spl$RL.smo["RL1.hut"],
                                 ef)
          polygon(c(gtab.M1.spl$RL.smo[, "age"], 
                  rev(gtab.M1.spl$RL.smo[, "age"])),
                c(pd1[ ,1], rev(pd1[ ,2])), 
                col=sexcolfill[isex], border=sexcolfill[isex]) 
          lines(gtab.M1.spl$RL.smo[, "age"], pd1[ ,1], col=sexcol[isex],
                lty=2)
          lines(gtab.M1.spl$RL.smo[, "age"], pd1[ ,2], col=sexcol[isex],
                lty=2)

          pd2 <- PermissibleDiff(gtab.M1.spl$RL.smo["RL1.hut"],
                                 gtab.M1.spl$RL.smo["RL2.hut"],
                                 gtab.M1.spl$RL.smo["RL2.hut"],
                                 ef)
          polygon(c(gtab.M1.spl$RL.smo[, "age"], 
                  rev(gtab.M1.spl$RL.smo[, "age"])),
                c(pd2[ ,1], rev(pd2[ ,2])), 
                col=sexcolfill[isex], border=sexcolfill[isex]) 
          lines(gtab.M1.spl$RL.smo[, "age"], pd2[ ,1], col=sexcol[isex],
                lty=2)
          lines(gtab.M1.spl$RL.smo[, "age"], pd2[ ,2], col=sexcol[isex],
                lty=2)
        }

        if (plot.details & gtab.M1.spl$plot.ok) 
        { # Put additional details into the plot (for WW)

          # @@@ replace by function for display of overlapping intervals
          # abline(v=age.limits,col="gray")  # von CCLM nicht erlaubt

          #  @@@ also step below 
          # Re-draw tickmarks at age limits (were overwritten by abline() )
          #par(tcl=par.tcl)
          #par(las=par.las)
          #axis(1,at=age.limits, labels=FALSE)

          if ((sexcode.list.n > 0) & plot.details)
          { legend("topleft",sexcode.list,col=sexcolFM,pch=sexpchFM,
                   lty=1,lwd=1,cex=0.8)
          }

          if (!is.na(xc.RL1.gen) ) { abline(h=xc.RL1.gen, col=gencol, lty=2) }
          if (!is.na(xc.RL2.gen) ) { abline(h=xc.RL2.gen, col=gencol, lty=2) }

          #  Write parameter information into the plot
          ok <- Infoblock(RunId, yplotmin150.010, yplotmax150.010, 
                     outname.file, x.n14, subsample.n,  oh.sym, dev.sym,
                     round.unit, smooth.hist1, smooth.hist2, 
                     n.per.bin.min, bins.n.min, bins.n.max,
                     lambda.min, lambda.max,  
                     l.fact, p.fact, r.fact, s.fact, w.fact, 
                     x.tr.prop.min, x.tr.prop.max, p.fit.min)
        }
      } else
      {
        #  No valid RL estimates - return NA
        gtab.M1.spl <- list(RL.emp=data.frame(age=NA, RL1=NA,RL2=NA,
                                  RL1.cilo=NA, RL2.cilo=NA,
                                  RL1.tilo=NA, RL2.tilo=NA),
                          RL.smo=data.frame(age=NA,
                                  RL1.hut=NA,
                                  RL2.hut=NA),
                            plot.ok=FALSE,
                            spline.ok=FALSE)
      }

      #  Table with smoothed values produced for this sexcode?
 
      if (gtab.M1.spl$spline.ok)
      { #  Ja
        temp <- data.frame(gtab.M1.spl[["RL.smo"]]$age,
                         gtab.M1.spl[["RL.smo"]]$RL1.hut,
                         gtab.M1.spl[["RL.smo"]]$RL2.hut)

        #  Eindeutige Namen vergeben
        colnames(temp) <- c("Age",
                           paste(sexcode,"_RL1",sep=""),  
                           paste(sexcode,"_RL2",sep="") )     
        if (isex == 1)
        { # 1. sexcode: Gesamt-Tabelle anlegen
          RL.spline <- temp
          #print(sexcode)
          #print(RL.spline)

        } else
        { # späterer sexcode, mergen - should never happen here  
          #print(sexcode)
          RL.spline <- merge(RL.spline,temp,all=TRUE)
          #print(temp)
          #print(RL.spline)
        }
      }  # if (nrow(gtab.M ...

    }    #   if (gtab.M1.spl$spline.ok) ...    
 
    if (!eval.rep) { savePlot(file=PlotFile,type=figtype) }

  }      #   if (gtab.M1.n > ...
}   #  if (plot ...
if (print.log.message) { cat("%%%   TMC_seg160_PlotRL_vs_Age.R  End\n") }
