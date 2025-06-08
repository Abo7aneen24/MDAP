# TMC_seg150_PlotRL_vs_Age_by_S.R 

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

#  10.05.2021 Segment expanded for plotting results of various methods,
#             method contained in variable "RL.age.meth". Plot window and 
#             plotfile name are supplied in "PlotId" and "PlotFile".
#  09.04.2021 oh.sym, dev.sym added
#  02.12.2002 File name changed. Vars x.RL1, x.RL2, ... renamed to RL1, 
#             RL2, ...
#  16.07.2020 Common par parameters for plotting
#  27.02.2020 Infoblock replaced by function
#  24.01.2017 No details into plot if no plot was generated
#  08.01.2020 sexcol, sexpch moved to seg030_Style
#             Var names modified for TMC4
#  18.12.2019 Common y scale fitting for all sexes
#  27.08.2019 graphics file name corrected
#  13.05.     Legend, gridlines controlled by switch
#  08.05.     Legende ausgeschaltet
#  01.05.     Check for no points in input data
#  17.03.2019 Start from WWW17_seg47  (directory ProgWW17_RH)

# ===========================================================================

if (print.log.message) { cat("%%%   TMC_seg150_PlotRL_vs_Age_by_s  Start\n") }

# ===========================================================================

#  sexcol, sexpch set in the styles file

# ----------------------------------------------------------------------------

if (PlotId)        #  plot window
{ 
  #  Reduce gtab to method requested and subset.type  
  gtab.M1 <- gtab[ (gtab[ ,"method"]==RL.age.meth) &
                   (gtab[ ,"subset.type"]==4)   , ]

  #  Weiter, wenn gtab.M1 nicht leer
  gtab.M1.n <- nrow(gtab.M1)

  if (gtab.M1.n > 1)
  { #  Is a matrix with > 1 rows
     
    #  RL.label == NA: leave construction to RL.age.spline()
    if (RL.age.meth == "tmc" & !is.na(RL.label.tmc) )
    { RL.label <- RL.label.tmc }   #  user-defined label

    #  Find existing sex categories
  
    sexcode.list0   <- sort(unique(gtab.M1[ ,"Sex"]))
    sexcode.list0.n <- length(sexcode.list0)

    #  Remove sex codes containing a '+', but keep a list of '+' codes   
    sexcode.list  <- "dummy"
    pluscode.list <- "dummy"

    for (i in 1:sexcode.list0.n)
    { mit.plus <- FALSE
      s.code <- sexcode.list0[i] 
      for (j in 1:nchar(s.code))
      {  
        mit.plus <- mit.plus |(substr(s.code,j,j) == "+") 
      }
      if ( mit.plus) { pluscode.list <- c(pluscode.list,s.code)}
      if (!mit.plus) { sexcode.list <- c(sexcode.list,s.code)}
    }
    sexcode.list   <- sexcode.list[-1]      #  get rid of dummy entry
    sexcode.list.n <- length(sexcode.list)

    pluscode.list   <- pluscode.list[-1]      #  get rid of dummy entry
    pluscode.list.n <- length(pluscode.list)

    #group.shift <- (age.limits[age.limits.n]-age.limits[1]) * 
    #                0.005 * sexcode.list.n
    group.shift <- (age.class[age.class.n, "hi"]-age.class[1, "lo"]) * 
                    0.005 * sexcode.list.n
    group.shift <- 0

    #  'ok' == T signals code with '+'
    ok <- rep(FALSE, times=sexcode.list0.n)
    if (pluscode.list.n == 1)
    { ok <- gtab.M1[ ,"Sex"] == pluscode.list }

    if (pluscode.list.n > 1)
    { for (i in 2:sexcode.list.n)
      { ok <- ok | (gtab.M1[ ,"Sex"] == pluscode.list[i] ) }
    }
    ok[is.na(ok)] <- FALSE

    #  Go through the non-'+' sex codes

    newplot  <- TRUE
    spline.n <- 0       #  counts the number of existing splines
    isex     <- 0
    for (sexcode in sexcode.list)
    { isex <- isex + 1
      ok <- gtab.M1[ ,"Sex"] == sexcode
      gtab.M1.subset <- gtab.M1[ok, ]

      ok <-  !is.na(as.numeric(gtab.M1.subset[ ,"x.RL1"]))

      #  Plot only if at least 2 support points (then plot points only)
      #  If >= 4 points, plot splines also
      if (sum(ok) > 1)
      {
        gtab.M1.subset <- gtab.M1.subset[ok, ] 

        gtab.M1.spl <- RL.age.spline(PlotId,
                      gtab.M1.subset,
                      sexcode,RL.age.meth,group.shift*isex,
                      age.class[1, "lo"], age.class[age.class.n, "hi"],
                      yplotmin150.010,yplotmax150.010,
                      yplotby1150.010,yplotby2150.010,   
                      age.label, RL.label, RL.age.meth,
                      age.clip.min, age.clip.max, age.clip.by1, age.clip.by2,
                      par.tcl, par.las,
                      datafile,sexcol[isex],sexpch[isex],
                      arrow.length=0.04,newplot=newplot)
        newplot <- FALSE
       
        #  Permissible difference is calculated here from the smoothed 
        #  RL-age relation. Not from estimated RL per age class, as in Tab2c

        if (plot.perm.diff.band & gtab.M1.spl$plot.ok)
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
          { legend("topleft",sexcode.list,col=sexcol,pch=sexpch,
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
      } else  #  if (sum(ok) > 1) ...
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
        spline.n <- spline.n + 1       

        #  temp construction depends on whether the permissible
        #  difference was constructed or not
        
        if (plot.perm.diff.band)
        {
          temp <- data.frame(gtab.M1.spl[["RL.smo"]]$age,
                           gtab.M1.spl[["RL.smo"]]$RL1.hut,
                           pd1[ ,1],
                           pd1[ ,2],
                           gtab.M1.spl[["RL.smo"]]$RL2.hut,
                           pd2[ ,1],
                           pd2[ ,2]  )

          #  Eindeutige Namen vergeben
          colnames(temp) <- c("Age",
                           paste(sexcode,"_RL1",sep=""),  
                           paste(sexcode,"_RL1.pdlo",sep=""),
                           paste(sexcode,"_RL1.pdhi",sep=""),
                           paste(sexcode,"_RL2",sep=""),
                           paste(sexcode,"_RL2.pdlo",sep=""),
                           paste(sexcode,"_RL2.pdhi",sep="")
                           )     
        } else
        { # no permissible difference
          temp <- data.frame(gtab.M1.spl[["RL.smo"]]$age,
                           gtab.M1.spl[["RL.smo"]]$RL1.hut,
                           gtab.M1.spl[["RL.smo"]]$RL2.hut )

          #  Eindeutige Namen vergeben
          colnames(temp) <- c("Age",
                           paste(sexcode,"_RL1",sep=""),  
                           paste(sexcode,"_RL2",sep="")  )     
        } 

        if (spline.n == 1)
        { # This is the first existing spline: Install table
          RL.spline <- temp
        } else
        { # This is the second or later spline: merge to previous table 
          RL.spline <- merge(RL.spline,temp,all=TRUE)
        }
      }  # if (gtab.M1.spl$spline.ok) ...

    }    #   for (sexcode in ...    
  
    if (!eval.rep) { savePlot(file=PlotFile,type=figtype) }
 
  }      #   if (gtab.M1.n > ...
}  #  if (plot...

if (print.log.message) { cat("%%%   TMC_seg150_PlotRL_vs_Age_by_s  End\n") }

