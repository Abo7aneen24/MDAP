#  TMC_seg100_Analysis.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Segment doing the analysis (all requested methods) per stratum

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================


# ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Start\n") }
# ==========================================================================


if (time.check)  { tmc.t1 <- Sys.time() }

# ---------------------------------------------------------------------------
#  From here on processing of a new data subset.
#  Prevent previous results from being treated as actual results

if (exists("KDE"))             { rm(KDE) }
if (exists("x.kde"))           { rm(x.kde) }
if (exists("x.kde.mode"))      { rm(x.kde.mode) }
if (exists("HIST"))            { rm(HIST) }
if (exists("x.hist"))          { rm(x.hist) }
if (exists("x.hist.mode.idx")) { rm(x.hist.mode.idx) }
if (exists("tab.strag"))       { rm(tab.strag) }

# ---------------------------------------------------------------------------
#  Initialise error messages

ErrTxt   <- c("Not enough values in stratum",
              "Not enough histogram bins of requested size in stratum",
              "Requested number of bins reduced")

ErrTxt.n <- length(ErrTxt)

ErrMsg <- data.frame(ErrNum=1:ErrTxt.n,
                     Occurred=rep(FALSE, times=ErrTxt.n),
                     ErrTxt,
                     stringsAsFactors=FALSE)

#  Set back execution control. go.on == FALSE means no further execution
#  in this stratum
go.on <- TRUE

# ---------------------------------------------------------------------------
#  Define stratum specific output file names

# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  100\n") }

source("TMC_seg110_NamesPerStratum.R")  

# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  101\n") }

if (print.log.message) 
{ 
  cat("\n", sep.e, 
    "\n   Start analysis by seg100 for", 
    "\n   ", outname.stra,
    "\n", sep.e, 
    "\n", sep="")
}

# ---------------------------------------------------------------------------
#  Graph windows

source("TMC_seg099_OpenWindows.R")

# ---------------------------------------------------------------------------
#  Create results tables
#  All results from this segment, which refer to one stratum,
#  go into 'tab.stra'. Column names are defined in IniTab()
#  Information on data generation was read by read.table and is kept in 
#  the relevant variables - see below

#  Create output table for the actual stratum

tab.template <- IniTab()      #  further use see below
tabc.template <- tab.template[["tabc"]]
tabn.template <- tab.template[["tabn"]]

tabc.stra     <- tabc.template
tabn.stra     <- tabn.template

# ============================================================================
#  Describe the actual stratum

# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 200\n") }

subtitle <- paste(RunId," / ",outname.stra," / n =",x.n)

#  Calculate descriptive quantities
x.n   <- length(x)
x.min <- min(x)
x.max <- max(x)

#  Nonparametric quantiles of the total dataset. Subsequent RL estimates 
#  must lie in [x.Q1, x.Q2]

x.Q1 <- Quantile(x,probs=RL1.p)
x.Q2 <- Quantile(x,probs=RL2.p)

#  Wie differenziert sind die Daten (wieviel verschiedene Werte) ?
x.table <- table(x)

x.val   <- as.numeric(names(x.table))
x.val.n <- length(x.table)


#  Adjust  n.per.bin.min  by dataset size
#  Deactivated 15.09.2022
# n.per.bin.min.adj <- max(Round(0.0025*x.n, 1), n.per.bin.min)

#  Activated 15.09.2022 (bins.n.max now refer to inner range of data)
#n.per.bin.min.adj <- max(Round(x.n/(2*bins.n.max), 1), n.per.bin.min)

# 24.09.22
n.per.bin.min.adj <- max(Round(x.n/(3*bins.n.max), 1), n.per.bin.min)

#  Wieviele Werte unter NWG?
if (!is.na(detect.limits.max))
{
  x.lt.detect.n <- sum(x <= detect.limits.max)
} else
{
  x.lt.detect.n <- 0
}
x.lt.detect.pct <- 100*x.lt.detect.n/x.n

if (exists("age")) 
{ Age.mea <- mean(age) } else
{ Age.mea <- NA }

# ----------------------------------------------------------------------------
#  Put basic information in tab.###, also for subsets that will not be 
#  analysed below due to too few cases

tabn.template["irep"]        <- irep
tabc.template["Sex"]         <- sexlabel
tabc.template["Age"]         <- agelabel
tabn.template["Age.mea"]     <- Age.mea
tabn.template["subset.type"] <- subset.type
tabn.template["n"]           <- x.n
tabn.template["pct.lt.DL"]   <- x.lt.detect.pct

if ("bha" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 250\n") }
  tabn.bha           <- tabn.template
  tabc.bha           <- tabc.template
  tabc.bha["method"] <- "bha"
}

if ("npa" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 300\n") }
  tabn.npa           <- tabn.template
  tabc.npa           <- tabc.template
  tabc.npa["method"] <- "npa"
}

if ("qqw" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 400\n") }
  #  This will be repeated later in the iti loop
  tabn.qqw           <- tabn.template
  tabc.qqw           <- tabc.template
  tabc.qqw["method"] <- "qqw"
}

if ("tmc" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 500\n") }
  #  This will be repeated later in the iti loop
  tabn.tmc           <- tabn.template
  tabc.tmc           <- tabc.template
  tabc.tmc["method"] <- "tmc"
}

if ("tmu" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 600\n") }
  tabn.tmu           <- tabn.template
  tabc.tmu           <- tabc.template
  tabc.tmu["method"] <- "tmu"
}

if ("tuk" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 700\n") }
  tabn.tuk           <- tabn.template
  tabc.tuk           <- tabc.template
  tabc.tuk["method"] <- "tuk"
}

# -------------------------------------------------------------------------

#  Analysis only if data subset has sufficient cases (absolute minimum)
#  Data distribution may additionally cause too few bins of required size - 
#  is checked below. Values below detection limit are checked below.

if (x.n < x.n.min)
{ 
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 800\n") }

  #  Not enough observations - no analysis
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
      "\nData subset for    ", sexlabel,"  ",agelabel,
      "\n  has only         ", x.n," values",
      "\n  Required minimum:", x.n.min,
      "\n                   ==> no analysis",
      "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
 
  #  Set switch and error message 
  ErrMsg[1,"Occurred"] <- TRUE
  go.on <- FALSE
}  

if (go.on && (x.lt.detect.pct > x.lt.detect.pct.max))
{ 
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 850\n") }

  #  Too many observations below DL
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
      "\nData subset for    ", sexlabel,"  ",agelabel,
      "\n  has              ", x.lt.detect.pct, "% of values below DL (", 
                                           detect.limits.max, ")",
      "\n  Allowed maximum :", x.lt.detect.pct.max, "%", 
      "\n                   ==> no analysis",
      "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
 
  #  Set switch and error message 
  ErrMsg[1,"Occurred"] <- TRUE
  go.on <- FALSE
}  

#  If no estimation, create output matrix with only basic information
if (!go.on)
{
  TAB <-   
  CompileTabcn(meth.list, tabc.stra, tabn.stra,    
               tabc.npa, tabc.qqw, tabc.bha, tabc.tmc, tabc.tmu, tabc.tuk,
               tabn.npa, tabn.qqw, tabn.bha, tabn.tmc, tabn.tmu, tabn.tuk)
  tabc.stra <- TAB[["tabc.stra"]]
  tabn.stra <- TAB[["tabn.stra"]]

  tabc.stra[ , "errcode"] <- ".."   # denotes no result

  #  Stitch tabc and tabn together
  tab.stra <- data.frame(tabc.stra, tabn.stra)

  #  New row names
  row.names(tab.stra) <- tab.stra[ ,"method"]    
} 

if (go.on)
{ 
  # ===========================================================================
  #  Required minimum number of observations in a stratum exists,
  #  number of values below DL is acceptable: do analysis

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  900\n") }

  # .........................................................................
  #  Empirical distribution function

  x.emp.cdf <- c(0, cumsum(x.table)/x.n)
  x.emp.x   <- c(0, x.val)

  # ===========================================================================
  #  Is there indication of problematic rounding (eg "go to the even value"
  #  or inconsistent rounding)?  Do diagnostic plot.

  if (plot.fig100.004)
  {
    dev.set(fig100.004)

    par(las=par.las)    
    par(tcl=par.tcl)    
    plot(as.numeric(names(x.table)), as.numeric(x.table),type="h",
         log="x",
         ylim=c(0,max(x.table)),
         xlab=xlabel, ylab="Count",
         main="Frequencies of observed values",
         sub=subtitle,cex.sub=0.7)
    # if (!eval.rep) savePlot(file=fig100.004.file,type=figtype)
  }


  # ==========================================================================
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1000\n") }

  #  Calculate kde
  KDE <- FindKDE(x, round.unit, kernel.n, 
                 fig100.001, subtitle, kdecol,
                 Q.lo=0.025, Q.hi=0.975, b.fact=1.05)
  x.kde      <- KDE[["x.kde"]]
  x.kde.mode <- KDE[["x.kde.mode"]]

  xsupp      <- x.kde$x
  xsupp.clip <- seq(x.clip.min, x.clip.max, length.out=101) # for plots

  x.lt.mode.n <- sum(x <  x.kde.mode)
  x.ge.mode.n <- sum(x >= x.kde.mode)

  # ==========================================================================
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1100\n") }

  #  Produce a reasonable histogram on the original scale, together with a 
  #  smoothed version
 
  HIST <- FindHisto(x, round.unit, detect.limits.max, x.val, x.kde, 
                    n.per.bin.min.adj, bins.n.min, bins.n.max,
                    smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
                    bordercol1, histcol1, bordercol2, histcol2, 
                    bordercol3, histcol3, kdecol, fig100.002, fastnull)

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1150\n") }
 
  x.hist         <- HIST[["x.hist.sc"]] # non-equidistant, smoothed, collapsed

  #  Histogram construction may have left too few bins (==> x.hist == NA).
  x.hist.exists <- (length(x.hist) > 1)

  if (x.hist.exists)
  {
    counts.n   <-  length(x.hist$counts)
    breaks.n   <-  length(x.hist$breaks)

   # #  Which bin in x.hist contains x.kde.mode?
   # x.kde.mode.idx <- FindModeIndex(x.hist$breaks, x.hist$counts, x.kde.mode)

    #  Which bin in x.hist contains the histogram density mode?
    x.hist.mode.idx <- which.max(x.hist$density)
    if (length(x.hist.mode.idx) > 1)
    { cat("\n [seq100] x histogram has several local maxima, first taken \n")
      x.hist.mode.idx <- x.hist.mode.idx[1]
    }

    # =========================================================================
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1200\n") }

    if (!eval.rep) { sink(file=outname.stra.txt, append=TRUE, split=TRUE) }

    cat("\n\n", sep.p,
    "\nStart of analysis                      ", format(time.start),
    "\nData file                              ",
    "\n   ",infile,
    "\nAnalyte                                ",xlabel,
    "\n",
    "\nAnalysis of the stratum defined by",
    "\n  Age                                  ",agelabel,
    "\n  Sex                                  ",sexlabel,    
    "\nAge class definition style             ",age.class.style,
    "\nMin. size for dynamic age classes      ",x.n.min,
    "\nRecords for analysis in this stratum   ",x.n,
    "\nNumber of distinct values              ",x.val.n,
    "\nAbsolute minimum count per bin         ",n.per.bin.min,
    "\nAdjusted minimum count per bin         ",n.per.bin.min.adj,
    "\nMinimum number of bins                 ",bins.n.min,
    "\nMaximum number of bins (central)       ",bins.n.max,
    "\nValues lie in the range                ",x.min, x.max,
    "\nSmoothing parameter for x              ",x.kde$bw,
    "\nKDE mode(x)                            ",x.kde.mode,
    "\n", sep.p, 
    "\n\n")

    if (!eval.rep) { sink() }

    # ---------------------------------------------------------------------------
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1300\n") }

  }  else    #  x.hist.exists
  {
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1400\n") }
    #  x.hist does not exist
    counts.n <- 0
    breaks.n <- 0
  } 

  # =========================================================================
  # =========================================================================
  #  Start RL calculation

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1500\n") }

  # =========================================================================
  #  If group membership is known (typically for simulated data):
  #  Calculate the ideal direct result, if requested  
  #  Ideal direct result: all nonpathological persons are correctly identified
  #  No pathological cases lie in the truncation interval (though no 
  #  truncation interval is needed for RLx.npa)
  # 
  
  if (!is.na(spalte.g) & ("npa" %in% meth.list))
  {
    #  npa is possible and requested
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1600\n") }

    # source("TMC_seg102_Analysis_npa.R")
    tabn.npa <- npa(x, grp, round.unit, x.kde.mode, tabn.npa, 
                    RL1.p, RL2.p, print.log.message=FALSE)

  }  else
  {
    #  npa not possible
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1700\n") }
  }

  # ==========================================================================
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1800\n") }

  if (x.hist.exists)
  { 
    #  There is a reasonable number of bins
    #  Find TI candidates with required properties (proportion of values /
    #  mode in TI)
  
    #  Proportions V8 13.08.2022. Contains proper sorting
    #  x.tr.prop.min, max: user input overrules other settings
    #  If not present, x.tr.prop.min, max depend on whether a DL exists
    #  Processing the input is done once for all segments in seg040 and
    #  produces the effective x.tr.prop.min, x.tr.prop.max

    TI.can <- Proportions(x.hist$counts, x.hist.mode.idx, 
                          x.tr.bins.min, x.tr.prop.min, x.tr.prop.max,
                          (!is.na(detect.limits.max) &  TI.left.only.if.DL) | 
                          TI.left.only, 
                          deltap.lo, deltap.hi, print.log.message=FALSE)

    TI.can.n <- nrow(TI.can)
    if (is.na(TI.can[1,1]) ) TI.can.n <- 0

    #  Minimal proportion and inclusion of range around mode already organized

    if (TI.can.n == 0)
    {
      cat("\n\n [seg100] +++ No truncation interval candidates +++\n\n")

      stop("+++ Stop in seg 100 - no truncation interval candidates +++")

    }  else
    {  
      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1900\n") }
      #  Define a subset for QQW if the dataset is too large  

      #  @@@ move this to Proportion()
      red.proportion <- red.limit/x.n
      if (red.proportion > 1)
      { 
        #  Data set has < red.limit values: take data directly
        x.red <- x
      } else
      {
        #  Data set has >= red.limit values: replace data by select quantiles 
        # x.red   <- c( unname(Quantile(x, probs=seq(0.01, 0.99, by=0.01))) )
        x.red   <- c( unname(Quantile(x, probs=seq(0.005, 0.995, by=0.005))) )
      }
      x.red.n   <- length(x.red)
      x.red.cdf <- (1:x.red.n)/x.red.n

      # ---------------------------------------------------------------------
      #  TMC does loop over acceptable TI candidates, going down from 
      #  good quality (3) to not good quality (1), starting within a quality
      #  group with the largest TIs.
      #  Initial value for nlm is always taken from previous Tukey-QQ loop
      #  Stop criterion: 
      #  sol.score == sol.score.max   and
      #  p.fit >= p.fit.min           and
      #  p.rt  >= p.rt.min 
      #  Consideration of further TIs stops with the first TI that fulfils
      #  these criteria. The result from this TI sets TmcSolutAcc <- TRUE.
      #  If the stop criterion is not reached, the next candidate from TI.can
      #  is used.
      #  From iti = 2 on, the actual result is compared with the so far best
      #  result. This "relative best" solution is reported as result if 
      #  the global optimum is not reached.
      #  Comparison of two solutions: 
      #  1) larger sol.score is better 2) if sol.scores are identical, (larger 
      #  p.fit) smaller opt.crit is better. p.fit may be numerically unstable
      #  with bad fits.

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2400\n") }

      if (time.check)  
      { tmc.t2 <- Sys.time()
        sink(file="../temp/UsedTime.txt", append=TRUE)
        cat("\n tmc; t2-t1;", 
        format(difftime(tmc.t2, tmc.t1)))
        sink()
      }

      # ---------------------------------------------------------------------
      #  Generate a matrix of transformed x.red for all lambda \in lambda.seq
      #  as preparation of initial value calculation 

      Y.red <- matrix(NA, nrow=length(x.red), ncol=lambda.seq.n)
      for (il in 1:lambda.seq.n)
      {
        Y.red[ ,il] <- BoxCox(x.red, lambda.seq[il]) 
      }

      # ---------------------------------------------------------------------
      #  Start loop over TI candidates
      #  Preparation: store partially filled result vectors for use below

      tabn.qqw0 <- tabn.qqw
      tabn.tmc0 <- tabn.tmc

      TmcSolutAcc   <- FALSE
      iti           <- 0      

      while ( (iti < TI.can.n) & !TmcSolutAcc)
      {
        iti <- iti + 1
        if (time.check)  { tmc.t3 <- Sys.time() }

        cat("\n [seg100] iti/TI.can.n=", iti,"/", TI.can.n, 
                 " TI.qual=", TI.can[iti, "TI.qual"],
                 sexlabel, agelabel,
                 " -------------") 

        x.tr.n  <- TI.can[iti, "x.tr.n"]
        prop.lo <- TI.can[iti, "prop.lo"]
        prop.hi <- TI.can[iti, "prop.hi"]
        subset  <- (prop.lo < x.red.cdf) & (x.red.cdf < (1-prop.hi))

        tabn.qqw              <- tabn.qqw0
        tabn.qqw["ilo"]       <- TI.can[iti, "ilo"]
        tabn.qqw["ihi"]       <- TI.can[iti, "ihi"]
        tabn.qqw["x.tr.lo"]   <- x.hist$breaks[TI.can[iti, "ilo"]]
        tabn.qqw["x.tr.hi"]   <- x.hist$breaks[TI.can[iti, "ihi"]+1]
        tabn.qqw["x.tr.n"]    <- TI.can[iti, "x.tr.n"]
        tabn.qqw["x.tr.bins"] <- tabn.qqw["ihi"] - tabn.qqw["ilo"] + 1
        tabn.qqw["x.lt.tr.n"] <- TI.can[iti, "x.lt.tr.n"]
        tabn.qqw["x.ge.tr.n"] <- TI.can[iti, "x.ge.tr.n"]
        tabn.qqw["x.tr.prop"] <- TI.can[iti, "x.tr.prop"]
        tabn.qqw["prop.lo"]   <- TI.can[iti, "prop.lo"]
        tabn.qqw["prop.hi"]   <- TI.can[iti, "prop.hi"]
        tabn.qqw["TI.qual"]   <- TI.can[iti, "TI.qual"]

        tabn.tmc              <- tabn.tmc0
        tabn.tmc["ilo"]       <- TI.can[iti, "ilo"]
        tabn.tmc["ihi"]       <- TI.can[iti, "ihi"]
        tabn.tmc["x.tr.lo"]   <- x.hist$breaks[TI.can[iti, "ilo"]]
        tabn.tmc["x.tr.hi"]   <- x.hist$breaks[TI.can[iti, "ihi"]+1]
        tabn.tmc["x.tr.bins"] <- tabn.tmc["ihi"] - tabn.tmc["ilo"] + 1
        tabn.tmc["x.tr.n"]    <- TI.can[iti, "x.tr.n"]
        tabn.tmc["x.lt.tr.n"] <- TI.can[iti, "x.lt.tr.n"]
        tabn.tmc["x.ge.tr.n"] <- TI.can[iti, "x.ge.tr.n"]
        tabn.tmc["x.tr.prop"] <- TI.can[iti, "x.tr.prop"]
        tabn.tmc["prop.lo"]   <- TI.can[iti, "prop.lo"]
        tabn.tmc["prop.hi"]   <- TI.can[iti, "prop.hi"]
        tabn.tmc["TI.qual"]   <- TI.can[iti, "TI.qual"]

        # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2420\n") }

        if ("npa" %in% meth.list)
        { 
          x.RL1.npa <- tabn.npa["x.RL1"]
          x.RL2.npa <- tabn.npa["x.RL2"]
        }  else
        {
          x.RL1.npa <- NA
          x.RL2.npa <- NA
        }  

        TM <- tmc.master(x, Y.red, subset, xsupp, 
                       round.unit, detect.limits.max, 
                       x.Q1, x.Q2, x.RL1.npa, x.RL2.npa, 
                       x.hist,
                       idx.fix, idx.est,
                       lambda.min, lambda.max, lambda.seq, 
                       l.fact, p.fact, r.fact, rt.fact, w.fact,
                       df.est, df.con, 
                       RL1.p, RL2.p,
                       prev.acc.lo, prev.acc.hi,
                       gencol, kdecol,tmccol, NA,
                       tabn.qqw, tabn.tmc,
                       fastnull, fastnull.chi,
                       print.summary=FALSE, print.log.message=FALSE)

        # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2440\n") }

        tabn.tmc <- TM[["tabn.tmc"]]  
        tabn.qqw <- TM[["tabn.qqw"]]  

        #  Store results for the case of no acceptable result
        if (iti == 1) 
        { 
          tabc.tmc.all <- tabc.tmc
          tabc.qqw.all <- tabc.qqw
          tabn.tmc.all <- tabn.tmc
          tabn.qqw.all <- tabn.qqw
        } else
        {
          tabc.tmc.all <- rbind(tabc.tmc.all, tabc.tmc)
          tabc.qqw.all <- rbind(tabc.qqw.all, tabc.qqw)
          tabn.tmc.all <- rbind(tabn.tmc.all, tabn.tmc)
          tabn.qqw.all <- rbind(tabn.qqw.all, tabn.qqw)
        }   

        # --------------------------------------------------------------------
        #  Is the actual solution acceptable?
        #
        #  New 03.09.2022 
        TmcSolutAcc <- (tabn.tmc["sol.score"] == sol.score.max) & 
                       (tabn.tmc["p.fit"] >= p.fit.min) &
                       (tabn.tmc["p.rt"]  >= p.rt.min)
        TmcSolutAcc <-  unname(TmcSolutAcc)
    
      }     # while (iti ...    

      cat("\n\n [seg100] End loop iti/TI.can.n  ----------------------------\n") 
      cat("\n")

      tabn.tmc.all.n <- iti

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3200\n") }

      #  Decide about final solution:
      #  If an acceptable solution was found, take it as final solution, 
      #  even if the corresponding opt.crit is not the global minimum
      #  If no acceptable solution was found, take the solution with minimal
      #  opt.crit as the final solution

      if (tabn.tmc.all.n > 1)
      {
        #  Replace actual row names by iti
        rownames(tabc.tmc.all) <- 1:tabn.tmc.all.n
        rownames(tabc.qqw.all) <- 1:tabn.tmc.all.n
        rownames(tabn.tmc.all) <- 1:tabn.tmc.all.n
        rownames(tabn.qqw.all) <- 1:tabn.tmc.all.n

        #  Sort order of results by solution score, opt.crit
        so <- order(-tabn.tmc.all[ ,"sol.score"], tabn.tmc.all[ ,"opt.crit"])
        tabc.tmc.all <- tabc.tmc.all[so, ]
        tabc.qqw.all <- tabc.qqw.all[so, ]
        tabn.tmc.all <- tabn.tmc.all[so, ]
        tabn.qqw.all <- tabn.qqw.all[so, ]
      }  else
      {
        tabc.tmc.all <- RowToMatrix(tabc.tmc.all)
        tabc.qqw.all <- RowToMatrix(tabc.qqw.all)
        tabn.tmc.all <- RowToMatrix(tabn.tmc.all)
        tabn.qqw.all <- RowToMatrix(tabn.qqw.all)

        rownames(tabc.tmc.all) <- 1
        rownames(tabc.qqw.all) <- 1
        rownames(tabn.tmc.all) <- 1
        rownames(tabn.qqw.all) <- 1
      }

      if ( TmcSolutAcc )
      { #  tabc.tmc, tabn.tmc and  tabc.qqw, tabn.qqw are the solutions
        cat("\n [seg100] TMC estimation finished ",
            "\n          Acceptable solution found for TI ", iti,
            "\n")

        #  Last TI giving tabn.qqw and tabn.tmc is optimal 
        iti.opt     <- iti
      } else
      { #  first line of tabc.tmc.all, tabn.tmc.all and  
        #  tabc.qqw.all, tabn.qqw.all are the solutions
       
        iti.opt <- rownames(tabn.tmc.all)[1]

        cat("\n [seg100] TMC estimation finished ",
            "\n          NO ACCEPTABLE SOLUTION FOUND",
            "\n          Solution for TI ", iti.opt, " is final solution",
            "\n")
        #  Find the best among the non acceptable solutions
        tabc.tmc <- tabc.tmc.all[iti.opt, ]  
        tabn.tmc <- tabn.tmc.all[iti.opt, ]  
        tabc.qqw <- tabc.qqw.all[iti.opt, ]  
        tabn.qqw <- tabn.qqw.all[iti.opt, ]  
      }

      # Show the solution found
      cat("\n [seg100] Solution found ",
          "\n")
      print(tabn.tmc[c("ilo", "ihi", "x.tr.prop", 
                       "lambda", "mue", "sigma",
                       "x.RL1", "x.RL2", "p.fit", "p.rt", "opt.crit", 
                       "sol.score", "rel.dist") ] )

      #  Show alternatives with good opt.crit
      #cat("\n [seg100] Alternative results, sorted by sol.score, opt.crit ",
      #    "\n")

      #print(tabn.tmc.all[1:min(5, nrow(tabn.tmc.all)),
      #                     c("ilo", "ihi", "x.tr.prop", 
      #                       "lambda", "mue", "sigma",
      #                       "x.RL1", "x.RL2", "p.fit", "p.rt", "opt.crit", 
      #                       "sol.score", "rel.dist") ] )

      # ----------------------------------------------------------------------
      # If requested, show properties of all solutions

       if (plot.sol.props)
       {
         ShowTiProps(tabn.tmc.all, p.fit.min, p.rt.min,
                     paste(xlabel, sexlabel, agelabel))
         cat("\n ***  Click into picture to proceed  *** \n")
         locator(1) 
       }

      # ----------------------------------------------------------------------
      # Fit criteria for tab.qqw not yet calculated, done below

      prop.qqw <- chi2trunc(x.hist,tabn.qqw["ilo"], tabn.qqw["ihi"],
                             tabn.qqw["x.lt.tr.n"], tabn.qqw["x.ge.tr.n"],
                             x.Q1, x.Q2, RL1.p, RL2.p,
                             tabn.qqw["lambda"],
                             tabn.qqw["mue"],
                             tabn.qqw["sigma"],
                             df.est, df.con,
                             l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                             rt.fact=rt.fact, w.fact=w.fact, 
                             opt.crit.only=FALSE,fastnull=fastnull)
  
      tabn.qqw["opt.crit"] <- unname(prop.qqw[["res"]]["opt.crit"])
      tabn.qqw["p.fit"]    <- unname(prop.qqw[["res"]]["chi2.trun.p"])
      tabn.qqw["p.rt"]     <- unname(prop.qqw[["res"]]["p.rt"])
      tabn.qqw["prev.l"]   <- prop.qqw$res["prev.l.tmc"]  # @@@ change in 
      tabn.qqw["prev.c"]   <- prop.qqw$res["prev.c.tmc"]  # @@@ chi2trunc
      tabn.qqw["prev.r"]   <- prop.qqw$res["prev.r.tmc"]

      # ----------------------------------------------------------------------
      # Asymptotic confidence intervals for RLs
      CI <-  CIQuant.PNV(tabn.qqw["n"], RL1.p, tabn.qqw["x.RL1"],
                         tabn.qqw["lambda"], tabn.qqw["mue"], tabn.qqw["sigma"],
                         alpha=0.05, fastnull=fastnull)
      tabn.qqw["x.RL1.cilo"] <- CI[1]
      tabn.qqw["x.RL1.cihi"] <- CI[2]
      CI <-  CIQuant.PNV(tabn.qqw["n"], RL2.p, tabn.qqw["x.RL2"],
                         tabn.qqw["lambda"], tabn.qqw["mue"], tabn.qqw["sigma"],
                         alpha=0.05, fastnull=fastnull)
      tabn.qqw["x.RL2.cilo"] <- CI[1]
      tabn.qqw["x.RL2.cihi"] <- CI[2]

      CI <-  CIQuant.PNV(tabn.tmc["n"], RL1.p, tabn.tmc["x.RL1"],
                         tabn.tmc["lambda"], tabn.tmc["mue"], tabn.tmc["sigma"],
                         alpha=0.05, fastnull=fastnull)
      tabn.tmc["x.RL1.cilo"] <- CI[1]
      tabn.tmc["x.RL1.cihi"] <- CI[2]
      CI <-  CIQuant.PNV(tabn.tmc["n"], RL2.p, tabn.tmc["x.RL2"],
                         tabn.tmc["lambda"], tabn.tmc["mue"], tabn.tmc["sigma"],
                         alpha=0.05, fastnull=fastnull)
      tabn.tmc["x.RL2.cilo"] <- CI[1]
      tabn.tmc["x.RL2.cihi"] <- CI[2]

      # ----------------------------------------------------------------------
      #  Collect extra information for histogram plotting
      #  Variables required in F_PlotHistResid.R
      #c("UG", "OG", "nobs", "nerw", "diff")

      prop.tmc <- chi2trunc(x.hist,
                             tabn.tmc["ilo"], tabn.tmc["ihi"],
                             tabn.tmc["x.lt.tr.n"], tabn.tmc["x.ge.tr.n"],                             
                             x.Q1, x.Q2, RL1.p, RL2.p,
                             tabn.tmc["lambda"],
                             tabn.tmc["mue"],
                             tabn.tmc["sigma"],
                             df.est, df.con,
                             l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                             rt.fact=rt.fact, w.fact=w.fact,
                             penalty.details=TRUE, 
                             opt.crit.only=FALSE,fastnull=fastnull)

      # Rename tables in prop.tmc to fit plot function names                   
      tt.names <- names(prop.tmc$tab)
      tt.names[1] <- "UG"
      tt.names[2] <- "OG"
      tt.names[3] <- "nobs"
      tt.names[10] <- "nerw"
      tt.names[11] <- "diff"
      names(prop.tmc$tab) <- tt.names
      if (length(prop.tmc$tab.lo) > 1 ) 
         { names(prop.tmc$tab.lo) <- tt.names }
      if (length(prop.tmc$tab.hi) > 1 ) 
         { names(prop.tmc$tab.hi) <- tt.names }

      if (time.check)  
      { tmc.t7 <- Sys.time() 

        sink(file="../temp/UsedTime.txt", append=TRUE)
        cat("\n tmc; t7-t3;", 
        format(difftime(tmc.t7, tmc.t3)))
        sink()
      }

      # =========================================================================
      #  TMC estimation is finished
      # =========================================================================

      #  Print and plot result   
      #  Replace pseudo-zero estimate by zero for nicer printout
      if (abs(tabn.tmc["lambda"]) < fastnull) 
      { 
        tabn.tmc["lambda"] <- 0 
      }
    
      TI <- tabn.tmc["ilo"] : tabn.tmc["ihi"]
      counts.range <- paste(formatC(min(x.hist$counts[TI]), format="f", digits=0),
                        "-", 
                        formatC(max(x.hist$counts[TI]), format="f", digits=0), 
                        sep="")
      x.tr.prop.range <-               
           paste(formatC(100*tabn.tmc["x.tr.prop"], format="f", width=5, 
                 digits=1), "%", sep="")

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3400\n") }

      if (plot.fig100.016)
      { 
        #  Plot result - basic part
        #  No plot file name, save is done here below
        fig100.016.lim <- PlotHistFit(fig100.016, NA, figtype, 
                x.hist, xlabel, "TMC result", subtitle, n.per.bin.min.adj,
                NA, x.val, counts.range, x.tr.prop.range,
                x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                x.clip.type, NA, par.las, par.tcl,
                bordercol, histcol, difcol, kdecol, denlty)

        #  Plot result - tmc specific part 
        #  Add estimated tmc solution, RLs, TI
        #  Calculate estimated TMC density (not yet scaled by prev.c)
        xc.pdf.tmc <- pdf.PN(xsupp.clip, tabn.tmc["lambda"], 
                                  tabn.tmc["mue"], tabn.tmc["sigma"])    
  
        PlotMetRes(fig100.016, fig100.016.file, figtype, 
                   fig100.016.lim["yplotmin"],
                   xsupp.clip, tabn.tmc["prev.c"], xc.pdf.tmc, tmccol,
                   x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol, 
                   tabn.tmc["x.RL1"], tabn.tmc["x.RL2"], 
                   0.85*fig100.016.lim["yplotmax"], 
                   tabn.tmc["x.tr.lo"], tabn.tmc["x.tr.hi"], 
                   0.90*fig100.016.lim["yplotmin"], tmccol,
                   prop.tmc)

        # -------------------------------------------------------------
        # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3500\n") }

        #  If the data is generated : show unaffected range
        if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
        { 
          lines(c(x.unaffected.lo, x.unaffected.hi), 
                    0.60*fig100.016.lim["yplotmin"]*c(1,1), col="chartreuse", 
                    lwd=2)
        }

        # -------------------------------------------------------------
        # Save tmc result figure
        savePlot(file=fig100.016.file, type=figtype)
 
      }      #  if (plot.fig100.016)

 
      # ================================================================
      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3600\n") }

      if ("tmu" %in% meth.list)
      { 
        #  Calculate TMU solution, using lambda from tmc
        source("TMC_seg106_Analysis_tmu.R") 
       }  

      # ================================================================
      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3700\n") }           

      if ("tuk" %in% meth.list)
      {       
        #  Calculate Tukey solution, using lambda from tmc
        source("TMC_seg108_Analysis_tuk.R")       
      }  

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3800\n") }           

      # ================================================================

    }                #  x.hist exists

    # ==========================================================================
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3900\n") }

  }    #  TI.can.n > 0?
  # ==========================================================================

  #  Combine essential results from all methods for this stratum
  #  tab.stra has a dummy line from definition in the beginning here

  TAB <-   
  CompileTabcn(meth.list, tabc.stra, tabn.stra,    
               tabc.npa, tabc.qqw, tabc.bha, tabc.tmc, tabc.tmu, tabc.tuk,
               tabn.npa, tabn.qqw, tabn.bha, tabn.tmc, tabn.tmu, tabn.tuk)
  tabc.stra <- TAB[["tabc.stra"]]
  tabn.stra <- TAB[["tabn.stra"]]

  # --------------------------------------------------------------------------
  #  Error code for result. Considers goodness of fit (p.fit)
  #  and further properties of the solution. Version 2022-08-05.
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 4000\n") }

  tabc.stra[ , "errcode"] <- "!!"   # should never appear

  ok1 <- (tabn.stra[ , "p.fit"] > p.fit.min) 
  ok2 <- (tabn.stra[ , "p.rt"]  > p.rt.min) 

  ok1[is.na(ok1)] <- FALSE
  ok2[is.na(ok2)] <- FALSE

  tabc.stra[ok1 & ok2, "errcode"]   <- "++"
  tabc.stra[ok1 & !ok2, "errcode"]  <- "+-"
  tabc.stra[!ok1 & ok2, "errcode"]  <- "-+"
  tabc.stra[!ok1 & !ok2, "errcode"] <- "--"

  # No mark, if p.fit or p.rt not determined
  tabc.stra[is.na(tabn.stra[ , "p.fit"]) |
            is.na(tabn.stra[ , "p.rt"])    , "errcode"] <- ".."

  # --------------------------------------------------------------------------
  #  Toleranzintervalle für die tmc-RLs per bootstrapping berechnen
  #  NOT YET ADAPTED TO THIS VERSION OF  tmc.master

  source("TMC_seg109_TMC_BS.R")

  # -------------------------------------------------------------------------

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 5000\n") }

  # -------------------------------------------------------------------------
  #  Stitch tabc and tabn together
  tab.stra <- data.frame(tabc.stra, tabn.stra)

  #  New row names
  row.names(tab.stra) <- tab.stra[ ,"method"]    

  #  ==================================================================
  #  Print results tables per stratum (not if replicates are analysed) 

  if (!eval.rep) { sink(file=outname.stra.txt, append=TRUE, split=TRUE) }

  cat("\n\n", sep.e, 
      "\n\nRL ESTIMATES, ALL REQUESTED METHODS, ONE age-sex STRATUM",
      "\nUser                                   ", user,
      "\nData information taken from            ", info.file,
      "\nAnalysis date                          ", date(),
      "\nCopy of this table in                  ",
      "\n    ",                                   outname.stra.txt,
      "\nData file name                         ",
      "\n    ",                                   infile,
      "\nFilter  outpatient/hospitalized        ", use.oh,
      "\nFilter  device                         ", use.dev,
      "\nFilter  analysis time of day           ", paste(ana.hour.min, ":",
                                                         ana.minu.min, " - ",
                                                         ana.hour.max, ":",
                                                         ana.minu.max, sep=""),
      "\nSex group                              ", sexlabel,    
      "\nAge interval                           ", agelabel,
      "\nInput values are multiplied by         ", scale.fact,
      "\nScaled values are limited by           ", x.lo.limit,x.hi.limit,
      "\nScaled values are rounded here to      ",round.unit,
      "\nTotal number of all values in file     ",x.n1,
      "\nTotal number of numerical values       ",x.n14,
      "\nNumber of values in this stratum       ",x.n,
      "\nNumber of distinct values in stratum   ",x.val.n,
      "\nMaximal detection/determination limit  ",detect.limits.max,
      "\n# of values / % below detection limit  ",x.lt.detect.n, " / ",
                                                  format(x.lt.detect.pct,digits=3),
                                                  "%",
      "\nValue range                            ",x.min, x.max,
      "\nEmpirical P025, P975, complete stratum ",x.Q1, x.Q2,
      "\n",
      "\nSmooth empirical histogram?            ",smooth.hist1, smooth.hist2,
      "\nMin. # of values per bin (adj.)        ",n.per.bin.min.adj,
      "\nMin. # of bins in stratum              ",bins.n.min,     
      "\nMin. # of bins in truncation interval  ",x.tr.bins.min,     
      "\nMin. # of values in stratum            ",x.n.min,     
      "\nActual # of bins in stratum            ",counts.n,
      "\nUser-specified min. prop. in TI        ", x.tr.prop.min, 
      "\nUser-specified max. prop. in TI        ", x.tr.prop.max, 
      "\nWeight of truncation interval length   ",l.fact,
      "\nPenalty factor 'wrong prevalence'      ",p.fact,
      "\nPenalty factor 'wrong prediction'      ",w.fact,
      "\nPenalty factor 'outside P025, P975'    ",r.fact,
      "\nPenalty factor 'runs test'             ",rt.fact,
      "\nMinimal lambda                         ",lambda.min,
      "\nMaximal lambda                         ",lambda.max,
      "\n",
      "\nTMC: Truncation interval               ",tabn.tmc["x.tr.lo"], 
                                                  tabn.tmc["x.tr.hi"],   
      "\nTMC: TI contains proportion ...        ",tabn.tmc["x.tr.prop"],
      "\nTMC: estimated lambda                  ",tabn.tmc["lambda"],
      "\nTMC: estimated mu                      ",tabn.tmc["mue"],
      "\nTMC: estimated sigma                   ",tabn.tmc["sigma"],
      "\nTMC: chi^2 from truncation interval    ",format(
                                                    tabn.tmc["chi2.trun"],
                                                    digits=4),
      "\nTMC: chi^2 from outside TI             ",format(
                                                   tabn.tmc["chi2.total"]-
                                                   tabn.tmc["chi2.trun"],
                                                   digits=4),
      "\n",
      "\nFalse negative probability  <RL1 / <RL2", RL1.p,"/ ",RL2.p,
      "\nTMC: estimated RLs                       RL1 = ", 
                                             format(tabn.tmc["x.RL1"],
                                                    digits=4),      
      "\n                                         RL2 = ", 
                                             format(tabn.tmc["x.RL2"],
                                                    digits=4), 
      "\nAssessment of TMC solution",
      "\np for chi^2 goodness of fit            ",format(
                                                    tabn.tmc["p.fit"],
                                                    digits=4), 
      "\np for runs.test                        ",format(
                                                    tabn.tmc["p.rt"],
                                                    digits=4), 
      "\nRL1 >= x.Q1, RL2 <= x.Q2 ?             ", tabn.tmc["crit1"],  
                                                   tabn.tmc["crit2"], 
      "\nprev.l ok / prev.c ok /prev.r ok?      ", tabn.tmc["crit3"],
                                                   tabn.tmc["crit4"], 
                                                   tabn.tmc["crit5"], 
      "\nTotal score for solution (max =",sol.score.max,")    ", 
                                                tabn.tmc["sol.score"],    
      "\nAcceptable solution found?             ", 
                                                TmcSolutAcc,    
      "\n", sep.e, 
      "\n\n")


  if (user == "WW")
  { #  Ausführung für WW

    #cat("\n All details \n")
    #print(tabn.stra[meth.list, ])
    #cat("\n")
  }

  if (!eval.rep) { sink() }
 
  # ==========================================================================
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6000\n") }
  # ==========================================================================

}  #  if (go.on) ...

# ---------------------------------------------------------------------------
#  Copy results in tab.stra (one stratum) to gtab

# ==========================================================================
# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6100\n") }
# ==========================================================================

#  Reduce tab.stra to those methods which are requested for gtab
meth.idx <- rep(NA, times=gtab.list.n)
imeth <- 0
for (meth in gtab.list)
{
  imeth <- imeth + 1
  meth.idx[imeth] <- which(tab.stra[ ,"method"] == meth)
}
meth.idx <- meth.idx[!is.na(meth.idx)]

tab.strag <- tab.stra[meth.idx, ]

# Add variables to tab.stra that are needed for printing
tab.strag <- data.frame(SexPrt=tab.strag[ , "Sex"],
                        AgePrt=tab.strag[ , "Age"],
                        tab.strag )
if (gtab.list.n > 1)
{ 
  tab.strag[2:gtab.list.n, "SexPrt"] <- " "
  tab.strag[2:gtab.list.n, "AgePrt"] <- " "
}

iblock <- iblock + 1

#  Append tab.strag to gtab
if (iblock == 1)
{
  #  First block at all, establish gtab
  gtab <- tab.strag
}  else
{
  #  Append block to gtab
  gtab <- rbind(gtab, tab.strag)
}

# ==========================================================================
  if (print.log.message) { cat("%%%   TMC_seg100_Analysis  End\n") }
# ==========================================================================

cat("\n", sep.e, 
    "\n   End analysis by seg100 for", 
    "\n   ", outname.stra,
    "\n", sep.e, 
    "\n", sep="")
   
#############################################################################