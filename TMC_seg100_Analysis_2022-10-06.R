#  TMC_seg100_Analysis.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Segment doing the analysis per stratum

#  #  CHANGE HISTORY

#  01.07.2022
#  13.12.2021 Change of optimization sequence: Tukey+QQ only per lambda
#  13.12.2021 Safety prior to change of optimization sequence
#  03.10.2021 This is the test version for seg 100, original see seg 100 backup
#             Calculation of initial values for tmc modified
#  21.08.2021 Decision sequenc for TI search streamlined
#  02.05.2021 Result plots per method moved to method-specific segments
#  25.04.2021 Tukey approach added
#  03.03.2021 Basic information per stratum always in output table
#  22.02.2021 Calculation of x.tr.prop.limits moved to 040_ProcessSettings 
#  15.02.2021 Subsegments outsourced (100_Analysis_QQW, ...)
#  26.01.2021 meth.list made effective
#  13.01.2021 r.fact introduced
#  11.01.2021 Search sequence for TI changed to large -> small with 
#             intermediate checks for acceptability
#  04.01.2021 Smoothing strategy changed, applied to uncollapsed equidistant 
#             data 
#  17.12.2020 Graph windows opened by seg046
#  07.12.2002 Safety actions if QQW produces no result
#  01.12.2020 Organisation of histogram calculation, smoothing and reduction
#             changed, now done by FindHisto()
#  16.11.2020 x.kde added to call QQW()
#  14.11.2020 bins.n.min added to call QQW()
#  09.11.2020 New approach of analysis started. TMC is now a function, 
#             qqw and Bhattacharya added. New approach for kde bandwidth 
#             and histogram bandwidth.

#  21.08.2020 Save plot100.040 only for regular analysis, not simulation 
#             Unaffected range taken from DataInfoFile, not calculated here
#  18.08.2020 x.clip.by1 == NA or x.clip.by2 == NA correctly processed 
#  21.07.2020 Pr.lo, Pr.hi set to 0.001, 0.999 (previously 0.01, 0.95)
#             Setting moved to TMC_seg015_DefaultSettings.R
#  20.07.2020 Detail output for startum in txt file re-organized
#  16.07.2020 Common par parameters for plotting
#  01.07.2020 Check generated data unified
#  30.05.2020 Logarithmic x axis in fig100.010
#  25.05.2020 RB3[2, ] expanded: rule for y[1]==0 added
#     05.2020 Plotting of density with new approach for density bw added
#  11.05.2020 RB3[2, ] geändert
#  29.04.2020 Plots newly arranged. kde plot still excluded.
#  21.04.2020 x.kde calculated from only an interior range of the data
#  29.03.2020 opt.crit als einheitliche Bezeichnung
#  11.03.2020 Strategie vorübergehend zurückgesetzt, neuen Namen beibehalten
#  10.03.2020 Strategie für Bandbreitenberechnung geändert 
#  19.02.2020 tmc.master4 replaced by tmc.master5 (hierarchical search of
#             truncation interval)
#  24.01.2020 New structure for stopping execution if not enough data
#  23.01.2020 Break values for histogram calculated by find.breaks,
#             construction keeps n.per.bin.min and x.bins.min, if possible
#  18.01.2020 QQ plot based calculation of initial values
#  05.01.2020 Plot numbers changed
#  04.01.2020 Notation unified: always a . after x, y in x.hist, x.breaks, ...
#             No x or y with x.counts: only 'counts'
#  03.01.2020 New organisation: only tmc is done, new strategy for finding
#             the optimal truncation interval, initial value selection 
#             simplified, results are stored in 'tab' only  
#  05.11.2019 Plots for bootstrapping updated
#  27.10.2019 More information for BS result in output table
#  23.10.2019 TMC estimation applied to smoothed histogram
#  12.10.2019 BS calculation of CI activated
#  10.10.2019 TMC part of the analysis moved to TMC_seg101_Analysis_TMC.R
#  28.08.2019 prop.tmc.seq calculation revised to chi2 contribution
#             p.fit.min introduced, see Expert settings
#  27.08.2019 prop.tmc.seq calculation revised
#  14.08.2019 prop.tmc.seq shifted to interval above detection limit
#  26.07.2019 In RL.est: xc.mode, xc.sig, xc.P50 added, also for
#             methods other than IFCC 
#             In RL.est: x.mode,  y.mode,  y.sig replaced by
#                        xc.mode, yc.mode, yc.sig 
#  24.07.2019 Calculation prop.qqw changed
#             New file naming convention, see seg050
#  22.07.2019 Typing error in optimality criterion removed
#  13.07.2019 Sequence for prop.tmc candidates defined in 
#             TMC_seg060_ExpertSettings.R
#  02.07.2019 names of TMC result variables corrected
#  01.07.2019 RL1.p, RL2.p  instead of fixed probabilities
#  08.06.2019 Empirische Quantile grundsätzlich berechnet
#  12.05.2019 errcode für tmc noch einmal geändert
#             prop.qqw and prop.tmc are dynamcally calculated, depending
#             on the number of values < detection limit
#  08.05.2019 errcode für tmc geändert
#  07.05.2019 qqw operates on x with ties replaced by increasing sequence
#  05.04.2019 y.mode hier berechnet
#  12.03.2019 Start from WW17_seg40  (directory ProgWW17_RH)

# ===========================================================================

# if (print.log.message) 
#{ cat("\n ==================================================================",
#      "\n   %%%   TMC_seg100_Analysis  Start",
#      "\n ==================================================================",
#      "\n") 
#}

# ===========================================================================

tmc.t1 <- Sys.time()

# ......................................................................
#  Clean up. From here on processing of a new data subset.

if (exists("x.kde"))  { rm(x.kde) }
if (exists("y.kde"))  { rm(y.kde) }
if (exists("x.hist")) { rm(x.hist) }
if (exists("y.hist")) { rm(y.hist) }

# ......................................................................
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

# ......................................................................
#  Define stratum specific output file names

source("TMC_seg110_NamesPerStratum.R")  

# ......................................................................
#  Graph windows

source("TMC_seg099_OpenWindows.R")

# ......................................................................
#  Create results tables
#  All results from this segment, which refer to one stratum,
#  go into 'tab.stra'. Column names are given in tab.names, defined
#  in seg050.
#  Information on data generation was read by read.table and is kept in 
#  the relevant variables - see below

#  Create output table for the actual stratum
tab.stra           <- matrix(NA,nrow=meth.list.n, ncol=length(tab.names))
dimnames(tab.stra) <- list(meth.list, tab.names)
tab.stra           <- data.frame(tab.stra)

# ......................................................................
#  Separator lines for output tables
sep.p <- str_dup(".", 78)
sep.m <- str_dup("-", 78)
sep.e <- str_dup("=", 78)

# ============================================================================
#  Decribe the actual stratum

# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 100\n") }

subtitle <- paste(RunId," / ",outname.stra," / n =",x.n)

x.n <- length(x)

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
#  Put basic information in tab.stra, also for subsets that will not be 
#  analysed below due to too few cases

if ("bha" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 200\n") }

  tab.stra["bha", "irep"]    <- irep
  tab.stra["bha", "method"]  <- "bha"

  tab.stra["bha", "sexlabel"]    <- sexlabel
  tab.stra["bha", "agelabel"]    <- agelabel
  tab.stra["bha", "Age.mea"]     <- Age.mea

  tab.stra["bha","subset.type"] <- subset.type
  tab.stra["bha", "n"]          <- x.n
  tab.stra["bha", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("npa" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 300\n") }

  tab.stra["npa", "irep"]    <- irep
  tab.stra["npa", "method"]  <- "npa"

  tab.stra["npa", "sexlabel"]    <- sexlabel
  tab.stra["npa", "agelabel"]    <- agelabel
  tab.stra["npa", "Age.mea"]     <- Age.mea

  tab.stra["npa","subset.type"] <- subset.type
  tab.stra["npa", "n"]          <- x.n
  tab.stra["npa", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("qqw" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 400\n") }

  tab.stra["qqw", "irep"]    <- irep
  tab.stra["qqw", "method"]  <- "qqw"

  tab.stra["qqw", "sexlabel"]    <- sexlabel
  tab.stra["qqw", "agelabel"]    <- agelabel
  tab.stra["qqw", "Age.mea"]     <- Age.mea

  tab.stra["qqw","subset.type"] <- subset.type
  tab.stra["qqw", "n"]          <- x.n
  tab.stra["qqw", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tmc" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 500\n") }

  tab.stra["tmc", "irep"]    <- irep
  tab.stra["tmc", "method"]  <- "tmc"

  tab.stra["tmc", "sexlabel"]    <- sexlabel
  tab.stra["tmc", "agelabel"]    <- agelabel
  tab.stra["tmc", "Age.mea"]     <- Age.mea

  tab.stra["tmc","subset.type"] <- subset.type
  tab.stra["tmc", "n"]          <- x.n
  tab.stra["tmc", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tmu" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 600\n") }

  tab.stra["tmu", "irep"]    <- irep
  tab.stra["tmu", "method"]  <- "tmu"

  tab.stra["tmu", "sexlabel"]    <- sexlabel
  tab.stra["tmu", "agelabel"]    <- agelabel
  tab.stra["tmu", "Age.mea"]     <- Age.mea

  tab.stra["tmu","subset.type"] <- subset.type
  tab.stra["tmu", "n"]          <- x.n
  tab.stra["tmu", "pct.lt.DL"]  <- x.lt.detect.pct
}

if ("tuk" %in% meth.list)
{
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 700\n") }

  tab.stra["tuk", "irep"]    <- irep
  tab.stra["tuk", "method"]  <- "tuk"

  tab.stra["tuk", "sexlabel"]    <- sexlabel
  tab.stra["tuk", "agelabel"]    <- agelabel
  tab.stra["tuk", "Age.mea"]     <- Age.mea

  tab.stra["tuk","subset.type"] <- subset.type
  tab.stra["tuk", "n"]          <- x.n
  tab.stra["tuk", "pct.lt.DL"]  <- x.lt.detect.pct
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

if (x.lt.detect.pct > x.lt.detect.pct.max)
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

if (go.on)
{ 
  # ===========================================================================
  #  Required minimum number of observations in a stratum exists,
  #  number of values below DL is acceptable: do analysis

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  900\n") }

  #  Calculate descriptive quantities
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

  # .........................................................................
  #  Empirical distribution function

  x.emp.cdf <- c(0, cumsum(x.table)/x.n)
  x.emp.x   <- c(0, x.val)

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
 
  #HIST <- FindHisto(x, round.unit, detect.limits.max, x.val, x.kde, 
  #                  n.per.bin.min.adj, bins.n.min, bins.n.max,
  #                  smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
  #                  bordercol1, histcol1, bordercol2, histcol2, 
  #                  bordercol3, histcol3, kdecol, fig100.002, fastnull,
  #                  bins.n.des=bins.n.max)
  HIST <- FindHisto(x, round.unit, detect.limits.max, x.val, x.kde, 
                    n.per.bin.min.adj, bins.n.min, bins.n.max,
                    smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
                    bordercol1, histcol1, bordercol2, histcol2, 
                    bordercol3, histcol3, kdecol, fig100.002, fastnull)

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1150\n") }

 # x.hist.eqd     <- HIST[["x.hist.eqd"]]      # equidistant, raw
 # x.hist.eqd.smo <- HIST[["x.hist.eqd.smo"]]  # equidistant, smoothed
 # x.hist         <- HIST[["x.hist"]]          # collapsed, possibly smoothed

  #  17.01.2022
  #x.hist         <- HIST[["x.hist.eqd.smo"]]  # equidistant, smoothed
  
  #  14.09.2022 
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
    "\nData file                              ",datafile, 
    "\nAnalyte                                ",xlabel, 
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

  }  else    #  x.hist.exists
  {
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1400\n") }
    #  x.hist does not exist
    counts.n <- 0
    breaks.n <- 0
  } 

  # =========================================================================
  #  Start RL calculation

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1500\n") }

  #  Initialize output tables 

  source("TMC_seg101_Analysis_Ini.R")

  # =========================================================================
  #  If group membership is known (typically for simulated data):
  #  Calculate the ideal direct result  
  #  Ideal direct result: all nonpathological persons are correctly identified
  #  No pathological cases lie in the truncation interval (though no 
  #  truncation interval is needed for RLx.npa)
  # 
  
  if (!is.na(spalte.g))
  {
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1600\n") }

    source("TMC_seg102_Analysis_npa.R")

  }  else
  {
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 1700\n") }

    #  Output table tab.npa initialized in seg_1010 
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
                          !is.na(detect.limits.max) | TI.left.only, 
                          deltap.lo, deltap.hi)

    TI.can.n <- nrow(TI.can)
    if (is.na(TI.can[1,1]) ) TI.can.n <- 0

    #  Minimal proportion and inclusion of range around mode already organized

    if (TI.can.n == 0)
    {
      cat("\n\n [seg100] +++ Not enough truncation interval candidates +++\n\n")

      #  @@@ stop processing here

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
      #  opt.crit < p.fit.min
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

      tmc.t2 <- Sys.time()

      sink(file="../temp/UsedTime.txt", append=TRUE)
      cat("\n tmc; t2-t1;", 
      format(difftime(tmc.t2, tmc.t1)))
      sink()

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

      TmcSolutAcc   <- FALSE
      iti           <- 0      

      while ( (iti < TI.can.n) & !TmcSolutAcc)
      {
        iti <- iti + 1
        tmc.t3 <- Sys.time()
    
        cat("\n [seg100] iti/TI.can.n=", iti,"/", TI.can.n, 
                 " TI.qual=", TI.can[iti, "TI.qual"],
                 sexlabel, agelabel,
                 " -------------") 
        x.tr.n <- TI.can[iti, "x.tr.n"]

        #  Calculate QQ estimate for the data in TI
        prop.lo <- TI.can[iti, "prop.lo"]
        prop.hi <- TI.can[iti, "prop.hi"]
        subset   <- (prop.lo < x.red.cdf) & (x.red.cdf < (1-prop.hi))

        #  Find and keep optimal rsq from QQ plot estimate
        rsq.opt.qqw <- 0

        #  qq loop over lambda candidates 

        #  apply should be faster
        r2.tab <- apply(Y.red, MARGIN=2, FUN=QQEstimateIni,
                        prop.lo=prop.lo, prop.hi=prop.hi, subset=subset, 
                        figA=NA)
        r2.tab <- t(r2.tab)

        #  Add lambda
        r2.tab <- cbind(lambda.seq, r2.tab)
        colnames(r2.tab) <- c("lambda", "r2", "mue", "sigma")

        #  Find optimum
        idx <- which.max(r2.tab[ ,"r2"])

        tmc.t4 <- Sys.time()
        sink(file="../temp/UsedTime.txt", append=TRUE)
        cat("\n tmc; t4-t3;", 
        format(difftime(tmc.t4, tmc.t3)))
        sink()

        # --------------------------------------------------------------------  
        #  Calculate TMC estimate for this interval
        ##  TI limits are implanted into tab.qqw
        tab.qqw[1, "bin.start"] <- TI.can[iti, "ilo"]
        tab.qqw[1, "bin.end"]   <- TI.can[iti, "ihi"]

        #  Uses tab.qqw containing inital values for theta
        tab.qqw[1, "lambda.qqw"] <- r2.tab[idx ,"lambda"]
        tab.qqw[1, "mue.qqw"]    <- r2.tab[idx ,"mue"]
        tab.qqw[1, "sigma.qqw"]  <- r2.tab[idx ,"sigma"]
        tab.qqw[1, "r2"]         <- r2.tab[idx ,"r2"]
        tab.qqw[1, "x.RL1.qqw"]  <- r2.tab[idx ,"sigma"] <- 
              q.PN(RL1.p, tab.qqw[1, "lambda.qqw"], tab.qqw[1, "mue.qqw"], 
                            tab.qqw[1, "sigma.qqw"])
        tab.qqw[1, "x.RL2.qqw"]  <- r2.tab[idx ,"sigma"] <- 
              q.PN(RL2.p, tab.qqw[1, "lambda.qqw"], tab.qqw[1, "mue.qqw"], 
                            tab.qqw[1, "sigma.qqw"])

        # seg105 produces estimate for the current TI. If this is better
        # than the previous best solution, the new solution is the best 
        # solution. 

        # if (print.log.message) { cat("\n%%%   TMC_seg100_Analysis 2500\n") }

        tmc.t5 <- Sys.time()

        source("TMC_seg105_Analysis_tmc.R")

        tmc.t6 <- Sys.time()
        sink(file="../temp/UsedTime.txt", append=TRUE)
        cat("\n tmc; t6-t5;", 
        format(difftime(tmc.t6, tmc.t5)))
        sink()


        # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2600\n") }

        #  returns tab.tmc0
        #  Save this result in tab.tmc (=the table of all results)
        if (iti == 1) { tab.tmc <- tab.tmc0 } else
                      { tab.tmc <- rbind(tab.tmc, tab.tmc0) }

   
        #  Is this result better than the actual optimum?
        if (iti == 1) 
        { 
          # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2700\n") }
          iti.opt      <- 1
          TI.qual.opt  <- TI.can[1, "TI.qual"]
          tab.tmc.opt  <- tab.tmc0 
          # temp.tmc.opt <- temp.tmc0  
        }  else
        { 
          #  iti > 1, compare actual solution with actual optimum

          if (tab.tmc0["sol.score"] == tab.tmc.opt["sol.score"])
          { 
            # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2750\n") }
            # Both scores are identical
            # Further decision by opt.crit
            if (tab.tmc0["opt.crit"] < tab.tmc.opt["opt.crit"])
            { # tmc0 is better
              # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2800\n") }
              tab.tmc.opt  <- tab.tmc0
              # temp.tmc.opt <- temp.tmc0  
                             #  @@@ Ort der Berechnung von temp.tmc prüfen!
              iti.opt      <- iti
              TI.qual.opt  <- TI.can[iti, "TI.qual"]

            cat("\n [seg100] (Improved opt.crit) iti.opt, sol.score, opt.crit", 
                  iti.opt, sol.score, tab.tmc.opt["opt.crit"] )
            }
          }  else   
          { 
            # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2900\n") }
            # both scores not identical
            #  Scores differ, solution with higher score is better
            if (tab.tmc0["sol.score"] > tab.tmc.opt["sol.score"])
            { # tab.tmc0 has higher score, is better
              # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3000\n") }
              tab.tmc.opt  <- tab.tmc0
              #temp.tmc.opt <- temp.tmc0  
              iti.opt      <- iti
              TI.qual.opt  <- TI.can[iti, "TI.qual"]

              cat("\n [seg100] (Improved sol.score) iti.opt, sol.score, opt.crit",
                   iti.opt, sol.score,  tab.tmc.opt["opt.crit"])
            }
            # No action needed, if score(tab.tmc0) <= score(tab.tmc.opt)
          }
        }        #  iti == 1

        # --------------------------------------------------------------------
        #  Is the actual solution acceptable?
        #
        #  New 03.09.2022 
        TmcSolutAcc <- (tab.tmc0["sol.score"] == sol.score.max) & 
                       (tab.tmc0["p.fit"] >= p.fit.min) &
                       (tab.tmc0["p.rt"]  >= p.rt.min)
        TmcSolutAcc <-  unname(TmcSolutAcc)

        #cat("\n [seg100] sol.score, p.fit, p.rt, TmcSolutAcc ",
        #tab.tmc0["sol.score"], 
        #formatC(tab.tmc0["p.fit"], format="f", width=6, digits=4),
        #formatC(tab.tmc0["p.rt"], format="f", width=6, digits=4), 
        #TmcSolutAcc, "\n")

        # --------------------------------------------------------------------

        #  Report assessment
        # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3100\n") }

        if (FALSE)
        { 
          #  Check: which decision taken?
          #  Report solution status
  cat("\n [100_tmc] Actual iti                 ", iti, tab.tmc0["ilo"], 
                                                       tab.tmc0["ihi"],                   
      "\n           initial lambda for TMC     ", tab.qqw[1, "lambda.qqw"],
      "\n           actual theta.tmc           ", tab.tmc0["lambda.tmc"],
                                                  tab.tmc0["mue.tmc"], 
                                                  tab.tmc0["sigma.tmc"], 
      "\n           opt.crit of actual optimum ", tab.tmc0["opt.crit"],
      "\n           score of actual solution   ", tab.tmc0["sol.score"],
      "\n           p.fit of actual solution   ", tab.tmc0["p.fit"],
      "\n           p.rt of actual solution    ", tab.tmc0["p.rt"],
      "\n           opt.crit of global optimum ", tab.tmc.opt["opt.crit"],
      "\n           score of global solution   ", tab.tmc.opt["sol.score"],
      "\n           iti.opt                    ", iti.opt,
      "\n           TmcSolutAcc                ", TmcSolutAcc,
      "\n")
        }
    
      } # while (iti ...    

      cat("\n\n [seg100] End loop iti/TI.can.n  ============================\n") 
      cat("\n")

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3200\n") }

      #  Decide final solution:
      #  If an acceptable solution was found, take it as final solution, 
      #  even if the corresponding opt.crit is not zhe global minimum
      #  If no acceptable solution was found, take the solution with minimal
      #  opt.crit as the final solution

      if ( TmcSolutAcc )
      { cat("\n [seg100] TMC estimation finished ",
            "\n          Acceptable solution found for TI ", iti,
            "\n")
        iti.opt     <- iti
        tab.tmc.opt <- tab.tmc0     
        tab.tmc.opt <- c(iti=iti.opt, TI.qual=TI.qual.opt, tab.tmc.opt)
      } else
      { 
        cat("\n [seg100] TMC estimation finished ",
            "\n          NO ACCEPTABLE SOLUTION FOUND",
            "\n          Solution for TI ", iti.opt, " is final solution",
            "\n")
        tab.tmc.opt <- c(iti=iti.opt, TI.qual=TI.qual.opt, tab.tmc.opt)
      }

      tab.tmc   <- RowToMatrix(tab.tmc)
      tab.tmc.n <- nrow(tab.tmc)

      rownames(tab.tmc) <- NULL

      #cat("\n [seg100] The first five solutions, in calculation order \n")
      #print(tab.tmc[1:min(5, tab.tmc.n), 
      #      c("ilo", "ihi", "x.tr.prop", 
      #        "lambda.tmc", "mue.tmc", "sigma.tmc",
      #        "x.RL1.tmc", "x.RL2.tmc", "p.fit", "p.rt", "opt.crit", 
      #        "sol.score") ] )

      if ( tab.tmc.n > 1)
      { 
        tab.tmc <- tab.tmc[order(-tab.tmc[ ,"sol.score"],
                                  tab.tmc[ ,"opt.crit"],
                                 -tab.tmc[ ,"p.fit"]),  ]
      }

      #cat("\n [seg100] The final solution", "\n")
      #print(tab.tmc.opt[  
      #      c("ilo", "ihi", "x.tr.prop", 
      #        "lambda.tmc", "mue.tmc", "sigma.tmc",
      #        "x.RL1.tmc", "x.RL2.tmc", "p.fit", "p.rt", "opt.crit", 
      #        "sol.score") ] )

      #cat("\n [seg100] The five best opt.crit solutions",
      #    "\n")

      #print(tab.tmc[1:min(5, tab.tmc.n), 
      #      c("ilo", "ihi", "x.tr.prop", 
      #        "lambda.tmc", "mue.tmc", "sigma.tmc",
      #        "x.RL1.tmc", "x.RL2.tmc", "p.fit", "p.rt", "opt.crit", 
      #        "sol.score") ] )

      #cat("\n [seg100] All solutions (first = optimal opt.crit)",
      #    "\n")

      #print(tab.tmc[ , 
      #      c("ilo", "ihi", "x.tr.prop", 
      #        "lambda.tmc", "mue.tmc", "sigma.tmc",
      #        "x.RL1.tmc", "x.RL2.tmc", "p.fit", "p.rt", "opt.crit", 
      #        "sol.score") ] )

      #reld <- sqrt((tab.tmc[ , "x.RL1.tmc"] - 10)^2 +  
      #             (tab.tmc[ , "x.RL2.tmc"] - 35)^2)

      #cat("\n [seg100] All solutions (first = min dist of RLs)",
      #    "\n")

      #print(tab.tmc[order(reld) , 
      #      c("ilo", "ihi", "x.tr.prop", 
      #        "lambda.tmc", "mue.tmc", "sigma.tmc",
      #        "x.RL1.tmc", "x.RL2.tmc", "p.fit", "p.rt", "opt.crit", 
      #        "sol.score") ] )

      # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3300\n") }

      # ----------------------------------------------------------------------
      #  Collect extra information for histogram plotting
      #  Variables required in F_PlotHistResid.R
      #c("UG", "OG", "nobs", "nerw", "diff")

      temp.tmc.opt <- chi2trunc(x.hist,
                             as.numeric(tab.tmc.opt["ilo"]),
                             as.numeric(tab.tmc.opt["ihi"]),
                             TI.can[iti.opt, "x.lt.tr.n"], 
                             TI.can[iti.opt, "x.ge.tr.n"],                             
                             x.Q1, x.Q2, RL1.p, RL2.p,
                             as.numeric(tab.tmc.opt["lambda.tmc"]),
                             as.numeric(tab.tmc.opt["mue.tmc"]),
                             as.numeric(tab.tmc.opt["sigma.tmc"]),
                             df.est, df.con,
                             l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                             rt.fact=rt.fact, w.fact=w.fact,
                             penalty.details=TRUE, 
                             opt.crit.only=FALSE,fastnull=fastnull)

      # Rename tables in temp.tmc.opt to fit plot function names                   
      tt.names <- names(temp.tmc.opt$tab)
      tt.names[1] <- "UG"
      tt.names[2] <- "OG"
      tt.names[3] <- "nobs"
      tt.names[10] <- "nerw"
      tt.names[11] <- "diff"
      names(temp.tmc.opt$tab) <- tt.names
      if (length(temp.tmc.opt$tab.lo) > 1 ) 
         { names(temp.tmc.opt$tab.lo) <- tt.names }
      if (length(temp.tmc.opt$tab.hi) > 1 ) 
         { names(temp.tmc.opt$tab.hi) <- tt.names }

      tmc.t7 <- Sys.time()

      sink(file="../temp/UsedTime.txt", append=TRUE)
      cat("\n tmc; t7-t3;", 
      format(difftime(tmc.t7, tmc.t3)))
      sink()

      # =========================================================================
      #  TMC Estimation is finished
      # =========================================================================

      #  Print and plot result   
      #  Replace pseudo-zero estimate by zero for nicer printout
      if (abs(tab.tmc.opt["lambda.tmc"]) < fastnull) 
      { 
        tab.tmc.opt["lambda.tmc"] <- 0 
      }
    
      TI <- tab.tmc.opt["ilo"] : tab.tmc.opt["ihi"]
      counts.range <- paste(formatC(min(x.hist$counts[TI]), format="f", digits=0),
                        "-", 
                        formatC(max(x.hist$counts[TI]), format="f", digits=0), 
                        sep="")
      x.tr.prop.range <-               
           paste(formatC(100*tab.tmc.opt["x.tr.prop"], format="f", width=5, 
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
        #  Calculate estimated TMC density (not yet scaled by prev.c.tmc)
        xc.pdf.tmc <- pdf.PN(xsupp.clip, tab.tmc.opt["lambda.tmc"], 
                                  tab.tmc.opt["mue.tmc"], tab.tmc.opt["sigma.tmc"])    
  
        PlotMetRes(fig100.016, fig100.016.file, figtype, 
                   fig100.016.lim["yplotmin"],
                   xsupp.clip, tab.tmc.opt["prev.c.tmc"], xc.pdf.tmc, tmccol,
                   x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol, 
                   tab.tmc.opt["x.RL1.tmc"], tab.tmc.opt["x.RL2.tmc"], 
                   0.85*fig100.016.lim["yplotmax"], 
                   tab.tmc.opt["x.tr.lo"], tab.tmc.opt["x.tr.hi"], 
                   0.90*fig100.016.lim["yplotmin"], tmccol,
                   temp.tmc.opt)

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
      # if (print.log.message) 
      # { cat("%%%   TMC_seg100_Analysis 3600", date(), "\n") }

      if ("tmu" %in% meth.list)
      { #  Calculate TMU solution, using lambda from tmc

      source("TMC_seg106_Analysis_tmu.R") 
 
      }  

      # ================================================================
      # if (print.log.message) 
      { cat("%%%   TMC_seg100_Analysis 3700", date(), "\n") }           

      #  Calculate Tukey solution, using lambda from tmc

      if ("tuk" %in% meth.list)
      { #  Calculate Tukey solution
        source("TMC_seg108_Analysis_tuk.R")       
      }  

      # if (print.log.message) 
      # { cat("%%%   TMC_seg100_Analysis 3800", date(), "\n") }           

      # ================================================================

    }                #  x.hist exists

    # ==========================================================================
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis 3900\n") }
    # ==========================================================================

  }    #  TI.can.n > 0?

  #  Combine essential results from all methods for this stratum
  #  tab.stra defined in the beginning

  if ("npa" %in% meth.list)
  {
    tab.stra["npa", "x.tr.lo"] <- tab.npa["x.tr.lo.npa"]
    tab.stra["npa", "x.tr.hi"] <- tab.npa["x.tr.hi.npa"]
    tab.stra["npa", "RL1"]     <- tab.npa["x.RL1.npa"]
    tab.stra["npa", "RL2"]     <- tab.npa["x.RL2.npa"]
    tab.stra["npa", "prop"]    <- tab.npa["x.tr.prop.npa"]
    tab.stra["npa", "prev.l"]  <- tab.npa["prev.l.npa"]
    tab.stra["npa", "prev.c"]  <- tab.npa["prev.c.npa"]
    tab.stra["npa", "prev.r"]  <- tab.npa["prev.r.npa"]
  }

  if ( "qqw" %in% meth.list )
  {
    tab.stra["qqw", "x.tr.lo"] <- tab.qqw[1, "x.tr.lo"]
    tab.stra["qqw", "x.tr.hi"] <- tab.qqw[1, "x.tr.hi"]
    tab.stra["qqw", "bins.n"]  <- tab.qqw[1, "bins.n"]
    tab.stra["qqw", "prop"]    <- tab.qqw[1, "prop.qqw"]
    tab.stra["qqw", "lambda"]  <- tab.qqw[1, "lambda.qqw"]
    tab.stra["qqw", "mue"]     <- tab.qqw[1, "mue.qqw"]
    tab.stra["qqw", "sigma"]   <- tab.qqw[1, "sigma.qqw"]
    tab.stra["qqw", "RL1"]     <- tab.qqw[1, "x.RL1.qqw"]
    tab.stra["qqw", "RL2"]     <- tab.qqw[1, "x.RL2.qqw"]
    tab.stra["qqw", "RL1.cilo"] <- tab.qqw[1, "x.RL1.cilo"]
    tab.stra["qqw", "RL1.cihi"] <- tab.qqw[1, "x.RL1.cihi"]
    tab.stra["qqw", "RL2.cilo"] <- tab.qqw[1, "x.RL2.cilo"]
    tab.stra["qqw", "RL2.cihi"] <- tab.qqw[1, "x.RL2.cihi"]
    tab.stra["qqw", "bins.n"]  <- tab.qqw[1, "bins.n"]

    if (x.n < red.limit)
    { #   Next 2 entries only if QQW used real data, not a reduced version
      tab.stra["qqw", "x.tr.n"]  <- tab.qqw[1, "x.tr.n"] 
      tab.stra["qqw", "prop"]    <- tab.qqw[1, "prop.qqw"]
    }

    # Fit criteria not yet calculated

    temp.qqw <- chi2trunc(x.hist,tab.qqw[1, "bin.start"],
                          tab.qqw[1, "bin.end"],
                             sum(x <  tab.qqw[1, "x.tr.lo"]), 
                             sum(x >= tab.qqw[1, "x.tr.hi"]),
                             x.Q1, x.Q2, RL1.p, RL2.p,
                             tab.qqw[1, "lambda.qqw"],
                             tab.qqw[1, "mue.qqw"],
                             tab.qqw[1, "sigma.qqw"],
                             df.est, df.con,
                             l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                             rt.fact=rt.fact, w.fact=w.fact, 
                             opt.crit.only=FALSE,fastnull=fastnull)

    tab.stra["qqw", "opt.crit"] <- unname(temp.qqw[["res"]]["opt.crit"])
    tab.stra["qqw", "p.fit"]   <- unname(temp.qqw[["res"]]["chi2.trun.p"])
    tab.stra["qqw", "prev.l"]  <- temp.qqw$res["prev.l.tmc"]
    tab.stra["qqw", "prev.c"]  <- temp.qqw$res["prev.c.tmc"]
    tab.stra["qqw", "prev.r"]  <- temp.qqw$res["prev.r.tmc"]
  }

  if ("bha" %in% meth.list)
  {
    #  @@@ subpopulation auswählen
    tab.stra["bha", "RL1"]    <- tab.bha[["subpop"]]["comp01", "RL1.2"]
    tab.stra["bha", "RL2"]     <- tab.bha[["subpop"]]["comp01", "RL2.2"]
  }

  if ("tmc" %in% meth.list)
  {
    tab.stra["tmc", "x.tr.lo"] <- tab.tmc.opt["x.tr.lo"]
    tab.stra["tmc", "x.tr.hi"] <- tab.tmc.opt["x.tr.hi"]
    tab.stra["tmc", "x.tr.n"]  <- tab.tmc.opt["x.tr.n"]
    tab.stra["tmc", "bins.n"]  <- tab.tmc.opt["ihi"] - tab.tmc.opt["ilo"] + 1
    tab.stra["tmc", "prop"]    <- tab.tmc.opt["x.tr.prop"]  
    tab.stra["tmc", "rc"]      <- tab.tmc.opt["rc"]  
    tab.stra["tmc", "iter"]    <- tab.tmc.opt["iter"]  
    tab.stra["tmc", "lambda"]  <- tab.tmc.opt["lambda.tmc"]
    tab.stra["tmc", "mue"]     <- tab.tmc.opt["mue.tmc"]
    tab.stra["tmc", "sigma"]   <- tab.tmc.opt["sigma.tmc"]
    tab.stra["tmc", "RL1"]     <- tab.tmc.opt["x.RL1.tmc"]
    tab.stra["tmc", "RL2"]     <- tab.tmc.opt["x.RL2.tmc"]

    tab.stra["tmc", "RL1.cilo"] <- tab.tmc.opt["x.RL1.cilo"]
    tab.stra["tmc", "RL1.cihi"] <- tab.tmc.opt["x.RL1.cihi"]
    tab.stra["tmc", "RL2.cilo"] <- tab.tmc.opt["x.RL2.cilo"]
    tab.stra["tmc", "RL2.cihi"] <- tab.tmc.opt["x.RL2.cihi"]

    tab.stra["tmc", "opt.crit"] <- tab.tmc.opt["opt.crit"]
    tab.stra["tmc", "p.fit"]   <- tab.tmc.opt["p.fit"]
    tab.stra["tmc", "p.rt"]    <- tab.tmc.opt["p.rt"]
    tab.stra["tmc", "bins.n"]  <- tab.tmc.opt["ihi"] - 
                                  tab.tmc.opt["ilo"] + 1
    tab.stra["tmc", "prop"]    <- tab.tmc.opt["x.tr.prop"]
    tab.stra["tmc", "prev.l"]  <- tab.tmc.opt["prev.l.tmc"]
    tab.stra["tmc", "prev.c"]  <- tab.tmc.opt["prev.c.tmc"]
    tab.stra["tmc", "prev.r"]  <- tab.tmc.opt["prev.r.tmc"]
  }

  if ("tmu" %in% meth.list)
  {
    tab.stra["tmu", "x.tr.lo"] <- tab.tmc.opt["x.tr.lo"]
    tab.stra["tmu", "x.tr.hi"] <- tab.tmc.opt["x.tr.hi"]
    tab.stra["tmu", "x.tr.n"]  <- tab.tmc.opt["x.tr.n"]
    tab.stra["tmu", "bins.n"]  <- tab.tmc.opt["ihi"] - tab.tmc.opt["ilo"] + 1
    tab.stra["tmu", "prop"]    <- tab.tmc.opt["x.tr.prop"]  
    tab.stra["tmu", "rc"]      <- tab.tmu["rc.tmu"]  
    tab.stra["tmu", "iter"]    <- tab.tmu["iter.tmu"]  
    tab.stra["tmu", "lambda"]  <- tab.tmc.opt["lambda.tmc"]
    tab.stra["tmu", "mue"]     <- tab.tmu["mue.tmu"]
    tab.stra["tmu", "sigma"]   <- tab.tmu["sigma.tmu"]
    tab.stra["tmu", "RL1"]     <- tab.tmu["x.RL1.tmu"]
    tab.stra["tmu", "RL2"]     <- tab.tmu["x.RL2.tmu"]

    tab.stra["tmu", "RL1.cilo"] <- tab.tmu["x.RL1.cilo"]
    tab.stra["tmu", "RL1.cihi"] <- tab.tmu["x.RL1.cihi"]
    tab.stra["tmu", "RL2.cilo"] <- tab.tmu["x.RL2.cilo"]
    tab.stra["tmu", "RL2.cihi"] <- tab.tmu["x.RL2.cihi"]
    tab.stra["tmu", "opt.crit"] <- tab.tmu["opt.crit.tmu"]
    tab.stra["tmu", "p.fit"]   <- tab.tmu["p.fit.tmu"]
    tab.stra["tmu", "prev.l"]  <- tab.tmu["prev.l.tmu"]
    tab.stra["tmu", "prev.c"]  <- tab.tmu["prev.c.tmu"]
    tab.stra["tmu", "prev.r"]  <- tab.tmu["prev.r.tmu"]
  }                #  do tmu

  if ("tuk" %in% meth.list)
  {
    tab.stra["tuk", "x.tr.lo"] <- tab.tuk["x.tr.lo"]
    tab.stra["tuk", "x.tr.hi"] <- tab.tuk["x.tr.hi"]
    tab.stra["tuk", "x.tr.n"]  <- tab.tuk["x.tr.n"]
    tab.stra["tuk", "bins.n"]  <- NA
    tab.stra["tuk", "prop"]    <- (x.n-tab.tuk["out.lo"]-tab.tuk["out.hi"])/ 
                                  x.n  
    tab.stra["tuk", "rc"]      <- NA 
    tab.stra["tuk", "iter"]    <- NA  
    tab.stra["tuk", "lambda"]  <- tab.tmc.opt["lambda.tmc"]
    tab.stra["tuk", "mue"]     <- tab.tuk["mean"]
    tab.stra["tuk", "sigma"]   <- tab.tuk["sd.adj"]
    tab.stra["tuk", "RL1"]     <- tab.tuk["x.RL1.tuk"]
    tab.stra["tuk", "RL2"]     <- tab.tuk["x.RL2.tuk"]

    tab.stra["tuk", "RL1.cilo"] <- NA
    tab.stra["tuk", "RL1.cihi"] <- NA
    tab.stra["tuk", "RL2.cilo"] <- NA
    tab.stra["tuk", "RL2.cihi"] <- NA
    tab.stra["tuk", "opt.crit"] <- tab.tuk["opt.crit.tuk"]
    tab.stra["tuk", "p.fit"]    <- tab.tuk["p.fit.tuk"]
    tab.stra["tuk", "prev.l"]   <- tab.tuk["prev.l.tuk"]
    tab.stra["tuk", "prev.c"]   <- tab.tuk["prev.c.tuk"]
    tab.stra["tuk", "prev.r"]   <- tab.tuk["prev.r.tuk"]
  }                #  do tuk

  # --------------------------------------------------------------------------
  #  Error code for result. Considers goodness of fit (p.fit)
  #  and further properties of the solution. Version 2022-08-05.
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 4000\n") }


  tab.stra[ , "errcode"] <- "!!"   # should never appear

  ok1 <- (tab.stra[ , "p.fit"] > p.fit.min) 
  ok2 <- (tab.stra[ , "p.rt"]  > p.rt.min) 

  ok1[is.na(ok1)] <- FALSE
  ok2[is.na(ok2)] <- FALSE

  tab.stra[ok1 & ok2, "errcode"]   <- "++"
  tab.stra[ok1 & !ok2, "errcode"]  <- "+-"
  tab.stra[!ok1 & ok2, "errcode"]  <- "-+"
  tab.stra[!ok1 & !ok2, "errcode"] <- "--"

  # No mark, if p.fit an p.rt not determined
  tab.stra[is.na(tab.stra[ , "p.fit"]) |
           is.na(tab.stra[ , "p.rt"])    , "errcode"] <- ".."
  tab.stra[ , "errcode"]

  # --------------------------------------------------------------------------
  #  Toleranzintervalle für die tmc-RLs per bootstrapping berechnen
  #  NOT YET ADAPTED TO NEW CALL OF tmc.master

  source("TMC_seg109_TMC_BS.R")

  # --------------------------------------------------------------------------

  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 5000\n") }

 
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
      "\nTMC: Truncation interval               ",tab.tmc.opt["x.tr.lo"], 
                                                  tab.tmc.opt["x.tr.hi"],   
      "\nTMC: TI contains proportion ...        ",tab.tmc.opt["x.tr.prop"],
      "\nTMC: estimated lambda                  ",tab.tmc.opt["lambda.tmc"],
      "\nTMC: estimated mu                      ",tab.tmc.opt["mue.tmc"],
      "\nTMC: estimated sigma                   ",tab.tmc.opt["sigma.tmc"],
      "\nTMC: chi^2 from truncation interval    ",format(
                                                    tab.tmc.opt["chi2.trun"],
                                                    digits=4),
      "\nTMC: chi^2 from outside TI             ",format(
                                                   tab.tmc.opt["chi2.total"]-
                                                   tab.tmc.opt["chi2.trun"],
                                                   digits=4),
      "\n",
      "\nFalse negative probability  <RL1 / <RL2", RL1.p,"/ ",RL2.p,
      "\nTMC: estimated RLs                       RL1 = ", 
                                             format(tab.tmc.opt["x.RL1.tmc"],
                                                    digits=4),      
      "\n                                         RL2 = ", 
                                             format(tab.tmc.opt["x.RL2.tmc"],
                                                    digits=4), 
      "\nAssessment of TMC solution",
      "\np for chi^2 goodness of fit            ",format(
                                                    tab.tmc.opt["p.fit"],
                                                    digits=4), 
      "\np for runs.test                        ",format(
                                                    tab.tmc.opt["p.rt"],
                                                    digits=4), 
      "\nRL1 >= x.Q1, RL2 <= x.Q2 ?             ", tab.tmc.opt["crit1"],  
                                                   tab.tmc.opt["crit2"], 
      "\nprev.l ok / prev.c ok /prev.r ok?      ", tab.tmc.opt["crit3"],
                                                   tab.tmc.opt["crit4"], 
                                                   tab.tmc.opt["crit5"], 
      "\nTotal score for solution (max =",sol.score.max,")    ", 
                                                tab.tmc.opt["sol.score"],    
      "\nAcceptable solution found?)            ", 
                                                TmcSolutAcc,    
      "\n", sep.e, 
      "\n\n")


  if (any(ErrMsg[ ,"Occurred"]))
  { cat("\n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
    print(ErrMsg[ErrMsg[ ,"Occurred"], "ErrTxt"])
    cat(   "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  } 

  if (user == "WW")
  { #  Ausführung für WW

    #cat("\n All details \n")
    #print(tab.stra[meth.list, ])
    #cat("\n")
  }

  if (!eval.rep) { sink() }
 
  # ==========================================================================
  # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6000\n") }
  # ==========================================================================

}  #  if (x.n < x.n.min)

# ---------------------------------------------------------------------------
#  Copy results in tab.stra (one stratum) to tab.f (one file)
#  Preliminary: gtab instead of tab.f

# ==========================================================================
# if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6100\n") }
# ==========================================================================

#  block.start controls the identification by sexlabel and agelabel in the 
#  printout

#  Put only rows with information in gtab

if (!is.na(tab.stra[1, "irep"]) )
{
  block.start          <- TRUE
  iblock               <- iblock + 1
  for (meth in gtab.list)
  { igtab              <- igtab + 1

    meth.idx <- tab.stra[ ,"method"] == meth
  
    #  Append tab.stra to gtab
    if (iblock == 1 & igtab == 1)
    {
      #  First block at all, establish gtab
      gtab <- tab.stra[meth.idx, ]
    }  else
    {
      #  Append block to gtab
      gtab <- rbind(gtab, tab.stra[meth.idx, ])
    }

    if (block.start)
    { 
      gtab[igtab,"Sex"]     <- sexlabel
      gtab[igtab,"Age"]     <- agelabel
    } else 
    { 
      gtab[igtab,"Sex"]     <- " "
      gtab[igtab,"Age"]     <- " "
    }
    block.start          <- FALSE

    # ==========================================================================
    # if (print.log.message) { cat("%%%   TMC_seg100_Analysis  Position 6400\n") }
    # ==========================================================================

  }    #  information for stratum is available
  cat("\n", sep.e, 
    "\n   End analysis by seg100 for", 
    "\n   ", outname.stra,
    "\n", sep.e, 
    "\n", sep="")
   
}      #  x.n > x.n.min

#############################################################################