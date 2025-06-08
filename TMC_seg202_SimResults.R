#  TMC_seg202_SimResults.R
#  (c) wwosniok@math.uni-bremen.de

#  Presentation of results if a sequence of data files is evaluated, 
#  typically when analysing simulated data

#  INPUT:
#  soutfile.tmc.csv   contains results from bha, npa, qqw, tmc, tmu
#                     usually produced by TMC
#  soutfile.tml.csv   contains results from tml. Produced by FAs
#                     RLE, modified for this output. May not exist.

#  OUTPUT
#  outname.scen.rep.csv   contains results from all methods, per replicate
#  outname.scen.sum.csv   summary over all replicates, per method


#  CHANGE HISTORY
#  17.04.2021 csv output files for a scenario new organized. New file names,
#             the *.rep file has only results specific to a replicate,
#             the *.sum file has all control values
#  11.12.2020 New organisation with figure production by dedicated segments/
#             functions. See TMC4 logbook, ch. 7.2 .
#
# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg202_SimResults  Start\n") }
# ===========================================================================

source("TMC_seg201_OpenWindows.R")

#  Results from TMC are still in 'stab' with a copy in soutfile.tmc.csv
#  Read this file for easier development and put in stab.TMC
#  'stab' will lateron be used for joining results from all methods
#  '.TMC'  indicates that the data frame is calculated by TMC (but usually
#  contains more than method=tmc).
 
stab.TMC <- read.table(file=soutfile.tmc.csv, header=TRUE, sep=";", dec=".",
                       stringsAsFactors=FALSE)
stab.TMC.nrow  <- nrow(stab.TMC)

#  Construct meth.list.202 from stab.TMC. TML not contained in stab.TMC.
meth.list.202   <- sort(unique(stab.TMC[ ,"method"]))
meth.list.202.n <- length(meth.list.202)

#  If requested: check if there is a TML result and join it, if present.

if (do.tml & file.exists(soutfile.tml.csv))  
{ # yes, read it
  stab.tml0 <- read.table(file=soutfile.tml.csv, header=TRUE, sep=";", dec=".",
                         stringsAsFactors=FALSE)

  #  Keep TML results only for replicates in stab.TMC
  rep.list.tmc <- sort(unique(stab.TMC[ , "irep"]))
  rep.list.tml <- stab.tml0[ , "irep"]
  rep.list     <- intersect(rep.list.tmc, rep.list.tml)
  
  #  Reduce rep.list.tmc, rep.list.tml to the common replicates in rep.list
  ok.tmc <- rep(FALSE, times=nrow(stab.TMC))
  ok.tml <- rep(FALSE, times=nrow(stab.tml0))
 
  for (irep in 1:length(rep.list))
  { 
    ok.tmc <- ok.tmc | (stab.TMC[ ,"irep"] == rep.list[irep])
    ok.tml <- ok.tml | (stab.tml0[ ,"irep"] == rep.list[irep])
  }
  stab.TMC <- stab.TMC[ok.tmc, ]
  stab.tml0 <- stab.tml0[ok.tml, ]

  stab.TMC.nrow  <- nrow(stab.TMC)
  stab.tml0.nrow <- nrow(stab.tml0)

  #  Join TMC and TML as stab, noting that tml has less variables than stab
  stab     <- stab.TMC
  stab.tml <- data.frame(matrix(NA, nrow=stab.tml0.nrow, ncol=ncol(stab)))
  colnames(stab.tml)     <- colnames(stab)
  stab.tml[ , "irep"]    <- stab.tml0[ ,"irep"]
  stab.tml[ , "method"]  <- "tml"
  stab.tml[ , "RL1"]     <- stab.tml0[ ,"RL1"]
  stab.tml[ , "RL2"]     <- stab.tml0[ ,"RL2"]
  stab.tml[ , "lambda"]  <- stab.tml0[ ,"lambda"]
  stab.tml[ , "mue"]     <- stab.tml0[ ,"mue"]
  stab.tml[ , "sigma"]   <- stab.tml0[ ,"sigma"]

  stab <- rbind(stab,stab.tml)
  stab <- stab[order(stab[ ,"irep"], stab[ ,"method"]), ]
 
  meth.list.202   <- c(meth.list.202, "tml")
  meth.list.202.n <- length(meth.list.202)
}

# ----------------------------------------------------------------------------
#  Select methods to display. All methods listed in gtab.list plus TML (if
#  requested and existent) will be shown in the summary plots for simulations. gtab.list 
#  is set in seg_015.


# ----------------------------------------------------------------------------
#  Create summarizing table over all replicates, one row per method

scen.sum.colnames <- c("n.per.bin.min","prop.qqw",
                     "prop.mw","prop.sd","prop.min","prop.max",
                     "RL1.mw","RL1.med","RL1.sd","RL1.min","RL1.max",
                     "RL2.mw","RL2.med","RL2.sd","RL2.min","RL2.max",
                     "cov.RL1RL2",
                     "prev.l.mw","prev.l.sd","prev.l.min","prev.l.max",
                     "prev.c.mw","prev.c.sd","prev.c.min","prev.c.max",
                     "prev.r.mw","prev.r.sd","prev.r.min","prev.r.max",
                     "lambda.mw","lambda.sd","lambda.min","lambda.max",
                     "mue.mw","mue.sd","mue.min","mue.max",
                     "sig.mw","sig.sd","sig.min","sig.max",
                     "cov.muesig",
                     "p.fit.mw","p.fit.sd","p.fit.min","p.fit.max")
scen.sum <- matrix(NA,nrow=meth.list.202.n,ncol=length(scen.sum.colnames))
rownames(scen.sum) <- meth.list.202
colnames(scen.sum) <- scen.sum.colnames

#  Find blocks with missing estimates. If too few data points or bins,
#  tmc result (and others) is missing.

is.ok <- rep(TRUE, times=nrep*meth.list.202.n)

for (i in 1+seq(0,(nrep-1)*meth.list.202.n, by=meth.list.202.n))
{ 
  idx <- i:(i+meth.list.202.n-1)
  if (any(is.na(stab[idx, "RL1"]))) 
  {
    is.ok[idx] <- FALSE
  }
}

#  Reduce stab accordingly
stab1 <- stab[is.ok, ]
nrep1 <- nrow(stab1)/meth.list.202.n

cat("\n Planned # of replicates ", nrep,
    "\n Successful estimates    ", nrep1,
    "\n")

subtitle.202 <- paste(RunId, " / ",datafile0, "_", spalte.w,
                     " (r",r.start," - r",r.ende,") / n=", x.n, sep="")

# ----------------------------------------------------------------------------
#  More additional measures describing the estimate
#  Calculate the False Positive Rates per method
#  FPR refers to the generating parameters. There is no way to use npa as
#  reference, because this method generates no distribution parameters.

#  yc.mode.gen is the location parameter of the PND = yc.mue.gen

FPR1               <- cdf.PN(stab1[ ,  "RL1"],
                                 lambda.c.gen, yc.mode.gen, yc.sig.gen)
FPR2               <- 1 - cdf.PN(stab1[ ,  "RL2"],
                                 lambda.c.gen, yc.mode.gen, yc.sig.gen)
FPR.sum            <- abs(FPR1) + abs(FPR2)
delta.FPR          <- abs(FPR1-alpha/2) + abs(FPR2-alpha/2)

FPR1      <- 100 * FPR1
FPR2      <- 100 * FPR2
FPR.sum   <- 100 * FPR.sum
delta.FPR <- 100 * delta.FPR

#  Switch to FPR in %
stab1[ , "FPR1"]      <- FPR1
stab1[ , "FPR2"]      <- FPR2
stab1[ , "FPR.sum"]   <- FPR.sum
stab1[ , "delta.FPR"] <- delta.FPR

# ============================================================================
#  Data is organized.

# ----------------------------------------------------------------------------
#  TMC_seg202_SimResultsA.R
#  Generates 
#  fig202.010 Joint distribution of RL1, RL2 in original units, 
#             no contour background
 
source("TMC_seg202_SimResultsA.R")

source("TMC_seg203_SimResultsA.R")

# ----------------------------------------------------------------------------
#  TMC_seg202_SimResultsD.R (B and C are below)
#  Generates 
#  fig202.011 (contour) Joint distribution of RL1, RL2 in original 
#                       units, contour background, plot limits according
#                       to estimated RL ranges 
#  fig202.012 legend for fig202.011

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.011))
{
  figA.dev <- fig202.011
  figB.dev <- fig202.012
  figA.file <- fig202.011.file
  figB.file <- fig202.012.file

  show.pD         <- FALSE
  show.extra.axes <- FALSE
  
  xplotmin <- min(as.numeric(stab1[ ,"RL1"]),na.rm=TRUE)
  xplotmax <- max(as.numeric(stab1[ ,"RL1"]),na.rm=TRUE)
  yplotmin <- min(as.numeric(stab1[ ,"RL2"]),na.rm=TRUE)
  yplotmax <- max(as.numeric(stab1[ ,"RL2"]),na.rm=TRUE)

  source("TMC_seg202_SimResultsD.R")
}

# ----------------------------------------------------------------------------
#  Fig202e: Joint distribution of RL1, RL2 in original units, contour 
#  background, plot limits according to estimated permissible limits.
#  Generates fig202.013 (contour) and fig202.014 (legend)

#  TMC_seg202_SimResultsD.R (B and C are below)
#  Generates 
#  fig202.013 (contour) Joint distribution of RL1, RL2 in original units,
#                       contour background, plot limits according to estimated
#                       permissible limits.
#  fig202.014 legend for fig202.013

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.013))
{
  figA.dev <- fig202.013
  figB.dev <- fig202.014
  figA.file <- fig202.013.file
  figB.file <- fig202.014.file

  show.pD         <- TRUE
  show.extra.axes <- FALSE

  pD1 <- PermissibleDiff(xc.RL1.gen, xc.RL2.gen, xc.RL1.gen, ef)
  pD2 <- PermissibleDiff(xc.RL1.gen, xc.RL2.gen, xc.RL2.gen, ef)

  xplotmin <- unname(as.numeric(Round(pD1[1] - round.unit, round.unit)))
  xplotmax <- xc.RL1.gen + (xc.RL1.gen - xplotmin) 
  yplotmax <- unname(as.numeric(Round(pD2[2] + round.unit, round.unit)))
  yplotmin <- xc.RL2.gen - (yplotmax - xc.RL2.gen)

  source("TMC_seg202_SimResultsD.R")
}

# ----------------------------------------------------------------------------
#  TMC_seg202_SimResultsD.R

#  Generates 
#  fig202.015 (contour) Joint distribution of RL1, RL2 in original units, 
#                       contour background, plot limits according to user 
#                       specification.
#  fig202.016 legend for fig202.015

#  All summaries only if there is > 1 replicate left

if ((nrep1 > 1) & !is.na(fig202.015) & (!is.na(RL1.clip.min)))
{
  figA.dev <- fig202.015
  figB.dev <- fig202.016
  figA.file <- fig202.015.file
  figB.file <- fig202.016.file

  show.pD         <- TRUE      # may not be possible
  show.extra.axes <- TRUE

  pD1 <- PermissibleDiff(xc.RL1.gen, xc.RL2.gen, xc.RL1.gen)
  pD2 <- PermissibleDiff(xc.RL1.gen, xc.RL2.gen, xc.RL2.gen)

  xplotmin <- RL1.clip.min     
  xplotmax <- RL1.clip.max     
  yplotmin <- RL2.clip.min     
  yplotmax <- RL2.clip.max     

  source("TMC_seg202_SimResultsD.R")
}
# ----------------------------------------------------------------------------
#  TMC_seg202_SimResultsB.R
#  Generates
#  fig202.020  Joint distribution of RL1, RL2 expressed as % of the npa 
#              results 
 
source("TMC_seg202_SimResultsB.R")

# ----------------------------------------------------------------------------
#  TMC_seg202_SimResultsC.R
#  Generates fig202.030 Absolute False Positive Rates per replicate
#            fig202.040 PDF of False Positive Rates
#            fig202.050 CDF of absolute sum of FPR errors

source("TMC_seg202_SimResultsC.R")

# ----------------------------------------------------------------------------
#  More approaches for graphical presentation: see seg201

# -------------------------------------------------------------------------
cat("\nZusammenfassung der Ergebnisse aus allen Replikationen",
    "\n   round.unit      ", round.unit,
    "\n   smooth.hist1    ", smooth.hist1, 
    "\n   smooth.hist2    ", smooth.hist2, 
    "\n   l.fact          ", l.fact    ,
    "\n   p.fact          ", p.fact    ,
    "\n   r.fact          ", r.fact    ,
    "\n   s.fact          ", s.fact    ,
    "\n   w.fact          ", w.fact    ,
    "\n   n.per.bin.min   ", n.per.bin.min,
    "\n   bins.n.min      ", bins.n.min,
    "\n   bins.n.max      ", bins.n.max,
    "\n   x.tr.prop.min   ", x.tr.prop.min,
    "\n   x.tr.prop.max   ", x.tr.prop.max,
    "\n   p.fit.min       ", p.fit.min,
    "\n   crit6.p         ", crit6.p,
    "\n\n") 
print(scen.sum)


#  Einzelergebnisse als csv ("method" ist bereits eine Spalte,
#  Zeilennamen sind laufende Nr., nicht rausschreiben))
write.table(stab1,file=outname.scen.rep.csv,quote=FALSE,sep=";",dec=".",  
            append=FALSE, row.names=FALSE)
cat("\nResults per replication data set written to",
    "\n",outname.scen.rep.csv,
    "\n")

#  Damit die csv-Datei die richtigen Spaltennamen bekommt:
#  Zeilennamen als erste Spalte mit Namen "method" einbauen, dann
#  ohne Zeilennamen rausschreiben (Zeilennamen erzeugen Spalte ohne
#  Spaltennamen)
scen.sum.csv <- scen.sum
scen.sum.csv <- data.frame("method"=row.names(scen.sum),scen.sum)
scen.sum.csv.nrow <- nrow(scen.sum.csv)
scen.sum.csv <- data.frame(time.start=rep(time.start,times=scen.sum.csv.nrow),
                           r.start=rep(r.start,times=scen.sum.csv.nrow),
                           r.ende=rep(r.ende,times=scen.sum.csv.nrow),
                           round.unit=rep(round.unit,times=scen.sum.csv.nrow),
                           l.fact=rep(l.fact,times=scen.sum.csv.nrow),
                           p.fact=rep(p.fact,times=scen.sum.csv.nrow),
                           w.fact=rep(w.fact,times=scen.sum.csv.nrow),
                           n.per.bin.min=rep(n.per.bin.min,
                                             times=scen.sum.csv.nrow),
                           scen.sum.csv) 

write.table(scen.sum.csv,file=outname.scen.sum.csv,append=TRUE,
              quote=FALSE,sep=";",dec=".",row.names=FALSE) 

cat("\nSummary over all replication datasets written to",
      "\n",outname.scen.sum.csv,
      "\n\n")

#  COmpare 

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg202_SimResults  End\n") }
# ===========================================================================

##############################################################################