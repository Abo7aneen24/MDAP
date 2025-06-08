#  TMC_seg052_MasterTables.R 

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Produce general result tables

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg057_MasterTables  Start\n") }
# ===========================================================================

  #  Detailed tables
  #  Variables related to bootstrapping deactivated

  # -------------------------------------------------------------------------
  #  Create table for printing. Here: keep both tables (gtab and gtab.print), 
  #  but remove empty rows
  #  Produce table containing only non-empty rows
  gtab         <- gtab[!is.na(gtab[ , "method"]), ] 
  gtab.print   <- gtab
  gtab.print.n <- nrow(gtab.print)

  # -------------------------------------------------------------------------

  #  Calculate variables that result from existing estimates
  if (!is.na(xc.RL1.gen))
  { #  Relative difference beween observed and generated
    gtab[ ,"delta.RL1"] <- (gtab[ ,"x.RL1"]- xc.RL1.gen)/xc.RL1.gen
    gtab[ ,"delta.RL2"] <- (gtab[ ,"x.RL2"]- xc.RL2.gen)/xc.RL2.gen
  }

  #  Alternative PD calculation 28.02.2022
  #  PD calculation based on estimated RLs. For all existing RLs.
  # PermissibleDiff <- function(RL1, RL2, xi)
  
  for (i in 1:gtab.print.n)
  { if(!is.na(gtab[i ,"x.RL1"]))
    { pd <- PermissibleDiff(gtab[i ,"x.RL1"], gtab[i ,"x.RL2"], 
                            gtab[i ,"x.RL1"], ef)
      gtab[i, "RL1.pdlo"] <- pd["xi.lo"]
      gtab[i, "RL1.pdhi"] <- pd["xi.hi"]

      pd <- PermissibleDiff(gtab[i ,"x.RL1"], gtab[i ,"x.RL2"], 
                            gtab[i ,"x.RL2"], ef)
      gtab[i, "RL2.pdlo"] <- pd["xi.lo"]
      gtab[i, "RL2.pdhi"] <- pd["xi.hi"]
    }
  }

  # -------------------------------------------------------------------------
  #  Format some numbers for nicer printing
  #  width: size of output field
  #  digits: number of digits to show. Overrides width, if necessary.
  #  flag="#": show trailing zeroes  
  #  FC is shortcut for formatC(x, format="f", width=w, digits=d, flag="#")
  #  including correct treatment of NA:
  #  FC <- function(x, w, d)

  gtab.print[ ,"Age.mea"]       <- FC(gtab[ ,"Age.mea"], w=NA, d=1)

  gtab.print[ ,"chi2.total"]    <- format(gtab[ ,"chi2.total"],digits=4)
  gtab.print[ ,"chi2.total.df"] <- format(gtab[ ,"chi2.total.df"],digits=4)
  gtab.print[ ,"chi2.trun"]     <- format(gtab[ ,"chi2.trun"],digits=4)
  gtab.print[ ,"chi2.trun.df"]  <- format(gtab[ ,"chi2.trun.df"],digits=4)
  gtab.print[ ,"chi2.path"]     <- format(gtab[ ,"chi2.path"],digits=4)
  gtab.print[ ,"opt.crit"]      <- FC(gtab[ ,"opt.crit"], w=NA, d=2) 

  gtab.print[ ,"p.fit"]         <- FC(gtab[ ,"p.fit"], w=6, d=4) 
  gtab.print[ ,"p.rt"]          <- FC(gtab[ ,"p.rt"], w=6, d=4) 
  gtab.print[ ,"pct.lt.DL"]     <- FC(gtab[ ,"pct.lt.DL"], w=6, d=1)
  gtab.print[ ,"prev.l"]        <- FC(gtab[ ,"prev.l"], w=6, d=3)
  gtab.print[ ,"prev.c"]        <- FC(gtab[ ,"prev.c"], w=6, d=3)
  gtab.print[ ,"prev.r"]        <- FC(gtab[ ,"prev.r"], w=6, d=3)
  gtab.print[ ,"x.RL1"]         <- FC(gtab[ ,"x.RL1"],w=NA, d=3)
  gtab.print[ ,"x.RL2"]         <- FC(gtab[ ,"x.RL2"],w=NA, d=3)
  gtab.print[ ,"x.tr.prop"]     <- FC(gtab[ ,"x.tr.prop"], w=6, d=3)
  gtab.print[ ,"delta.RL1"]     <- FC(gtab[ ,"delta.RL1"], w=NA, d=3)
  gtab.print[ ,"delta.RL2"]     <- FC(gtab[ ,"delta.RL2"], w=NA, d=3)

  gtab.print[ ,"RL1.cilo"]      <- FC(gtab[ ,"x.RL1.cilo"], w=NA, d=3)
  gtab.print[ ,"RL1.cihi"]      <- FC(gtab[ ,"x.RL1.cihi"], w=NA, d=3)
  gtab.print[ ,"RL2.cilo"]      <- FC(gtab[ ,"x.RL2.cilo"], w=NA, d=3)
  gtab.print[ ,"RL2.cihi"]      <- FC(gtab[ ,"x.RL2.cihi"], w=NA, d=3)
   
  #  Do formatting 
  gtab.print[ ,"RL1.pdlo"]      <- FC(gtab[ ,"RL1.pdlo"], w=NA, d=3)
  gtab.print[ ,"RL1.pdhi"]      <- FC(gtab[ ,"RL1.pdhi"], w=NA, d=3)
  gtab.print[ ,"RL2.pdlo"]      <- FC(gtab[ ,"RL2.pdlo"], w=NA, d=3)
  gtab.print[ ,"RL2.pdhi"]      <- FC(gtab[ ,"RL2.pdhi"], w=NA, d=3)

  #  Bootstrapping deactivated
  #gtab.print[ ,"RL1.bs.mea"] <- FC(gtab[ ,"RL1.bs.mea"], w=NA, d=4)
  #gtab.print[ ,"RL1.tilo"]   <- FC(gtab[ ,"RL1.tilo"], w=NA, d=4)
  #gtab.print[ ,"RL1.tihi"]   <- FC(gtab[ ,"RL1.tihi"], w=NA, d=4)
  #gtab.print[ ,"RL2.bs.mea"] <- FC(gtab[ ,"RL2.bs.mea"], w=NA, d=4)
  #gtab.print[ ,"RL2.tilo"]   <- FC(gtab[ ,"RL2.tilo"], w=NA, d=4)
  #gtab.print[ ,"RL2.tihi"]   <- FC(gtab[ ,"RL2.tihi"], w=NA, d=4)

  gtab.print[ ,"lambda"]   <- FC(gtab[ ,"lambda"], w=6, d=3)
  gtab.print[ ,"mue"]      <- FC(gtab[ ,"mue"], w=6, d=3)
  gtab.print[ ,"sigma"]    <- FC(gtab[ ,"sigma"], w=6, d=3)

  #  Rename col names for printing
  gtab.names       <- colnames(gtab)
  gtab.print.names <- gtab.names

  gtab.print.names[which(gtab.names=="SexPrt")]      <- "Sex" 
  gtab.print.names[which(gtab.names=="AgePrt")]      <- "Age" 
  gtab.print.names[which(gtab.names=="subset.type")] <- "sst" 

  gtab.print.names[which(gtab.names=="errcode")]    <- "ok?"
  gtab.print.names[which(gtab.names=="pct.lt.DL")]  <- "pct.DL"
  gtab.print.names[which(gtab.names=="method")]     <- "meth"
  gtab.print.names[which(gtab.names=="prev.l")]     <- "prevL"
  gtab.print.names[which(gtab.names=="prev.c")]     <- "prevC"
  gtab.print.names[which(gtab.names=="prev.r")]     <- "prevR"

  gtab.print.names[which(gtab.names=="x.RL1")]      <- "RL1"
  gtab.print.names[which(gtab.names=="x.RL2")]      <- "RL2"
  gtab.print.names[which(gtab.names=="x.tr.prop")]  <- "prop"

  gtab.print.names[which(gtab.names=="x.RL1.cilo")] <- "RL1.cilo"
  gtab.print.names[which(gtab.names=="x.RL1.cihi")] <- "RL1.cihi"
  gtab.print.names[which(gtab.names=="x.RL2.cilo")] <- "RL2.cilo"
  gtab.print.names[which(gtab.names=="x.RL2.cihi")] <- "RL2.cihi"

  gtab.print.names[which(gtab.names=="chi2.total")]    <- "c2.to"
  gtab.print.names[which(gtab.names=="chi2.total.df")] <- "to.df"
  gtab.print.names[which(gtab.names=="chi2.trun")]     <- "c2.tr"
  gtab.print.names[which(gtab.names=="chi2.trun.df")]  <- "tr.df"
  gtab.print.names[which(gtab.names=="chi2.path")]     <- "c2.pa"

  gtab.print.names[which(gtab.names=="x.tr.bins")]        <- "xtbn"

  names(gtab.print)    <- gtab.print.names
  rownames(gtab.print) <- NULL
  
  # ...........................................................................

  tab057.file <- paste(path.tab.file,outname,"-Tab13.txt", sep="")   
  sink(file=tab057.file, append=FALSE, split=TRUE)

  #  Part 1
  cat("\n--------------------------------------------------------------------------",
      "\n       Table 1: Main results",
      "\n--------------------------------------------------------------------------",
      "\n")
  print(gtab.print[ ,c("Sex","Age","n", "meth",
                       "RL1","RL2", 
                       "prevC",  "p.fit", "p.rt", "opt.crit", "ok?") ])
  cat("--------------------------------------------------------------------------",
      "\n")
 
    #  Part 2a
    cat("\n--------------------------------------------------------------------------",
      "\n       Table 2a: Asymptotic CIs for RILs",
      "\n--------------------------------------------------------------------------",
      "\n")
    print(gtab.print[ ,c("Sex","Age", "n", "meth", "pct.DL",
                 "RL1.cilo",  "RL1.cihi", "RL2.cilo",  "RL2.cihi")])
    cat("-------------------------------------------------------------------------",
      "\n")

#    if (nbs > 0)
#    {
#    cat("\n--------------------------------------------------------------------------",
#      "\n       Table 2b: Bootstrap CIs for RILs (nbs =", nbs, ")",
#      "\n--------------------------------------------------------------------------",
#      "\n")
#    print(gtab.print[ ,c("Sex","Alter", "n", "meth",
#                 "RL1.bs.mea", "RL1.tilo",  "RL1.tihi", 
#                 "RL2.bs.mea", "RL2.tilo",  "RL2.tihi")])
#    cat("--------------------------------------------------------------------------",
#      "\n")
#    }

    if ("tmc" %in% gtab.list)
    {
    cat("\n--------------------------------------------------------------------------",
      "\n       Table 2c: Permissible difference for TMC estimates",
      "\n                 Equivalence factor: ", ef,
      "\n--------------------------------------------------------------------------",
      "\n")
    print(gtab.print[ ,c("Sex","Age","Age.mea", "n", "meth",
                 "RL1.pdlo", "RL1.pdhi", "RL2.pdlo", "RL2.pdhi")])
    cat("--------------------------------------------------------------------------",
      "\n")

    #  Print more details
    #  Part 3
    cat("\n--------------------------------------------------------------------------",
      "\n       Table 3: Estimated PND parameters",
      "\n--------------------------------------------------------------------------",
      "\n")
    print(gtab.print[ ,c("Sex","Age","n",  "prop", "meth",
                 "rc", "iter", "lambda","mue","sigma") ])
    cat("--------------------------------------------------------------------------",
      "\n")
    sink()

    }

  if (print.6tables)
  {
    sink(file=tab057.file, append=TRUE, split=TRUE)

    #  Part 4
    #cat("\n--------------------------------------------------------------------------",
    #  "\n       Table 4: Properties of the estimated PND distribution",
    #  "\n--------------------------------------------------------------------------",
    #  "\n")
    #print(gtab.print[ ,c("Sex","Age", "n", "meth",
    #             "xc.mode","xc.mean","sigma","xc.P50") ])
    #cat("----------------------------------------------------------------------------",
    #  "\n")

    #  Part 5
    cat("\n--------------------------------------------------------------------------",
      "\n       Table 5: Diagnostic quantities: fit statistics",
      "\n--------------------------------------------------------------------------",
      "\n")
    print(gtab.print[ ,c("Sex","Age","n",  "meth",
                 "c2.to", "to.df", "c2.tr", "tr.df", "c2.pa",
                 "p.fit", "opt.crit") ])
    cat("----------------------------------------------------------------------------",
      "\n")

    sink()
  }              #  print.6tables ....

  #  Kopie von gtab (!) als csv-Datei 
  write.table(gtab,file=outname.file.csv,row.names=FALSE,col.names=TRUE,
              sep=";",dec=".")

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg057_MasterTables  End\n") }
# ===========================================================================
