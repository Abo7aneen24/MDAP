#  TMC_seg050_Master.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Master segment. Organises the analysis:
#  - Interpret user input
#  - Set up output files
#  - Do data analysis per subgroup
#  - Produce summary output

#  Naming convention

#  Output files are either plot files (directory figs) or
#  tables/ text files (directory tabs).
#  Output files refer either 
#  - to a stratum (=combination of sex and age) - these go in subdirectory 
#    1_stra under figs or tabs, respectively
#  - to a file (=complete data file) - these go in subdirectory 
#    2_file under figs or tabs, respectively
#  - to a sequence of files (a scenario, typically from a simulation) - these go in 
#    subdirectory 3_sce under figs or tabs, respectively

#  File names have the form outname_suffix, where
#  outname <- paste(file name, column number of measurand in the data file, 
#                   subset categories (if given))
#  This gives a unique identification also for data files which contain data
#  on several measurands. 

#  The suffix for output referring to a stratum contains of a stratum 
#  description (e.g. F_18-29) for females, 18-29 years plus a further
#  description of the contents (e.g. F100.010 for plot 010 from segment 100).

#  The suffix for output referring to a file is a further
#  description of the contents (e.g. F095.010 for plot 010 from segment 095).

#  The suffix for output referring to a scenario describes some conditions
#  of the simulation (e.g. -Z-D2-w 0.0 summary, distance 2, w.fact=0).

#  Data file names are defined:
#  outname  (the name basis)  - in seg050
#  outname files for strata   - in seg110
#  outname files for files    - in seg080
#  outname files for scenario - in seg211

#  Plot file names are defined in the ###_OpenWindows segments

#  TODO       

#  CHANGE HISTORY
#  15.11.2022 Last revision (editorial changes)
#  08.11.2022 Last essential modification
# ===========================================================================


# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg50_master  Start\n") }
# ===========================================================================

if (time.check) { master.t1 <- Sys.time() }

# ----------------------------------------------------------------------------
#  Version identification

VersionNr   <- 13
RevisionDate <- "2022-11-15"

cat("\n*********************************************************************",
    "\n  TMC: Truncated Minimum Chi-Square Estimation of Reference limits",
    "\n             Version ", VersionNr, "  Revision ", RevisionDate,
    "\n*********************************************************************",
    "\n\n")

# ----------------------------------------------------------------------------
#  Process parameter settings
source("TMC_seg040_ProcessSettings.R")  # Process parameters that may have 
                                        # been changed by the user

#  ..........................................................................
#  Read details of data file from info file
info <- read.table(file=info.file,sep=";",dec=".",
                   header=TRUE,stringsAsFactors=FALSE)

info.names <- colnames(info)

# Belo the full set from development. 
# [1] "FileNo"          "FileName"        "Path"            "Source"         
# [5] "DateReceived"    "Label"           "Value"           "Rounding"       
# [9] "Sex"             "SexCodeM"        "SexCodeF"        "Age"            
#[13] "DateTime"        "OH"              "Device"          "DecChar"        
#[17] "Group"           "lambda.l.gen"    "xl.RL1.gen"      "xl.RL2.gen"     
#[21] "xl.prev.gen"     "lambda.c.gen"    "xc.RL1.gen"      "xc.RL2.gen"     
#[25] "xc.prev.gen"     "xc.mue.gen"      "xc.sig.gen"      "yc.mue.gen"     
#[29] "yc.sig.gen"      "lambda.r.gen"    "xr.RL1.gen"      "xr.RL2.gen"     
#[33] "xr.prev.gen"     "drift.amp"       "diurn.amp"       "ntotal"         
#[37] "x.unaffected.lo" "x.unaffected.hi" "xc.unaffected.n" "xc.unaffected.p"
#[41] "x.unaffected.p"  "FP1.change"      "FP2.change"      "FP.change" 

if (time.check) { master.t2 <- Sys.time() }

#  Check for presence of mandatory columns
go.on <- TRUE

if (!("FileNo" %in% info.names))
{ cat("+++ Column 'FileNo' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("FileName" %in% info.names))
{ cat("+++ Column 'FileName' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("Path" %in% info.names))
{ cat("+++ Column  'Path'  missing in ",info.file,"  - Programme stops +++\n")
    go.on <- FALSE
} 

if (!("Value" %in% info.names))
{ cat("+++ Column 'Value' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("Rounding" %in% info.names))
{ cat("+++ Column 'Rounding' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("Sex" %in% info.names))
{ cat("+++ Column 'Sex' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("SexCodeF" %in% info.names))
{ cat("+++ Column 'SexCodeF' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("SexCodeM" %in% info.names))
{ cat("+++ Column 'SexCodeM' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

if (!("Age" %in% info.names))
{ cat("+++ Column 'Sex' missing in ",info.file,"  - programme stops +++\n")
  go.on <- FALSE
} 

#  ..........................................................................
#  Does the requested file number exist in the data information file?
#  Index of the data file to analyse in the data information file 
idx <- which(info[ ,"FileNo"] == FileNo)

if ( length(idx) == 0 )
{ cat("\n\n ++++++  Requested file no FileNo =",FileNo,
        "\n         not in                   ",info.file,
        "\n         ==> Programme stopped \n\n")
  stop("++++++  Programme stopped because file information missing  ++++")
}

if ( length(idx) > 1  )
{ cat("\n\n ++++++  Requested file no FileNo   =",FileNo,
        "\n         appears more than once in  ",info.file,
        "\n         ==> Programme stopped \n\n")
  stop("++++++  Programme stopped because of ambiguous file information ++++")
}

#  ..........................................................................

#  Interpret user input 

if (!go.on) 
{ 
  stop("+++  Programme stopped because of incomplete file information +++")
}

#  Collect information for reading the data

datafile    <- info[idx,"FileName"]

#  Die Datei mit Dateninformationen enthält für echte Daten den vollen 
#  Dateinamen einschl. ".csv" am Ende. Suffix abtrennen, wenn vorhanden
datafile.nchar <- nchar(datafile)
datafile.suff  <- substr(datafile,datafile.nchar-3,datafile.nchar)
if (tolower(datafile.suff) == ".csv") 
{ 
  datafile <- substr(datafile,1,datafile.nchar-4)
}
datafile0   <- datafile

#  datafile and datafile0 are used below to construct the file name to read,
#  thereby accounting for the construction of names when running a 
#  simulation.

inpath      <- info[idx, "Path"]
xlabel      <- info[idx, "Label"]  # will be extended below by OH category,
                                   # if spalte.OH and use.oh are given
spalte.w    <- as.numeric(info[idx,"Value"])

#  If no label is given, use value column as label 
if (xlabel == "") { xlabel <- paste("Column", spalte.w) }

spalte.s    <- as.numeric(info[idx,"Sex"])
spalte.a    <- as.numeric(info[idx,"Age"])

dec         <- info[idx,"DecChar"]

#  round unit in info.file refers to the rounding in the original scale
#  Use this, if there was no user specification in the start file.
if (is.na(round.unit)) 
{
  round.unit  <- as.numeric(info[idx,"Rounding"])
}

#  For further use, multiply by scale.fact
round.unit  <- round.unit * scale.fact

#  Columns of non-mandatory variables

if ("OH" %in% info.names)
{ spalte.oh   <- as.numeric(info[idx,"OH"]) 

  #  Complete xlabel by the OH category, if present 
  #  (OH = outpatient/hospitalized)
  #  no - not for publication
  #if (is.numeric(spalte.oh) & !is.na(use.oh)) 
  #{ xlabel <- paste(xlabel, "-", use.oh, sep="") 
  #} 
} else
{ spalte.oh   <- NA }

if ("Device" %in% info.names)
{ spalte.dev   <- as.numeric(info[idx,"Device"]) } else
{ spalte.dev   <- NA }

if ("DateTime" %in% info.names)
{ spalte.dz   <- as.numeric(info[idx,"DateTime"]) } else
{ spalte.dz   <- NA }

#  'Group' contains group membership, typically (but not exclusive)
#  available for generated data

if ("Group" %in% info.names)
{ spalte.g   <- as.numeric(info[idx,"Group"]) } else
{ spalte.g   <- NA }

#  Parameters of data generation. It is assumed that if generated data is used,
#  all parameters below are given in the file information

if ("lambda.c.gen" %in% info.names)
{ 
  lambda.c.gen <- as.numeric(info[idx,"lambda.c.gen"])
  xc.RL1.gen  <- as.numeric(info[idx,"xc.RL1.gen"])
  xc.RL2.gen  <- as.numeric(info[idx,"xc.RL2.gen"])
  xc.mode.gen <- as.numeric(info[idx,"xc.mode.gen"])
  yc.mode.gen <- as.numeric(info[idx,"yc.mue.gen"])
  yc.sig.gen  <- as.numeric(info[idx,"yc.sig.gen"])

  lambda.r.gen <- as.numeric(info[idx,"lambda.r.gen"])
  xr.RL1.gen  <- as.numeric(info[idx,"xr.RL1.gen"])
  xr.RL2.gen  <- as.numeric(info[idx,"xr.RL2.gen"])
  xr.prev.gen <- as.numeric(info[idx,"xr.prev.gen"])

  lambda.l.gen <- as.numeric(info[idx,"lambda.l.gen"])
  xl.RL1.gen  <- as.numeric(info[idx,"xl.RL1.gen"])
  xl.RL2.gen  <- as.numeric(info[idx,"xl.RL2.gen"])
  xl.prev.gen <- as.numeric(info[idx,"xl.prev.gen"])

  x.unaffected.lo <- as.numeric(info[idx,"x.unaffected.lo"])
  x.unaffected.hi <- as.numeric(info[idx,"x.unaffected.hi"])

  #  Safety action, if DataFileInfo has not all variables associated with 
  #  data generation

  if (is.null(xc.RL1.gen))      { xc.RL1.gen <- NA }  
  if (is.null(xc.RL2.gen))      { xc.RL2.gen <- NA } 
  if (is.null(xc.mode.gen))     { xc.mode.gen <- NA }
  if (is.null(yc.mode.gen))     { yc.mode.gen <- NA }
  if (is.null(yc.sig.gen))      { yc.sig.gen <- NA } 
  if (is.null(xr.RL1.gen))      { xr.RL1.gen <- NA } 
  if (is.null(xr.RL2.gen))      { xr.RL2.gen <- NA } 
  if (is.null(xr.prev.gen))     { xr.prev.gen <- NA }
  if (is.null(xl.RL1.gen))      { xl.RL1.gen <- NA } 
  if (is.null(xl.RL2.gen))      { xl.RL2.gen <- NA } 
  if (is.null(xl.prev.gen))     {xl.prev.gen  <- NA } 
  if (is.null(x.unaffected.lo)) { x.unaffected.lo <- NA } 
  if (is.null(x.unaffected.hi)) { x.unaffected.hi <- NA } 

  # Action if no right or no left pathological data
  if (is.null(xl.prev.gen)) { xl.prev.gen <- 0 }
  if (is.null(xr.prev.gen)) { xr.prev.gen <- 0 }

} else
{ 
  lambda.c.gen  <- NA
  xc.RL1.gen  <- NA
  xc.RL2.gen  <- NA
  xc.mode.gen <- NA
  yc.mode.gen <- NA
  yc.sig.gen  <- NA

  lambda.r.gen  <- NA
  xr.RL1.gen  <- NA
  xr.RL2.gen  <- NA
  xr.prev.gen <- NA

  lambda.l.gen  <- NA
  xl.RL1.gen  <- NA
  xl.RL2.gen  <- NA
  xl.prev.gen <- NA
  x.unaffected.lo <- NA
  x.unaffected.hi <- NA
}

theta.c.gen <- c(lambda.c=lambda.c.gen, 
                 yc.mode.gen=yc.mode.gen, 
                 yc.sig.gen=yc.sig.gen)

# ...........................................................................
# If the data is not generated data, then there is no information on
# the *.gen parameters. Remove "npa" from the list of methods to execute
# or display.

if (is.na(spalte.g))
{ #  *.gen data is not present
  #  "npa" in gtab.list?
  npa.idx <- which(gtab.list=="npa")
  if (length(npa.idx) > 0)
  { # Is present - remove it
    gtab.list <- gtab.list[-npa.idx]
    gtab.list.n  <- length(gtab.list)
    #cat("\n [seg050_Master]  'npa' removed from gtab.list",
    #    "\n                  Remaining methods to display in gtab:",
    #    "\n")
    #print(gtab.list)
  }

  #  "npa" in plot.list?
  npa.idx <- which(plot.list=="npa")
  if (length(npa.idx) > 0)
  { # Is present - remove it
    plot.list <- plot.list[-npa.idx]
    plot.list.n  <- length(plot.list)
    #cat("\n [seg050_Master]  'npa' removed from plot.list",
    #    "\n                  Remaining methods to display in gtab:",
    #    "\n")
    #print(gtab.list)
  }
} 

#  List of methods to execute. If tmc is present, qqw must also be present,
#  even if qqw is not plotted and does not appear in the final table.

meth.list   <- union(gtab.list, plot.list)
if ( ("tmc" %in% meth.list) & !("qqw" %in% meth.list) )
{ meth.list <- c(meth.list, "qqw") }
meth.list.n <- length(meth.list)


if (gtab.list.n == 0 | meth.list.n == 0)
{
  cat("\n ++++++ No methods to execute:\n")
  print(gtab.list)
  print(plot.list)
  print(meth.list)

  go.on = FALSE

  stop("+++  Stop because of missing entries in gtab.list / meth.list +++ ")
}

# ---------------------------------------------------------------------------

#  Separator lines for output tables
sep.p <- str_dup(".", 78)
sep.m <- str_dup("-", 78)
sep.e <- str_dup("=", 78)

# ---------------------------------------------------------------------------

if (eval.rep)
{ #  Evaluate test data (> 1 replicates): no separation by sex
  use.sex <- NA
} else
{ #  Evaluate real data (1 replicate): separate by sex
  #  Eigenmächtigkeit von Excel beheben: Eintrag "F" wird zu FALSE
  if (!is.na(info[idx,"SexCodeF"]) &&
      info[idx,"SexCodeF"] == "FALSE") {info[idx,"SexCodeF"] <- "F" }

  use.sex <- c(info[idx,"SexCodeF"], info[idx,"SexCodeM"])
}

# =================================================================
#  Widersprüche in der Eingabe beheben
#  Wenn x.lo.limit und break.min.user angegeben: 
#  Es muss gelten     break.min.user <= x.lo.limit 
if ( (!is.na(x.lo.limit)) & (!is.na(break.min.user)) &  
     (break.min.user > x.lo.limit) )
{
  break.min.user <- x.lo.limit
  cat("\n+++  Interval start for histograms corrected to ", break.min.user,
      "\n")
}
#  Weiterer Widerspruch: wenn x.lo.limit  nicht angegeben,  
#  break.min.user angegeben und break.min.user < x.min: Korrektur 
#  erst nach dem Lesen der Daten möglich, siehe dort

#  Für die Dimensionierung interner Tabellen

#  sex  ................................................................
#  NA in use.sex sind sinnlos, rauswerfen
use.sex[use.sex == ""] <- NA
use.sex      <- use.sex[!is.na(use.sex)]
use.sex.n    <- length(use.sex)

#  Wenn eval.rep, dann nicht nach sex trennen
if (eval.rep) { use.sex.n <- 0 }

if (use.sex.n > 0)
{ # Mindestens ein echter Code - nehmen
  sex.val   <- use.sex
  sex.val.n <- length(sex.val)
} else
{ # Kein Code angegeben oder Auswertung von Replikationen -  alles nehmen 
  #  @@@ This piece of code is overruled in seg070_ReadData
  sex.val   <- "All"
  sex.val.n <- 1
} 

#  age  ................................................................
if (!is.na(spalte.a))
{ #  Data contains age. 
  #  Age classes may be defined by age limits (cutpoints, non-overlapping
  #  classes) or as age class matrix (limits are arbitrary). 
  #  Use age.class, if given, even if age.limits exists as well.
  #  If age.limits is given, transform age.limits to age.class
  #  If neither age.limits nor age.class givven (both set to NA in seg030),
  #  create 1 age.class with min age, max age.

  if (is.na(age.class[1]) & is.na(age.limits[1]) )
  { #  No limits / classes given
    #  min age, max age are available only after data reading. 
    #  Set surrogate values (large enough to deal with age in days)
    age.class  <- AgeLimitsToClass(c(0, 36500))
  }   

  if (is.na(age.class[1]))
  {
    #  age.class is not given, create age.class from age.limits
    #  ...............................................................
    #  age.limits may contain NA in first or last position
    #  Replace by 0 and 99999. Latter will be set to actual age maximum  
    #  after data reading.
    age.limits.n <- length(age.limits)
 
    if (is.na(age.limits[1])) 
    {   
      age.limits[1] <- 0
    }
    if (is.na(age.limits[age.limits.n])) 
    { 
      age.limits[age.limits.n] <- 99999
    }

    #  Create age classes 
    age.class  <- AgeLimitsToClass(age.limits)
    if (length(age.class) == 1)
    { # No valid age limits
      cat("\n\n +++  No valid age limits, 'age' ignored +++ \n\n")
    } else
    {
      age.class.n  <- nrow(age.class)
    }
  } else      #  age.class[1] is NA  @@ else check of correct age.class

  { # age.class exists
    age.class.n  <- nrow(age.class)
  }

  #  Variable age is given, age.class exists: generate age class labels
  if (age.class.n > 0)
  { 
    agelabel.vec <- rep(NA, times=age.class.n)

    for (iage.class in 1:age.class.n)
    { #  Age classes
      agelabel.vec[iage.class] <- paste(age.class[iage.class, "lo"],
                                        age.class[iage.class, "hi"],sep="-")
    }
  } else
  { # keine Altersklassen
    agelabel <- "All"
  } 
#
}  else
{ #  Alter ist nicht angegeben - age.limits intern setzen
  age.limits   <- c(NA,NA)
  age.limits.n <- length(age.limits)
  age.class    <- RowToMatrix(c(lo=NA, hi=NA))
  age.class.n  <- 1
}

#  Wenn eval.rep, dann keine Altersunterscheidung
if (eval.rep)
{ #  age.limits intern setzen
  age.limits   <- c(NA,NA)
  age.limits.n <- length(age.limits)
  age.class    <- RowToMatrix(c(lo=NA, hi=NA))
  age.class.n  <- 1
}

#  outpatient / hospitalized ..............................................
#  NULL in use.oh ist sinnlos, rauswerfen

use.oh[use.oh == ""] <- NA
use.oh      <- use.oh[!is.na(use.oh)]
use.oh.n    <- length(use.oh)
if (use.oh.n == 0) 
{ 
  use.oh   <- NA
  oh.val.n <- 1 
}

if (use.oh.n > 0)
{ # Mindestens ein echter Code - nehmen
  oh.val   <- use.oh
  oh.val.n <- length(oh.val)
} else
{ # Kein Code angegeben -  alles nehmen 
  oh.val   <- "All"
  oh.val.n <- 1
} 

#  device ..................................................................
#  NULL in use.dev ist sinnlos, rauswerfen
use.dev[use.dev == ""] <- NA
use.dev      <- use.dev[!is.na(use.dev)]
use.dev.n    <- length(use.dev)

if (use.dev.n == 0) 
{ 
  use.dev   <- NA
  oh.dev.n  <- 1 
}

if (use.dev.n > 0)
{ # Mindestens ein echter Code - nehmen
  dev.val   <- use.dev
  dev.val.n <- length(dev.val)
} else
{ # Kein Code angegeben -  alles nehmen 
  dev.val   <- "All"
  dev.val.n <- 1
} 

#  ---------------------------------------------------------------
#  Schleife über Replikationen vorbereiten

if (!eval.rep) 
{ 
  #  Keine Replikationen auswerten 
  r.start <- 1
  r.ende  <- 1
  nrep    <- 1 
} 

#  -----------------------------------------------------------------
#  Create output tables

#  tab.s Contains all results (for all methods) for one strata (age/sex 
#        combinations) of a dataset.
#        Is overwritten by the analysis of the next stratum, if more than
#        on replicate (nrep > 1) is analysed.
#        Is created in seg100_Analysis.
#  

#  Ausgabetabelle für Ergebnisse der Auswertung von Replikationen 
#  anlegen
#  stab enthält komplette Ergebnisse für die Auswertung aller Replikationen
#  Nur, wenn dieses ein Lauf über Replikationen von simulierten Daten ist
#  Dann: keine Auswertung nach Sex und Alter getrennt, d.h. nur 1 Ergebnis

#  "Sex" is used for printing, "Age" as well
#  changed to SexPrt, AgePrt
#  "Sex.num" removed 22.03.2022
#  "Age.num" removed 22.03.2022

tab.names <- c("irep", "Sex", "sexlabel", "Age","agelabel",
               "Age.mea", "method",              
               "subset.type","n", "pct.lt.DL", "x.tr.n", "prop",
               "n.per.bin.min",
               "n.per.bin.min.eff", "n.per.bin.max.eff", "bins.n",  
               "x.tr.bin.n", 
               "x.tr.lo", "x.tr.hi", "x.unaffected.lo", "x.unaffected.hi",
               "xl.max", "xr.min",           
               "l.fact", "w.fact", "p.fact",
               "rc", "iter", "errcode",
               "lambda","mue","sigma",
               "RL1","RL2","RL1.CI","RL2.CI",
               "xc.mode","xc.mean","xc.sigma","xc.P50",
               "prev.l", "prev.c", "prev.r", "n.l", "n.c", "n.r",
               "FPR1", "FPR2", "FPR.sum",
               "opt.crit", "p.fit", "p.rt",
               "chi2.total" ,"chi2.total.df",
               "chi2.trun", "chi2.trun.df",  "chi2.trun.p", "chi2.path", 
               "prev.l.tmc.pen", "prev.r.tmc.pen", "delta.RL1", "delta.RL2",
               "RL1.cilo","RL1.cihi","RL2.cilo","RL2.cihi",
               "RL1.pdlo","RL1.pdhi","RL2.pdlo","RL2.pdhi",
               "RL1.bs.mea", "RL2.bs.mea",
               "RL1.tilo","RL1.tihi","RL2.tilo","RL2.tihi")

if (nrep > 1)
{ 
  #  Namensmenge in stab muss Teilmenge von tab.names sein (??)
  #  Bei Auswertung von Replikationen wird nur der gesamte Datensatz
  #   ausgewertet, keine Trennung nach Kategorien von 
  #  Sex oder Alter
  #  Change 17.04.2021: 
  #  'stab' contains only information that is replicate-specific. 
  #  Settings per file (and automatically also per scenario) are stored in
  #  settings / settings.csv
  
  stab.names <- c("irep", "n", "pct.lt.DL", "method", 
                  "n.per.bin.min", "n.per.bin.min.eff",
                  "rc", "iter", "errcode",
                  "prop", "x.tr.lo", "x.tr.hi",
                  "lambda","mue","sigma","p.fit","RL1","RL2",
                  "prev.l", "prev.c", "prev.r",
                  "FPR1", "FPR2", "FPR.sum")

  stab       <- matrix(NA,nrow=meth.list.n*nrep,
                       ncol=length(stab.names))
  colnames(stab)  <- stab.names
  stab            <- data.frame(stab)
  istab           <- 1              
}

#  ==========================================================================
#  Input is checked for formal errors

if (time.check) {  master.t3 <- Sys.time() }

#  Process plot requests (plot level is defined in 060)
source("TMC_seg045_PlotRequest.R")

#  ==========================================================================
#  Show programme versions and all variables defined so far

source("TMC_seg042_SettingsSummary.R")  

#  ==========================================================================
#  Loop over replicates 

oh.sym <- "All"
if (!is.na(use.oh)) { oh.sym <- use.oh }

dev.sym <- "All"
if (!is.na(use.dev)) { dev.sym <- use.dev }

for (irep in r.start:r.ende)
{
  #  Construct input data file name.
  #  For real data (not simulated) only add the sufffix, otherwise
  #  add replicate number and suffix

  if (!eval.rep)
  { #  Names for single-file analysis   
    infile   <- paste(inpath,datafile0,".csv",sep="")
    
    #  Construct outname, the basic component for all output files
    outname <- paste(datafile0,
                     "_c",spalte.w,
                     "_o", oh.sym, 
                     "_d", dev.sym,    sep="")

  } else
  { #  Names for file sequence analysis
    #  This is the analysis of (simulated) test data. Add 3-digit replicate
    #  number to the base name.
    #  irep == 0: reserved for expected values of the data, no random error
    datafile <- paste(datafile0,formatC(irep,flag="0",width=3),sep="") 
    infile   <- paste(inpath,datafile,".csv",sep="")

    #  Construct outname, the basic component for all output files
    outname <- paste(datafile,
                     "_c",spalte.w,
                     "_o", oh.sym, 
                     "_d", dev.sym,    sep="")
    
    #  Construct more names for output files
    source("TMC_seg211_NamesPerScena.R")
  }

  #  ...............................................................
  #  If first replicate: 
  #  - define file names
  #  - remove all output files beginning with outname
  #  - organize graph windows

  if (irep == 1)
  { 

    #  Clean up
    source("TMC_seg060_RemoveOldFiles.R")
  }

if (time.check) { master.t4 <- Sys.time() }

  #  -----------------------------------------------------------------
  #  Read input data

  cat("\n\n",  
      "\n############################################################################",
      "\n#  Start of analysis for data file                                         #",
      "\n   ",infile,
      "\n############################################################################",
      "\n\n")
  
  source("TMC_seg070_ReadData.R")

if (time.check) { master.t5 <- Sys.time() }

  #  Now available:
  #  Input file name 
  #  Input data (in dataset)
  #  Various control parameters

  #  -----------------------------------------------------------------
  #  If dynamic age classes are requested: determine age classes per 
  #  sex group, minimum size = x.n.min, set in 030
  
  # names(dataset)
  # [1] "PatId"    "Value"    "Age"      "Sex"      "DateTime" "OH"      
  # [7] "Device"   "Grp"      "ValueRaw" "for.npa"  "ana.day"  "ana.mon" 
  #[13] "ana.yea"  "ana.hour" "ana.minu" "ana.wday" "DM"  

  #dataset[1, ]
  #    PatId Value Age Sex            DateTime OH Device Grp ValueRaw for.npa
  #687   687   133  53   F 22.01.2017 15:31:00  s    D02   C 132.9755   FALSE
  #    ana.day ana.mon ana.yea ana.hour ana.minu ana.wday         DM
  #687      22       1    2017       15       31        6 2017-01-15

  
  #  -----------------------------------------------------------------
  #  Set up plot windows for plots on file level

  source("TMC_seg046_OpenWindows.R")

 #  ...............................................................
  #  For plotting: if x.clip.min or x.clip.max are NA, replace by 
  #  x.min or x.max, resp.

  if (is.na(x.clip.min)) x.clip.min <- x.min
  if (is.na(x.clip.max)) x.clip.max <- x.max

  # ==========================================================================
  #  Start analysis 
  # ==========================================================================

  #  -----------------------------------------------------------------
  #  Process acceptable proportions in truncation intervals

if (time.check) { master.t6 <- Sys.time() }

  prop.lim <- ProcessPropLimits(x.tr.prop.min, x.tr.prop.max,
                                x.tr.prop.min.DL, x.tr.prop.max.DL,
                                x.tr.prop.min.noDL, x.tr.prop.max.noDL,
                                detect.limits.max, 
                                logfile=logfile, print.log.message=FALSE)
  x.tr.prop.min <- prop.lim[1]
  x.tr.prop.max <- prop.lim[2]

  #  -----------------------------------------------------------------
  #  From here on: use only age.class, not age.limits

  #  -----------------------------------------------------------------
  #  Auswertung 1: ganzer Datensatz (mit Ausnahme der Bereiche, die global
  #  weggefiltert wurden)

  iblock          <- 0       # Index for row blocks in gtab

  verlaufsplot    <- FALSE   # Ist folgende Auswertung Teil einer Auswertung 
                             # über Kategorien von Sex oder Alter?

  #  x, age, sex, grp für weitere Verarbeitung aus dem data frame ziehen

  x   <- dataset[ ,spalte.w]
  if (!is.na(spalte.a)) { age <- dataset[ ,spalte.a] }
  if (!is.na(spalte.s)) { sex <- toupper(dataset[ ,spalte.s]) }
  if (!is.na(spalte.g)) { grp <- toupper(dataset[ ,spalte.g]) }

  sexlabel <- sex.val[1]
  isex     <- 1
  while(isex < use.sex.n)
  { isex <- isex + 1
    sexlabel   <- paste(sexlabel,use.sex[isex],sep="+")
  }

  age.min      <- age.class[1, "lo"]
  age.max      <- age.class[age.class.n, "hi"]
  agelabel     <- paste(age.min,"-",age.max,sep="")
  iage         <- NA

  #  Darstellung der vollständigen Verteilung (alle Alter), getrennt nach Sex,
  #  nicht bei Simulationen
  if (!eval.rep)
  {
    sink(file=outname.file.txt,split=TRUE, append=TRUE)
    if (plot.fig090.010 & plot.fig090.020)
    { 
      source("TMC_seg090_Plot_Age_Val.R")
    }
    sink()

    #  Darstellung der Werteverteilung über Wochentage und Tageszeit
    sink(file=outname.file.txt,split=TRUE, append=TRUE)
    if (plot.fig095.010 & plot.fig095.020)
    { 
      source("TMC_seg095_Plot_DateTime_Val.R") 
    }
    sink()
  }

  # ------------------------------------------------------------------------
  #  Time check
  if (time.check) 
  {  master.t7 <- Sys.time()

    sink(file="../temp/UsedTime.txt", append=TRUE)
    cat("\n master; t2-t1;", 
      format(difftime(master.t2, master.t1)))
    cat("\n master; t3-t2;", 
      format(difftime(master.t3, master.t2)))
    cat("\n master; t4-t3;", 
      format(difftime(master.t4, master.t3)))
    cat("\n master; t5-t4;", 
      format(difftime(master.t5, master.t4)))
    cat("\n master; t6-t5;", 
      format(difftime(master.t6, master.t5)))
    cat("\n master; t7-t6;", 
      format(difftime(master.t7, master.t6)))  
    cat("\n master; t7-t1;", 
      format(difftime(master.t7, master.t1)))
    sink()
  }

  # ------------------------------------------------------------------------
  #  Prepare the pattern of subgroup analyses

  A <- matrix(FALSE, nrow=sex.val.n+1, ncol=age.class.n+1)
  rnamesA <- sex.val
  if (sex.val.n > 1) 
  { 
    rnamesA <- c(rnamesA, paste(sex.val[1], sex.val[2], sep="+"))
  }

  cnamesA <- age.class
  cnamesA <- c(paste(age.class[ , "lo"], age.class[ , "hi"], sep="-"))
  cnamesA <- c(cnamesA, paste(age.class[1, "lo"], 
                              age.class[age.class.n, "hi"], sep="-"))

  dimnames(A) <- list(rnamesA, cnamesA)

  # ------------------------------------------------------------------------
  #  Analyse complete data set, all sexes, all ages jointly
  subset.type <- 1
 
  if ( ana.type1 ) 
  { 
    source("TMC_seg100_Analysis.R") 

    #  Bookkeeping
    A[sex.val.n+1, age.class.n+1] <- TRUE 
    if (age.class.n == 1) { A[sex.val.n+1, 1] <- TRUE }
    if (sex.val.n == 1)   { A[1, 1:age.class.n] <- TRUE }
  }

  # ------------------------------------------------------------------------
  #  Teildatensätze via Sex und Altersklasse definieren,
  #  per Indikatorvektor die betreffenden Sätze aus dataset extrahieren,
  #  dabei bisheriges x und y überschreiben
  #  Gesamtdaten bleiben im data frame "dataset" 

  # ------------------------------------------------------------------------
  #  If sex given in data, analyse per sex group, regardless of age
  subset.type <- 2

  if (!is.na(spalte.s) && ana.type2) 
  {
    agelabel  <- paste(age.class[1,"lo"],"-",age.class[age.class.n, "hi"],
                       sep="")
      
    for (isex in 1:sex.val.n)
    { 
      #  Not, if this job is already done
      if (!A[isex, age.class.n+1])
      {
        sexlabel <- sex.val[isex]
        subset <- dataset[ ,spalte.s] == sex.val[isex]
        x   <- dataset[subset,spalte.w]
        x.n <- length(x)

        if (!is.na(spalte.a)) { age <- dataset[subset,spalte.a] }
        if (!is.na(spalte.s)) { sex <- dataset[subset,spalte.s] }
        if (!is.na(spalte.g)) { grp <- dataset[subset,spalte.g] }

        # Auswertung von x, Sex = Ausprägung isex
        source("TMC_seg100_Analysis.R")

        #  Bookkeeping
        A[isex, age.class.n+1] <- TRUE 
        if (age.class.n == 1) { A[isex, 1] <- TRUE }
        if (sex.val.n == 1)   { A[1, 1:age.class.n] <- TRUE }

      } #  if (!A[isex, age.class.n+1] ...
    }   #  for (sex ...
  }     #  if ... ana.type2?

  # -------------------------------------------------------------------------
  #  If age given in data and age classes are defined, analysis by age class,
  #  jointly for all sexes  
  #  Different organisation depending on age.class.style ("fix" or "dyn")
  subset.type <- 3
  
  if (!is.na(spalte.a) && ana.type3 ) 
  {
    if (age.class.style == "fix")
    { #  age.limits or age.class must have been given

      source("TMC_seg051_AgeClassFix3.R")
    }

    #  next loop not yet prepared for ana.type@@
    if (age.class.style == "dyn")
    { #  age.class is determined from data, depending on x.n.min

      source("TMC_seg052_AgeClassDyn3.R")
    }

  }   #  !is.na(spalte.a) &&

  # -------------------------------------------------------------------------
  #  If sex and age available, analysis per sex and age class
  subset.type    <- 4

  if (!is.na(spalte.a) && !is.na(spalte.s) && ana.type4 ) 
  {
    if (age.class.style == "fix")
    { #  age.limits or age.class are (must be) given 

      if (print.log.message) { cat("%%%   TMC_seg50_master  150\n") }
      source("TMC_seg051_AgeClassFix4.R")
      if (print.log.message) { cat("%%%   TMC_seg50_master  160\n") }
    }

    if (age.class.style == "dyn")
    { #  age.class is determined from data, depending on x.n.min
  
      source("TMC_seg052_AgeClassDyn4.R")
    }
  }     #   if (ana.type4 ...

  # ---------------------------------------------------------------------------
  #  Output for a complete file, all Sex*Age combinations

  # ...........................................................................
  #  Overview

  sink(file=outname.file.txt,split=TRUE, append=FALSE)

  cat("\n\n===========================================================================",
      "\n\nEstimation of Reference Limits ",
      "\nSUMMARY OF RESULTS FOR ALL DATA SUBSETS (STRATA)",
      "\nIdentification of this run             ", RunId,
      "\nComplete list of parameters in         ", logfile,
      "\n",
      "\nUser                                   ", user,
      "\nData file information in               ", info.file,
      "\nStart of this run                      ", format(time.start),
      "\nData file                              ",
      "\n   ",infile,
      "\nAnalyte                                ",xlabel, 
      "\n",
      "\nResults in file                        ",
      "\n   ",outname.stra.txt,
      "\nFilter  outpatient/hospitalized        ", use.oh,
      "\nFilter  device                         ", use.dev,
      "\nFilter  analysis time of day           ", paste(ana.hour.min, ":",
                                                         ana.minu.min, " - ",
                                                         ana.hour.max, ":",
                                                         ana.minu.max, sep=""),
      "\nFilter  date between (limits included) ", start.date, end.date, 
      "\nInput values multiplied by             ", scale.fact,
      "\nScaled values are limited by           ", x.lo.limit,x.hi.limit,
      "\nScaled values are rounded here to      ",round.unit,
      "\nSurrogate factor for values < DL       ", s.fact, 
      "\nTotal # of values in data file         ",x.n1,
      "\n# of values after filtering, sampling  ",x.n,
      "\nAge class definition style             ",age.class.style,
      "\nFor dynamic age classes only: min. size",x.n.min,
      "\n",
      "\nSmooth empirical histogram?            ",smooth.hist1, smooth.hist2,
      "\nMin. # of values per bin               ",n.per.bin.min,
      "\nMin. # of bins per stratum             ",bins.n.min,
      "\nMax. # of central bins per stratum     ",bins.n.max,
      "\nMin. # of bins in truncation interval  ",x.tr.bins.min,     
      "\nMin. # of values in a stratum          ",x.n.min,
      "\nUser-specified min. prop. in TI        ", x.tr.prop.min, 
      "\nUser-specified max. prop. in TI        ", x.tr.prop.max,     
      "\nWeight of truncation interval length   ",l.fact,
      "\nPenalty factor 'wrong prevalence'      ",p.fact,
      "\nPenalty factor 'wrong prediction'      ",w.fact,
      "\nMinimal/ maximal lambda                ",lambda.min, lambda.max,
      "\nMinimal p for goodness of fit in TI    ",p.fit.min,
      "\nMinimal p residuals runs test in TI    ",p.rt.min,
      "\n",
      "\nProbabilities < RL1 / < RL2            ",RL1.p,"/ ",RL2.p,   
      "\n===========================================================================",
      "\n\n")

  # ...........................................................................
  #  Before printing detailed result tables: repeat warnings

  # if (print.log.message) { cat("%%%   TMC_seg50_master  220\n") }

  # ...........................................................................
  #  Consistent rounding?
  #  Calculations were made in seg070

  cat("\n", rep("=", times=76), "\n", sep="")
  cat("\n Check for inconsistent rounding in the input data\n")

  if (ldt.p < 0.05)
  {
    cat("\n", rep("+", times=76), "\n", sep="")
    cat("\n p value for H0: 'Last digits are equally likely': ", ldt.p,"\n")
    cat("\n +++ Consistent rounding for the whole data set is questionable +++\n")
    cat("\n    Frequencies of the last digit", 
           "('Observed' should be close to 'Expected') \n\n")
    print(ldt)
    cat("\n", rep("+", times=76), "\n", sep="")
  }  else      
  {
    cat("\n p value for H0: 'Last digits are equally likely': ", ldt.p,"\n")
    cat("\n *** Rounding seems consistent ***\n")
  }
  cat("\n", rep("=", times=76), "\n\n", sep="")

  # if (print.log.message) { cat("%%%   TMC_seg50_master  230\n") }

  # ...........................................................................
  #  Drift over time? 
  #  Only if requested and if date of measurement is provided. Then 'dataset'
  #  contains ana.mon, ana.yea
  if (plot.fig050.010 && 
      (("ana.mon" %in% colnames(dataset)) & ("ana.yea" %in% colnames(dataset))))
  { 
    dev.set(fig050.010)
    par(mfrow=c(2,1))

    #  ana.mon, ana.yea are present. ana.day not needed as we are looking for
    #  monthly medians

    dataset <- data.frame(dataset, 
                          DM=as.Date(paste("15.", dataset[ , "ana.mon"], ".", 
                                                  dataset[ , "ana.yea"],
                                           sep=""), format = "%d.%m.%Y") )
    dataset[1:10, ]
    as.numeric(dataset[1:10, "DM"])

    dataset <- dataset[order(dataset[ ,"DM"]), ]

    #  Check drift for F
    ok <- dataset[ , spalte.s] == "F"
    subt <- paste(RunId, " / " ,outname, "_F / n=", sum(ok), sep="")

    #  Take empirical quantiles of all data as RLs. No age stratification.
    RL1.F.emp <- Quantile(dataset[ok, spalte.w], probs=RL1.p)
    RL2.F.emp <- Quantile(dataset[ok, spalte.w], probs=RL2.p)

    drift.F <- Drift(dataset[ok, "DM"], dataset[ok, spalte.w], 
                     RL1.F.emp, RL2.F.emp, fig050.010, xlabel, 
                     "Drift for sex=F", ef,
                     subt=subt, gamcol=blue.32)

    ok <- dataset[ , spalte.s] == "M"
    subt <- paste(RunId, " / " ,outname, "_M / n=", sum(ok), sep="")

    #  Take empirical quantiles of all data as RLs. No age stratification.
    RL1.M.emp <- Quantile(dataset[ok, spalte.w], probs=RL1.p)
    RL2.M.emp <- Quantile(dataset[ok, spalte.w], probs=RL2.p)

    drift.M <- Drift(dataset[ ,"DM"], dataset[ , spalte.w], 
                     RL1.M.emp, RL2.M.emp, fig050.010, xlabel, 
                     "Drift for sex=M", ef,
                     subt=subt, gamcol=blue.32)

    cat("\n", rep("=", times=76), "\n", sep="")

    cat("\n Check for drift of values over time\n")
    cat("\n The blue confidence band should overlap",
            "the gray range at all times\n\n")
    cat("\n See window ", fig050.010,
        "\n and        ", fig050.010.file,
        "\n") 

    drifttext <- "No drift"
    if (drift.F) { drifttext <- "Drift detected" }
    cat("\n    Sex = F: ", drifttext) 

    drifttext <- "No drift"
    if (drift.M) { drifttext <- "Drift detected" }
    cat("\n    Sex = M: ", drifttext) 

    savePlot(file=fig050.010.file, type=figtype)

    cat("\n\n", rep("=", times=76), "\n\n", sep="")
  }
  
  # ===========================================================================
  #  If this analysis is not an analysis of several replicates (usually
  #  simulated datasets): plot estimated RLs by age, if the data allows
  #  (ana.type >) 1 
  #  Calculate splines and permissible differences
  #  Prepare entries in gtab
  
  if ( ana.type3 | ana.type4 )
  {
    source("TMC_seg056_MasterPlots.R")
  }

  # ===========================================================================
  #  If this analysis is not an analysis of several replicates (usually
  #  simulated datasets): produce general output table per file
  
  source("TMC_seg057_MasterTables.R")
  

  # ===========================================================================
  #  Falls dieses die Auswertung einer Replikation ist:
  #  Ergebnis dieser Analyse in stab einfüllen
  if (nrep > 1)
  {
    stab[istab:(istab+meth.list.n-1), stab.names] <- gtab[ ,stab.names]
    istab <- istab + meth.list.n
  }

  # ===========================================================================
  time.end   <- Sys.time()

  cat("\n\n",  
      "\n#####################################################################",
      "\n     End of analysis for data file                                   ", 
      "\n    ",infile,
      "\n",
      "\n    Execution time: ",format(difftime(time.end,time.start)),
      "\n#####################################################################",
      "\n\n")

  cat("\n*********************************************************************",
      "\n  TMC: Truncated Minimum Chi-Square Estimation of Reference limits",
      "\n             Version ", VersionNr, "  Revision ", RevisionDate,
      "\n*********************************************************************",
      "\n\n")

  #  End of results output to total table for this data file
  sink()

}  #  for (irep in 1:nrep) 

# ===========================================================================
#  If this analysis *is* an analysis of several replicates (usually
#  simulated datasets): plot estimated RLs over all replicates

if (nrep > 1)
{
  #  Organize filenames and plot windows 

  write.table(stab, file=soutfile.tmc.csv, 
              row.names=FALSE, col.names=TRUE, 
              quote=FALSE, sep=";", dec=".")

  source("TMC_seg202_SimResults.R")
}

# ===========================================================================
#  End of execution

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg50_master  End\n") }
# ===========================================================================
