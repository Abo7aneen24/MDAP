#  TMC_seg015_DefaultSettings.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment containing default settings. Can be modified in the start segment.

#  TODO
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg15_DefaultSettings  Start\n") }
# ==========================================================================

time.check <- FALSE          #  If TRUE, write time stamps to 
                             #  ../temp/UsedTime.txt

#  Record start time for total time
time.start <- Sys.time()     #  "2020-02-26 11:55:02 CET"
                             #   12345678901234567890123
                             #            11111111112222 

if (time.check) { DefaultSettings.t1 <- time.start }

# ===========================================================================
#  Packages

#  Set repository from wich packages are to be loaded
repos <- "https://ftp.fau.de/cran/"

# ...........................................................................
#  Try to load package date
date.loaded <- require(date)

# if load was not successful, then date.loaded == FALSE
if (!date.loaded)
{ 
  install.packages("date", dependencies=TRUE, repos=repos)
}
library(date)

# ...........................................................................
#  Try to load package mgcv  (for gam according to Simon Wood)
mgcv.loaded <- require(mgcv)

# if load was not successful, then mgcv.loaded == FALSE
if (!mgcv.loaded)
{ 
  install.packages("mgcv", dependencies=TRUE, repos=repos) 
}
library(mgcv)

# ...........................................................................
#  Try to load package RColorBrewer
RColorBrewer.loaded <- require(RColorBrewer)

# if load was not successful, then RColorBrewer.loaded == FALSE
if (!RColorBrewer.loaded)
{ 
  install.packages("RColorBrewer", dependencies=TRUE, repos=repos)
}
library(RColorBrewer)

# ...........................................................................
#  Try to load package stringr
stringr.loaded <- require(stringr)

# if load was successful, then stringr.loaded == FALSE
if (!stringr.loaded)
{ 
  install.packages("stringr", dependencies=TRUE, repos=repos) 
}
library(stringr)

if (time.check) { DefaultSettings.t2 <- Sys.time() }

# ---------------------------------------------------------------------------

#  Load styles and functions

source("TMC_seg022_Functions.R")        # Non-standard functions
source("TMC_seg030_Style.R")            # Define style (colors, lines, ...)

if (time.check) { DefaultSettings.t3 <- Sys.time() }

# ---------------------------------------------------------------------------

#  Selected functions (new or modified versions in ../func or /dev)

# source("../func/F_ShowTiProps.R")

if (time.check) { DefaultSettings.t4 <- Sys.time() }

# ---------------------------------------------------------------------------

#  Create identification of this analysis
RunId <- paste(substr(time.start,  1, 10), "_", 
               substr(time.start, 12, 13), 
               substr(time.start, 15, 16), 
               substr(time.start, 18, 19), sep="")

#  Logfile containing (nearly) all run parameters
logfile      <- paste("../log/", RunId,".txt", sep="")      #  txt version
logfile.rda  <- paste("../log/", RunId,".RData", sep="")    #  RData
                                                            #  only simulation

#  --------------------------------------------------------------------
#  General constants
# ----------------------------------------------------------------------       

fastnull     <- 1.e-10
fastnull.chi <- 0.1

# ---------------------------------------------------------------------------
#  Produce a new data set with first values? If so, provide
#  ... the column containing the patient identification
spalte.i <- NA

#  ... the column containing the date of analysis (or the patient age
#  at the time of analysis)
spalte.d <- NA

#   ... the name of the new file with first values 
DataFileFirstValue <- NA 

#  If spalte.i or spalte.d or DataFileFirstValue == NA:
#  no file with first values is produced

# -------------------------------------------------------------------------

#  Values for non-mandatory data columns 
spalte.oh   <- NA      #  data column outpatient/hospitalized
spalte.dev  <- NA      #  data column device
spalte.dz   <- NA      #  data column date-time
use.oh      <- NA      #  select outpatient/hospitalized
use.dev     <- NA      #  select device

# -------------------------------------------------------------------------
#  Default for user-defined rounding: NA <=> take from DataFileInfo
#  Specification in start file overrides this
round.unit <- NA

# -------------------------------------------------------------------------
#  Factor for constructing surrogate values for values < detection/
#  determination limit
s.fact <- 0.67

#  Option for dealing wit "< DL" in the data
DL.option <- 1          # replace values of the form "< DL" with 
                        # max(DL) * s.fact, return detect.limits.max
                        # More options available, but not recommended

# -------------------------------------------------------------------------
#  Draw a subsample of size ...
subsample.n <- NA

# -------------------------------------------------------------------------
#       Time range selection: use only data obtained before / after 
split.date <- NA                # no split 
split.time <- NA                # no split
# split.date <- "04.02.2019"    # "dd.mm.yyyy"
# split.time <- "13:25"         # "hh:mm"

# split.date.n <- as.numeric(mdy.date(02, 04, 2019)) #  split example

use.split      <- NA        # split to use

#       Date selection: use only data from a data range
                           #  Select a date range for analysis.
                           # Format is "dd.mm.yyyy"
                           # Limits belong to interval 
                           # Both values must be given
                           # Both values == NA means to use all
# -------------------------------------------------------------------------
#       Alternative time range selection: use only data between start.date,
#       end.date, limits included
#       format: "dd.mm.yyyy"
start.date <- NA 
end.date   <- NA

# -------------------------------------------------------------------------
#   (6) Daytime selection: use only data obtained between  ... 
#      (no selection: set all 4 values to NA)
ana.hour.min <- 0
ana.minu.min <- 0 

#  (7) and ...
ana.hour.max <- 24
ana.minu.max <- 00

# -------------------------------------------------------------------------
#  (7a) Weekday selection:
use.wday <- "all"       #  "all"   or a subset of 
                        #  c("mo", "tu", "we", "th", "fr", "sa", "su")
                        #  use.wday <- c("sa", "su")

# -------------------------------------------------------------------------
#  PND properties
#  If lambda.min = lambda.max, this value will be used as fixed value 
#  of lambda. Influences df.est.
lambda.min <-    0.0 # Minimal allowed lambda
lambda.max <-    0.0 # Maximal allowed lambda  21.12.2021

# -------------------------------------------------------------------------
#  Degrees of freedom for estimates of all 3 PND parameters

#  df to pay for estimation  = # of parameters
#  df.est    <- 2 + (abs(lambda.min - lambda.max) > 1.e-6) @@@ in process
#                                               settings überarbeiten
df.est <- 3

#  df to pay for constraints
df.con    <- 1

#  Df to pay must be subtracted from number of bins when assessing the chi^2 
#  value 

# -------------------------------------------------------------------------
#     Histogram breaks start and width, maximimum will always be calculated.
#     Not usually needed, histogram is by default automatically created,
#     following rules in the next section
#     NA in both variables: automatic calculation
break.min.user   <-  NA       # NA: automatic selection
break.delta.user <-  NA

# -------------------------------------------------------------------------

#  Standard minimum number of bins in histogram. Default value for data which
#  has enough bins. Will automatically be reduced if not enough bins available.
#  Former variable x.bins.n removed. Now only bins.n.min effective. 15.04.2021
bins.n.min <- df.est + df.con + 4        #  Standard  df.est + df.con + 4

#  Standard maximum number of bins in histogram. Is used to prevent 
#  excessively  many bins that require too much analysis time.
bins.n.max <- 50                # 29.09.2022 reduced because of new meaning,
                                # refers to inner part of the data oly 

#  Minimum number of values per bin
n.per.bin.min <- 50          #  28.11.2021
                             #  15.09.2022 new meaning: absolute minimum.
                             #  Effective value is calculated in seg 100
                             #  (n.per.bin.min.adj) depending on x.n and 
                             #  bins.n.max
#  Truncation intervals used for estimation must contain at least 
#  x.tr.bins.min   bins
x.tr.bins.min <- df.est + df.con + 2    #  23.04.2021

#  Truncation intervals used for TMC estimation must contain at least 
#  a proportion of ... of all data
#  These limits are now set in seg040 and depend on whether a DL exists.
#  Defaults values for both cases are given here
#  User input of x.tr.prop.min, x.tr.prop.max overrides all defaults
x.tr.prop.min      <- NA
x.tr.prop.min.noDL <- 0.55   #  29.09.2022
x.tr.prop.min.DL   <- 0.40   #  29.09.2022

#  Truncation intervals used for TMC estimation may contain at most 
#  a proportion of ... of all data
x.tr.prop.max      <- NA     #  29.09.2022
x.tr.prop.max.noDL <- 1.00   #  29.09.2022
x.tr.prop.max.DL   <- 1.00   #  29.09.2022

#  A truncation interval must contain a proportion of ... left of the mode
deltap.lo <- 0.15

#  An optimal truncation interval must contain a proportion of ... right of the mode
deltap.hi <- 0.15

#  Required minimum number of observations in a stratum
x.n.min <- n.per.bin.min*x.tr.bins.min

#  Use only TIs starting with leftmost bin? 
TI.left.only       <- FALSE

#  Use only TIs starting with leftmost bin if a decision limit (DL) exists 
#  in the data? 
TI.left.only.if.DL <- FALSE

#  Maximum allowed percentage of values below detection limit
x.lt.detect.pct.max <- 65         #  in % !!  

# -------------------------------------------------------------------------

#  Limit for using a reduced dataset in QQW: 
red.limit <- 1000     #  activated 30.09.2022

# -------------------------------------------------------------------------
#
age.limits   <- NA
age.limits.n <- NA

age.class   <- NA
age.class.n <- NA

age.class.style <- "fix"        #  "fix" | "dyn"

# -------------------------------------------------------------------------
#  Select probability levels for RLs

RL1.p <- 0.025
RL2.p <- 0.975

# -------------------------------------------------------------------------
#  Age/ sex stratification --- needs correction @@@@

ana.type1 <- TRUE            # whole dataset
ana.type2 <- TRUE            # by sex group, all ages
ana.type3 <- TRUE            # all sexes, by age group
ana.type4 <- TRUE            # by sex*age group  

# -------------------------------------------------------------------------
#  Plot layout
#
par.tcl <- 0.5            #  0.5: tickmarks  inside the graph area
                          # -0.5: tickmarks outside the graph area
par.las <- 1              # 0: all labels parallel to the axis
                          # 1: all labels horizontal
                          # For plots 150, 151, 160 
age.label    <- "Age"     # default if not changed by user 
RL.label     <- NA
RL.label.tmc <- NA        # automatically generated if not changed by user 
# -------------------------------------------------------------------------
#  When plotting histograms with x range limited by x.clip: use
#  the limits as specified or move them to the next value in x.breaks?

x.clip.type <- "as.specified"     #  "coincide"

# -------------------------------------------------------------------------
# Smooth the empirical histogram? Useful if data is rounded by "go to the 
# even number". At most 1 option may be selected.

smooth.hist1 <- TRUE       #  5 point moving average / 23.04.2021
smooth.hist2 <- FALSE      #  use kernel density estimate for smoothing
      
# -------------------------------------------------------------------------
#  Parameters for QQW
#  This will be changed in 040_ProcessSettings if the user restricts lambda 
#  such that this sequence makes no sense

lambda.seq <- c(-0.10, 0.00, 0.25, 0.50, 0.75, 1.00) # activated 30.09.2022

# -------------------------------------------------------------------------
#  Parameters for TMC

#  Compensation factor for length of truncation interval (worse chi2.per.df
#  is compensated if interval is larger)
#  lt = 0 means no compensation
#  
l.fact <- 0.0              #  Default in last published version: ---
                           #  LEAVE AT 0!!

#  Weighting factor for for penalty term "negative estimated prevalence"  
#  Controls degree of penalisation, p.fact <- 0 means no penalisation 
p.fact <- 0.0            

#  Weighting factor for penalty term "RL not in empirical [x.Q1, x.Q2]"
#  Controls degree of penalisation, r.fact <- 0 means no penalisation 
#  
r.fact <-  0              #  Default set by experimentation

#  Weighting factor for penalty term "Non-random fluctuation in TI"
#  Controls degree of penalisation, rt.fact <- 0 means no penalisation 
#  
rt.fact <- 1               # 18.09.2022 
                           # Default set by experimentation

#  Weighting factor for penalty term "too large prediction"  
#  Controls degree of penalisation, w.fact <- 0 means no penalisation 
#  
w.fact <-  1.0             # 1.0 set because of simulation study 

# -------------------------------------------------------------------------
#  TMC: Fit citeria for the choice of of a truncation interval
p.fit.min     <- 0.025   #  05.08.2022: 0.05
                         #  08.10.2022: 0.025

p.rt.min      <- 0.025   #  05.08.2022: 0.05
                         #  08.10.2022: 0.025

#  Accepted range for estimated prevalence
prev.acc.lo <- -0.05     #  -0.02
prev.acc.hi <-  1.05

sol.score.max <- 5

# -------------------------------------------------------------------------
#  Parameters for TUK
IQR.fact     <- 1.5

iter.max.tuk <- 5   

# -------------------------------------------------------------------------
#  False positive rate for confidence intervals 
#  (1-alpha = coverage probability)
alpha <- 0.05

# -------------------------------------------------------------------------
#  When analysing a sequence of replicates: first replicate? 
r.start <-  NA     # (NA: no replicates)

#  When analysing a sequence of replicates: last replicate?
r.ende  <-  NA    # (NA: no replicates)

# -------------------------------------------------------------------------
#  Kernel density estimation
kernel.n     <- 4096   # Number of support points  Pr.lo <- 0.001 
Pr.lo        <- 0.001  # Lower quantile of the interval approximated by kde 
Pr.hi        <- 0.999  # Upper quantile 

qtype        <- 4      # Type of quantile to calculate,
                       # see RLib\Quantile.R

# -------------------------------------------------------------------------
#  Bootstrappimg
nbs          <-  0     # Number of bootstrap replicates for tolerance intervalsle
                       # 0: no TI (instead asymptotic CI)                   
                       # temporarily deactivated!

# -------------------------------------------------------------------------
#  Descriptive output
# Print total age distribution?

print.daytime.distr    <- FALSE  # value distribution over daytime and weekday
print.age.contours     <- FALSE  # contours of value by age
print.daytime.contours <- FALSE  # course of values over taytime and weekday
print.age.distr        <- FALSE  # age distribution

# -------------------------------------------------------------------------
#  Results from which methods wanted in the summmary table per data 
#  file?
#  "npa" will be removed for data not containing the grp variable 
#  (typically only present in simulated data)
#  Do not include "tml", because it will not be executed from TMC, but 
#  externally. TML results are processed only by seg202, but the inclusion of 
#  TML results in seg202 presentations is controlled by do.tml below. 

# gtab.list <- c("bha", "npa", "qqw", "tmc", "tmu", "tuk")   # full list
# gtab.list   <- c("npa", "qqw", "tmc")                      # WW
gtab.list   <- c("tmc")                                      # Usr
gtab.list.n <- length(gtab.list)

#  Plot results for which methods?
# plot.list  <- c("bha", "npa", "qqw", "tmc", "tmu", "tuk")   # full list
# plot.list  <- c("npa", "qqw", "tmc")                      # WW
plot.list    <- c("tmc")                                    # Usr
plot.list.n  <- length(plot.list)

#  Which methods have to be executed? In principle 
# meth.list   <- union(gtab.list, plot.list)
# meth.list.n <- length(meth.list)
#  but "npa" must be removed from all lists, if the grp variable is not 
#  present in the data. This is done in seg050 after reading the data 
#  description.

# -------------------------------------------------------------------------
#  Plot bands of permissible differences in plots of RL vs age?
plot.perm.diff.band <- TRUE

#  If bands of permissible differences are plotted, which equivalence factor?
ef <-  1.28       #  default according to publication

# -------------------------------------------------------------------------
#  Which estimation methods shall be executed?

do.tml <- TRUE    #  external calculation!!!

# -------------------------------------------------------------------------
#  Methods for which to plot the age-RL dependency  --- OBSOLETE
Verlauf.Meth1 <- "tmc"
if ("tmu" %in% gtab.list) 
{ Verlauf.Meth2 <- "tmu" } else
{ Verlauf.Meth2 <- NA }

#  Which plots? Keyword is interpreted in TMC_seg045_PlotRequest.R
 
#plot.level   <- "all"    # Do all plots - not recommended for routine
plot.level <- "auto"      # Most relevant plots. Automatic distinction
                          # between the analysis of single files and 
                          # the analysis of file sequences (simulation).

#  Plot details in certain plots?
plot.details <- FALSE      # TRUE: legend, gridlines, infoblock,
                           # especially in plots 150.10, 160.010
if (user == "WW") { plot.details <- TRUE } 

plot.sol.props <- FALSE    # TRUE: properties of tmc solutions

# -------------------------------------------------------------------------
#   If no age axis parameters for the age axis in plots 150.10, 160.10, set
#   corresponding variables to NA
if (!exists("age.clip.min"))
{ age.clip.min   <- NA
  age.clip.max   <- NA
  age.clip.by1   <- NA
  age.clip.by2   <- NA
}

# -------------------------------------------------------------------------
#   If no axis parameters for the RL axis in plots 150.10, 160.10, set
#   corresponding variables to NA
if (!exists("RL.clip.min"))
{ RL.clip.min   <- NA
  RL.clip.max   <- NA
  RL.clip.by1   <- NA
  RL.clip.by2   <- NA
}

# -------------------------------------------------------------------------
#  For validation only: axis limits in fig 202.010

RL1.clip.min <- NA
RL1.clip.max <- NA
RL1.clip.by1 <- NA
RL1.clip.by2 <- NA

RL2.clip.min <- NA
RL2.clip.max <- NA
RL2.clip.by1 <- NA
RL2.clip.by2 <- NA

# -------------------------------------------------------------------------
#  Graphics file format? (emf, bmp, ...)
figtype <- "bmp"

#  End of default settings
# =================================================================

if (time.check)
{
  DefaultSettings.t5 <- Sys.time()

  if (file.exists(file="../temp/UsedTime.txt")) 
  { file.remove("../temp/UsedTime.txt") }

  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n DefaultSettings; t2-t1;", 
      format(difftime(DefaultSettings.t2, DefaultSettings.t1)))
  cat("\n DefaultSettings; t3-t2;", 
      format(difftime(DefaultSettings.t3, DefaultSettings.t2)))
  cat("\n DefaultSettings; t4-t3;", 
      format(difftime(DefaultSettings.t4, DefaultSettings.t3)))
  cat("\n DefaultSettings; t5-t4;", 
      format(difftime(DefaultSettings.t5, DefaultSettings.t4)))

  cat("\n DefaultSettings; t5-t1;", 
      format(difftime(DefaultSettings.t5, DefaultSettings.t1)))
  sink()
}

# ==========================================================================
if (print.log.message) { cat("%%%   TMC_seg15_DefaultSettings  End\n") }
# ===========================================================================
