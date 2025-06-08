#  TMC_seg000_Start_Demo.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  For method description see:
#  W. Wosniok, R. Haeckel: A new indirect estimation of reference intervals:
#  truncated minimum chi-square (TMC) approach. 
#  https://doi.org/10.1515/cclm-2018-1341

#  For the programme manual: see the text in the 'text' directory
#  Template starting programme for the analysis of single data files.

#  CHANGE HISTORY
#  08.11.2022 Start
# ===========================================================================

# Read external files

source("TMC_seg010_User.R")            # clean up & define user characteristics
source("TMC_seg015_DefaultSettings.R") # Default parameters

#info.file <- "/data/DataFileInfo.csv"        # default 

                                              # Demo program uses its
                                              # own info.file 
info.file <- "../data/DataFileInfoDemo.csv"   # describes demo data file

FileNo     <-  9999        # Test data, similar to AST

scale.fact <-     1        # scaling factor for measurand

x.lo.limit <-    NA        # Exclude values < x.lo.limit from all analysis
x.hi.limit <-    NA        # Exclude values > x.hi.limit

use.oh     <-    NA        #  Select outpatient /  hospitalized
use.dev    <-    NA        #  Select analytical device 

                           #  Select time range within a day 
ana.hour.min <- NA         #  Example: 06
ana.minu.min <- NA         #  Example: 30
ana.hour.max <- NA         #  Example: 18
ana.minu.max <- NA         #  Example: 45 

                           #  Select weekdays:
use.wday <- "all"          #  "all"   or a subset of 
                           #  c("mo", "tu", "we", "th", "fr", "sa", "su")
                           #  e.g.   use.wday <- c("sa", "su")

                           #  Split data be date-time
split.date <- NA           # Example: "04.02.2019", format is "dd.mm.yyyy"
split.time <- NA           # Example: "13:25",      format is "hh:mm"
use.split  <- NA           # Select par to use
                           # 1: before /  2: after given date-time
                        
                           #  Alternative time selection by dates:
                           #  use only data between start.date, end.date, 
                           #  limits included
                           #       format: "dd.mm.yyyy"
start.date <- NA           #  "01.01.2022" 
end.date   <- NA           #  "30.11.2022"

age.limits <- c(18, 30, 40, 50, 60, 70, 80, 90, 100)

#   Alternative age definition. Classes may overlap.
#age.class  <- matrix(c(18, 29,
#                       30, 49,
#                       50, 69,
#                       70, 100), byrow=TRUE, ncol=2)
#colnames(age.class) <- c("lo", "hi")

#  Select subgroups (strata) to analyze
ana.type1 <- TRUE             # full dataset            
ana.type2 <- TRUE            # by sex group, all ages  
ana.type3 <- TRUE            # all sexes, by age group 
ana.type4 <- TRUE            # by sex*age group  

#  Clip range of measurand in plots. Does not restrict analysis. 
#  Must be given in scaled units
x.clip.min   <-       0    # left limit
x.clip.max   <-      80    # right limit
x.clip.by1   <-      10    # horizontal tickmarks with labels
x.clip.by2   <-       2    # horizontal tickmarks without labels

#  Define axis for plots of RL vs age
age.clip.min   <-       10  # left limit
age.clip.max   <-      100  # right limit
age.clip.by1   <-       10  # horizontal tickmarks with labels
age.clip.by2   <-        5  # horizontal tickmarks without labels
 
#  Define RL axis for plots of RL vs age. Must be given in scaled units.
RL.clip.min   <-       0    # left limit
RL.clip.max   <-      50    # right limit
RL.clip.by1   <-      10    # horizontal tickmarks with labels
RL.clip.by2   <-       2    # horizontal tickmarks without labels

# --------------------------------------------------------------------------
####  Default settings are in      TMC_seg015_DefaultSettings.R
####  (Nearly) everything can be changed by giving a new value here.
####  Put changes of default settings after this line 


# =================================================================
#  Start execution
# =================================================================

source("TMC_seg050_Master.R")      # Master file, controls execution

# =================================================================
# =====  Programme end  ===========================================
# ================================================================
