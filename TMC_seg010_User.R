info.file <- "../Data/DataFileInfoLab.csv"    #  routine for real data

#  TMC_seg010_User.R
#  
#  (c) wwosniok@math.uni-bremen.de
#  Truncated minimum chi-square estimation
#
#  Clean up, define user-specific settings, particularly the amount of 
#  control printout and paths for input / output files

#  TODO       

#  CHANGE HISTORY

#  08.11.2022 Last modification (clean up), installation for user Usr only
# =================================================================


#  Remove (nearly) all objects, all graphics windows
#  3 variables must be kept as loop controls when analysing file sequences

idx <- NA          # to get it removed below
varlist <- ls()

#  Commands below do nonsense if a searched names does not exist
#varlist <- varlist[-which(varlist=="x.tr.prop.min0")]
#varlist <- varlist[-which(varlist=="smooth.hist1.0")]
#varlist <- varlist[-which(varlist=="subsample.n0")]

idx <- which(varlist=="x.tr.prop.min0")
if (length(idx) > 0) { varlist <- varlist[-idx] }
idx <- which(varlist=="smooth.hist1.0")
if (length(idx) > 0) { varlist <- varlist[-idx] }
idx <- which(varlist=="subsample.n0")
if (length(idx) > 0) { varlist <- varlist[-idx] }

rm(list=varlist)

graphics.off()

# =================================================================
#  remove all output diversions

nsink <- sink.number()
if (nsink > 0)
{ for (isink in 1:nsink) { sink()  } }

# =================================================================
#  Which user? Controls path settings, the amount of output etc.

#  Upper case / lower case matters in the user name
#  For distribution:  Usr  (the user may create new user names)
#  For development:   WW
# 
#  If user != "WW": if this programme runs on the WW computer, paths are
#  redirected such that the programme can be tested on WWs machine 

user <- "Usr"

#  Ask for machine settings

sysinfo <- Sys.info()

#  Arrange user-specific paths etc.

if (user == "Usr")
{
  #  Path and name of the file containing data file information
  info.file <- "DataFileInfo.csv"
 
  #  Directories for plot files
  path.fig.stra   <- "../figs/1_stra/"
  path.fig.file   <- "../figs/2_file/"
  path.fig.scen   <- "../figs/3_scen/"

  #  Directories for text and data output
  path.tab.stra   <- "../tabs/1_stra/"
  path.tab.file   <- "../tabs/2_file/"
  path.tab.scen   <- "../tabs/3_scen/"

  #  Print log messages during execution? 
  print.log.message <- FALSE

  #  Print all 6 result tables?
  print.6tables <- FALSE            #  TRUE activates long output
}

# =================================================================
