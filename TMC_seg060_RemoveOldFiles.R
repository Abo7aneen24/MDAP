#  TMC_seg060_RemoveOldFiles.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment removing old output file names having the the same initial part
#  as outname

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg060_MasterTables  Start\n") }
# ===========================================================================

#  Text files under ../tabs

#  Remove strata text files referring to outname
files <- list.files(path = path.tab.stra, pattern=outname)

if (length(files) > 0)
{for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.tab.stra, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_055   ",paste(path.tab.stra, files[i],sep=""), 
          "  could not be removed\n")
    } 
  }
}

#  Remove file text files referring to outname 
files <- list.files(path = path.tab.file, pattern=outname)

if (length(files) > 0)
{ for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.tab.file, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_055   ",paste(path.tab.file, files[i],sep=""), 
            "  could not be removed\n")
    } 
  }
}

#  Remove scenario text files referring to outname 
files <- list.files(path = path.tab.scen, pattern=outname)
if (length(files) > 0)
{ for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.tab.scen, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_055   ",paste(path.tab.scen, files[i],sep=""), 
            "  could not be removed\n")
    } 
  }
}

#  Plot files under ../figs

#  Remove strata plot files referring to outname
files <- list.files(path = path.fig.stra, pattern=outname)

if (length(files) > 0)
{ for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.fig.stra, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_080   ",paste(path.fig.stra, files[i],sep=""), 
            "  could not be removed\n")
    } 
  }
}
#  Remove file plot files referring to outname 
files <- list.files(path = path.fig.file, pattern=outname)

if (length(files) > 0)
{ for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.fig.file, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_055   ",paste(path.fig.file, files[i],sep=""), 
            "  could not be removed\n")
    } 
  }
}
#  Remove scenario plot files referring to outname 
files <- list.files(path = path.fig.scen, pattern=outname)
if (length(files) > 0)
{ for (i in 1:length(files))
  { 
    ok <- file.remove(paste(path.fig.scen, files[i],sep=""))
    if (!ok) 
    { 
      cat("\nseg_088   ",paste(path.fig.scen, files[i],sep=""), 
            "  could not be removed\n")
    } 
  }
}

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg060_MasterTables  End\n") }
# ===========================================================================


