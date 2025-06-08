#  TMC_seg042_SettingsSummary.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Produce a summary of all input parameters before reading the datafile 
#  information, write to log

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg042_SettingsSummary  Start\n") }
# ===========================================================================

#  ..........................................................................
#  Show programm segment versions

#  Next sections takes everything from folder /prog
seg.list <- dir()
seg.list.n <- length(seg.list)

j <- 0
if (seg.list.n > 0)
{ for (i in 1:seg.list.n) 
  { 
    temp <- seg.list[i]
    if (substr(temp, 1, 7) == "TMC_seg")
    { #  We are looking at a program segement
      j <- j + 1  
      # cat("\nInspecting", seg.list[i])
      if (j == 1)
      {
        finfo <- file.info(temp)["mtime"]
      }  else
      {
        finfo <- rbind(finfo, file.info(temp)["mtime"])
      }
    }     
  }
}

sink(file=logfile, split=FALSE, append=FALSE)

cat("\n =====================================================================",
    "\n  Truncated minimum chi-square estimation of reference interval limits",
    "\n                               (TMC)",
    "\n =====================================================================",
    "\n",
    "\n   Identification of this run: ", RunId,  
    "\n   Start at                    ", format(time.start),
    "\n   Execution in directory      ", getwd(),
    "\n",
    "\nProgramme segments versions\n")

print(finfo)

cat("\nInitial parameter settings\n")
sink()

if (file.exists("../temp/ShowValue1.txt")) 
{ 
  ok <- file.remove("../temp/ShowValue1.txt") 
}
if (file.exists("../temp/ShowValue2.txt")) 
{ 
  ok <- file.remove("../temp/ShowValue2.txt") 
}

obj1.list <- ls()

for (obj in obj1.list)
{
  # cat(obj, class(get(obj)), length(get(obj)), "\n")
  if ( (length(get(obj)) == 1) &&
       (class(get(obj)) %in% c("character", "numeric", "logical")))
  {
    cat('cat("', obj[1], '", ', '";",   ', obj[1], ', fill=256)', 
        fill=256, sep="", file="../temp/ShowValue1.txt", append=TRUE)
  }
}

sink("../temp/ShowValue2.txt", split=FALSE, append=FALSE)
source("../temp/ShowValue1.txt")
sink()

if(exists("A1")) { rm(list="A1") }

A1 <- read.table(file="../temp/ShowValue2.txt", header=FALSE, sep=";", 
                stringsAsFactors=FALSE)
names(A1) <- c("VarName", "Value1") 

sink(file=logfile,split=FALSE, append=TRUE)

print(A1)

cat("\n",
    "\n =====================================================================",
    "\n")

sink()

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg042_SettingsSummary  End\n") }
# ===========================================================================
