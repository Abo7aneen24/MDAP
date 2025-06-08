#  TMC_seg075_FinalParSettings.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Produce a summary of all input parameters including those that depend
#  on the data file information, write to log

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg075_FinalParSettings  Start\n") }
# ===========================================================================

if (file.exists("../temp/ShowValue1.txt")) 
{ 
  ok <- file.remove("../temp/ShowValue1.txt") 
}
if (file.exists("../temp/ShowValue2.txt")) 
{ 
  ok <- file.remove("../temp/ShowValue2.txt") 
}

obj2.list <- ls()
obj2.list

for (obj in obj2.list)
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

if(exists("A2")) { rm(list="A2") }

A2 <- read.table(file="../temp/ShowValue2.txt", header=FALSE, sep=";", 
                stringsAsFactors=FALSE)
names(A2) <- c("VarName", "Value2") 

A  <- merge(A1, A2, by="VarName", all.y=TRUE) 

#  If this is a simulation: Put A in tab/3_scen/<folder>
if (eval.rep) { save(A, file=outname.scen.log.rda) }

#  Show only new or modified values in txt file
new.or.mod <-      (A[ ,"Value1"] !=  A[ ,"Value2"] ) |
              is.na(A[ ,"Value1"] !=  A[ ,"Value2"] )

sink(file=logfile,split=FALSE, append=TRUE)

cat("\nFinal parameter settings, new or different from previous list\n")

print(A[new.or.mod, ])

cat("\n",
    "\n =====================================================================",
    "\n")

sink()

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg075_FinalParSettings  End\n") }
# ===========================================================================
