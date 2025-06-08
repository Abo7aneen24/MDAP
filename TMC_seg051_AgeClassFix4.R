#   TMC_seg051_AgeClassFix3.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Organize TMC analysis by age class * sex
#  Age class is given fixed in age.limits or age.class 

#  TODO
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassFix4  Start\n") }
# ==========================================================================

#   Data subsets are sex*age groups, age limits fixed by user

  if (sex.val.n > 0 &  age.class.n > 0)  #  doubled conditions

  { for (isex in 1:sex.val.n)
    { subset.s <- dataset[ ,spalte.s] == sex.val[isex]
      sexlabel <- sex.val[isex]

      for (iage.class in 1:age.class.n)
      { 
        #  Only if job is not yet done
        if (!A[isex, iage.class] )
        { 
          subset <- subset.s &
                   (age.class[iage.class, "lo"] <= dataset[ ,spalte.a]) & 
                   (dataset[ ,spalte.a]   <= age.class[iage.class, "hi"])

          agelabel <- agelabel.vec[iage.class]  # defined in in 050

          x   <- dataset[subset,spalte.w]
          x.n <- length(x)

          if (!is.na(spalte.a)) { age <- dataset[subset,spalte.a] }
          if (!is.na(spalte.s)) { sex <- toupper(dataset[subset,spalte.s]) }
          if (!is.na(spalte.g)) { grp <- toupper(dataset[subset,spalte.g]) }

          # Auswertung von x, Sex = Ausprägung isex und Alter = iage
          verlaufsplot   <- TRUE

          source("TMC_seg100_Analysis.R")

          #  Bookkeeping
          A[isex, iage.class] <- TRUE 

        }   # if (!A[isex, iage.clas ...
      }     #  iage in ...
    }       #  isex in ...
  }         #  if (sex.val.n > 1 &  age.class.n > 1) ...

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassFix4  End\n") }
# ==========================================================================


