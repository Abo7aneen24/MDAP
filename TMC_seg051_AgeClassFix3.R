#   TMC_seg051_AgeClassFix3.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Organize TMC analysis by age class, sexes combined
#  Age class is given fixed in age.limits or age.class 

#  TODO
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassFix3  Start\n") }
# ==========================================================================

#   Data subsets are age groups, limits fixed by user
 
    #  If > 1  age classes, analysis per age class, jointly for all sexes
    if (age.class.n > 1)        #  this is a doubled condition
    { #  Create label for combined sexes
      sexlabel <- use.sex[1]
      isex <- 1
      while(isex < use.sex.n)
      { isex <- isex + 1
        sexlabel   <- paste(sexlabel,use.sex[isex],sep="+")
      }
       
      for (iage.class in 1:age.class.n)
      { 
        #  Only if job is not yet done
        if (!A[sex.val.n+1, iage.class] )
        { subset <- (age.class[iage.class, "lo"] <= dataset[ ,spalte.a]) & 
                    (dataset[ ,spalte.a]   <=  age.class[iage.class, "hi"])

          agelabel <- agelabel.vec[iage.class]
          x   <- dataset[subset,spalte.w]
          x.n <- length(x)

          if (!is.na(spalte.a)) { age <- dataset[subset,spalte.a] }
          if (!is.na(spalte.s)) { sex <- toupper(dataset[subset,spalte.s]) }
          if (!is.na(spalte.g)) { grp <- toupper(dataset[subset,spalte.g]) }

          # Auswertung von x, Alter = Ausprägung iage
          verlaufsplot <- TRUE

          source("TMC_seg100_Analysis.R")

          #  Bookkeeping
          A[sex.val.n+1, iage.class] <- TRUE 

        }  #  (!A[sex.val.n+1,  ...
      } #  iage.class in 1:age.class.n
    }   #  age.class.n > 1  

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassFix3  End\n") }
# ==========================================================================

