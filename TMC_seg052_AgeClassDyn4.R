#  TMC_seg052_AgeClassDyn4.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Organize TMC analysis by age class * sex
#  Data subsets are sex*age groups, age group limits derived from data, 
#  based on minimum stratum size = x.n.min, set in 030

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ==============================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassDyn4  Start\n") }
# ==========================================================================

  if (sex.val.n > 1 &  age.class.n > 1)
  { for (isex in 1:sex.val.n)
    { subset.s <- dataset[ ,spalte.s] == sex.val[isex]
      sexlabel <- sex.val[isex]

      #  Find age groups of size x.n.min in this sex group

      age.class <- DynAgeLimits(dataset[subset.s, spalte.a], x.n.min)
      age.class.n <- nrow(age.class)

      if (age.class.n > 1)
      {
        agelabel.vec <- rep(NA, times=age.class.n)

        for (iage.class in 1:age.class.n)
        { #  Age classes
          agelabel.vec[iage.class] <- paste(age.class[iage.class, "lo"],
                                        age.class[iage.class, "hi"],sep="-")
        }
      } else
      { # Only 1 age class 
        agelabel <- "All"
      } 

      if (age.class.n > 1)
      {
        for (iage.class in 1:age.class.n)
        { 
          subset <- subset.s &
                   (age.class[iage.class, "lo"] <= dataset[ ,spalte.a]) & 
                   (dataset[ ,spalte.a]   <= age.class[iage.class, "hi"])

          agelabel <- agelabel.vec[iage.class]  # defined in
                                              # TMC_seg051_AgeClassFix3.R

          x   <- dataset[subset,spalte.w]
          x.n <- length(x)

          if (!is.na(spalte.a)) { age <- dataset[subset,spalte.a] }
          if (!is.na(spalte.s)) { sex <- toupper(dataset[subset,spalte.s]) }
          if (!is.na(spalte.g)) { grp <- toupper(dataset[subset,spalte.g]) }

          # Auswertung von x, Sex = Ausprägung isex und Alter = iage
          verlaufsplot   <- TRUE

          source("TMC_seg100_Analysis.R")

        } #  age.class.n > 1 ...
      }   #  iage in ...
    }     #  isex in ...
  }       #  if (sex.val.n > 1 &  age.class.n > 1) ...

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassDyn4  End\n") }
# ==========================================================================


