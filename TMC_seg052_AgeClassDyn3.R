#  TMC_seg052_AgeClassDyn3.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Organize TMC analysis by age class, sexes combined
#  Data subsets are age groups, group limits derived from data, 
#  based on minimum stratum size = x.n.min, set in 030

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ==============================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassDyn3  Start\n") }
# ==========================================================================

#  Find age groups of size x.n.min

age.class   <- DynAgeLimits(dataset[spalte.a], x.n.min)
age.class.n <- nrow(age.class)

#  From here on: nearly the same as TMC_seg051_AgeClassFix3.R 

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

    #  If > 1  age classes, analysis per age class, jointly for all sexes
    if (age.class.n > 1)
    { sexlabel <- use.sex[1]
      isex <- 1
      while(isex < use.sex.n)
      { isex <- isex + 1
        sexlabel   <- paste(sexlabel,use.sex[isex],sep="+")
      }
      for (iage.class in 1:age.class.n)
      { 
        subset <- (age.class[iage.class, "lo"] <= dataset[ ,spalte.a]) & 
                  (dataset[ ,spalte.a]   <=  age.class[iage.class, "hi"])

        agelabel <- agelabel.vec[iage.class]
        x   <- dataset[subset,spalte.w]
        x.n <- length(x)

        if (!is.na(spalte.a)) { age <- dataset[subset,spalte.a] }
        if (!is.na(spalte.s)) { sex <- toupper(dataset[subset,spalte.s]) }
        if (!is.na(spalte.g)) { grp <- toupper(dataset[subset,spalte.g]) }

        # Auswertung von x, Alter = Ausprägung iage
        verlaufsplot <- TRUE
        subset.type  <- 3
        source("TMC_seg100_Analysis.R")

      } #  iage.class in 1:age.class.n
    }   #  age.class.n > 1  

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg051_AgeClassDyn3  End\n") }
# ==========================================================================

