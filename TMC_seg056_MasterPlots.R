#  TMC_seg051_MasterPlots.R 

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Organize plots referring to the entire data file
#  Calculate splines and permissible differences
#  Prepare entries in gtab

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ==========================================================================
if (print.log.message) { cat("%%%   TMC_seg056_MasterPlots  Start\n") }
# ===========================================================================

  #  Find plot scales for the methods involved
  #  Plot scale depends only on the range of the methods to plot
  
  Verlauf.liste <- c(Verlauf.Meth1,Verlauf.Meth2)
  Verlauf.liste <- Verlauf.liste[!is.na(Verlauf.liste)]
  ok <- gtab[ ,"method"] %in% Verlauf.liste

  #  If there are user definitions for the RL axis use these, else
  #  calculate own values
  if (!is.na(RL.clip.min))
  {
    yplotmin150.010 <- RL.clip.min
    yplotmax150.010 <- RL.clip.max
    yplotby1150.010 <- RL.clip.by1
    yplotby2150.010 <- RL.clip.by2
  } else
  {
    yplotmin150.010 <- min(as.numeric(gtab[ok,"RL1.cilo"]),na.rm=TRUE)
    yplotmin150.010 <- 0
    yplotmax150.010 <- max(as.numeric(gtab[ok,"RL2.cihi"]),na.rm=TRUE)
    yplotby1150.010 <- NA
    yplotby2150.010 <- NA
  }

  dx <- 0.1       #  Versatz method 1 gegen method 2

  # --------------------------------------------------------------------------
  #  More than  > 1 Sex and  > 1 age class?

  if (ana.type4 & sex.val.n > 1 & age.class.n > 1 & 
      (!is.na(Verlauf.Meth1) | !is.na(Verlauf.Meth2) ))
  {
    if (print.log.message) { cat("%%%   TMC_seg056_MasterPlots  1000\n") }

    #  Verlauf über Altersklassen plotten, 
    #  getrennt nach Sex, 
    #  2 Methoden in ein Bild
    #  Altersangabe = Grenzen der Altersklassen

    #  Methode zum Plotten als Verlauf wurden bei den speziellen 
    #  Benutzereingaben ausgewählt (Verlauf.Meth1, Verlauf.Meth2)
    #  Farben aus dem style file holen

    meth1col <- NA
    meth1pch <- NA
    meth2col <- NA
    meth2pch <- NA

    if (!is.na(Verlauf.Meth1))
    { meth1col <- meth.dsg[which(meth.dsg[ ,1] == Verlauf.Meth1),2] 
      meth1pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1] == 
                                                 Verlauf.Meth1),3])
    }
    if (!is.na(Verlauf.Meth2))
    { meth2col <- meth.dsg[which(meth.dsg[ ,1] == Verlauf.Meth2),2] 
      meth2pch <- as.numeric(meth.dsg[which(meth.dsg[ ,1] == 
                                                  Verlauf.Meth2),3])
    }

  }   # if (sex.val.n > 1 & age.class.n > 1)

  # ----------------------------------------------------------------------
  #  Verlauf mit Glättung, Kurven für F/ M/ F+M, 1 method
  if (ana.type4 & age.class.n > 1)
  {
    #  Verlauf über Altersklassen plotten, 
    #  getrennt nach Sex (F / M / F+M), 
    #  1 methode in ein Bild
    #  Altersangabe = Grenzen der Altersklassen
    #  Konfidenz/Toleranzintervall nach Verfügbarkeit

    #  Bilder mit reduzierter Information, F und M getrennt, 
    #  für Veröffentlichung

    #  Produce plots of RL vs age for all requested methods 
    #  Graph windows

    # ------------------------------------------------------------------------
    #  Plots for F and M separately 

    source("TMC_seg149_OpenWindows.R")

    if ("npa" %in% plot.list)
    {
      RL.age.meth <- "npa"
      PlotId      <- fig150.010
      PlotFile    <- fig150.010.file
      source("TMC_seg150_PlotRL_vs_Age_by_S.R")
    }
    if ("qqw" %in% plot.list)
    {
      RL.age.meth <- "qqw"
      PlotId      <- fig150.020
      PlotFile    <- fig150.020.file
      source("TMC_seg150_PlotRL_vs_Age_by_S.R")
    }
    if ("bha" %in% plot.list)
    {
      RL.age.meth <- "bha"
      PlotId      <- fig150.030
      PlotFile    <- fig150.030.file
      source("TMC_seg150_PlotRL_vs_Age_by_S.R")
    }
    if ("tuk" %in% plot.list)
    {
      RL.age.meth <- "tuk"
      PlotId      <- fig150.040
      PlotFile    <- fig150.040.file
      source("TMC_seg150_PlotRL_vs_Age_by_S.R")
    }
    if ("tmc" %in% plot.list)
    {
      RL.age.meth <- "tmc"

      if (!plot.perm.diff.band)
      { PlotId      <- fig150.050
        PlotFile    <- fig150.050.file
        source("TMC_seg150_PlotRL_vs_Age_by_S.R")
      }

      if (plot.perm.diff.band)
      { PlotId      <- fig150.051
        PlotFile    <- fig150.051.file
        source("TMC_seg150_PlotRL_vs_Age_by_S.R")
      }

      RL.spline.tmc <- NA
      if (gtab.M1.spl$spline.ok) RL.spline.tmc <- RL.spline
    }

    if ("tmu" %in% plot.list)
    {
      RL.age.meth <- "tmu"
      PlotId      <- fig150.060
      PlotFile    <- fig150.060.file
      source("TMC_seg150_PlotRL_vs_Age_by_S.R")
    }
  }   
  
  # ----------------------------------------------------------------------

  if ((ana.type3 | ana.type4) & age.class.n > 1 & 
      (!is.na(Verlauf.Meth1) | !is.na(Verlauf.Meth2) ))
  {
    if (print.log.message) { cat("%%%   TMC_seg056_MasterPlots  2000\n") }
    #  Plots for F and M jointly 

    source("TMC_seg159_OpenWindows.R")

    if ("npa" %in% plot.list)
    {
      RL.age.meth <- "npa"
      PlotId      <- fig160.010
      PlotFile    <- fig160.010.file
      source("TMC_seg160_PlotRL_vs_Age.R")
    }
    if ("qqw" %in% plot.list)
    {
      RL.age.meth <- "qqw"
      PlotId      <- fig160.020
      PlotFile    <- fig160.020.file
      source("TMC_seg160_PlotRL_vs_Age.R")
    }
    if ("bha" %in% plot.list)
    {
      RL.age.meth <- "bha"
      PlotId      <- fig160.030
      PlotFile    <- fig160.030.file
      source("TMC_seg160_PlotRL_vs_Age.R")
    }
    if ("tuk" %in% plot.list)
    {
      RL.age.meth <- "tuk"
      PlotId      <- fig160.040
      PlotFile    <- fig160.040.file
      source("TMC_seg160_PlotRL_vs_Age.R")
    }
    if ("tmc" %in% plot.list)
    {
      RL.age.meth <- "tmc"

      if (!plot.perm.diff.band)
      { PlotId      <- fig160.050
        PlotFile    <- fig160.050.file
        source("TMC_seg160_PlotRL_vs_Age.R")
      }

      if (plot.perm.diff.band)
      { PlotId      <- fig160.051
        PlotFile    <- fig160.051.file
        source("TMC_seg160_PlotRL_vs_Age.R")
      }

      RL.spline.tmc <- NA
      if (gtab.M1.spl$spline.ok) RL.spline.tmc <- RL.spline
    }

    if ("tmu" %in% plot.list)
    {
      RL.age.meth <- "tmu"
      PlotId      <- fig160.060
      PlotFile    <- fig160.060.file
      source("TMC_seg160_PlotRL_vs_Age.R")
    }
  }   # if ((ana.type3 | ana.type4) & age.class.n > 1 ...

# ==========================================================================
if (print.log.message) { cat("%%%   TMC_seg056_MasterPlots  End\n") }
# ===========================================================================

