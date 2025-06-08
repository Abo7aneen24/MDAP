#  TMC_seg094_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg094_OpenWindows.R  Start\n") }
# ============================================================================


# ============================================================================
# Open.095.010.R
# ==============================================================================
  if (plot.fig095.010)
  { #  Plot is wanted
    if (!exists("fig095.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig095.010 <- dev.cur()
    }
    fig095.010.file <- paste(path.fig.file,outname.file,"-F095.010.",
                             figtype,sep="")
    if (file.exists(fig095.010.file)) { file.remove(fig095.010.file) }
  } else
  { #  Plot is not wanted
    fig095.010 <- NA
  }


# ============================================================================
# Open.095.020.R
# ==============================================================================
  if (plot.fig095.020)
  { #  Plot is wanted
    if (!exists("fig095.020"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig095.020 <- dev.cur()
    }
    fig095.020.file <- paste(path.fig.file,outname.file,"-F095.020.",
                             figtype,sep="")
    if (file.exists(fig095.020.file)) { file.remove(fig095.020.file) }
  } else
  { #  Plot is not wanted
    fig095.020 <- NA
  }


# ============================================================================
# Open.095.030.R
# ==============================================================================
  if (plot.fig095.030)
  { #  Plot is wanted
    if (!exists("fig095.030"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig095.030 <- dev.cur()
    }
    fig095.030.file <- paste(path.fig.file,outname.file,"-F095.030.",
                             figtype,sep="")
    if (file.exists(fig095.030.file)) { file.remove(fig095.030.file) }
  } else
  { #  Plot is not wanted
    fig095.030 <- NA
  }


# ============================================================================
# Open.095.040.R
# ==============================================================================
  if (plot.fig095.040)
  { #  Plot is wanted
    if (!exists("fig095.040"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig095.040 <- dev.cur()
    }
    fig095.040.file <- paste(path.fig.file,outname.file,"-F095.040.",
                             figtype,sep="")
    if (file.exists(fig095.040.file)) { file.remove(fig095.040.file) }
  } else
  { #  Plot is not wanted
    fig095.040 <- NA
  }

# ==============================================================================
if (print.log.message) { cat("%%%   TMC_seg094_OpenWindows.R  End\n") }
# ==============================================================================


