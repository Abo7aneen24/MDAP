#  TMC_seg149_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  05.06.2021 Graphics devices changed to X11()
#  03.01.2021 Target directory corrected
#  17.12.2020 Start
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg149_OpenWindows.R  Start\n") }
# ============================================================================

# ============================================================================
# Open.150.010.R
# ==============================================================================
  if (plot.fig150.010)
  { #  Plot is wanted
    if (!exists("fig150.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.010 <- dev.cur()
    }
    fig150.010.file <- paste(path.fig.file,outname,"-F150.npa.",
                             figtype,sep="")
    if (file.exists(fig150.010.file)) { file.remove(fig150.010.file) }
  } else
  { #  Plot is not wanted
    fig150.010 <- NA
  }

# ============================================================================
# Open.150.020.R
# ==============================================================================
  if (plot.fig150.020)
  { #  Plot is wanted
    if (!exists("fig150.020"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.020 <- dev.cur()
    }
    fig150.020.file <- paste(path.fig.file,outname,"-F150.qqw.",
                             figtype,sep="")
    if (file.exists(fig150.020.file)) { file.remove(fig150.020.file) }
  } else
  { #  Plot is not wanted
    fig150.020 <- NA
  }

# ============================================================================
# Open.150.030.R
# ==============================================================================
  if (plot.fig150.030)
  { #  Plot is wanted
    if (!exists("fig150.030"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.030 <- dev.cur()
    }
    fig150.030.file <- paste(path.fig.file,outname,"-F150.bha.",
                             figtype,sep="")
    if (file.exists(fig150.030.file)) { file.remove(fig150.030.file) }
  } else
  { #  Plot is not wanted
    fig150.030 <- NA
  }

# ============================================================================
# Open.150.040.R
# ==============================================================================
  if (plot.fig150.040)
  { #  Plot is wanted
    if (!exists("fig150.040"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.040 <- dev.cur()
    }
    fig150.040.file <- paste(path.fig.file,outname,"-F150.tuk.",
                             figtype,sep="")
    if (file.exists(fig150.040.file)) { file.remove(fig150.040.file) }
  } else
  { #  Plot is not wanted
    fig150.040 <- NA
  }

# ============================================================================
# Open.150.050.R
# ==============================================================================
  if (plot.fig150.050)
  { #  Plot is wanted
    if (!exists("fig150.050"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.050 <- dev.cur()
    }
    fig150.050.file <- paste(path.fig.file,outname,"-F150.tmc.",
                             figtype,sep="")
    if (file.exists(fig150.050.file)) { file.remove(fig150.050.file) }
  } else
  { #  Plot is not wanted
    fig150.050 <- NA
  }

# ============================================================================
# Open.150.051.R
# ==============================================================================
  if (plot.fig150.051)
  { #  Plot is wanted
    if (!exists("fig150.051"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.051 <- dev.cur()
    }
    fig150.051.file <- paste(path.fig.file,outname,"-F151.tmc.",
                             figtype,sep="")
    if (file.exists(fig150.051.file)) { file.remove(fig150.051.file) }
  } else
  { #  Plot is not wanted
    fig150.051 <- NA
  }

# ============================================================================
# Open.150.060.R
# ==============================================================================
  if (plot.fig150.060)
  { #  Plot is wanted
    if (!exists("fig150.060"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig150.060 <- dev.cur()
    }
    fig150.060.file <- paste(path.fig.file,outname,"-F150.tmu.",
                             figtype,sep="")
    if (file.exists(fig150.060.file)) { file.remove(fig150.060.file) }
  } else
  { #  Plot is not wanted
    fig150.060 <- NA
  }


# ============================================================================
# ============================================================================
