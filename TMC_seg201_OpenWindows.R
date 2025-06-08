#  TMC_seg201_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  01.08.2021 Additional devices for plots with contour background
#             Note special definition for contour legends (012, 014, 016) 
#  05.06.2021 Graphics devices changed to X11()
#  18.04.2020 File names corrected, 202.040, 202.050 added
#  17.12.2020 Start
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg201_OpenWindows.R  Start\n") }
# ============================================================================

# ============================================================================
# Open.202.010.R
# ============================================================================
  if (plot.fig202.010)
  { #  Plot is wanted
    if (!exists("fig202.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.010 <- dev.cur()
    }
    fig202.010.file <- paste(path.fig.scen.Id, outname.stra,"-F202.010.",
                             figtype,sep="")
    if (file.exists(fig202.010.file)) { file.remove(fig202.010.file) }
  } else
  { #  Plot is not wanted
    fig202.010 <- NA
  }


# ============================================================================
# Open.202.011.R
# ============================================================================
  if (plot.fig202.011)
  { #  Plot is wanted
    if (!exists("fig202.011"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.011 <- dev.cur()
    }
    fig202.011.file <- paste(path.fig.scen.Id, outname.stra,"-F202.011.",
                             figtype,sep="")
    if (file.exists(fig202.011.file)) { file.remove(fig202.011.file) }
  } else
  { #  Plot is not wanted
    fig202.011 <- NA
  }


# ============================================================================
# Open.202.012.R
# ============================================================================
  if (plot.fig202.012)
  { #  Plot is wanted
    if (!exists("fig202.012"))
    { #  Open graph new window
      X11(width=3.5, height=7)    #  win.graph()
      fig202.012 <- dev.cur()
    }
    fig202.012.file <- paste(path.fig.scen.Id, outname.stra,"-F202.012.",
                             figtype,sep="")
    if (file.exists(fig202.012.file)) { file.remove(fig202.012.file) }
  } else
  { #  Plot is not wanted
    fig202.012 <- NA
  }


# ============================================================================
# Open.202.013.R
# ============================================================================
  if (plot.fig202.013)
  { #  Plot is wanted
    if (!exists("fig202.013"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.013 <- dev.cur()
    }
    fig202.013.file <- paste(path.fig.scen.Id, outname.stra,"-F202.013.",
                             figtype,sep="")
    if (file.exists(fig202.013.file)) { file.remove(fig202.013.file) }
  } else
  { #  Plot is not wanted
    fig202.013 <- NA
  }


# ============================================================================
# Open.202.014.R
# ============================================================================
  if (plot.fig202.014)
  { #  Plot is wanted
    if (!exists("fig202.014"))
    { #  Open graph new window
      X11(width=3.5, height=7)    #  win.graph()
      fig202.014 <- dev.cur()
    }
    fig202.014.file <- paste(path.fig.scen.Id, outname.stra,"-F202.014.",
                             figtype,sep="")
    if (file.exists(fig202.014.file)) { file.remove(fig202.014.file) }
  } else
  { #  Plot is not wanted
    fig202.014 <- NA
  }


# ============================================================================
# Open.202.015.R
# ============================================================================
  if (plot.fig202.015)
  { #  Plot is wanted
    if (!exists("fig202.015"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.015 <- dev.cur()
    }
    fig202.015.file <- paste(path.fig.scen.Id, outname.stra,"-F202.015.",
                             figtype,sep="")
    if (file.exists(fig202.015.file)) { file.remove(fig202.015.file) }
  } else
  { #  Plot is not wanted
    fig202.015 <- NA
  }


# ============================================================================
# Open.202.016.R
# ============================================================================
  if (plot.fig202.016)
  { #  Plot is wanted
    if (!exists("fig202.016"))
    { #  Open graph new window
      X11(width=3.5, height=7)    #  win.graph()
      fig202.016 <- dev.cur()
    }
    fig202.016.file <- paste(path.fig.scen.Id, outname.stra,"-F202.016.",
                             figtype,sep="")
    if (file.exists(fig202.016.file)) { file.remove(fig202.016.file) }
  } else
  { #  Plot is not wanted
    fig202.016 <- NA
  }


# ============================================================================
# Open.202.020.R
# ============================================================================
  if (plot.fig202.020)
  { #  Plot is wanted
    if (!exists("fig202.020"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.020 <- dev.cur()
    }
    fig202.020.file <- paste(path.fig.scen.Id, outname.stra,"-F202.020.",
                             figtype,sep="")
    if (file.exists(fig202.020.file)) { file.remove(fig202.020.file) }
  } else
  { #  Plot is not wanted
    fig202.020 <- NA
  }


# ============================================================================
# Open.202.030.R
# ============================================================================
  if (plot.fig202.030)
  { #  Plot is wanted
    if (!exists("fig202.030"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.030 <- dev.cur()
    }
    fig202.030.file <- paste(path.fig.scen.Id, outname.stra,"-F202.030.",
                             figtype,sep="")
    if (file.exists(fig202.030.file)) { file.remove(fig202.030.file) }
  } else
  { #  Plot is not wanted
    fig202.030 <- NA
  }

# ============================================================================
# Open.202.040.R
# ============================================================================
  if (plot.fig202.040)
  { #  Plot is wanted
    if (!exists("fig202.040"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.040 <- dev.cur()
    }
    fig202.040.file <- paste(path.fig.scen.Id, outname.stra,"-F202.040.",
                             figtype,sep="")
    if (file.exists(fig202.040.file)) { file.remove(fig202.040.file) }
  } else
  { #  Plot is not wanted
    fig202.040 <- NA
  }

# ============================================================================
# Open.202.050.R
# ============================================================================
  if (plot.fig202.050)
  { #  Plot is wanted
    if (!exists("fig202.050"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig202.050 <- dev.cur()
    }
    fig202.050.file <- paste(path.fig.scen.Id, outname.stra,"-F202.050.",
                             figtype,sep="")
    if (file.exists(fig202.050.file)) { file.remove(fig202.050.file) }
  } else
  { #  Plot is not wanted
    fig202.050 <- NA
  }

# ============================================================================
# Open.203.010.R
# ============================================================================
  if (plot.fig203.010)
  { #  Plot is wanted
    if (!exists("fig203.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig203.010 <- dev.cur()
    }
    fig203.010.file <- paste(path.fig.scen.Id, outname.stra,"-F203.010.",
                             figtype,sep="")
    if (file.exists(fig203.010.file)) { file.remove(fig203.010.file) }
  } else
  { #  Plot is not wanted
    fig203.010 <- NA
  }

# ============================================================================
# Open.203.020.R
# ============================================================================
  if (plot.fig203.020)
  { #  Plot is wanted
    if (!exists("fig203.020"))
    { #  Open graph new window
      X11(width=3.5, height=7)    #  win.graph()
      fig203.020 <- dev.cur()
    }
    fig203.020.file <- paste(path.fig.scen.Id, outname.stra,"-F203.020.",
                             figtype,sep="")
    if (file.exists(fig203.020.file)) { file.remove(fig203.020.file) }
  } else
  { #  Plot is not wanted
    fig203.020 <- NA
  }

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg201_OpenWindows.R  End\n") }
# ============================================================================
