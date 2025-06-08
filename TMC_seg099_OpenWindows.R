#  TMC_seg099_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg099_OpenWindows.R  Start\n") }
# ============================================================================

# ============================================================================
# Open.100.001.R
# ==============================================================================
  if (plot.fig100.001)
  { #  Plot is wanted
    if (!exists("fig100.001"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.001 <- dev.cur()
    }
    fig100.001.file <- paste(path.fig.stra,outname.stra,"-F100.001.",
                             figtype,sep="")
    if (file.exists(fig100.001.file)) { file.remove(fig100.001.file) }
  } else
  { #  Plot is not wanted
    fig100.001 <- NA
  }

# ============================================================================
# Open.100.002.R
# ==============================================================================
  if (plot.fig100.002)
  { #  Plot is wanted
    if (!exists("fig100.002"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.002 <- dev.cur()
    }
    fig100.002.file <- paste(path.fig.stra,outname.stra,"-F100.002.",
                             figtype,sep="")
    if (file.exists(fig100.002.file)) { file.remove(fig100.002.file) }
  } else
  { #  Plot is not wanted
    fig100.002 <- NA
  }

# ============================================================================
# Open.100.003.R
# ==============================================================================
  if (plot.fig100.003)
  { #  Plot is wanted
    if (!exists("fig100.003"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.003 <- dev.cur()
    }
    fig100.003.file <- paste(path.fig.stra,outname.stra,"-F100.003.",
                             figtype,sep="")
    if (file.exists(fig100.003.file)) { file.remove(fig100.003.file) }
  } else
  { #  Plot is not wanted
    fig100.003 <- NA
  }

# ============================================================================
# Open.100.004.R
# ==============================================================================
  if (plot.fig100.004)
  { #  Plot is wanted
    if (!exists("fig100.004"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.004 <- dev.cur()
    }
    fig100.004.file <- paste(path.fig.stra,outname.stra,"-F100.004.",
                             figtype,sep="")
    if (file.exists(fig100.004.file)) { file.remove(fig100.004.file) }
  } else
  { #  Plot is not wanted
    fig100.004 <- NA
  }


# ============================================================================
# Open.100.005.R
# ============================================================================
  if (plot.fig100.005)
  { #  Plot is wanted
    if (!exists("fig100.005"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.005 <- dev.cur()
    }
    fig100.005.file <- paste(path.fig.stra,outname.stra,"-F100.005.",
                             figtype,sep="")
    if (file.exists(fig100.005.file)) { file.remove(fig100.005.file) }
  } else
  { #  Plot is not wanted
    fig100.005 <- NA
  }

# ============================================================================
# Open.100.006.R
# ============================================================================
  if (plot.fig100.006)
  { #  Plot is wanted
    if (!exists("fig100.006"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.006 <- dev.cur()
    }
    fig100.006.file <- paste(path.fig.stra,outname.stra,"-F100.006.",
                             figtype,sep="")
    if (file.exists(fig100.006.file)) { file.remove(fig100.006.file) }
  } else
  { #  Plot is not wanted
    fig100.006 <- NA
  }


# ============================================================================
# Open.100.007.R
# ============================================================================
  if (plot.fig100.007)
  { #  Plot is wanted
    if (!exists("fig100.007"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.007 <- dev.cur()
    }
    fig100.007.file <- paste(path.fig.stra,outname.stra,"-F100.007.",
                             figtype,sep="")
    if (file.exists(fig100.007.file)) { file.remove(fig100.007.file) }
  } else
  { #  Plot is not wanted
    fig100.007 <- NA
  }


# ============================================================================
# Open.100.008.R
# ============================================================================
  if (plot.fig100.008)
  { #  Plot is wanted
    if (!exists("fig100.008"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.008 <- dev.cur()
    }
    fig100.008.file <- paste(path.fig.stra,outname.stra,"-F100.008.",
                             figtype,sep="")
    if (file.exists(fig100.008.file)) { file.remove(fig100.008.file) }
  } else
  { #  Plot is not wanted
    fig100.008 <- NA
  }


# ============================================================================
# Open.100.009.R
# ============================================================================
  if (plot.fig100.009)
  { #  Plot is wanted
    if (!exists("fig100.009"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.009 <- dev.cur()
    }
    fig100.009.file <- paste(path.fig.stra,outname.stra,"-F100.009.",
                             figtype,sep="")
    if (file.exists(fig100.009.file)) { file.remove(fig100.009.file) }
  } else
  { #  Plot is not wanted
    fig100.009 <- NA
  }


# ============================================================================
# Open.100.010.R
# ============================================================================
  if (plot.fig100.010)
  { #  Plot is wanted
    if (!exists("fig100.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.010 <- dev.cur()
    }
    fig100.010.file <- paste(path.fig.stra,outname.stra,"-F100.010.",
                             figtype,sep="")
    if (file.exists(fig100.010.file)) { file.remove(fig100.010.file) }
  } else
  { #  Plot is not wanted
    fig100.010 <- NA
  }


# ============================================================================
# Open.100.011.R
# ============================================================================
  if (plot.fig100.011)
  { #  Plot is wanted
    if (!exists("fig100.011"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.011 <- dev.cur()
    }
    fig100.011.file <- paste(path.fig.stra,outname.stra,"-F100.011.",
                             figtype,sep="")
    if (file.exists(fig100.011.file)) { file.remove(fig100.011.file) }
  } else
  { #  Plot is not wanted
    fig100.011 <- NA
  }


# ============================================================================
# Open.100.012.R
# ============================================================================
  if (plot.fig100.012)
  { #  Plot is wanted
    if (!exists("fig100.012"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.012 <- dev.cur()
    }
    fig100.012.file <- paste(path.fig.stra,outname.stra,"-F100.012.",
                             figtype,sep="")
    if (file.exists(fig100.012.file)) { file.remove(fig100.012.file) }
  } else
  { #  Plot is not wanted
    fig100.012 <- NA
  }


# ============================================================================
# Open.100.013.R
# ============================================================================
  if (plot.fig100.013)
  { #  Plot is wanted
    if (!exists("fig100.013"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.013 <- dev.cur()
    }
    fig100.013.file <- paste(path.fig.stra,outname.stra,"-F100.013.",
                             figtype,sep="")
    if (file.exists(fig100.013.file)) { file.remove(fig100.013.file) }
  } else
  { #  Plot is not wanted
    fig100.013 <- NA
  }


# ============================================================================
# Open.100.014.R
# ============================================================================
  if (plot.fig100.014)
  { #  Plot is wanted
    if (!exists("fig100.014"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.014 <- dev.cur()
    }
    fig100.014.file <- paste(path.fig.stra,outname.stra,"-F100.014.",
                             figtype,sep="")
    if (file.exists(fig100.014.file)) { file.remove(fig100.014.file) }
  } else
  { #  Plot is not wanted
    fig100.014 <- NA
  }


# ============================================================================
# Open.100.015.R
# ============================================================================
  if (plot.fig100.015)
  { #  Plot is wanted
    if (!exists("fig100.015"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.015 <- dev.cur()
    }
    fig100.015.file <- paste(path.fig.stra,outname.stra,"-F100.015.",
                             figtype,sep="")
    if (file.exists(fig100.015.file)) { file.remove(fig100.015.file) }
  } else
  { #  Plot is not wanted
    fig100.015 <- NA
  }


# ============================================================================
# Open.100.016.R
# ============================================================================
  if (plot.fig100.016)
  { #  Plot is wanted
    if (!exists("fig100.016"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.016 <- dev.cur()
    }
    fig100.016.file <- paste(path.fig.stra,outname.stra,"-F100.016.",
                             figtype,sep="")
    if (file.exists(fig100.016.file)) { file.remove(fig100.016.file) }
  } else
  { #  Plot is not wanted
    fig100.016 <- NA
  }


# ============================================================================
# Open.100.017.R
# ============================================================================
  if (plot.fig100.017)
  { #  Plot is wanted
    if (!exists("fig100.017"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig100.017 <- dev.cur()
    }
    fig100.017.file <- paste(path.fig.stra,outname.stra,"-F100.017.",
                             figtype,sep="")
    if (file.exists(fig100.017.file)) { file.remove(fig100.017.file) }
  } else
  { #  Plot is not wanted
    fig100.017 <- NA
  }

# ============================================================================
# Open.108.010.R
# ============================================================================
  if (plot.fig108.010)
  { #  Plot is wanted
    if (!exists("fig108.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig108.010 <- dev.cur()
    }
    fig108.010.file <- paste(path.fig.stra,outname.stra,"-F108.010.",
                             figtype,sep="")
    if (file.exists(fig108.010.file)) { file.remove(fig108.010.file) }
  } else
  { #  Plot is not wanted
    fig108.010 <- NA
  }

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg099_OpenWindows.R  End\n") }
# ============================================================================
