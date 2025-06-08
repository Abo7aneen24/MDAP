#  TMC_seg089_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg089_OpenWindows.R  Start\n") }
# ============================================================================


# ============================================================================
# Open.090.010.R
# ==============================================================================
  if (plot.fig090.010)
  { #  Plot is wanted
    if (!exists("fig090.010"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig090.010 <- dev.cur()
    }
    fig090.010.file <- paste(path.fig.file,outname.file,"-F090.010.",
                             figtype,sep="")
    if (file.exists(fig090.010.file)) { file.remove(fig090.010.file) }
  } else
  { #  Plot is not wanted
    fig090.010 <- NA
  }


# ============================================================================
# Open.090.020.R
# ==============================================================================
  if (plot.fig090.020)
  { #  Plot is wanted
    if (!exists("fig090.020"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig090.020 <- dev.cur()
    }
    fig090.020.file <- paste(path.fig.file,outname.file,"-F090.020.",
                             figtype,sep="")
    if (file.exists(fig090.020.file)) { file.remove(fig090.020.file) }
  } else
  { #  Plot is not wanted
    fig090.020 <- NA
  }

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg089_OpenWindows.R  End\n") }
# ============================================================================

