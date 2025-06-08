#  TMC_seg159_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  05.06.2021 Graphics devices changed to X11()
#  03.01.2021 Target directory corrected
#  17.12.2020 Start
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg159_OpenWindows.R  Start\n") }
# ============================================================================


# ============================================================================
# Open.160.050.R
# ============================================================================
  if (plot.fig160.050)
  { #  Plot is wanted
    if (!exists("fig160.050"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig160.050 <- dev.cur()
    }
    fig160.050.file <- paste(path.fig.file,outname.file,"-F160.050.",
                             figtype,sep="")
    if (file.exists(fig160.050.file)) { file.remove(fig160.050.file) }
  } else
  { #  Plot is not wanted
    fig160.050 <- NA
  }

# ============================================================================
# Open.160.051.R
# ============================================================================
  if (plot.fig160.051)
  { #  Plot is wanted
    if (!exists("fig160.051"))
    { #  Open graph new window
      X11()    #  win.graph()
      fig160.051 <- dev.cur()
    }
    fig160.051.file <- paste(path.fig.file,outname.file,"-F160.051.",
                             figtype,sep="")
    if (file.exists(fig160.051.file)) { file.remove(fig160.051.file) }
  } else
  { #  Plot is not wanted
    fig160.051 <- NA
  }


# ============================================================================
# ============================================================================
