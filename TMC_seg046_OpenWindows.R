#  TMC_seg046_OpenWindows.R
#  (c) wwosniok@math.uni-bremen.de

#  Open graph windows and define plot file names for plots requested in
#  TMC_seg045_PlotRequest.R

#  05.06.2021 Graphics devices changed to X11()
#  23.04.2021 New organisation: open only requested windows on file level.
#             Windows for involving date or time of day are opened in the
#             segments that do these plots  
#  17.12.2020 Start
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg046_OpenWindows.R  Start\n") }
# ============================================================================

# ============================================================================
# Open.050.010.R
# ==============================================================================
  if (plot.fig050.010)
  { #  Plot is wanted
    if (!exists("fig050.010"))
    { #  Open graph new window
      X11()      #  win.graph()
      fig050.010 <- dev.cur()
    }
    fig050.010.file <- paste(path.fig.file,outname,"-F050.010.",
                             figtype,sep="")
    if (file.exists(fig050.010.file)) { file.remove(fig050.010.file) }
  } else
  { #  Plot is not wanted
    fig050.010 <- NA
  }

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg046_OpenWindows.R  End\n") }
# ============================================================================

