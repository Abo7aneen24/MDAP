#  TMC_seg045_PlotRequest.R
#  (c) wwosniok@math.uni-bremen.de

#  Define the plots to produce
#  Graph windows are generated (later) only for requested plots

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ============================================================================

# ============================================================================
if (print.log.message) { cat("%%%   TMC_seg045_PlotRequest.R  Start\n") }
# ============================================================================
# ----------------------------------------------------------------------------
#  Automatic selection of plots to produce depends on variable plot.level 
#  and nrep

# ---------------------------------------------------------------------------
if (plot.level == "all")
{ # Plot all

plot.fig050.010 <- TRUE  #  Drift plot
plot.fig090.010 <- TRUE  #  Countour plot measurand vs age, males
plot.fig090.020 <- TRUE  #  Countour plot measurand vs age, females
plot.fig095.010 <- TRUE  #  Plot of measurand vs date-time in week, means
plot.fig095.020 <- plot.fig095.010  
                         #  Plot of measurand vs date-time in week, medians
plot.fig095.030 <- TRUE  #  Plot of measurand vs time of day, means
plot.fig095.040 <- plot.fig095.030  
                         #  Plot of measurand vs time of day, medians
plot.fig100.002 <- TRUE  #  Show kde bandwidth determination
plot.fig100.005 <- TRUE  #  Compare empirical distribution with kde
plot.fig100.010 <- TRUE  #  frequency of each value 
plot.fig100.020 <- TRUE  #  Histogram, full linear scale, breakpoints: TMC
plot.fig100.030 <- TRUE  #  Histogram, clipped linear scale, bin size (min/eff.) =
plot.fig100.040 <- TRUE  #  QQ plot for determination of initial values
plot.fig100.050 <- TRUE  #  QQ plot for TMC result
plot.fig100.125 <- TRUE  #  Estimated RLs vs prop.tmc when bootstrapping
plot.fig100.130 <- TRUE  #  Convergence of RL estimates when bootstrapping
plot.fig100.140 <- TRUE  #  Bootstrapped histograms and densities
plot.fig100.145 <- TRUE  #  Densities of bootstrapped results
plot.fig100.160 <- TRUE  #  Empirical and fitted cdfs
plot.fig102.010 <- TRUE  #  Proportion in TI vs optimality and RLs
plot.fig108.010 <- TRUE  #  Tukey plot

plot.fig120.010 <- TRUE  #  TMC result, linear, clipped x axis, delta hist, ex 16
plot.fig130.010 <- TRUE  #  TMC result, linear, clipped x axis, delta hist, ex 16
plot.fig135.010 <- TRUE  #  TMC result, linear, clipped x axis, delta hist, ex 16
plot.fig140.010 <- TRUE  #  Plot RL vs age class, smoothed, by M / F / M+F, ex 18

plot.fig150.010 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , npa
plot.fig150.020 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , qqw
plot.fig150.030 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , bha
plot.fig150.040 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , tuk
plot.fig150.050 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , tmc
plot.fig150.060 <- TRUE  #  Plot RL vs age class, smoothed, by M / F , tmu

plot.fig160.010 <- TRUE  #  Plot RL vs age class, smoothed, by M + F , tmc
plot.fig202.010 <- TRUE  #  RL1 vs RL2 from simulation
plot.fig202.011 <- TRUE  #  RL1 vs RL2 from simulation, contour background 1
plot.fig202.012 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 1
plot.fig202.013 <- TRUE  #  RL1 vs RL2 from simulation, contour background 2
plot.fig202.014 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 2
plot.fig202.015 <- TRUE  #  RL1 vs RL2 from simulation, contour background 3
plot.fig202.016 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 3
plot.fig202.020 <- TRUE  #  Relative differences between RLx, RL.npa
plot.fig202.030 <- TRUE  #  FPRx from  simulation
plot.fig202.040 <- TRUE  #  Density of FPR1+FPR2
plot.fig202.050 <- TRUE  #  Density of delta.FPR

plot.fig203.010 <- TRUE  # SimResult
plot.fig203.020 <- TRUE  # SimResult

}  #  plot.level == "all"

# ---------------------------------------------------------------------------
if (plot.level == "auto" & nrep == 1)
{ # Single file, standard plots

plot.fig050.010 <- TRUE    # 2 Drift plot
plot.fig090.010 <- TRUE    # 3 Contour plot measurand vs age, females
plot.fig090.020 <- TRUE    # 4 Contour plot measurand vs age, males
plot.fig095.010 <- TRUE    # 5 Plot of measurand vs date-time in week, all data,
                           #  means
plot.fig095.020 <- plot.fig095.010  
                           # 6 Plot of measurand vs date-time in week, all data,
                           #  medians
plot.fig095.030 <- TRUE    # 7 Plot of measurand vs time of day, Mo-Fr, means
plot.fig095.040 <- plot.fig095.030  
                           # 8 Plot of measurand vs time of day, Mo-Fr, medians
plot.fig100.001 <- FALSE   #  FindKDE
plot.fig100.002 <- FALSE   #  FindHisto
plot.fig100.003 <- FALSE   #  FindHisto
plot.fig100.004 <- TRUE    # 9 Frequencies of observed values, Rounding check
plot.fig100.005 <- FALSE   #  QQW, collapsed and reduced histo
plot.fig100.006 <- FALSE   #  QQW
plot.fig100.007 <- FALSE   #  QQW  optimal QQW solutions   @@@ leer
plot.fig100.008 <- FALSE   #  QQW
plot.fig100.009 <- FALSE   #  QQW
plot.fig100.010 <- FALSE   #  QQW
plot.fig100.011 <- FALSE   #  QQW                          @@@ leer
plot.fig100.012 <- FALSE   #  tmc.master, contour plot     @@@ leer
plot.fig100.013 <- FALSE   #  Bhattacharya
plot.fig100.014 <- FALSE   #  Bhattacharya
plot.fig100.015 <- FALSE   #  TMC, from function
plot.fig100.016 <- TRUE    # 10 TMC, from seg100
plot.fig100.017 <- FALSE   #  TMU, from seg100
plot.fig108.010 <- FALSE   #  Tukey plot

plot.fig150.010 <- ("npa" %in% meth.list)
plot.fig150.020 <- ("qqw" %in% meth.list)
plot.fig150.030 <- ("bha" %in% meth.list)
plot.fig150.040 <- ("tuk" %in% meth.list)
plot.fig150.050 <- ("tmc" %in% meth.list) & (!plot.perm.diff.band)
plot.fig150.051 <- ("tmc" %in% meth.list) & (plot.perm.diff.band)  #  11
plot.fig150.060 <- FALSE

plot.fig160.010 <- FALSE  #  Plot RL vs age class, smoothed, by M + F
plot.fig160.020 <- FALSE  #  Plot RL vs age class, smoothed, by M + F
plot.fig160.030 <- FALSE  #  Plot RL vs age class, smoothed, by M + F
plot.fig160.040 <- FALSE  #  Plot RL vs age class, smoothed, by M + F
plot.fig160.050 <- ("tmc" %in% meth.list) & (!plot.perm.diff.band)
plot.fig160.051 <- ("tmc" %in% meth.list) & (plot.perm.diff.band)  #  12
plot.fig160.060 <- FALSE  #  Plot RL vs age class, smoothed, by M + F

plot.fig202.010 <- FALSE  #  RL1 vs RL2 from simulation
plot.fig202.011 <- FALSE  #  RL1 vs RL2 from simulation, contour background 1
plot.fig202.012 <- FALSE  #  RL1 vs RL2 from simulation, contour legend 1
plot.fig202.013 <- FALSE  #  RL1 vs RL2 from simulation, contour background 2
plot.fig202.014 <- FALSE  #  RL1 vs RL2 from simulation, contour legend 2
plot.fig202.015 <- FALSE  #  RL1 vs RL2 from simulation, contour background 3
plot.fig202.016 <- FALSE  #  RL1 vs RL2 from simulation, contour legend 3
plot.fig202.020 <- FALSE  #  Relative differences between RLx, RL.npa
plot.fig202.030 <- FALSE  #  FPRx from  simulation
plot.fig202.040 <- FALSE  #  Density of FPR1+FPR2
plot.fig202.050 <- FALSE  # Density of delta.FPR

plot.fig203.010 <- FALSE  # SimResult
plot.fig203.020 <- FALSE  # SimResult


}  #  plot.level == "auto" & nrep == 1


# ---------------------------------------------------------------------------
if (plot.level == "auto" & nrep > 1)
{ # Analysis of a file squence (simulation data)

plot.fig050.010 <- FALSE   #  Drift plot
plot.fig090.010 <- FALSE   #  Contour plot measurand vs age, males
plot.fig090.020 <- FALSE   #  Contour plot measurand vs age, females
plot.fig095.010 <- FALSE   #  Plot of measurand vs date-time in week, means
plot.fig095.020 <- plot.fig095.010  
                           #  Plot of measurand vs date-time in week, medians
plot.fig095.030 <- FALSE   #  Plot of measurand vs time of day, means
plot.fig095.040 <- plot.fig095.030  
                           #  Plot of measurand vs time of day, medians
plot.fig100.001 <- FALSE   #  FindKDE
plot.fig100.002 <- FALSE   #  FindHisto
plot.fig100.003 <- FALSE   #  FindHisto
plot.fig100.004 <- FALSE    #  Rounding problem
plot.fig100.005 <- FALSE   #  QQW, collapsed and reduced histo
plot.fig100.006 <- FALSE   #  QQW
plot.fig100.007 <- FALSE   #  QQW
plot.fig100.008 <- FALSE   #  QQW
plot.fig100.009 <- FALSE   #  QQW
plot.fig100.010 <- FALSE   #  QQW
plot.fig100.011 <- FALSE   #  QQW
plot.fig100.012 <- FALSE    #  QQW
plot.fig100.013 <- FALSE   #  Bhattacharya
plot.fig100.014 <- FALSE   #  Bhattacharya
plot.fig100.015 <- FALSE   #  TMC, from function
plot.fig100.016 <- TRUE   #  TMC, from seg100
plot.fig100.017 <- FALSE   #  TMU, from seg100
plot.fig108.010 <- TRUE  #  Tukey plot

plot.fig150.010 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , npa
plot.fig150.020 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , qqw
plot.fig150.030 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , bha
plot.fig150.040 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , tuk
plot.fig150.050 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , tmc
plot.fig150.060 <- FALSE  #  Plot RL vs age class, smoothed, by M / F , tmu

plot.fig160.010 <- FALSE  #  Plot RL vs age class, smoothed, by M + F , ex 20
plot.fig202.010 <- TRUE  #  RL1 vs RL2 from simulation
plot.fig202.011 <- TRUE  #  RL1 vs RL2 from simulation, contour background 1
plot.fig202.012 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 1
plot.fig202.013 <- TRUE  #  RL1 vs RL2 from simulation, contour background 2
plot.fig202.014 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 2
plot.fig202.015 <- TRUE  #  RL1 vs RL2 from simulation, contour background 3
plot.fig202.016 <- TRUE  #  RL1 vs RL2 from simulation, contour legend 3
plot.fig202.020 <- TRUE  #  Relative differences between RLx, RL.npa
plot.fig202.030 <- TRUE  #  FPRx from  simulation
plot.fig202.040 <- TRUE  #  Density of FPR1+FPR2
plot.fig202.050 <- TRUE  #  Density of delta.FPR

plot.fig203.010 <- TRUE  # SimResult
plot.fig203.020 <- TRUE  # SimResult


}  #  plot.level == "auto" & nrep> 1


# ----------------------------------------------------------------------------
#  This is the command sequence to set up graph window and plot file
#  Care for datafile / outfile.base and the directory!
#
#  if (plot.fig100.010)
#  { if (!exists("fig100.010"))
#    { #  Open graph new window
#      win.graph()
#      fig100.010 <- dev.cur()
#    }
#    fig100.010.file <- paste(path.fig.stra,outfile.base,"-F100.010",
#                             figtype,sep="")
#    if (file.exists(fig100.010.file)) { file.remove(fig100.010.file) }
# 
#    dev.set(fig100.010)
#
#   @@@@@@@@ plot
#
#
#    savePlot(file=fig100.010.file,type=figtype)
#  }

# =============================================================================
if (print.log.message) { cat("%%%   TMC_seg045_PlotRequest.R  End\n") }
