#  TMC_seg211_NamesPerScena.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment defining output file names for csv files that are scenario specific.
#  Files for graphics are defined in the OPEN scripts. 

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  21.04.2021 log file in RData format for simulations added  
#  17.04.2021 csv output files for a scenario new organized. New file names,
#             the *.rep file has only results specific to a replicate,
#             the *.sum file has all control values
#  15.12.2020 Paths for numerical and graphical output more clearly separated
#  13.12.2020 Path structure for scenario changed
#  07.07.2019 Start from TMC_seg210_NamesPerScena.R
#             Files are reorganized acccording to TMC_seg201_SimResults.R

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg211_NamesPerScena  Start\n") }
# ===========================================================================

# ---------------------------------------------------------------------------
#  NUMERICAL results

#  Files containing results from all methods, either per replicate or 
#  summarizing over replicates, go into folder 
#  tabs/3_scen/<RunId>/

path.tab.scen.Id <- paste(path.tab.scen, RunId, "/", sep="")
if (!dir.exists(path.tab.scen.Id))
{ 
  # Create the folder, if it does not exist
  dir.ok <- dir.create(path.tab.scen.Id, showWarnings = TRUE, 
                       recursive = FALSE, mode = "0777")
}

#  Results produced by the methods in TMC for each file (replicate) of a 
#  scenario
soutfile.tmc.csv <- paste(path.tab.scen.Id, datafile0,"-S-TMC.csv", sep="")
if (file.exists(soutfile.tmc.csv)) { file.remove(soutfile.tmc.csv) }

#  Results obtained by preceding call of RLe / TML, per replicate
#  now stored in a subfolder of RLe
#  Example for TML summary result file: 106_Na_010000_r-S-TML
#  datafile0 is "106_Na_010000_r" = data.file.base in TML

# soutfile.tml.csv <- paste(path.tab.scen.Id, datafile0,"-S-TML.csv", sep="")
soutfile.tml.csv <- paste("../../RLE/RLs1/summary/", 
                           datafile0,"-S-TML.csv", sep="")

#  Results obtained by all applied methods, per replicate. Joins the preceding
#  two files.
#soutfile.all.csv <- paste(path.tab.scen.Id, datafile0,"-S-ALL.csv", sep="")
#  renamed 17.04.2021
outname.scen.rep.csv <- paste(path.tab.scen.Id, datafile0,"-scen-rep.csv", 
                              sep="")

#  Summary over all replicates, by method.
# zoutfile.csv <- paste(path.tab.scen.Id, datafile0, "-Z.csv", sep="")
#  renamed 17.04.2021
outname.scen.sum.csv <- paste(path.tab.scen.Id, datafile0, "-scen-sum.csv",
                              sep="")

#  Method assessment
outname.scen.asm.csv <- paste(path.tab.scen.Id, datafile0, "-scen-asm.csv", 
                              sep="")

#  Log file containing (nearly) all parameters
outname.scen.log.rda <- paste(path.tab.scen.Id, datafile0, "-scen-log.RData", 
                              sep="")

# Do not delete!
# if (file.exists(zoutfile.csv)) { file.remove(zoutfile.csv) }

# ---------------------------------------------------------------------------
#  GRAPHICAL results

#  Files containing results from all methods go into folder 
#  figs/3_scen/<RunId>/

path.fig.scen.Id <- paste(path.fig.scen, RunId, "/", sep="")
if (!dir.exists(path.fig.scen.Id))
{ 
  # Create the folder, if it does not exist
  dir.ok <- dir.create(path.fig.scen.Id, showWarnings = TRUE, 
                       recursive = FALSE, mode = "0777")
}

#  Full file names are constructed in the Open.xxx.yyy files 

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg211_NamesPerScena  End\n") }
# ===========================================================================

