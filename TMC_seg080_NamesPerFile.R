#  TMC_seg080_NamesPerFile.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment defining output file names that are file specific.
#  Stratum-specific files are defined in the segment 110

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg080_NamesPerFile  Start\n") }
# ===========================================================================

# ===========================================================================

#  Base name for results per file 
outname.file  <- outname

#  Ergebnis für alle Strata eines Datensatzes (gtab), 
#  txt-Datei mit Zeilenumbruch
outname.file.txt <- paste(path.tab.file,outname,"-g.txt",sep="")
if (file.exists(outname.file.txt)) { file.remove(outname.file.txt) }

#  Ergebnis für alle Strata eines Datensatzes (gtab), 
#  csv-Datei ohne Zeilenumbruch
outname.file.csv <- paste(path.tab.file,outname,"-g.csv",sep="")
if (file.exists(outname.file.csv)) { file.remove(outname.file.csv) }

#  Ausgabedatei RL1, RL2, empirisch, spaltenweise getrennt nach sex,
RLempfile <- paste(path.tab.file,outname,"_RLemp.csv",sep="")
if (file.exists(RLempfile)) { file.remove(RLempfile) }

#  Ausgabedatei RL1, RL2, geglättet, spaltenweise getrennt nach sex,
RLsmofile <- paste(path.tab.file,outname,"_RLsmo.csv",sep="")
if (file.exists(RLsmofile)) { file.remove(RLsmofile) }

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg080_NamesPerFile  End\n") }
# ===========================================================================

