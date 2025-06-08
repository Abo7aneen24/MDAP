#  TMC_seg110_NamesPerStratum.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Segment defining output file names that are stratum-specific.
#  See head of seg050 for the naming convention
#  File-specific names are defined in segment 080
#  Scenario-specific names are defined in segment 211

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  20.07.2020 Sex code and age class added to outfile.stra 
#  24.07.2019 File names changed according to new convention
#             outfile.base -> outname.stra
#             outfile.txt  <- outname.stra.txt
#  12.03.2019 Start from WWW17_seg40  (directory ProgWW17_RH)

# ===========================================================================

# ......................................................................
outname.stra <- paste(outname,"_",sexlabel,"_",agelabel,sep="") 

#  Text file with results per stratum (= separated by sex and age), 
#  no smoothing
outname.stra.txt  <- paste(path.tab.stra,outname.stra, ".txt",sep="")
if (file.exists(outname.stra.txt)) { file.remove(outname.stra.txt) }

#  Histogram data from Fig02, equidistant bins 
outhist.rdata <- paste(path.tab.stra,outname.stra,".RData",sep="")
if (file.exists(outhist.rdata)) { file.remove(outhist.rdata) }

#  Histogram data for estimated pathological data, as in fig 16.1 
outhistpath.rdata <- paste(path.tab.stra,outname.stra,"_path.RData",sep="")
if (file.exists(outhistpath.rdata)) { file.remove(outhistpath.rdata) }

# ......................................................................
