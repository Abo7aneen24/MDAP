# ==============================================================================
#  TMC ANALYSIS START FILE
# ==============================================================================
#
# PURPOSE:  This script runs the Truncated Minimum Chi-Square (TMC) analysis
#           on the prepared dataset to calculate age- and sex-specific
#           reference intervals for TSH and Free T4.
#
# AUTHOR:   Abdullah M. Algarni, MBBS (adapted from TMC template)
# DATE:     June 8, 2025
#
# ==============================================================================
#  SECTION 1: USER-DEFINED PARAMETERS
# ==============================================================================
#
# --- File to Analyze ---
# This number points to the row in your DataFileInfo.csv file.
# Set to 1 for TSH analysis.
# Set to 2 for Free T4 analysis.
FileNo <- 1

# --- Define Age Groups ---
# These are the lower limits of each age bin.
# e.g., c(18, 30, 40) creates intervals: 18-29, 30-39, and 40+.
# We will use 10-year intervals as discussed.
age.limits <- c(18, 30, 40, 50, 60, 70, 80, 90)

# --- Define Reference Interval Percentiles ---
# These are the standard values for a 95% reference interval.
RL1.p <- 0.025 # Lower limit (2.5th percentile)
RL2.p <- 0.975 # Upper limit (97.5th percentile)

# --- Optional Data Filtering (Not needed for this analysis) ---
# Multiplies all values by this factor. Set to 1 for no change.
scale.fact <- 1
# Exclude values outside this range before analysis. Set to NA to disable.
x.lo.limit <- NA
x.hi.limit <- NA
# Select a subset based on a specific column (e.g., outpatient status). Set to NA to disable.
Use.oh <- NA
# Select a subset based on the lab device. Set to NA to disable.
Use.dev <- NA

# --- Optional Plotting Adjustments (Not needed initially) ---
# You can adjust plot axes here if needed later. Set to NA to use defaults.
x.clip.min <- NA
x.clip.max <- NA
x.clip.by1 <- NA
x.clip.by2 <- NA

age.clip.min <- NA
age.clip.max <- NA
age.clip.by1 <- NA
age.clip.by2 <- NA

# ==============================================================================
#  SECTION 2: RUN THE TMC PROGRAM (DO NOT MODIFY BELOW THIS LINE)
# ==============================================================================
# This section loads the TMC program files and executes the analysis
# based on the parameters you set above.

#### Do not change the following lines ####
source("TMC_seg005_Initialise.R")
source("TMC_seg010_User.R")
source("TMC_seg015_DefaultSettings.R")

#### Put changes of default settings after this line ####

#### Do not change the following lines ####
source("TMC_seg020_Libraries.R")
source("TMC_seg025_Functions.R")
source("TMC_seg030_Style.R")
source("TMC_seg035_ReadData.R")
source("TMC_seg040_PrepareData.R")
source("TMC_seg045_PlotRequest.R")
source("TMC_seg050_BasicPlots.R")
source("TMC_seg055_Stratify.R")
source("TMC_seg060_FinalPlots.R")
source("TMC_seg065_FinalTables.R")
# ==============================================================================
#  END OF SCRIPT
# ==============================================================================