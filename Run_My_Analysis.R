# ==============================================================================
#  TMC ANALYSIS START FILE (FINAL CORRECTED VERSION)
# ==============================================================================
#
# PURPOSE:  This script runs the Truncated Minimum Chi-Square (TMC) analysis.
#           All configuration issues have been resolved.
#
# AUTHOR:   Abdullah M. Algarni, MBBS (adapted from TMC template)
# DATE:     June 8, 2025
#
# ==============================================================================
#  SECTION 1: LOAD TMC SCRIPTS FIRST
# ==============================================================================

# Load user settings and default program settings FIRST
# (These may clear the environment, so we load them before setting our variables)
source("TMC_seg010_User.R")
source("TMC_seg015_DefaultSettings.R")

# ==============================================================================
#  SECTION 2: SET PARAMETERS (AFTER TMC SCRIPTS ARE LOADED)
# ==============================================================================

# --- File to Analyze (from your DataFileInfo.csv) ---
# Set to 1 for TSH analysis.
# Set to 2 for Free T4 analysis.
FileNo <- 1

# Verify FileNo is set
cat("FileNo is set to:", FileNo, "\n")

# --- Define Age Groups ---
# These are the lower limits of each age bin (e.g., 18-29, 30-39, etc.)
age.limits <- c(18, 30, 40, 50, 60, 70, 80, 90)

# --- Define Reference Interval Percentiles ---
# Standard values for a 95% reference interval.
RL1.p <- 0.025
RL2.p <- 0.975

# --- Optional Data Filtering ---
scale.fact <- 1
x.lo.limit <- NA
x.hi.limit <- NA
Use.oh <- NA
Use.dev <- NA
ana.hour.min <- NA; ana.minu.min <- NA; ana.hour.max <- NA; ana.minu.max <- NA
use.wday <- "all"
split.date <- NA; split.time <- NA; use.split <- NA

# --- Optional Plotting Adjustments ---
x.clip.min <- NA; x.clip.max <- NA; x.clip.by1 <- NA; x.clip.by2 <- NA
age.clip.min <- NA; age.clip.max <- NA; age.clip.by1 <- NA; age.clip.by2 <- NA
RL.clip.min <- NA; RL.clip.max <- NA; RL.clip.by1 <- NA; RL.clip.by2 <- NA

# --- Simulation Parameters (Not applicable for this analysis) ---
r.start <- NA
r.ende <- NA

# ==============================================================================
#  SECTION 3: OVERRIDE DEFAULT SETTINGS
# ==============================================================================

# Put changes of default settings after this line, if any are ever needed
# This line will override the default setting from TMC_seg010_User.R
# to point to your specific data file information.
info.file <- "../data/DataFileInfo.csv"

# ==============================================================================
#  SECTION 4: RUN THE TMC ANALYSIS
# ==============================================================================

# Process settings and run the main analysis
source("TMC_seg040_ProcessSettings.R")
source("TMC_seg050_Master.R")

# ==============================================================================
#  END OF SCRIPT
# ==============================================================================