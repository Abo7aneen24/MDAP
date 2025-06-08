#  TMC_seg040_ProcessSettings.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Process actual parameters (default or user defined)

#  TODO       
#  -  

#  #  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

if (print.log.message) { cat("%%%   TMC_seg040_ProcessSettings.R  Start\n") }

# ==========================================================================

# --------------------------------------------------------------------------
#  Check for reasonable smoothing options

if (smooth.hist1 & smooth.hist2)
{ if (print.log.message) { cat("%%%   TMC_seg040_ProcessSettings 100\n") }
  # Should not happen
    cat("\n\n ++++++  Please decide for either smooth.hist1 OR smooth.hist2",
        "\n ++++++  EXECUTION STOPPED")
    stop("++++++  EXECUTION STOPPED")
    cat("\n\n")
    sink()
} 

# --------------------------------------------------------------------------
#  Process options for lambda

if (lambda.min > lambda.max)
{
  cat("\n ++++++ ",
      "\n ++++++ lambda.min > lambda.max  --- execution stops",
      "\n ++++++ ",
      "\n\n")
  stop("++++++ Execution stops because invalid lambda interval ++++++")
}

#  lambda interval is ok

#  Reduce lamba.seq according to lambda.min, lambda.max

lambda.seq <- lambda.seq[(lambda.min <= lambda.seq) & 
                         (lambda.seq <= lambda.max) ]

lambda.seq.n <- length(lambda.seq)

# ----------------------------------------------------------------------------

#  Treatment of parameters during estimation:
#  theta    is the full parameter
#  test.fix is the fixed component. May be NULL.
#  test.est is the component that is estimated
#  So far, only the possibility of fixing lambda is realized!
 
#  Functions used during estimation:
#  theta.est <- Extract(theta, idx.est)  
#  theta     <- Compose(theta.fix, theta.est, idx.fix, idx.est)

#  Indicate the parameters to estimate
#  Default: all components are estimated
idx.fix <- NULL
idx.est <- c(1, 2, 3)

if (abs(lambda.min - lambda.max) < fastnull)
{
  #  lambda is set to the fixed value lambda.min
  lambda.seq       <- lambda.min

  idx.fix <- 1
  idx.est <- c(2, 3)
}

# --------------------------------------------------------------------------
#  Generate sequence of truncation interval limits (used in Proportions())
#  Unused??
#x.tr.prop.limits   <- seq(x.tr.prop.min, x.tr.prop.max, length.out=5)
#x.tr.prop.ints.n   <- length(x.tr.prop.limits) - 1

#  Increase the highest limit by a small amount such that the query 
#  lower limit <= x.tr.prop < upper limit captures the highest limit
#x.tr.prop.limits[x.tr.prop.ints.n+1] <- 
#                      x.tr.prop.limits[x.tr.prop.ints.n+1] + 0.01
                      

#  Formatted version of l.fact, for appending in file names
l.fact.fmt <- formatC(l.fact, width=4, digits=2, flag="#") 

#  Formatted version of p.fact, for appending in file names
p.fact.fmt <- formatC(p.fact, width=4, digits=2, flag="#") 

#  Formatted version of w.fact, for appending in file names
w.fact.fmt <- formatC(w.fact, width=4, digits=2, flag="#") 

#  Number of replicate data sets to analyze 
nrep     <- r.ende - r.start + 1

#  Is this an evaluation of replicates?
eval.rep <- !is.na(nrep)

# ===========================================================================

if (print.log.message) { cat("%%%   TMC_seg040_ProcessSettings.R  End\n") }

# ==========================================================================
