# TMC_seg102_Analysis_npa.R
#

#  11.10.2022 New variable system: all data per method in tabn.npa without
#             suffix for method 
#  18.10.2021 Treatment of no L or no R corrected
#  15.02.2021
# ==============================================================================

    if (print.log.message) { cat("%%%   TMC_seg102_Analysis_npa Start\n") }

    # xl.n.npa  <- sum(grp=="L")     # also under perfect identification unknown
    tabn.npa["xc.n"] <- sum(grp=="C") # this number would be known under perfect
                                 # identification of non-diseased subject
    # xr.n.npa  <- sum(grp=="R") # also under perfect identification unknown
    tabn.npa["prev.c"] <- tabn.npa["xc.n"] / tabn.npa["n"]

    #  Calculation of prev.l.npa and prev.r.npa here is different from
    #  the other methods, because prev.c.n is exactly known
    tabn.npa["xl.n"]   <- sum(x < x.kde.mode) - 
                          sum((x < x.kde.mode) & (grp=="C")) 
    tabn.npa["prev.l"] <- tabn.npa["xl.n"] / tabn.npa["n"]
    tabn.npa["xr.n"]   <- sum(x >= x.kde.mode) - 
                          sum((x >= x.kde.mode) & (grp=="C")) 
    tabn.npa["prev.r"] <- tabn.npa["xr.n"] / tabn.npa["n"]

    #  Optimal truncation interval (left limits are contained, right not)
    tabn.npa["x.tr.lo"] <- min(x[grp=="C"]) - round.unit/2 
    tabn.npa["x.tr.hi"] <- max(x[grp=="C"]) + round.unit/2

    tabn.npa["x.tr.n"]    <- length(x[(tabn.npa["x.tr.lo"] <= x) & 
                                      (x < tabn.npa["x.tr.hi"])])
    tabn.npa["x.tr.prop"] <- tabn.npa["x.tr.n"] / tabn.npa["n"]  

    tabn.npa["x.RL1"] <- Quantile(x[grp=="C"], probs=RL1.p)
    tabn.npa["x.RL2"] <- Quantile(x[grp=="C"], probs=RL2.p)

    if (print.log.message) { cat("%%%   TMC_seg102_Analysis_npa End\n") }
