#  TMC_seg105_Analysis_tmc.R
#
#  Do tmc iteration for one TI, calculate solution assessment criteria
#  Return tab.tmc containing iteration result including 
#  opt.crit, p.fit, sol.score
#
#  13.08.2022 No decision about acceptability, moved to seg100, only
#             calculation of criteria 
#  01.07.2022
#  07.04.2021 crit1 not used if values < DL exist 
#  15.02.2021 Start
# ==============================================================================

  # if (print.log.message) { cat("%%%   TMC_seg105_Analysis  100\n") }

  seg105.t1 <- Sys.time()

  #  tmc.master does grid search around theta.ini0 and subsequent 
  #  nlm optimization, for one TI

  theta.ini0 <-  c(tab.qqw[1, "lambda.qqw"], tab.qqw[1, "mue.qqw"], 
                   tab.qqw[1, "sigma.qqw"])
  names(theta.ini0) <- c("lambda", "mue", "sigma")

  tab.tmc0 <- tmc.master(x, x.hist, TI.can[iti, "ilo"], TI.can[iti, "ihi"],
                         TI.can[iti, "prop"], x.tr.n, 
                         TI.can[iti, "x.lt.tr.n"], TI.can[iti, "x.ge.tr.n"],
                         theta.ini0, 
                         theta.fix, theta.est, idx.fix, idx.est,
                         lambda.min, lambda.max, 
                         df.est, df.con,
                         l.fact, p.fact, r.fact, rt.fact, w.fact, 
                         x.Q1, x.Q2, RL1.p, RL2.p, fastnull, fastnull.chi,                     
                         print.log.message,
                         unlist(tab.npa["x.RL1.npa"]), 
                         unlist(tab.npa["x.RL2.npa"]), 
                         xsupp, gencol, kdecol, tmccol, 
                         NA, fig100.015)    # fig100.012  to show contour plot  

  # if (print.log.message) { cat("%%%   TMC_seg105_Analysis  150\n") }

  #  Calculate assessment criteria
  #  Estimated RLs in [x.Q1, x.Q2]?
  #  Consider crit 1 only if no values < DL, because x.Q1 is
  #  very unreliable if values < DL exist
  if (is.na(detect.limits.max))
  { 
    crit[1] <- unname(x.Q1-round.unit <= tab.tmc0["x.RL1.tmc"]) # 01.03.2022
  } else
  { crit[1] <- TRUE } 

  #  Test 28.02.2022
  crit[1] <- TRUE

  crit[2] <- unname(x.Q2+round.unit >= tab.tmc0["x.RL2.tmc"])  # 01.03.2022

  #  Estimated prevalences in [prev.acc.lo, prev.acc.hi]?
  crit[3] <- unname((prev.acc.lo <= tab.tmc0["prev.l.tmc"])) &
                   unname((tab.tmc0["prev.l.tmc"] <= prev.acc.hi))  
  crit[4] <- unname((prev.acc.lo <= tab.tmc0["prev.c.tmc"])) &
                   unname((tab.tmc0["prev.c.tmc"] <= prev.acc.hi))  
  crit[5] <- unname((prev.acc.lo <= tab.tmc0["prev.r.tmc"])) &
                   unname((tab.tmc0["prev.r.tmc"] <= prev.acc.hi))

  #  Basic solution score
  sol.score <- sum(crit[1:5])

  #  Former crit 6 now contained as penalty term in opt.crit (16.08.2022)
  
  tab.tmc0["crit1"] <- crit[1]
  tab.tmc0["crit2"] <- crit[2]
  tab.tmc0["crit3"] <- crit[3]
  tab.tmc0["crit4"] <- crit[4]
  tab.tmc0["crit5"] <- crit[5]

  tab.tmc0["sol.score"] <- sol.score
        
  ## # if (print.log.message) { cat("%%%   TMC_seg105_Analysis 200\n") }

  #  Result from this segment:
  #  tab.tmc         contains result
  #  temp.tmc        plot information

  # # if (print.log.message) { cat("%%%   TMC_seg105_Analysis  900\n") }

  seg105.t2 <- Sys.time()
  # cat("\n seg105 Gesamtzeit", format(difftime(seg105.t2, seg105.t1)), "\n")

##############################################################################