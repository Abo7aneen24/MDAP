#  TMC_seg106_Analysis_tmu.R
#
#  02.05.2021  Plot of residual histogram instead of residual density
#  15.02.2021
# ==============================================================================

y.tr <- BoxCox(x[(tab.tmc["x.tr.lo"] <= x) &
                 (x < tab.tmc["x.tr.hi"]) ],
               tab.tmc["lambda.tmc"] )
y.tr.lo <- BoxCox(tab.tmc["x.tr.lo"], 
                  tab.tmc["lambda.tmc"] )
y.tr.hi <- BoxCox(tab.tmc["x.tr.hi"], 
                  tab.tmc["lambda.tmc"] )

tab.tmu <- EstParentDist(y.tr.lo, y.tr.hi, 
                         y.tr, c(tab.tmc["mue.tmc"], tab.tmc["sigma.tmc"]),
                         alpha, tab.tmc["lambda.tmc"], RL1.p, RL2.p)

temp.tmu <- chi2trunc(x.hist,tab.tmc["x.tr.lo"],
                      tab.tmc["x.tr.hi"],
                      sum(x <  tab.tmc["x.tr.lo"]), 
                      sum(x >= tab.tmc["x.tr.hi"]),
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      tab.tmc["lambda.tmc"],
                      tab.tmu["mue.tmu"],
                      tab.tmu["sigma.tmu"],
                      df.est, df.con,
                      l.fact, p.fact, r.fact, w.fact, 
                      opt.crit.only=FALSE,fastnull=fastnull)
tab.tmu <- c(tab.tmu, 
             opt.crit.tmu=unname(temp.tmu[["res"]]["opt.crit"]),
             p.fit.tmu=unname(temp.tmu[["res"]]["chi2.trun.p"]),
             prev.l.tmu=unname(temp.tmu[["res"]]["prev.l.tmc"]),
             prev.c.tmu=unname(temp.tmu[["res"]]["prev.c.tmc"]),
             prev.r.tmu=unname(temp.tmu[["res"]]["prev.r.tmc"]) )

# ------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg100_Analysis 2000\n") }
  
if (plot.fig100.017)
{ 
  #  Plot result - basic part 
  #  No plot file name, save is done here below
  fig100.017.lim <- PlotHistFit(fig100.017, NA, figtype, 
                       x.hist, xlabel, "TMU result", subtitle, n.per.bin.min,
                       NA, x.val, counts.range, 
                       paste(formatC(100*tab.tmc["x.tr.prop"], format="f", width=5, 
                             digits=1), "%", sep=""),
                       x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                       x.clip.type, NA, par.las, par.tcl,
                       bordercol, histcol, difcol, kdecol, denlty)

  #  Plot result - tmu specific part 
  #  Add estimated tmu solution, RLs, TI
  #  Calculate estimated TMC density (not yet scaled by prev.c.tmu)
  xc.pdf.tmu <- pdf.PN(xsupp, tab.tmc["lambda.tmc"], 
                              tab.tmu["mue.tmu"],
                              tab.tmu["sigma.tmu"]) 

  #  Calculate difference between x.kde and xc.pdf.tmu
  #  Truncate difference below to yplotmin
  #  Switched off 01.05.2021
  # x.kde.tmu <- x.kde$y -  tab.tmu["prev.c.tmu"] * xc.pdf.tmu
  # x.kde.tmu[x.kde.tmu < fig100.017.lim["yplotmin"]] <- 
  #                                                 fig100.017.lim["yplotmin"]
             
  PlotMetRes(fig100.017, fig100.017.file, figtype, fig100.017.lim["yplotmin"],
             xsupp, tab.tmu["prev.c.tmu"], xc.pdf.tmu, tmucol,
             x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
             tab.tmu["x.RL1.tmu"], tab.tmu["x.RL2.tmu"], 
             0.85*fig100.017.lim["yplotmax"], 
             tab.tmc["x.tr.lo"], tab.tmc["x.tr.hi"], 
             0.90*fig100.017.lim["yplotmin"], tmccol,
             temp.tmu)

  # -----------------------------------------------------------
  #  If the data is generated : show unaffected range
  if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
  { 
    lines(c(x.unaffected.lo, x.unaffected.hi), 
            0.60*fig100.017.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
  }

  # -------------------------------------------------------------
  # Save tmu result figure
  savePlot(file=fig100.017.file, type=figtype)
}

