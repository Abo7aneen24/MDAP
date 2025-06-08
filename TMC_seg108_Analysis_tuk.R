#  TMC_seg108_Analysis_tuk.R 
#
#  Tukey analysis, separate solution with standard settings,
#  not to be confused with the Tukey-QQplot procedure used to get initial 
#  values for TMC

#  Transformation paramter lambda is taken from TMC solution

#
#  25.04.2021
# ==============================================================================

if (print.log.message) { cat("%%%   TMC_seg108_Analysis  100\n") }

#  Transform data according to tmc result
y <- BoxCox(x, tab.tmc["lambda.tmc"])

#  Run Tukey analysis

#tab.tuk <- TukeyWW(y, xlabel, IQR.fact=IQR.fact, RL1.p=RL1.p, RL2.p=RL2.p, 
#                    figA=NA)
tab.tuk <- TukeyWW(y, xlabel, IQR.fact=IQR.fact, iter.max=iter.max.tuk, 
                    x.tr.prop.min=x.tr.prop.min, RL1.p=RL1.p, RL2.p=RL2.p, 
                    print.details=FALSE, figA=NA)

#print(tab.tuk)

#return(c(mean=mean.hat, sd.raw=sigma.raw, sd.adj=sigma.adj, 
#           x.tr.lo=x.lo, x.tr.hi=x.hi, 
#           out.lo=out.lo, out.hi=out.hi, x.RL1=x.RL1, x.RL2=x.RL2, xc.n=xc.n)) 

#  Back transformation to x scale
x.tr.lo.tuk <- BoxCoxInv(tab.tuk["x.tr.lo"], tab.tmc["lambda.tmc"]) 
x.tr.hi.tuk <- BoxCoxInv(tab.tuk["x.tr.hi"], tab.tmc["lambda.tmc"]) 
x.RL1.tuk <- BoxCoxInv(tab.tuk["x.RL1"], tab.tmc["lambda.tmc"]) 
x.RL2.tuk <- BoxCoxInv(tab.tuk["x.RL2"], tab.tmc["lambda.tmc"]) 

tab.tuk["x.tr.lo"] <- x.tr.lo.tuk
tab.tuk["x.tr.hi"] <- x.tr.hi.tuk
tab.tuk["x.RL1.tuk"] <- x.RL1.tuk
tab.tuk["x.RL2.tuk"] <- x.RL2.tuk

temp.tuk <- chi2trunc(x.hist, x.tr.lo.tuk, x.tr.hi.tuk,
                      sum(x < x.tr.lo.tuk), 
                      sum(x >= x.tr.hi.tuk),
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      tab.tmc["lambda.tmc"],
                      tab.tuk["mean"],
                      tab.tuk["sd.adj"],
                      df.est, df.con,
                      l.fact, p.fact, r.fact, w.fact, 
                      opt.crit.only=FALSE,fastnull=fastnull)

#save(temp.tuk, file="../dev/temp.tuk")   # used for development of HistoResid


#  Estimated prevalence cannot be taken from chi2trunc because the Tukey 
#  approach decomposes the data into non-overlapping components, while tmc, 
#  tmu assume overlaps.
              
tab.tuk <- c(tab.tuk, 
             opt.crit.tuk=unname(temp.tuk[["res"]]["opt.crit"]),
             p.fit.tuk=unname(temp.tuk[["res"]]["chi2.trun.p"]),
             prev.l.tuk=unname(tab.tuk["out.lo"]/x.n),
             prev.c.tuk=unname((x.n-tab.tuk["out.lo"]-tab.tuk["out.hi"])/
                               x.n), 
             prev.r.tuk=unname(tab.tuk["out.hi"]/x.n) )

# ------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg108_Analysis  200\n") }
  
if (plot.fig108.010)
{ 
  #  Plot result - basic part 
  #  No plot file name, save is done here below
  fig108.010.lim <- PlotHistFit(fig108.010, NA, figtype,
                    x.hist, xlabel, paste("Tukey result, IQR.fact", IQR.fact), 
                    subtitle, n.per.bin.min,
                    NA, x.val, counts.range, 
                 paste(formatC(100*tab.tuk["xc.n"]/x.n, format="f", width=5,
                   digits=1), "%", sep=""),
                   x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                   x.clip.type, NA, par.las, par.tcl,
                   bordercol, histcol, difcol, kdecol, denlty)

  #  Plot result - tuk specific part 
  #  Add estimated tuk solution, RLs, TI
  #  Calculate estimated tuk density (not yet scaled by prev.c.tuk)
  xc.pdf.tuk <- pdf.PN(xsupp, tab.tmc["lambda.tmc"], 
                              tab.tuk["mean"],
                              tab.tuk["sd.adj"]) 

  #  Calculate difference between x.kde and xc.pdf.tuk
  #  Truncate difference below to yplotmin
  #  Switched off 01.05.2021
  #x.kde.tuk <- x.kde$y -  tab.tuk["prev.c.tuk"] * xc.pdf.tuk
  #x.kde.tuk[x.kde.tuk < fig108.010.lim["yplotmin"]] <- 
  #                                             fig108.010.lim["yplotmin"]
           
  PlotMetRes(fig108.010, fig108.010.file, figtype, fig108.010.lim["yplotmin"], 
             xsupp, tab.tuk["prev.c.tuk"], xc.pdf.tuk, tukcol,
             x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
             tab.tuk["x.RL1.tuk"], tab.tuk["x.RL2.tuk"], 
             0.85*fig108.010.lim["yplotmax"], 
             tab.tuk["x.tr.lo"], tab.tuk["x.tr.hi"], 
             0.90*fig108.010.lim["yplotmin"], tukcol,
             temp.tuk)

  # -----------------------------------------------------------
  #  If the data is generated : show unaffected range
  if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
  { 
    lines(c(x.unaffected.lo, x.unaffected.hi), 
          0.60*fig108.010.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
  }

  # -----------------------------------------------------------
  #  Save the tuk result plot

  savePlot(file=fig108.010.file, type=figtype) 
}
