#  TMC_seg104_Analysis_qqw.R
#
#  15.02.2021
# ==============================================================================

if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw  100\n") }

#  In next call: bins.n.min replaced by x.tr.bins.min
tab.qqw <- QQW(x, x.reso.red, round.unit, lambda.seq, 
                   x.hist, x.kde, x.kde.mode, x.kde.mode.idx, x.val, 
                   x.unaffected.lo, x.unaffected.hi,  
                   n.per.bin.min, bins.n.min, x.tr.bins.min, 
                   counts.range, x.tr.prop.range,
                   x.tr.prop.limits[iti], x.tr.prop.limits[iti+1], 
                   prev.acc.lo, prev.acc.hi, 
                   xsupp,
                   x.Q1, x.Q2, RL1.p, RL2.p,
                   lambda.c.gen, yc.mode.gen, yc.sig.gen,
                   df.est,df.con,
                   l.fact, p.fact, r.fact, w.fact, 
                   xlabel, subtitle,
                   fig100.006, fig100.007,
                   fig100.008, fig100.009,
                   fig100.010, 
                   fig100.011, fig100.011.file,
                   fig100.012, fig100.012.file, figtype,
                   x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                   x.clip.type, par.las, par.tcl,
                   bordercol2, histcol2, kdecol, qqrawcol, polycol, gencol, 
                   difcol, denlty,
                   print.log.message)

# --------------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 200\n") }

#  QQW may have produced no result. Then go to next range.
#  tab.qqw may be an NA vector

QqwSolution <- (length(tab.qqw) > 1) && !all(is.na(tab.qqw))
  
if (!QqwSolution)
{ #  No result
  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 300\n") }

  cat("\n +++ [seg104] No QQW result for proportion range ",
              x.tr.prop.range,
      "\n +++", subtitle,
      "\n\n")

  #  For safety
  temp.qqw <- NA 

}  else
{
  # QQW solution exists, add more details 

  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 400\n") }

  temp.qqw <- chi2trunc(x.hist, tab.qqw[1, "x.tr.lo"],
                        tab.qqw[1, "x.tr.hi"],
                        sum(x < tab.qqw[1, "x.tr.lo"]), 
                        sum(x >= tab.qqw[1, "x.tr.hi"]),
                        x.Q1, x.Q2, RL1.p, RL2.p,
                        tab.qqw[1, "lambda.qqw"],
                        tab.qqw[1, "mue.qqw"],
                        tab.qqw[1, "sigma.qqw"],
                        df.est, df.con,
                        l.fact, p.fact, r.fact, w.fact, 
                        opt.crit.only=FALSE,fastnull=fastnull)

  # -------------------------------------------------------------------------
  #  Plot the QQW result
  if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 420\n") }

  if (!is.na(fig100.012))
  { dev.set(fig100.012)
    bringToTop(fig100.012)

    #  Plot result - basic part x.kde set to NA
    fig100.012.lim <- PlotHistFit(fig100.012, fig100.012.file, figtype, 
                  x.hist, xlabel, "QQW result", subtitle, n.per.bin.min,
                  NA, x.val, counts.range, x.tr.prop.range,
                  x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                  x.clip.type, NA, par.las, par.tcl,
                  bordercol2, histcol2, difcol, kdecol, denlty)

    if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 440\n") }

    #  Plot result - qqw specific part 
    #  Add estimated QQW solution, RLs, TI
    #  Calculate estimated QQW density (not yet scaled by prev.c.qqw)
    xc.pdf.qqw <- pdf.PN(xsupp, unlist(tab.qqw[1, "lambda.qqw"]), 
                                unlist(tab.qqw[1, "mue.qqw"]), 
                                unlist(tab.qqw[1, "sigma.qqw"])) 

    #  Calculate difference between x.kde and xc.pdf.qqw
    #  Truncate difference below to yplotmin
    #  Switched off 01.05.2021
    #x.kde.qqw <- x.kde$y -  tab.qqw["prev.c.qqw"] * xc.pdf.qqw
    #x.kde.qqw[x.kde.qqw < fig100.012.lim["yplotmin"]] <- 
    #                                              fig100.012.lim["yplotmin"]

    if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 460\n") }

    #  No plot file name, save is done farther below
    PlotMetRes(fig100.012, NA, figtype, 
               fig100.012.lim["yplotmin"], 
               xsupp, tab.qqw[1, "prev.c.qqw"], xc.pdf.qqw, qqwcol,
               x.kde.mode, kdecol, NA, difcol, difbordercol, difhistcol,
               unlist(tab.qqw[1, "x.RL1.qqw"]), unlist(tab.qqw[1, "x.RL2.qqw"]), 
               0.85*fig100.012.lim["yplotmax"], 
               unlist(tab.qqw[1, "x.tr.lo"]), unlist(tab.qqw[1, "x.tr.hi"]), 
               0.90*fig100.012.lim["yplotmin"], qqwcol,
               temp.qqw)

    # --------------------------------------------------------------------------
    #  If the data is generated : show unaffected range
    if (!is.na(x.unaffected.lo) & !is.na(x.unaffected.hi))
    { 
      lines(c(x.unaffected.lo, x.unaffected.hi), 
      0.60*fig100.012.lim["yplotmin"]*c(1,1), col="chartreuse", lwd=2)
    }

    # -------------------------------------------------------------
    # Save qqw result figure

    savePlot(file=fig100.012.file, type=figtype) 
  }  
} 

if (print.log.message) { cat("%%%   TMC_seg104_Analysis_qqw 500\n") }

