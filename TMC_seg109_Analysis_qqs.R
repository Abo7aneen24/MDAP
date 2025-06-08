#  TMC_seg109_Analysis_qqs.R 
#
#  QQ-Plot solution, separate solution with standard settings,
#  not to be confused with the Tukey-QQplot procedure used to get initial 
#  values for TMC

#  Transformation paramter lambda is taken from TMC solution

#
#  30.05.2022
# ==============================================================================

if (print.log.message) { cat("%%%   TMC_seg109_Analysis  100\n") }

#  Transform data according to tmc result
y <- BoxCox(x, tab.tmc["lambda.tmc"])

#  Run QQ plot analysis  
#  Output is
#y.mue.qq0=y.mue.qq0,
#           y.sig.qq0=y.sig.qq0,
#           y.mue.qq1=unname(y.qq.lm1$coefficients[1]),
#           y.sig.qq1=unname(y.qq.lm1$coefficients[2]),
#           rsq=unname(y.qq.lm1.sum$adj.r.squared)) )

tab.qqs.names <- c("mean", "sd.adj", "x.tr.lo", "x.tr.hi", "x.RL1.qqs", 
                   "x.RL2.qqs", "xc.n")
tab.qqs <- rep(NA, times=length(tab.qqs.names))
names(tab.qqs) <- tab.qqs.names

temp.qqs <- QQPlotIMSS(y, RL1.p, RL2.p, NA, NA, fastnull=1.e-10)
tab.qqs["mean"]   <- temp["y.mue.qq1"]
tab.qqs["sd.adj"] <- temp["y.sig.qq1"]
tab.qqs["x.tr.lo"] <- tab.tuk["x.tr.lo"]
tab.qqs["x.tr.hi"] <- tab.tuk["x.tr.hi"]
tab.qqs["x.RL1.qqs"]   <- q.PN(RL1.p, tab.tmc["lambda.tmc"],
                           temp["y.mue.qq1"], temp["y.mue.qq1"])
tab.qqs["x.RL2.qqs"]   <- q.PN(RL2.p, tab.tmc["lambda.tmc"],
                           temp["y.mue.qq1"], temp["y.mue.qq1"])


#  output frum tukey is
#  return(c(mean=mean.hat, 
#           sd.raw=sigma.raw, sd.adj=unname(sigma.adj), 
#           x.tr.lo=unname(x.lo), x.tr.hi=unname(x.hi), 
#           z.tr.lo=unname(z.lo), z.tr.hi=unname(z.hi), 
#           out.lo=out.lo, out.hi=out.hi,            
#           x.RL1=x.RL1, x.RL2=x.RL2, 
#           xc.n=xc.n)) 


#print(tab.qqs)


temp.qqs <- chi2trunc(x.hist, tab.qqs["x.tr.lo"], tab.qqs["x.tr.hi"],
                      sum(x < tab.qqs["x.tr.lo"]), 
                      sum(x >= tab.qqs["x.tr.hi"]),
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      tab.tmc["lambda.tmc"],
                      tab.qqs["mean"],
                      tab.qqs["sd.adj"],
                      df.est, df.con,
                      l.fact, p.fact, r.fact, w.fact, 
                      opt.crit.only=FALSE,fastnull=fastnull)

#  Estimated prevalence cannot be taken from chi2trunc because the Tukey 
#  approach decomposes the data into non-overlapping components, while tmc, 
#  tmu assume overlaps.
              
tab.qqs <- c(tab.qqs, 
             opt.crit.qqs=unname(temp.qqs[["res"]]["opt.crit"]),
             p.fit.qqs=unname(temp.qqs[["res"]]["chi2.trun.p"]),
             prev.l.qqs=unname(tab.tuk["out.lo"]/x.n),
             prev.c.qqs=unname((x.n-tab.tuk["out.lo"]-tab.tuk["out.hi"])/
                               x.n), 
             prev.r.qqs=unname(tab.tuk["out.hi"]/x.n) )

# ------------------------------------------------------------
if (print.log.message) { cat("%%%   TMC_seg109_Analysis  200\n") }
  
#  for the rest from here see seg 108
