#  TMC_seg101_Ini_Analysis.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Initialize output tables for TMC_seg100

#  CHANGE HISTORY

#  15.02.2021 Start 
# ===========================================================================


tab.npa <- data.frame(method="npa",  
                      x.tr.lo.npa=NA, x.tr.hi.npa=NA,
                      x.tr.n.npa=NA, x.tr.n.npa=NA,
                      xl.n.npa=NA, xc.n.npa=NA, xr.n.npa=NA,
                      prev.l.npa=NA,prev.c.npa=NA, prev.r.npa=NA,
                      x.RL1.npa=NA, x.RL2.npa=NA)   

tab.qqw.names <- c(
 "lambda.qqw" ,     "bin.start"       ,"bin.end"    ,     "bins.n"    ,     
 "prop.qqw"   ,     "x.tr.lo"         ,"x.tr.hi"    ,     "y.tr.lo"   ,     
 "y.tr.hi"    ,     "mue.qqw"         ,"sigma.qqw"  ,     "r2"        ,     
 "opt.crit"   ,     "x.tr.n"          ,"rank.r2"    ,     "y.RL1.qqw" ,     
 "y.RL2.qqw"  ,     "x.RL1.qqw"       ,"x.RL2.qqw"  ,     "n.l.qqw"   ,     
 "prev.l.qqw" ,     "n.c.qqw"         ,"prev.c.qqw" ,     "n.r.qqw"   ,     
 "prev.r.qqw" ,     "neg.prev.sum.qqw","x.RL1.cilo" ,     "x.RL1.cihi",     
 "x.RL2.cilo" ,     "x.RL2.cihi")

tab.qqw <- matrix(NA, nrow=1, ncol=length(tab.qqw.names))
colnames(tab.qqw) <- tab.qqw.names

tab.tmc.names <- c(
 "ilo"           ,"ihi"         ,  "x.tr.lo"    ,   "x.tr.hi"       ,
 "x.tr.n"        ,"x.tr.prop"   ,  "rc"         ,   "iter"          ,
 "lambda.tmc"    ,"mue.tmc"     ,  "sigma.tmc"  ,   "x.RL1.tmc"     ,
 "x.RL2.tmc"     ,"x.RL1.cilo"  ,  "x.RL1.cihi" ,   "x.RL2.cilo"    ,
 "x.RL2.cihi"    ,"reldist.tmc" ,  "p.fit"      ,   "opt.crit"      ,
 "prev.l.tmc"    ,"prev.c.tmc"  ,  "prev.r.tmc" ,   "prev.l.tmc.pen",
 "prev.r.tmc.pen","xc.n.tmc"    ,  "chi2.total" ,   "chi2.total.df" ,
 "chi2.trun"     ,"chi2.trun.df",  "chi2.path"  ,   "sol.score" )

tab.tmc <- matrix(NA, nrow=1, ncol=length(tab.tmc.names))
colnames(tab.tmc) <- tab.tmc.names

crit <- rep(NA, times=5)

# colnames(tab.tmu)
tab.tmu.names <- c(
 "mue.tmu"    , "sigma.tmu"   ,"distance.tmu","rc.tmu"     , "iter.tmu"  , 
 "x.RL1.tmu"  , "x.RL2.tmu"   ,"x.RL1.cilo"  ,"x.RL1.cihi" , "x.RL2.cilo", 
 "x.RL2.cihi" , "opt.crit.tmu","p.fit.tmu"   ,"prev.l.tmu" , "prev.c.tmu", 
 "prev.r.tmu") 

tab.tmu <- matrix(NA, nrow=1, ncol=length(tab.tmu.names))
colnames(tab.tmu) <- tab.tmu.names

