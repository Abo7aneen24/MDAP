
#  TMC_seg022_Functions.R 
#  (c) wwosniok@math.uni-bremen.de 
 
#  Truncated minimum chi-square estimation 
 
#  Compilation of all functions required by TMC 
#  see ../func/CompileFunctions.R for the generation of this file 
 
#  #  CHANGE HISTORY 
#  10.04.2021 All actual experimental functions included here 
#  16.02.2021 Actual functions collected in TMC_seg022_Functions.R 
#  26.01.2021 Start as update of TMC_seg020_Functions.R 
# ============================================================================== 

#  Compilation date   Mon Nov 14 17:32:38 2022

#  Functions compiled in this file:

# [1] "F_AgeLimitsToClass.R"     "F_Bhatta_V2.R"           
# [3] "F_BoxCox.R"               "F_BoxCoxInv.R"           
# [5] "F_CalcChi2.R"             "F_CalcParDistance.R"     
# [7] "F_CalcPrev_V2.R"          "F_CalcR2InSubInt.R"      
# [9] "F_cdf.PN.R"               "F_Ceiling.R"             
#[11] "F_CheckRounding.R"        "F_chi2.PNV.tr.lms_V2.R"  
#[13] "F_chi2trunc_V14.R"        "F_CIQuant.LNV.R"         
#[15] "F_CIQuant.PNV.R"          "F_CollapseBins_V5.R"     
#[17] "F_color.palette.R"        "F_CompileTabcn.R"        
#[19] "F_Compose.R"              "F_Drift.R"               
#[21] "F_DynAgeLimits.R"         "F_EstParentDist.R"       
#[23] "F_Expand.R"               "F_Extract.R"             
#[25] "F_FC.R"                   "F_FilterWDay.R"          
#[27] "F_FindBw_V4.R"            "F_FindHisto_V13.R"       
#[29] "F_FindKDE.R"              "F_FindModeIndex.R"       
#[31] "F_Floor.R"                "F_FprContour.R"          
#[33] "F_FprLegend.R"            "F_Hist.table.R"          
#[35] "F_Infoblock.R"            "F_IniMatrix.R"           
#[37] "F_IniTab.R"               "F_Interpolx.R"           
#[39] "F_Interpoly.R"            "F_modified.tukey.R"      
#[41] "F_npa.R"                  "F_o2R.R"                 
#[43] "F_opt.crit2.R"            "F_pdf.PN.R"              
#[45] "F_PermissibleDiff.R"      "F_plot.res.R"            
#[47] "F_PlotConfElli.R"         "F_PlotDS.R"              
#[49] "F_PlotHistFit.R"          "F_PlotHistResid.R"       
#[51] "F_PlotMetRes.R"           "F_PlotRL.R"              
#[53] "F_ProcessDL.R"            "F_ProcessPropLimits_V1.R"
#[55] "F_Proportions_V9.R"       "F_q.PN.R"                
#[57] "F_QQEstimateIni_V1.R"     "F_QQPlotIMSS.R"          
#[59] "F_QQW.R"                  "F_Quantile.R"            
#[61] "F_R2o.R"                  "F_reldist.R"             
#[63] "F_RL.age.spline.R"        "F_Round.R"               
#[65] "F_RowToMatrix.R"          "F_ShowTiProps.R"         
#[67] "F_SmoothHist1_V7.R"       "F_SmoothHist2.R"         
#[69] "F_tmc.gridsearch.R"       "F_tmc.master_V8.R"       
#[71] "F_TransformKDE.R"         "F_TukeyWW.R"             
#[73] "F_WaWoTest.R"             "F_x.mean.PN.R"           
#[75] "F_x.mode.PN.R"            "F_x.pdf.PN.R"            
#[77] "F_x.Var.PN.R"             "F_xx.pdf.PN.R"           

# **************************************************************** 
# **************************************************************** 

# F_AgeLimitsToClass.R
#  25.02.2022
#
#  (c) wwosniok@math.uni-bremen.de
# ============================================================================

AgeLimitsToClass <- function(age.limits)
{ #  Change age.limits to age.class table
  #  
  #  INPUT 
  #  age.limits  vector aof age limit borders
  #
  #  OUTPUT
  #  age.class   Matrix of age class limits. All limits belong to the 
  #              interval
  # =========================================================================

  age.class     <- c(lo=NA, hi=NA)
  age.limits.u  <- sort(unique(age.limits))
  age.limits.n  <- length(age.limits.u)
  age.class.n   <- length(age.limits.u) - 1
  if (age.class.n > 0)
  {
    age.class <- matrix(NA, nrow=age.class.n, ncol=2)
    colnames(age.class) <- c("lo", "hi")
    age.class[ , "lo"] <- age.limits.u[1:(age.limits.n-1)]
    age.class[ , "hi"] <- age.limits.u[2:age.limits.n]-1
  }
  age.class[age.class.n, "hi"] <- age.class[age.class.n, "hi"] + 1
  age.class <- RowToMatrix(age.class)
  
  return(age.class)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# source("../func/F_RowToMatrix.R")

# age.limits <- c(18, 25, 35)  
# age.class  <- AgeLimitsToClass(age.limits) 
# print(age.class)

# age.limits <- c(18, 25, 35, 45, 50, 35, 23)  
# age.class  <- AgeLimitsToClass(age.limits) 
# print(age.class)

# age.limits <- 99
# age.class  <- AgeLimitsToClass(age.limits) 
# print(age.class)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_Bhatta_V2.R
#
#  Bhattacharya-Methode nach CG Bhattacharya, Biometrics 23 (1967) 115-135
#
#  (c) wwosniok@math.uni-bremen.de
#
#  To do
#  - Warum die getroffene Auswahl unter den 3 Lösungen für s - siehe @@@@ 
#    unten?
#  - Bhatta-Transformation glätten, siehe Oosterhuis 1990
#    Könnte leicht geschehen durch Verwendung von x.hist.smo

#  11.10.2020 xhist changed to x.hist in accordance with TMC use
#             x.hist must be provided, condition h not < round.unit must be 
#             ensured during histogram generation
#  10.10.2020 Histogram generation as in tmc
#  06.10.2020 Input revised
#  28.11.2018 Start
# 
# ============================================================================

Bhatta <- function(x.hist, x.kde, data.text, xlabel, subtitle, RL1.p,RL2.p, 
                   figD, figE,
                   kdecol, bhacol, bordercol1, histcol1)
{ #
  #  INPUT
  #
  # x.hist     data as equidistant histogram, bin width >= round.unit 
  #            in the original raw data
  # x.kde      kernel density estimate, for plotting only, may  be NA
  # data.text  description of the data
  # RL1.p      level for RL1
  # RL2.p      level for RL2
  # figD       figure for data histogram + kde + fitted distributions 
  # figE       figure for Bhattacharya plot

  #  OUTPUT
  #
  # =========================================================================

  #  Extract necessary information
  xcounts <- x.hist$counts
  xbreaks <- x.hist$breaks
  xmids   <- x.hist$mids
  xcounts.n <- length(xcounts)
  xbreaks.n <- length(xbreaks)
  x.n       <- sum(xcounts)
  h         <- xmids[2] - xmids[1]

  xsupp <- seq(x.kde$x[1], tail(x.kde$x, 1), length.out=101) 

  # ..........................................................................
  # Plot histogram

  if (!is.na(figD))
  {
    dev.set(figD)
    plot(x.hist,freq=FALSE,
         main="Histogram for Bhattacharya method",
         xlab=xlabel, 
         border=bordercol1, col=histcol1, sub=subtitle,cex.sub=0.7)
    if (length(x.kde) > 0) { lines(x.kde,col=kdecol) }
  }

  # ..........................................................................
  #  Calculate delta.log.y, the Bhattacharya transformation. 
  #  To avoid problems with empty cells: replace 0 by fastnull

  fastnull <- 0.5/x.n
  xcounts.1 <- xcounts
  xcounts.1[xcounts==0] <- fastnull

  delta.log.y   <- log(xcounts.1[2:xcounts.n]) - 
                   log(xcounts.1[1:(xcounts.n-1)])

  # ..........................................................................
  # Plot delta.log.y

  if (!is.na(figE))
  {  
    dev.set(figE)
    plot(xmids[1:(xcounts.n-1)], delta.log.y, type="o", col=bhacol,
         main="Bhattacharya transformation",
         xlab=xlabel, ylab="log(f_i+1 / f_i)",
         sub=subtitle,cex.sub=0.7)  
    abline(h=0)
  }

  # ..........................................................................
  #  Select subsets (truncation areas, at least 2 distributions are suspected)
  #  Automatic selection
  #  Find minima / maxima
  delta.log.y.n <- length(delta.log.y)
  M.names <- c("xc","xcounts","yl","yc","yr","is.min","is.max","is.ext")
  M <- data.frame(matrix(NA,nrow=delta.log.y.n,ncol=length(M.names)))
  colnames(M) <- M.names

  M[ ,"xc"]      <- xmids[1:delta.log.y.n]
  M[ ,"xcounts"] <- xcounts[1:delta.log.y.n]
  M[2:delta.log.y.n, "yl"] <- delta.log.y[1:(delta.log.y.n-1)]
  M[ ,"yc"] <- delta.log.y
  M[1:(delta.log.y.n-1), "yr"] <- delta.log.y[2:(delta.log.y.n)]

  innen <- 2:(delta.log.y.n-1)

  M[1,"is.min"]     <- (M[1,"yc"] < M[1,"yr"])
  M[innen,"is.min"] <- (M[innen,"yl"] > M[innen,"yc"]) & 
                       (M[innen,"yc"] < M[innen,"yr"]) 
  M[delta.log.y.n,"is.min"] <- (M[delta.log.y.n,"yl"] > M[delta.log.y.n,"yc"])

  M[1,"is.max"]     <- (M[1,"yc"] > M[1,"yr"])
  M[innen,"is.max"] <- (M[innen,"yl"] < M[innen,"yc"]) & 
                       (M[innen,"yc"] > M[innen,"yr"])
  M[delta.log.y.n,"is.max"] <- (M[delta.log.y.n,"yl"] < M[delta.log.y.n,"yc"])

  # is.min and is.max may be NA due to zeroes in the data. Set to FALSE
  M[is.na(M[ ,"is.min"]), "is.min"] <- FALSE
  M[is.na(M[ ,"is.max"]), "is.max"] <- FALSE

  #  Identical values of delta.log.y are possible. Check for regular min/max
  #  alternating in the extrema sequence.
  is.ext <- rep(NA,times=delta.log.y.n)
  M[M[ ,"is.min"],"is.ext"]  <-  "min"
  M[M[ ,"is.max"],"is.ext"]  <-  "max"

  #cat("\n[Bhatta] Matrix M\n")
  #print(M)

  #  @@@@ hier fehlt noch was (min/max-Reihenfolge überpüfen)

  #  Interesting pieces start with max
  #  Indices are indices in xmids
  max.idx <- which(M[ ,"is.max"])
  max.idx.n <- length(max.idx)

  min.idx <- which(M[ ,"is.min"])

  res.names <- c("istart","iend","nobs","adjrsq","mue.hut","sig.hut",
                 "RL1","RL2","RL1.1","RL2.1","RL1.2","RL2.2","RL1.3","RL2.3",
                 "p.hut","n.hut","prop.hut")
  res <- matrix(NA,nrow=length(max.idx),ncol=length(res.names))
  colnames(res) <- res.names

  i <- 0
  for (istart in max.idx)
  { i <- i + 1
    
    #cat("\n Starting search for component ", i, "------------------------ \n")

    res[i,"istart"] <- istart
    #  Find next minimum
    #iend <- min(min.idx[min.idx>istart])  #### @@@ iend könnte inf sein!
    
    # Avoid warnings:
    min.idx.up <- min.idx>istart
    if (sum(min.idx.up)== 0 )
    { # no more minima
      iend <- -Inf
    } else
    { # more minima
      iend <- min(min.idx[min.idx.up])
    }

    res[i,"iend"] <- iend

    #cat("\n[Bhatta] istart, iend ",istart, iend,"\n")

    if (is.finite(iend) & ((iend-istart)>1) )
    {
      if (!is.na(figE))
      {
        dev.set(figE)
        abline(v=M[istart,"xc"],col="red")
        abline(v=M[iend,"xc"],col="green4")
      }

      # subset: indices of points in the truncation interval
      subset <- istart:iend
      x.sub  <- M[subset,"xc"]
      y.sub  <- M[subset,"yc"]

      delta.log.y.lm <- lm(y.sub ~ x.sub)
      delta.log.y.lm.sum <- summary(delta.log.y.lm)   
      #print(delta.log.y.lm.sum)

      if (!is.na(figE))
      {
        dev.set(figE)
        abline(delta.log.y.lm$coefficients,col="green3",lty=2)
      } 

      beta0 <- delta.log.y.lm$coefficients[1]
      beta1 <- delta.log.y.lm$coefficients[2]
      res[i,"adjrsq"] <- unlist(delta.log.y.lm.sum["adj.r.squared"])

      lambda    <- -beta0 / beta1   #  this is not the PND lambda!
      mue.hut   <-  lambda + h/2
      s.hut.1   <-  sqrt(-h/beta1) - h^2/12   #  Bhattacharya, (2)
      s.hut.2   <-  sqrt(-h/beta1)            #  Sparre & Venema

      #  Lösung der quadratischen Gleichung für sigma^2
      discrim <- 1/(4*beta1^2) + h/(12*beta1)

      sigma2.1 <- NA
      sigma2.2 <- NA
      s.hut.3  <- 0
      if (discrim >= 0)
      {
        sigma2.1 <- -h/(2*beta1) + h*sqrt(discrim)
        sigma2.2 <- -h/(2*beta1) - h*sqrt(discrim)
        s.hut.3  <- sqrt(sigma2.1)
      }

      res[i,"mue.hut"] <- mue.hut
      res[i,"sig.hut"] <- s.hut.2  # @@@@@ warum diese Lösung als Resultat?

      #  Estimate the size of this subpopulation
      #  Probability in this truncation interval 
      p.tr <- pnorm(tail(x.sub,1),mean=mue.hut,sd=s.hut.2) -
              pnorm(x.sub[1],mean=mue.hut,sd=s.hut.2)
      res[i,"p.hut"] <- p.tr

      #  Counts in the truncation interval
      xcounts.tr.n  <- sum(M[subset,"xcounts"])
      res[i,"nobs"] <- xcounts.tr.n

      #  Bezeichnung xc ist nicht gut
      xc.hut <- xcounts.tr.n/p.tr
      res[i,"n.hut"] <- xc.hut
      res[i,"prop.hut"]  <- xc.hut/x.n
      res[i,"RL1"]   <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=res[i,"sig.hut"] )
      res[i,"RL2"]   <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=res[i,"sig.hut"] )
      res[i,"RL1.1"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.1 )
      res[i,"RL2.1"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.1 )
      res[i,"RL1.2"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.2 )
      res[i,"RL2.2"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.2 )
      res[i,"RL1.3"] <- qnorm(RL1.p,mean=res[i,"mue.hut"], sd=s.hut.3 )
      res[i,"RL2.3"] <- qnorm(RL2.p,mean=res[i,"mue.hut"], sd=s.hut.3 )

      #cat("\nBHATTACHARYA method",
      #"\nVersion 1   Mean  ", mue.hut,"  sd  ",s.hut.1,
      #"\nVersion 2   Mean  ", mue.hut,"  sd  ",s.hut.2,
      #"\nVersion 3   Mean  ", mue.hut,"  sd  ",s.hut.3,
      #"\nEstimated subpopulation size            ",round(xc.hut),
      #"\nCorresponding proportion of total data  ",xc.hut/x.n,
      #"\n")

      if (!is.na(figD))
      {
        dev.set(figD)
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.1),col="blue")
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.2),col="red",
            lwd=2)
        lines(xsupp,(xc.hut/x.n)*dnorm(xsupp,mean=mue.hut,sd=s.hut.3),
            col="orange")

        legend("topleft",c("Bhattacharya (2)","Sparre & Venema","Quad. equation"),
             col=c("blue","red","orange"),lwd=c(1,2,1),cex=0.7)
      } 

    }   # is.finite ...

    #cat("\n Search for component ", i, "  finished ---------------------- \n")

  }     # for (istart ...

  # Make sure that res is a matrix
  if (!is.matrix(res))
  { res <- matrix(res,nrow=1) 
    #  no sorting need
  } else
  { 
    res <- res[!is.na(res[ ,"n.hut"]), ]
    if (!is.matrix(res))
    { res <- matrix(res,nrow=1) 
      #  no sorting need
    } else
    {   
      res <- res[order(-res[ ,"n.hut"]), ]
    }
  }
  colnames(res) <- res.names
  rownames(res) <- paste("comp", 
                         formatC(1:nrow(res), format="f", digits=0, width=2, 
                                 flag="0"), sep="")

  return(list(x.hist.bha=x.hist,subpop=res))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_BoxCox.R
#
#  Box-Cox-Transformation mit Parameter lambda
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY
#
#  18.07.2018 WW   xgt0   <- (!is.na(x)) & (     x > fastnull)
#                  geändert in 
#                  xgt0   <- (!is.na(x)) & (     x >= fastnull)
#
# --------------------------------------------------------------------------
# Eingabe
# x         zu transformierender Vektor ("Originalskala" = beobachtet)
# lambda    Box-Cox-Parameter,
#           lambda = 0: entspricht log(x)
#           lambda = 1: entspricht x-1
# fastnull  Absicherung gegen kleine Abbruchfehler

# Ausgabe
# Tx        transformierter Vektor 
#                           x < 0  Tx = NA   (obwohl sinnvolle Werte für 
#                                             ganzzahliges lambda möglich)
#           lambda <= 0 und x = 0  Tx = -Inf
#           lambda =  0 und x > 0  Tx = log(x)
#           lambda >  0 und x = 0  Tx = 0
#           lambda <> 0 und x > 0  Tx = (x^lambda-1)/lambda
# ============================================================================
             
BoxCox <- function(x,lambda,fastnull=1.e-10)
{ 
  xgt0   <- (!is.na(x)) & (     x >= fastnull)

  xeq0   <- (!is.na(x)) & (abs(x) < fastnull)

  Tx     <- rep(NA,times=length(x))

  if (lambda < -fastnull)
  { Tx[xeq0] <- -Inf
    Tx[xgt0] <- (x[xgt0]^lambda-1)/lambda 
  }
  if (abs(lambda) <= fastnull) 
  { 
    Tx[xeq0] <- -Inf
    Tx[xgt0] <- log(x[xgt0])
  }
  if (lambda > fastnull)
  { 
    # Tx[xeq0] <- (x[xeq0]^lambda-1)/lambda 
    # Tx[xgt0] <- (x[xgt0]^lambda-1)/lambda 
    Tx         <- (x^lambda-1)/lambda 
  }

  return(Tx)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#BCx <- BoxCox(10.5,1.29247e-25)
#BoxCox(0,1.001) 

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_BoxCoxInv.R

#  Inverse Box-Cox-Transformation withrameter lambda
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY

#  27.06.2016 Start
# --------------------------------------------------------------------------

# INPUT
# Tx       vector to transform 
# lambda   Box Cox parameter
#          lambda = 0: same as exp(x)
#          lambda = 1: produces x+1
# fastnull value replacing zero 

# Ausgabe
# x        transformed vector, x = (1 + lambda*Tx)^(1/lambda)
# ============================================================================

BoxCoxInv <- function(Tx,lambda,fastnull=1.e-10)
{ 
  x <- rep(NA,times=length(Tx))
 
  if (lambda < -fastnull)
  { x <- (Tx*lambda+1)^(1/lambda) }

  if (lambda > fastnull)
  { x <- (Tx*lambda+1)^(1/lambda) }
  
  if (abs(lambda) <= fastnull)
  { xneInf       <- is.finite(Tx)
    x[xneInf]    <- exp(Tx[xneInf])
    if (any(!xneInf)) 
    { x[!xneInf] <- 0  }
  }
  
  return(x)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#x <- BoxCoxInv(5, 1)   # 6
#x <- BoxCoxInv(5, 0)   # 148.4132 = exp(5)
#x <- 5
#y <- BoxCox(x, 0)      # 1.609438 = log(5)
#xx <- BoxCoxInv(y, 0)  # 5
#x <- 5
#y <- BoxCox(x, 1)      # 4
#xx <- BoxCoxInv(y, 1)  # 5

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CalcChi2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate chi squared from multinomial data, given the parameters
#  of a lognormal distribution 
 
#  #  CHANGE HISTORY
#  28.08.2020 Start
#
# ==============================================================================

CalcChi2 <- function(breaks, counts, mue, sig)
{
  #  INPUT 
  #  breaks   
  #  counts
  #  mue
  #  sig

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------

  x.n        <- sum(counts)
  counts.n   <- length(counts)
  prop.emp   <- counts/x.n
  cdf        <- plnorm(breaks, meanlog=mue, sdlog=sig)
  prop.the   <- cdf[2:(counts.n+1)]-cdf[1:counts.n]
  counts.the <- x.n * prop.the
  tab <- data.frame(x.lo=breaks[1:counts.n],
                    x.hi=breaks[2:(counts.n+1)],
                    count=counts,
                    prop.emp=prop.emp,
                    prop.the=prop.the,
                    chi2.bin=((counts-counts.the)^2)/counts.the )
  return(tab)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CalcParDistance.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculates distance between observed truncated parameters
#  and estimated truncated parameters resulting
#  from the presumed parent distribution, given the truncation limits
 
#  #  CHANGE HISTORY
#  09.09.2020 Operation restricted to normally distributed data
#  19.08.2020 Start
#
# ==============================================================================

CalcParDistance <- function(theta.parent, theta.trunc.obs, a, b)
{
  #  INPUT
  #  theta.parent      Vector with components 
  #          mue.parent   (presumed) mean of the parent normal distribution
  #          sig.parent   (presumed) sd of the parent normal distribution
  #  theta.trunc.obs   Vector with components 
  #          mue.trunc    observed mean of the truncated distribution 
  #          sig.trunc    observed sd of the truncated distribution 
  #  a             left truncation limit
  #  b             right truncation limit

  #  OUTPUT 
  #  SS     squared distance between 
  #         (mue.trunc.presumed, sig.trunc.presumed) and (mue, sig)
  # ---------------------------------------------------------------------------

  #  Auxiliary terms
  mue.parent <- theta.parent[1] 
  sig.parent <- theta.parent[2] 
  mue.trunc.obs  <- theta.trunc.obs[1]
  sig.trunc.obs  <- theta.trunc.obs[2]

  alpha <- (a - mue.parent) / sig.parent
  beta  <- (b - mue.parent) / sig.parent

  A  <-  dnorm(alpha, mean=0, sd=1)
  B  <-  dnorm(beta, mean=0, sd=1)
  AA <-  pnorm(alpha, mean=0, sd=1)
  BB <-  pnorm(beta, mean=0, sd=1)
  C  <-  (B - A) / (BB - AA)

  #  Calculate the mean of the truncated distribution, given theta.parent
  mue.trunc.est <- mue.parent - sig.parent * C

  #  Calculate the sd of the truncated distribution, given theta.parent
  sig.trunc.est <- sig.parent * sqrt((1 - (beta*B - alpha*A) / (BB - AA) - C^2 ))

  R <- matrix(c(a, b,
                mue.parent, sig.parent,
                mue.trunc.est, sig.trunc.est,
                mue.trunc.obs, sig.trunc.obs), byrow=TRUE, ncol=2)
  d  <- c(mue.trunc.obs - mue.trunc.est, sig.trunc.obs-sig.trunc.est)
  #cat("\n [CarlcParDist] R\n")
  #print(R)              
  #cat("\n [CarlcParDist] d\n")
  #print(d)

  SS <- sqrt(d %*% d)

  return(SS)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#theta.parent <- c(140.011879, 2.671882)
#theta.trunc  <- c(139.933013, 1.326838)
#a <- 138
#b <- 142

#HW <- CalcParDistance(theta.parent, theta.trunc, a, b)
#print(HW)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_CalcPrev_V2.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate the sizes of xl, xc, xr using the truncation limits and the 
#  estimated parameters of xc. Values inside the truncation interval 
#  are assumed to be completely from xc.

#  Calculation on the untransformed basis using cdf.PN(x, lambda, mue, sigma)
 
#  #  CHANGE HISTORY
#
# 02.08.2022 Too large estimate of c.n corrected
# 06.12.2020 Start
# ============================================================================

CalcPrev <- function(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
                     x.tr.lo, x.tr.hi, lambda, mue, sigma)
{
  #  Central component of the mixture

  xc.tr.p <- (cdf.PN(x.tr.hi, lambda, mue, sigma) - 
              cdf.PN(x.tr.lo, lambda, mue, sigma))

  # May get zero if provided parameters are too wrong. Prevent overflow.
  xc.tr.p   <- max(1.e-6, xc.tr.p) 

  c.n    <- x.tr.n/xc.tr.p

  #  c.n may get larger than x.n if provided parameters are too wrong. 
  #  Set appropriate limit. 02.08.2022

  # c.n    <- min(c.n, x.n-x.lt.tr.n-x.ge.tr.n)
  c.n    <- min(c.n, x.n)
  prev.c <- c.n/x.n

  #  Size of the left component of the mixture. This is assumed to exist only 
  #  < x.tr.lo, because the truncation interval is assumed to contain only
  #  values from xc.
 
  l.n    <- x.lt.tr.n - c.n * cdf.PN(x.tr.lo, lambda, mue, sigma)
  prev.l <- l.n/x.n

  #  Size of the right component of the mixture. This is assumed to exist only 
  #  >= x.tr.hi, because the truncation interval is assumed to contain only
  #  values from xc.
  r.n    <- x.ge.tr.n - c.n * (1-cdf.PN(x.tr.hi, lambda, mue, sigma))
  prev.r <- r.n/x.n

  neg.prev.sum <- sum(min(prev.l, 0), min(prev.r, 0)) 

  l.n <- unname(l.n)
  c.n <- unname(c.n)
  r.n <- unname(r.n)
  prev.l <- unname(prev.l)
  prev.c <- unname(prev.c)
  prev.r <- unname(prev.r)
  neg.prev.sum <- unname(neg.prev.sum)

  return(c(l.n=l.n, prev.l=prev.l, c.n=c.n, prev.c=prev.c,  
           r.n=r.n, prev.r=prev.r, neg.prev.sum=neg.prev.sum))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#
#  Test for version 2. Data used is FileNo 10021, M only, x.n = 42400

#x.n <- 42400
#lambda.c <- 0
#mue.c    <- 2.928666
#sigma.c  <- 0.3201523

#x.kde.mode  <- 17.04995
#x.lt.kde.mode.n <- 16874    # sum(x < x.kde.mode) 
#x.ge.kde.mode.n <- 25526    # sum(x >= x.kde.mode) 

#x.tmc.mode <- exp(mue.c-sigma.c^2) # 16.88066
#x.lt.tmc.mode.n <- 13971    # sum(x < x.tmc.mode) 
#x.ge.tmc.mode.n <- 28429    # sum(x >= x.tmc.mode) 

#x.tr.lo <-  7.5
#x.tr.hi <- 27.5

#x.tr.n <- 35615     # sum((x.tr.lo <= x) & (x < x.tr.hi))
#x.lt.tr.n <- 85     # sum(x < x.tr.lo)
#x.ge.tr.n <- 6700   # sum(x >= x.tr.hi)

#res <- CalcPrev(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
#                x.tr.lo, x.tr.hi, lambda.c, mue.c, sigma.c)
#res
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CalcR2InSubInt.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate r^2 for subintervals of the QQ plot with sufficiently many bins
#  and cases, and also contining the mode (if requested)
#  Subintervals are defined by the breaks of the collapsed histogram.
 
#  #  CHANGE HISTORY
#  22.08.2021 Return indices in y.hist, not those from ytab 
#  21.08.2021 remove y.hist.mode.idx from call, is not reasonable here
#             replace by y.kde.mode
#  09.08.2021 y.poly.qq renamed to y.red.qq (as in calling function)
#  08.08.2021 bins.n.min.act in call replaced by x.tr.bins.min
#  02.12.2020 x.tr.n added to output
#  24.11.2020 Estimation of xc.n, xl.n, xr.n removed
#  16.11.2020 data.frames constructed with stringsAsFactors=FALSE
#  26.09.2020 Call changed
#  31.08.2020 More criteria for subinterval selection added
#  29.08.2020 Start
#
# ==============================================================================

CalcR2InSubInt <- function(x.hist, y.hist, y.red.qq, lambda, 
                           x.tr.bins.min,  x.tr.prop.min, x.tr.prop.max, 
                           y.kde.mode,
                           df.est,df.con,
                           x.Q1, x.Q2, RL1.p, RL2.p,
                           l.fact, p.fact, r.fact, w.fact,
                           print.log.message)
{
  #  INPUT 
  #  y.hist
  #  y.red.qq
  #  xlabel   label for data
  #  subtitle
  #  x.tr.bins.min
  #  figA

  #  OUTPUT 
  #  
  # ===========================================================================

  # --------------------------------------------------------------------------
  #if (print.log.message) { cat("%%%   CalcR2InSubInt   Start\n") }
  # --------------------------------------------------------------------------

  y0.tab <- hist.table(y.hist$breaks, y.hist$counts) 
 
  #  Due to possibly dealing with a reduced dataset, some lines in y.tab might 
  #  be empty 

  bins.n <- nrow(y0.tab)
  y.tab <- matrix(NA, nrow=bins.n, ncol=ncol(y0.tab))
  colnames(y.tab) <- colnames(y0.tab)

  y.tab[1, "x.lo"] <- y0.tab[1, "x.lo"]

  j <- 0
  for (i in 1:bins.n)
  {  
    if (y0.tab[i, "count"] > 0 )
    {
      j <- j + 1

      y.tab[j, "x.hi"]  <-  y0.tab[i, "x.hi"]
      y.tab[j, "count"] <-  y0.tab[i, "count"]
      y.tab[j, "prop"]  <-  y0.tab[i, "prop"]

      if (j > 1) { y.tab[j, "x.lo"] <- y.tab[j-1, "x.hi"] } 
    }
  }  

  #  Remove empty lines
  y.tab  <- y.tab[!is.na(y.tab[ ,"count"]), ]
  bins.n <- nrow(y.tab)

  #  As we are (possibly) working with a reduced table, the former 
  #  y.hist.mode.idx does not work anymore. Find the interval in y.tab that
  #  contains y.kde.mode and use that for checking "mode.ok"
 
  y.tab.mode.idx <- which( (y.tab[ , "x.lo"] <= y.kde.mode) &
                           (y.kde.mode < y.tab[ , "x.hi"]) )

  #  Calculate partial r2 subintervals having
  #  - at least x.tr.bins.min bins
  #  - at least a proportion of x.tr.prop.min values
  #  - at most  a proportion of x.tr.prop.max values
  #  - contains the mode
  #  Note: y.hist may be a reduced version of the true histogram.
  #  Proportions like x.tr.prop.min are used here, but applied
  #  possibly to the reduced histogram. Absolute sizes of xl, xc, xr
  #  are not computed here. Must be done in the calling programme.

  bins.n <- nrow(y.tab)
  x.n    <- sum(y.hist$counts)

  Pr1.names <- c("bin.start", "bin.end", "bins.n", "prop",  
                 "x.tr.lo", "x.tr.hi", "y.tr.lo", "y.tr.hi", 
                 "mue", "sig", "r2", "opt.crit", 
                 "bins.n.ok", "prop.ok", "mode.ok")
  Pr1 <- matrix(NA, nrow=bins.n * (bins.n-1) / 2, 
               ncol=length(Pr1.names))
  colnames(Pr1) <- Pr1.names
  Pr1 <-  data.frame(Pr1, stringsAsFactors=FALSE)
  Pr1.n <- 0

  #  Intervals to consider have limits now given in y.tab!
  #  y.tab has y.lo, y.hi in one line, different from y.breaks
  #  There may be no interval in y.tab with the required proportion range.
  #  Therefore calculate all descriptions and take intervals with the required
  #  description, if they exist, and take the closest otherwise


  for (bin.start in 1:(bins.n-x.tr.bins.min+1) )
  { y.tr.lo <- y.tab[bin.start, "x.lo"]

    for (bin.end in (bin.start+x.tr.bins.min-1):bins.n)
    { 
      Pr1.n                   <- Pr1.n + 1
      Pr1[Pr1.n, "bin.start"] <- bin.start 
      Pr1[Pr1.n, "bin.end"]   <- bin.end
      Pr1[Pr1.n, "y.tr.lo"]   <- y.tr.lo 
      Pr1[Pr1.n, "y.tr.hi"]   <- y.tab[bin.end, "x.hi"]
      ok <- (y.tab[bin.start, "x.lo"] <= y.red.qq$y) &
            (y.red.qq$y < y.tab[bin.end, "x.hi"] )
      Pr1[Pr1.n, "x.tr.n"]    <- sum(ok)   # may refer to a reduced histogram
      Pr1[Pr1.n, "prop"]      <- Pr1[Pr1.n, "x.tr.n"]/x.n
    }
  }
  Pr1 <- Pr1[1:Pr1.n, ]

  # Check subintervals for having required properties
  Pr1[ , "bins.n.ok"] <- (Pr1[ , "bin.end"]  - Pr1[ , "bin.start"] + 1) >= 
                         x.tr.bins.min
  Pr1[ , "prop.ok"]   <- (x.tr.prop.min <= Pr1[ , "prop"]) & 
                         (Pr1[ , "prop"] <= x.tr.prop.max)
#  Pr1[ , "mode.ok"]   <- (Pr1[ , "bin.start"] <= y.tab.mode.idx) &
#                         (y.tab.mode.idx <= Pr1[ , "bin.end"] )

#  Test: criterion removed
  Pr1[ , "mode.ok"]   <- rep(TRUE, times=Pr1.n)


  #  Remove intervals with unacceptable # of bins or not containing the mode
  Pr1 <- Pr1[Pr1[ , "bins.n.ok"] & Pr1[ , "mode.ok"], ]
 
cat("\n [CalcR2InSubInt] Filtering TIs by prop in ", x.tr.prop.min, 
x.tr.prop.max, "\n")
#print(Pr1[Pr1[ , "prop.ok"] ,] )

  #  If intervals fulfilling the prop requirement exist, use these
  #  If not, take the interval closest to the requirement

  if (sum(Pr1[ , "prop.ok"]) == 0)
  {
    # No fitting intervals
    cat("\n [CalcR2InSubInt] Pr1, insufficient prop, prior to selection \n")
    print(Pr1)

    prop.dist <- min(abs(Pr1[ , "prop"] - x.tr.prop.min),
                     abs(Pr1[ , "prop"] - x.tr.prop.max) )
    idx <- which.min(prop.dist)
    Pr1[idx, "prop.ok"] <- TRUE

    cat("\n [CalcR2InSubInt] Pr1, insufficient prop, after selection \n")
    print(Pr1)
  } 

  #  Filter by prop.ok
  Pr1 <- Pr1[Pr1[ , "prop.ok"], ]
  Pr1 <- RowToMatrix(Pr1) 

  Pr1.n <- nrow(Pr1)

cat("\n [CalcR2InSubInt] Pr1 used for regression \n")
print(Pr1[  , ] )
  
  for (ip in 1:Pr1.n) 
  { 
    ok <- (Pr1[ip,"y.tr.lo"] <= y.red.qq$y) &
          (y.red.qq$y < Pr1[ip ,"y.tr.hi"] )

    #  Calculate regression per accepted interval

    y.lm <- lm(y.red.qq$y[ok] ~ y.red.qq$x[ok])
    y.lm.sum <- summary(y.lm)
    Pr1[ip, "mue"]      <- y.lm$coefficients[1]
    Pr1[ip, "sig"]      <- y.lm$coefficients[2]
    Pr1[ip, "r2"]       <- y.lm.sum$adj.r.squared
  }       

  #  Check for reasonable regression lines
  #  @ tbd

  #  Sort by descending r2
  if (Pr1.n > 1) { Pr1 <- Pr1[order(-Pr1[ ,"r2"]), ] }

  #  Sort by ascending opt.crit
  #Pr1 <- Pr1[order(Pr1[ ,"opt.crit"]), ]

  Pr1[ , "bins.n"] <- Pr1[ , "bin.end"] - Pr1[ , "bin.start"] + 1
  Pr1 <- data.frame(Pr1, rank.r2=1:nrow(Pr1), stringsAsFactors=FALSE)
  Pr1 <- data.frame(Pr1, y.RL1=qnorm(RL1.p, mean=Pr1[ ,"mue"], sd=Pr1[ ,"sig"]),
                    stringsAsFactors=FALSE)
  Pr1 <- data.frame(Pr1, y.RL2=qnorm(RL2.p, mean=Pr1[ ,"mue"], sd=Pr1[ ,"sig"]),
                    stringsAsFactors=FALSE)

cat("\n [CalcR2InSubInt] Pr1 with regression \n")
print(Pr1[  , ] )

  #  Index values in Pr1 refer to ytab. Generate Pr2 with index values with 
  #  regard to y.hist
  
  Pr2 <- Pr1
  y.hist.bin.start <- rep(NA, times=Pr1.n)
  y.hist.bin.end   <- rep(NA, times=Pr1.n)

  for (ipr in 1:Pr1.n)
  { 
    y.hist.bin.start[ipr] <- which(abs(Pr1[ipr, "y.tr.lo"] - y.hist$breaks) <
                                                             1.e-5) 
    y.hist.bin.end[ipr]   <- which(abs(Pr1[ipr, "y.tr.hi"] - y.hist$breaks) <
                                                             1.e-5)
  }
  Pr2[ , "bin.start"] <- y.hist.bin.start 
  Pr2[ , "bin.end"]   <- y.hist.bin.end - 1

  # --------------------------------------------------------------------------
  #if (print.log.message) { cat("%%%   CalcR2InSubInt   End\n") }
  # --------------------------------------------------------------------------

  return(list(tab=y.tab, Pr1=Pr1, Pr2=Pr2))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#
# F_cdf.PN.R

# ----------------------------------------------------------------------       
cdf.PN <- function(x,lambda,mue,sigma,fastnull=1.e-10)
{ 
  #  Verteilungsfunktion der Power-Normalverteilung. 
  #  mue, sigma gelten auf der Y-Skala
  #  siehe Freeman, S. 767

  #  23Jul2018 WW K = 0 abgefangen

  cdf  <- rep(0, times=length(x))
  xgt0 <- x > fastnull
  T <- 1/(lambda*sigma) + mue/sigma
  Z <- (BoxCox(x[xgt0],lambda)-mue)/sigma
 
  if (lambda > fastnull)
  { 
    K <- pnorm(T,mean=0,sd=1)

    #  Overflow protection
    if (K > fastnull)
    {
      cdf[xgt0]  <- (1/K) * (pnorm(Z) - pnorm(-T))
    } else
    {
      cdf[xgt0]  <- 1
    }
  }
  if (abs(lambda) <= fastnull)
  { # LNV
    
    cdf <- plnorm(x,meanlog=mue,sdlog=sigma)
  }
  if (lambda < -fastnull)
  { 
    K <- pnorm(-T,mean=0,sd=1)

    #  Overflow protection
    if (K > fastnull)
    {
      cdf[xgt0]  <- (1/K) * pnorm(Z)
    } else
    {
      cdf[xgt0] <- 0 
    }
  }
  return(cdf)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#lambda <- 0
#lambda <- -0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#x.P50  <- q.PN(0.50, lambda, mue, sigma)
#x      <-  seq(1,20)
#x      <-  seq(77600, 77700, by=5)
#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf, x.cdf")
#print(cbind(x, x.pdf, x.cdf))

#cdf.PN(135.8146, 0.9, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.99, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.999, 140.0014, 2.646370)  
#cdf.PN(135.8146, 0.9999, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.00, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.0001, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.001, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.01, 140.0014, 2.646370)  
#cdf.PN(135.8146, 1.1, 140.0014, 2.646370)  

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Ceiling.R

# 08.01.2020

# =============================================================================
Ceiling <- function(x,unit)
{ # Rundet nach oben auf die nächsthöhere unit, sofern x nicht 
  # ganzzahlig Vielfaches der unit ist
  vielfaches <- trunc(x/unit)
  x.ceiling <- (vielfaches+(abs(x-vielfaches*unit) > 1e-15)) * unit
  return(x.ceiling)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CheckRounding.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Check data for inconsistent rounding. 
#  Do all possible values of the last digit after rounding have similar
#  frequencies? Produce a frequency table of the last digit. The rounding unit
#  is determined from the smallest difference between different values.
#  Rounding to powers of ten is assumed (not e.g. to units of 0.2)
 
#  #  CHANGE HISTORY
#  15.10.2020 Problems due to small numerical errors removed
#  08.06.2020 Start
# ============================================================================

CheckRounding <- function(x, nearlyzero=1.e-10)
{
  #  INPUT 
  #  x                input data (vector)
  #  nearlyzero       absolute smaller values are treated as zero

  #  OUTPUT 
  #  last.dig.table   frequency table of the last digit
  # ==========================================================================

  #  Find rounding unit

  x             <- sort(x)
  x.diff        <- diff(x, lag=1)
  x.diff        <- x.diff[x.diff > nearlyzero]
  x.diff.min    <- min(x.diff)
  x.diff.min.lg <- log10(x.diff.min)
  round.uni.data    <- 10^Round(x.diff.min.lg,1)

  #  Determine last digit 

  last.dig0      <- x - (round.uni.data*10) * 
                    Floor((x+nearlyzero)/(round.uni.data*10), 1)

  #  last.dig values may be different due to representation problems
  #  (differences after the tenth position)
  #  Reduce to a reasonable number of digits
  last.dig1 <- signif(last.dig0, digits=3)
  last.dig  <- Round(last.dig0, round.uni.data)

  D <- data.frame(x, last.dig0, last.dig1, last.dig)

  D <- D[order(D[ ,"last.dig"]), ]

  last.dig.table.values <- (0:9) * round.uni.data 

  temp.exp              <- data.frame(val=last.dig.table.values,
                                      exp=rep(length(x)/10, times=10))

  last.dig.table.obs    <- table(last.dig)

  #  There may be no rounding in the data
  if (all(last.dig.table.obs == 1))
  { 
    cat("\n  +++ No rounding detected  +++ \n") 
    last.dig.table <- NA
  } else
  {
    val.obs               <- as.numeric(names(last.dig.table.obs))
    val.obs               <- Round(val.obs, round.uni.data)
    temp.obs              <- data.frame(val=val.obs,
                                        obs=c(last.dig.table.obs))

    #  Merge observed and expected counts
    last.dig.table        <- merge(temp.exp, temp.obs, all.x=TRUE)
    last.dig.table[is.na(last.dig.table[ ,"obs"]) ,"obs"] <- 0

    #  Calculate chi square
    last.dig.table <- data.frame(last.dig.table, 
                                 chi2=(last.dig.table[ ,"obs"] - 
                                       last.dig.table[ ,"exp"])^2 / 
                                       last.dig.table[ ,"exp"])

    colnames(last.dig.table) <- c("Value", "Expected count", "Observed count", 
                                  "chi2")
    last.dig.table <- last.dig.table[ ,
                                c("Value", "Observed count", "Expected count",
                                  "chi2")]

    return(last.dig.table)
  } 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test

#x <- c(21, 23, 25, 25, 28)
#x <- x / 100

#ldt <- CheckRounding(x) 
#ldt
#ldt <- CheckRounding(seq(79.75, 79.85, by=0.01)) 

#  roughly glucose
#x0   <- rlnorm(1000, meanlog=4.94, sdlog=0.34)
#x0   <- seq(79.75, 79.85, by=0.01)
#x0   <- sort(x0)

#x1   <- round(x0, digits=2)
#x2   <- Round(x0, 0.01)
#xxx  <- data.frame(x0, x1, x2, diff=abs(x1-x2)) 
#print(xxx)
#print(max(xxx[ ,"diff"]))

#ldt <- CheckRounding(x1) 
#ldt
#1-pchisq(sum(ldt[ ,"chi2"]), df=9)

#ldt <- CheckRounding(x2) 
#ldt
#1-pchisq(sum(ldt[ ,"chi2"]), df=9)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ***************************************************************************
#  F_chi2.PNV.tr.lms_V2.R

# ===========================================================================
#  Zielfunktion für gestutzte gruppierte Schätzung von lambda, MW und SD
#  einer PNV 

#  23.08.2022 rt.fact added to call
#  14.05.2021 Call changed, processing of fixed theta components added
#  13.01.2012 Call changed, parameters for RL@.pen added, q.fact now 
#             alphabetically
#  06.12.2020 Call changed, x.kde.mode removed, x.lt.mode, x.ge.mode changed 
#  05.12.2020 Call changed, x.kde.mode added
#  24.11.2020 Call changed, x.lt.mode, x.ge.mode added
#  17.09.2019 Faktor für Strafterm "negative Prävalenz" eingeführt
#  28.01.2018
# ===========================================================================

chi2.PNV.tr.lms <- function(ttheta.ini.est, RB, 
                            ttheta.ini.fix, idx.fix, idx.est,
                            x.hist, ilo, ihi,
                            x.lt.tr.n, x.ge.tr.n,
                            l.fact, p.fact, r.fact, rt.fact, w.fact,
                            x.Q1, x.Q2, RL1.p, RL2.p,
                            df.est,df.con, opt.crit.only, fastnull 
                            )
{
  #  LF auf Grundlage der bedingten Dichte einer gestutzten PNV
  #  Diskrete Version
  #  Berechnet wird -log-LF (wegen Nutzung durch nlm() )

  #  ttheta    auf R transformierte Parameter der PNV
  #  RB        Matrix der Randbedingungen an theta
  #  idx.fix   contains indices of fixed components in theta
  #  idx.est   contains indices of theta components tat are being estimated
  #  x.hist    Histogramm mit Komponenten $breaks, $counts, lineare Skala
  #  x.tr.lo   unterer truncation point, lineare Skala
  #  x.tr.hi   oberer truncation point, lineare Skala
  #            truncation points müssen breaks sein!
  #  df.est    Freiheitsgrade für geschätzte Parameter
  #  opt.crit.only Steuert die Ausgabe, nur opt.crit oder mehr, siehe chi2trunc 
  #  fastnull  Ersatzwert für kleine Erwartungswerte
  #  l.fact    Faktor für Gewichtung "Länge des truncation Intervalls"
  #  p.fact    Faktor für Strafterm "negative Prävalenz" in chi2-Funktion
  #  r.fact    Faktor für Strafterm ...
  #  w.fact    Faktor für Strafterm "zu hohe Schätzung" in chi2-Funktion

  # Compose the full parameter
  ttheta <- Compose(ttheta.ini.fix, ttheta.ini.est, idx.fix, idx.est)

  theta  <- R2o(ttheta,RB)
  lambda <- theta[1] 
  mue    <- theta[2]
  sigma  <- theta[3]

  #  chi^2-Anteil aus der gestutzten Verteilung berechnen
  
  #  Old version before 02.01.2022
  #res <- chi2trunc(x.hist,x.tr.lo,x.tr.hi,
  #                 x.lt.tr.n, x.ge.tr.n,
  #                 x.Q1, x.Q2, RL1.p, RL2.p,
  #                 lambda,mue,sigma,df.est, df.con, 
  #                 l.fact, p.fact, r.fact, w.fact, 
  #                 opt.crit.only=opt.crit.only,fastnull=fastnull)

  #  New version since 02.01.2022
  res <- chi2trunc(x.hist,ilo, ihi,
                   x.lt.tr.n, x.ge.tr.n,
                   x.Q1, x.Q2, RL1.p, RL2.p,
                   lambda,mue,sigma,df.est, df.con, 
                   l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                   rt.fact=rt.fact, w.fact=w.fact, 
                   opt.crit.only=opt.crit.only,fastnull=fastnull)

  return(res)  
}

# ****************************************************************************
#  F_chi2trunc_V14.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  CHANGE HISTORY
#  25.09.2022  Runs.test around 0, not around median
#  07.09.2022  Goodness of fit calculated conditional on TI count
#  02.09.2022  Runs test by WaWo()
#  15.08.2022  Runs test for residuals (ex crit6) included in opt.crit
#  05.08.2022  w.fact does no more act on the central contribution
#  01.03.2022  Default values for penalty terms taken from seg030
#  06.02.2022  Action against too small expected values deactivated
#  23.01.2022  Input x.hist can contain too small sized bins, create
#              x.histc with collapsed bins in lo and hi respecting the 
#              TI boundaries 
#  02.01.2022  Call changed, new organization
#  18.12.2021  Function oc() calculates modified optimality criterion
#  17.08.2021  Safety action against numerical instability when looking
#              for TI limit indices
#  13.01.2021  x.Q1, x.Q2, RL1.p, RL2.p, r.fact added to call for penalty 
#              calculation. Sequence of @.fact now alphabetically.
#  06.12.2020  Calculation of xl.n, xr.n and prevalence changed back to old
#              approach using the TI limits, because usage of mode, though 
#              theoretically correct, gives sometimes nonsense, because
#              of rounding effects when counting empirical frequencies
#              below / above the mode
#  04.12.2020  New calculation of xl.n, xr.n and prevalence new, now separation
#              by mode 
#              Strategy x.tr.eff changed
#              'check' removed
#  19.11.2020  Optimality criterion set back to chi2/df, because standard test
#              data creates p values = 0 for initial values from qqw,
#              therefore no progress in nlm. Attempt to "normalize" chi2 to get
#              more independent from df? 
#  29.03.2020  Optimality criterion renamed to opt.crit
#              chi2.per.df removed from output
#  06.01.2020  Prevalences l, c, r introduced
#  01.01.2020  prev.tmc.pen new defined
#              chi2.per.df  new defined
# ============================================================================

chi2trunc <- function(x.hist,ilo,ihi,
                      x.lt.tr.n, x.ge.tr.n,
                      x.Q1, x.Q2, RL1.p, RL2.p,
                      lambda,mue,sigma,
                      df.est,df.con,
                      l.fact=0.0, p.fact=0.0, r.fact=0, rt.fact=1,
                      w.fact=6,
                      penalty.details=TRUE, 
                      opt.crit.only=FALSE,fastnull=1.e-10)
{  
   #  Calculates ch-square goodness-of-fit and runs test plus penalty terms
   #  for a truncated Power Normal Distribution
   #  V14 / 07.09.2022
   #
   #  EINGABE
   #  x.hist     Histogramm von x (gesamter Datensatz)
   #  ilo        leftmost bin in TI
   #  ihi        rightmost bin in TI  
   #  x.lt.tr.n  number of values <   lower TI limit
   #  x.ge.tr.n  number of values >=  upper TI limit
   #  lambda )
   #  mue    )   Parameter der angepassten PNV 
   #  sigma  )
   #  df.est     Freiheitsgrade für geschätzte Parameter (sind nicht 
   #             zwangsläufig 3!)
   #  l.fact     Gewichtungsfaktor "Länge des truncation intervals" - deactivated
   #  p.fact     Gewichtungsfaktor für Strafterm negative Prävalenz
   #  r.fact     Gewichtungsfaktor für Strafterm RLs outside total RLs
   #  rt.fact    Gewichtungsfaktor für Strafterm runs test
   #  w.fact     Gewichtungsfaktor für Strafterm zu hohe Vorhersage
   #  opt.crit.only  TRUE: nur der chi2-Wert für Optimierung wird zurückgegeben
   #             FALSE: vollständige Information, siehe AUSGABE
   #
   #  AUSGABE
   #  If opt.crit.only=TRUE
   #  Wert des Optimierungskriteriums (kleiner ist besser)

   #  If opt.crit.only=FALSE
   #  Full information, see return statement at end of function

  # ========================================================================= 

  #  Truncation interval has limits from breaks. This makes earlier variables
  #  superfluous.
  #  The TI contains bins ilo:ihi in counts.
  #  The corresponding interval is [breaks[ilo], breaks[ihi+1]) except for
  #  ihi == counts.n, for which    [breaks[ilo], breaks[ihi+1]] holds.
  #  Earlier variables have been replaced
  #  x.tr.eff.lo         now  x.tr.lo
  #  x.tr.eff.lo.idx     now  ilo
  #  x.tr.eff.hi         now  x.tr.hi
  #  x.tr.eff.hi.idx     now  ihi+1
  #  x.tr.hist.counts.n  now tr.counts.n

  #  Gesamtzahl Intervalle
  counts.n <- length(x.hist$counts)

  #  Gesamtzahl breaks
  breaks.n <- counts.n + 1

  #  Gesamtzahl Fälle
  x.n <- sum(x.hist$counts)

  #  Counts in the truncation interval 
  x.tr.n <- sum(x.hist$counts[ilo:ihi])

  #  Potential problem: ilo and ihi define too few bins -  then do not execute
  #  Moved to construction of TI.can!
  tr.counts.n <- ihi - ilo + 1

# {  indentation no more needed
    #  TI part of the histogram
    #  Corresponding breaks
    x.tr.hist.breaks <- x.hist$breaks[ilo:((ihi+1))]

    #  Corresponding counts
    x.tr.hist.counts <- x.hist$counts[ilo:((ihi+1)-1)]

    #  Check counts

    x.n2 <- x.tr.n + x.lt.tr.n + x.ge.tr.n
    if ( abs(x.n2 - x.n) > 1.e-5)
    {
      cat("\n [chi2trunc_V6]  Count partition incorrect",
          "\n              x.lt.tr.n =", x.lt.tr.n,         
          "\n              x.tr.n    =", x.tr.n,         
          "\n              x.ge.tr.n =", x.ge.tr.n,         
          "\n              x.n       =", x.n,
          "\n        +++++ Execution stops +++++",
          "\n")

      stop("++++ Execution stopped in chi2trunc +++++")                 
    }

    cp <- CalcPrev(x.n, x.tr.n, x.lt.tr.n, x.ge.tr.n, 
                   x.hist$breaks[ilo], x.hist$breaks[ihi+1], 
                   lambda, mue, sigma)
    prev.l.hat <- cp["prev.l"]
    prev.c.hat <- cp["prev.c"]
    prev.r.hat <- cp["prev.r"]
    xl.n.hat   <- cp["l.n"]
    xc.n.hat   <- cp["c.n"]
    xr.n.hat   <- cp["r.n"]

    #  Unconditional probabilities for xc per bin, 
    #  given theta=c(lambda, mue, sigma)

    xc.hist.cdf <- cdf.PN(x.hist$breaks,lambda,mue,sigma)
    xc.hist.pdf <- diff(xc.hist.cdf)

    #  If theta is extremely wrong or the domain extremely large, 
    #  components of xc.hist.pdf may be
    #  extremely small and thereby cause an overflow later.  Prevent that.   

    too.small <- (xc.hist.pdf<1.e-10)

    xc.hist.pdf[too.small] <- 1.e-10     #  before 09.09.2022
    #  xc.hist.pdf[too.small] <- 0          #  test 09.09.2022

    #  xc.hist.pdf is needed only in the TI, much less outside

    #  Change against V13: expected counts are conditional on the breaks
    #  range. No expansion of first and last bins.

    xc.hist.pdf            <- xc.hist.pdf/sum(xc.hist.pdf)

    #  Check sum of probabilities
    if ( abs(sum(xc.hist.pdf) - 1) > 1.e-5)
    {
      cat("\n [chi2trunc]  Calculation of bin probabilities incorrect",
          "\n        +++++ Execution stops +++++",
          "\n")

      stop("++++ Execution stopped in chi2trunc +++++")                 
    }

    #  Table of empirical probabilities
    tab <- Hist.table(x.hist$breaks, x.hist$counts)

    #  Add indicators for below / in / above TI
    tab <- data.frame(tab, is.lo=rep(FALSE, times=counts.n))
    tab <- data.frame(tab, is.TI=rep(FALSE, times=counts.n))
    tab <- data.frame(tab, is.hi=rep(FALSE, times=counts.n))
    if (ilo > 1) 
    { 
      tab[1:(ilo-1), "is.lo"]  <- TRUE  
      is.lo                    <- tab[ , "is.lo"]
    }

    tab[ilo:ihi, "is.TI"]  <- TRUE
    is.TI                  <- tab[ , "is.TI"]

    if (ihi < counts.n) 
    { 
      tab[(ihi+1):counts.n, "is.hi"] <- TRUE  
      is.hi                          <- tab[ , "is.hi"]
    }
   
    #  Add theoretical probs for xc (the presumed nonpath distribution),
    #  conditional on xc.n.hat. This is the proportion of xc expected to 
    #  appear in a bin.   
    #  Calculate expected count from xc for each bin
    tab <- data.frame(tab, xc.hist.n.hat=xc.n.hat * xc.hist.pdf)

    #  Calculate observed-expected for each bin. 
    diff.i <- tab[ , "count"] - tab[ , "xc.hist.n.hat"]

    tab    <- data.frame(tab, diff.i=diff.i )

    #  Calculate chi2 contribution for each bin. 
    #  Bin-wise values only used in central part.
    #  Modification 3 01.03.2022:
    #  p.fit    is the chi2 goodness of fit, calculated from unweighted chi2
    #           contributions chi2.i in TI (larger is better)
    #  p.rt     is the p value of the runs test (Wald-Wolfowitz test) in the
    #           TI  (larger is better)
    #  opt.crit is calculated from unweighted chi2 contributions in TI plus
    #           the weighted chi2 contributions from below/ above the TI
    #           The oc is calculated from these by oc2(), which accounts for
    #           the number of bins in the TI 
    chi2.i <- tab[ , "diff.i"]^2/tab[ , "xc.hist.n.hat"]
    tab    <- data.frame(tab, chi2.i)

    #  Additional to *central* table: runs test for random sequence of 
    #  differences. Ex crit6 computed in seg105. Now computation here, 
    #  result goes into component $res.

    #  Deactivated 25.09.2022
    #diff.median <- median(tab[is.TI, "diff.i"])
    #signs     <- 1.0 * (tab[is.TI, "diff.i"] > diff.median) # > more plausible
    
    #  Activated 25.09.2022
    signs     <- 1.0 * (tab[is.TI, "diff.i"] > 0) # > more plausible

    # runstest <- runs.test(signs, exact=FALSE, alternative="two.sided")

    wwt <- WaWo(signs, alternative="less")

    #  No test is done here, runstest.p gets part of opt.crit

    #  Summarizing statistics from the TI

    chi2.TI     <- sum(tab[is.TI, "chi2.i"]) 
    chi2.TI.df  <- tr.counts.n - df.est - df.con 
    
    #  Next term (unweighted!) is used as goodness of fit (p.fit)
    chi2.TI.p  <- 1-pchisq(chi2.TI, chi2.TI.df) #  = p.fit
    rt.runs    <- wwt$runs
    rt.p       <- wwt$p.value            # rt.p! name later changed to p.rt

    # Express (asymptotic) runs test test statistic as chi2 quantity (1 df)
    rt.z       <- wwt$z
    #  Only negative z are a problem, therefore component only for negative z
    rt.chi2    <- rt.z^2 * (rt.z < 0)
    rt.pen     <- rt.fact * rt.chi2

    #  Central table completed, calculate (sub-) sums

    # ------------------------------------------------------------------------ 
    #  From below the TI
    if (ilo == 1)
    { 
      # No data below TI
      chi2.lo     <- NA
      chi2.lo.w   <- NA
      chi2.lo.df  <- NA
    }  else
    { 
      #  Change 24.01.2022: chi^2 contribution not bin-wise but over 
      #  complete lo range
      xc.hist.n.hat.lo <- sum(tab[is.lo, "xc.hist.n.hat"])
      dif.lo           <- x.lt.tr.n - xc.hist.n.hat.lo

      chi2.lo          <- (dif.lo)^2 / xc.hist.n.hat.lo
      chi2.lo.w        <- w.fact * pchisq(chi2.lo, 1)
      if (dif.lo >= 0) {  chi2.lo.w <- 0 }
      chi2.lo.w        <- chi2.lo * chi2.lo.w
    }

    #  From above the TI
    if (ihi == counts.n)
    { 
      chi2.hi     <- NA
      chi2.hi.w   <- NA
      chi2.hi.df  <- NA
    }  else
    {
      #  Change 24.01.2022: chi^2 contribution not bin-wise but over 
      #  complete hi range
      xc.hist.n.hat.hi <- sum(tab[is.hi, "xc.hist.n.hat"])
      dif.hi           <- x.ge.tr.n - xc.hist.n.hat.hi

      chi2.hi          <- (dif.hi)^2 / xc.hist.n.hat.hi
      chi2.hi.w        <- w.fact * pchisq(chi2.hi, 1)
      if (dif.hi >= 0) {  chi2.hi.w <- 0 }
      chi2.hi.w        <- chi2.hi * chi2.hi.w
    }

    #  Optimization criterion. Uses goodness of fit in TI plus further
    #  penalty terms.
    
    #  All weighted chi2 contributions are used to calculate  opt.crit
    #  chi2.path.w  penalty from bad fit outside TI (below and above)
    #  chi2.TI      gof in the TI
    #  prev.l.hat.pen  penalty for negative prevalence left
    #  prev.r.hat.pen  penalty for negative prevalence right
    #  RL1.pen         penalty for RL1 estimated left of total RL1
    #  RL2.pen         penalty for RL1 estimated right of total RL2

    chi2.path.w     <- sum(c(chi2.lo.w, chi2.hi.w), na.rm=TRUE)   

    #  Penalty for estimated pathological prevalence < 0
    #  Approach 01.01.2020
    prev.l.hat.pen <- p.fact * (prev.l.hat < 0) * 
                      abs(prev.l.hat) * chi2.TI

    prev.r.hat.pen <- p.fact * (prev.r.hat < 0) * 
                        abs(prev.r.hat) * chi2.TI

    #  Penalty for RLs outside [x.Q1, x.Q2], where x.Q@ are the RL1.p and
    #  RL2.p quantiles of the complete stratum

    x.RL1 <- q.PN(RL1.p, lambda, mue, sigma)
    x.RL2 <- q.PN(RL2.p, lambda, mue, sigma)
    RL.diff  <- x.RL2 - x.RL1
  
    #  If theta is extremely wrong, RL.diff may become zero. Prevent overflow
    RL.diff  <- max(RL.diff, fastnull)
    RL1.pen  <- r.fact * max(0, (x.Q1-x.RL1)/RL.diff )
    RL2.pen  <- r.fact * max(0, (x.RL2-x.Q2)/RL.diff )

    #  Df for chi2.total.w: only bins in the TI are relevant.
    #  In the case of a correct fit no contribution arises outside the TI

    #  Penalty for non-random fluctuation of residuals (diff.i). Needs
    #  chi2.total.w.df = chi2.TI.df
    # rt.pen <- rt.fact * chi2.TI.df*(1-runstest.p)   # until 23.08.22
    # rt.pen <- rt.fact * (1-runstest.p)   # since 24.08.22
    # rt.pen defined above

    #  Weighting factor for TI length
    #  Deactivated 02.09.2022
    #l.pen <- l.fact *(x.n/x.tr.n)^2

    #  Combine fit contributions and penalties as input to oc / oc2
    #  runstest.p added 16.08.2022

    #  Weighting factors for penalty terms applied where penalties are defined

    #cat("\n chi2.TI       ", chi2.TI,
    # "\n chi2.path.w   ", chi2.path.w,
    # "\n prev.l.hat.pen", prev.l.hat.pen,
    # "\n prev.r.hat.pen", prev.r.hat.pen,
    # "\n RL1.pen       ", RL1.pen,
    # "\n rt.pen        ", rt.pen,
    # "\n")

    chi2.total.w  <- chi2.TI + chi2.path.w + 
                     prev.l.hat.pen + prev.r.hat.pen +
                     RL1.pen + RL2.pen +
                     rt.pen

    #  opt.crit uses basically p.fit. But if the fit of the curve is too bad
    #  the iteration may stop prematurely because p.fit is too close to zero
    #  and does not change any more for other theta. Therefore:

    #  New definition of opt.crit (18.12.2021)
    #  Uses essentially p.fit, but switches *smoothly* to 
    #  a non-horizontal straight line, if p.fit is too close to zero
    #  1 df added for runs test

    opt.crit <- oc2(chi2.total.w, chi2.TI.df+1)

    #  Apply weight for size of truncation interval, longer is better
    #  opt.crit gets larger for smaller TI size
    #  Deactivated 02.09.2022

    # opt.crit <- opt.crit * (1 + l.pen)

    if (FALSE)
    { 
      cat("\n [chi2trunc_V14] Penalty details",
          "\n   Actual parameters              ", lambda, mue, sigma,
          "\n   x.RL1, x.RL2                   ", x.RL1, x.RL2,
          "\n   chi2, df, p.fit from TI        ", chi2.TI, chi2.TI.df, chi2.TI.p,
          "\n   chi2, df, w.fact from below TI ", chi2.lo.w, 1, w.fact,
          "\n   chi2, df, w.fact from above TI ", chi2.hi.w, 1, w.fact,
          "\n   chi2.path.w                    ", chi2.path.w,
          "\n   p.fact, prev.l.hat.pen, prev.r.hat.pen ", 
                p.fact, prev.l.hat.pen, prev.r.hat.pen,
          "\n   r.fact, RL1.pen, RL2.pen       ", r.fact, RL1.pen, RL2.pen,
          "\n   rt.runs, rt.p                  ", rt.runs, rt.p,
          "\n   rt.fact, rt.chi2, rt.pen       ", rt.fact, rt.chi2, rt.pen,
          "\n   l.fact  (l.pen deactivated)    ", l.fact, 
          "\n   chi2.total.w, chi2.total.w.df  ", chi2.total.w, chi2.TI.df,
          "\n   opt.crit                       ", opt.crit,
          "\n") 
    }

    #  Report result
    if (opt.crit.only)
    { 
      return(unname(opt.crit)) 
    }  else
    { 
      #  use old names on the right to avoid changes in the rest of TMC

      opt.crit       <- unname(opt.crit) 
      chi2.total     <- unname(chi2.total.w) 
      chi2.total.df  <- unname(chi2.TI.df) 
      chi2.lo.w      <- unname(chi2.lo.w)
      chi2.hi.w      <- unname(chi2.hi.w)
      chi2.trun      <- unname(chi2.TI) 
      chi2.trun.df   <- unname(chi2.TI.df) 
      chi2.trun.p    <- unname(chi2.TI.p) 
      chi2.path.w    <- unname(chi2.path.w) 
      prev.l.hat     <- unname(prev.l.hat) 
      prev.c.hat     <- unname(prev.c.hat) 
      prev.r.hat     <- unname(prev.r.hat) 
      prev.l.hat.pen <- unname(prev.l.hat.pen) 
      prev.r.hat.pen <- unname(prev.r.hat.pen) 
      xl.n.hat       <- unname(xl.n.hat) 
      xc.n.hat       <- unname(xc.n.hat) 
      xr.n.hat       <- unname(xr.n.hat)
      p.rt           <- unname(rt.p)        # new name!
      rt.pen         <- unname(rt.pen)
 

      tab.lo <- NA
      if (ilo > 1) { tab.lo <- tab[is.lo, ] }

      tab.hi <- NA
      if (ihi < counts.n) { tab.hi <- tab[is.hi, ] }

      return(list(res=c(opt.crit=opt.crit,
                        chi2.total=chi2.total, 
                        chi2.total.df=chi2.total.df, 
                        chi2.lo.w=chi2.lo.w,                    
                        chi2.trun=chi2.trun, 
                        chi2.trun.df=chi2.trun.df,
                        chi2.trun.p=chi2.trun.p,
                        chi2.hi.w=chi2.hi.w,
                        chi2.path.w=chi2.path.w,
                        prev.l.tmc=prev.l.hat,
                        prev.c.tmc=prev.c.hat,
                        prev.r.tmc=prev.r.hat,
                        prev.l.tmc.pen=prev.l.hat.pen,
                        prev.r.tmc.pen=prev.r.hat.pen,
                        xl.n.tmc=xl.n.hat,
                        xc.n.tmc=xc.n.hat,
                        xr.n.tmc=xr.n.hat,
                        p.rt=p.rt,
                        rt.pen=rt.pen
                        ),
                  tab=tab[is.TI, ],
                  tab.lo=tab.lo,
                  tab.hi=tab.hi)
            )
    }  
  # } this indentation no more needed
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#p.lo <- cdf.PN( 67.5, 0.0648018, 8.787397, 5.20248)
#p.hi <- cdf.PN(108.5, 0.0648018, 8.787397, 5.20248)
#x.n * (p.hi - p.lo)

#p.lo <- cdf.PN(188.5, 0.0648018, 8.787397, 5.20248)
#p.hi <- cdf.PN(193.5, 0.0648018, 8.787397, 5.20248)
#x.n * (p.hi - p.lo)
# ****************************************************************************
#  F_CIQuant.LNV.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Asymptotics 1-alpha confidence interval for quantiles of the 
#  log-normal distribution

#  See R.J.Serfling, Approximation theorems of mathematical statistics,
#  Wiley 1980, S. 77, corollary B
 
#  #  CHANGE HISTORY
#  15.05.2017 Start
# ============================================================================

CIQuant.LNV <- function(n,p,Q,mue,sigma,alpha=0.05)
{ 
  #  INPUT
  #  n       Size of the sample from which to calculate the quantiles
  #  p       Probability associated with the quantile
  #  Q       Quantile
  #  mue     )  Parameters of the log-normal distribution
  #  sigma   )
  #  alpha   1 - confidence probability

  #  OUTPUT
  #  Q.ci.lo )
  #  Q.ci.hi ) confidence interval limits
  # ==========================================================================

  #  pdf at the quantile
  f.Q <- dlnorm(Q,meanlog=mue,sdlog=sigma)

  #  Asymptotic variance of the quantile
  Q.var <- p*(1-p) / (f.Q^2 * n)

  #  Asymptotic standard deviation of the quantile
  Q.sd  <- sqrt(Q.var)

  #  Asymptotic confidene interval
  z       <- qnorm(1-alpha/2)
  Q.ci.lo <- Q - z*Q.sd
  Q.ci.hi <- Q + z*Q.sd

  return(c(Q.ci.lo,Q.ci.hi))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_CIQuant.PNV.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Asymptotics 1-alpha confidence interval for quantiles of the 
#  Power normal distribution

#  See R.J.Serfling, Approximation theorems of mathematical statistics,
#  Wiley 1980, S. 77, corollary B
 
#  #  CHANGE HISTORY
#  09.08.2017 Start
# ============================================================================

CIQuant.PNV <- function(n,p,Q,lambda,mue,sigma,alpha=0.05, fastnull=1.e-10)
{ 
  #  INPUT
  #  n       Size of the sample from which to calculate the quantiles
  #  p       Probability associated with the quantile
  #  Q       Quantile
  #  lambda  )
  #  mue     )  Parameters of the power normal distribution
  #  sigma   )
  #  alpha   1 - confidence probability

  #  OUTPUT
  #  Q.ci.lo )
  #  Q.ci.hi ) confidence interval limits
  # ==========================================================================

  #  pdf at the quantile
  f.Q <- pdf.PN(Q,lambda,mue,sigma,fastnull=fastnull)

  #  Asymptotic variance of the quantile
  Q.var <- p*(1-p) / (f.Q^2 * n)

  #  Asymptotic standard deviation of the quantile
  Q.sd  <- sqrt(Q.var)

  #  Asymptotic confidene interval
  z       <- qnorm(1-alpha/2)
  Q.ci.lo <- Q - z*Q.sd
  Q.ci.hi <- Q + z*Q.sd

  return(c(Q.ci.lo,Q.ci.hi))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_CollapseBins_V5.R
# 
#  (c) wwosniok@math.uni-bremen.de
#
# =========================================================================

CollapseBins <- function(breaks, counts, P.lo.idx, P.hi.idx, n.per.bin.min)
{
  #  Collapse Histogram bins such that all bins have at least a size of
  #  n.per.bin.min. Inner bins (indices P.lo.idx:P.hi.idx) remain unchanged.
  #
  #  INPUT
  #  breaks    breaks from initial histogram
  #            Should not produce zero counts at the beginning or the end
  #            of the data range. Such cells will not be removed. 
  #  counts    counts from initial histogram 
  #  P.lo.idx  index (in counts) of the leftmost inner bin 
  #  P.hi.idx  index (in counts) of the rightmost inner bin 
  #  n.per.bin.min mininumim required count in a collapsed histogram

  #  OUTPUT
  #  hist.coll  the collapsed histogram object (class = "histogram")

  #  CHANGE HISTORY 
  #  17.09.2022  Start  V5
  # =========================================================================

  # cat("\n %%%% [CollapseBins V5]  Start\n")

  x.n    <- sum(counts)
  ht.raw <- Hist.table(breaks, counts)

  #  Reduce to relevant part
  ht.raw   <- ht.raw[ ,c("x.lo", "x.hi", "count")]
  ht.raw.n <- nrow(ht.raw)

  #  Remove outer zero counts to avoid unneccessary large outer bins
  count.eq.0 <- TRUE
  while (count.eq.0)
  {
    count.eq.0 <- (ht.raw[1 , "count"] == 0)
    if (count.eq.0)  ht.raw <- ht.raw[-1, ]
  }
  
  # Correct P.lo.idx / P.hi.idx
  P.lo.idx <- P.lo.idx - (ht.raw.n - nrow(ht.raw))
  P.hi.idx <- P.hi.idx - (ht.raw.n - nrow(ht.raw))

  ht.raw.n <- nrow(ht.raw)

  count.eq.0 <- TRUE
  while (count.eq.0)
  {
    count.eq.0 <- (ht.raw[nrow(ht.raw) , "count"] == 0)
    if (count.eq.0)  ht.raw <- ht.raw[-nrow(ht.raw), ]
  }
  ht.raw.n <- nrow(ht.raw)

  # -------------------------------------------------------------------------
  #  Collapsed tables starts with unchanged inner range of ht.raw

  ht.col <- ht.raw[P.lo.idx: P.hi.idx, ]

  #  Rows in collapsed table
  ht.col.n <- nrow(ht.col)

  # -------------------------------------------------------------------------
  #  Collapse bins left of inner range

  #  Nothing to collapse if P.lo.idx == 1

  if (P.lo.idx > 1)
  {
    #  We are going downwards!
    istart <- P.lo.idx-1            # istart left, iend right end of interval
    iend   <- istart

    go.on    <- TRUE

    #  Count in the  collapsed intervals
    count.col <- sum(ht.raw[istart:iend, "count"])
    size.ok   <- count.col >= n.per.bin.min 

    while (go.on)
    {
      while ( (!size.ok) & (istart > 1))
      { 
        # Increase input range downwards if possible
        istart <- istart - 1
        count.col <- count.col + ht.raw[istart, "count"]
        size.ok   <- (count.col >= n.per.bin.min) 
      }

      if (size.ok)
      { 
        #  Requested count reached
        ht0 <- c(x.lo=ht.raw[istart, "x.lo"], 
                 x.hi=ht.raw[iend,   "x.hi"],
                 count=count.col)

        # ht.col.n == 0 cannot occur 
        ht.col   <- rbind(ht0, ht.col)
        ht.col.n <- ht.col.n + 1

        #  Prepare next search, if possible
        if (istart > 1)
        { #  There are more rows in ht.raw
          istart    <- istart - 1
          iend      <- istart
          count.col <- sum(ht.raw[istart:iend, "count"])
          size.ok   <- (count.col >= n.per.bin.min)
          go.on     <- TRUE
        }  else

        { #  No more rows in ht.raw, finished on left
          go.on <- FALSE
        }
      }  else      # size.ok?
      { 
        #  At begin of ht.raw, but requested count not reached, 
        #  enlarge last complete new row

        ht0 <- c(x.lo=ht.raw[istart, "x.lo"], 
                 x.hi=ht.raw[iend+1, "x.hi"],
                 count=ht.col[1, "count"] + count.col)

        # Replace first row in ht.col 
        ht.col[1, ] <-  ht0

        go.on <- FALSE
      }   # if (!last.row
    }     # while(go.on
  }       # if (P.hi.idx < ht.raw.n


  # -------------------------------------------------------------------------
  #  Collapse bins right of inner range

  #  Nothing to collapse if P.hi.idx == ht.raw.n

  if (P.hi.idx < ht.raw.n)
  {
    istart <- P.hi.idx + 1
    iend   <- istart

    go.on <- TRUE

    count.col <- sum(ht.raw[istart:iend, "count"])
    size.ok   <- count.col >= n.per.bin.min 

    while (go.on)
    {
      while ( (!size.ok) & (iend < ht.raw.n) )
      { 
        # Increase input range if possible
        iend <- iend + 1
        count.col <- count.col + ht.raw[iend, "count"]
        size.ok   <- count.col >= n.per.bin.min 
      }

      if (size.ok)
      { 
        #  Requested count reached
        ht0 <- c(x.lo=ht.raw[istart, "x.lo"], 
                 x.hi=ht.raw[iend,   "x.hi"],
                 count=count.col)

        ht.col <- rbind(ht.col, ht0)
        ht.col.n <- ht.col.n + 1

        #  Prepare next search, if possible
        if (iend < ht.raw.n)
        { istart    <- iend + 1
          iend      <- istart
          count.col <- sum(ht.raw[istart:iend, "count"])
          size.ok   <- (count.col >= n.per.bin.min)
          go.on     <- TRUE

        } else    # size.ok?
        { #  No more rows in ht.raw, finished on right
          go.on <- FALSE
        }
      }  else  # size.ok?
      { 
        #  Requested count not reached, enlarge last complete new row
        ht0 <- c(x.lo=ht.col[ht.col.n, "x.lo"], 
                 x.hi=ht.raw[iend,   "x.hi"],
                 count=ht.col[ht.col.n, "count"] + count.col)

        # Replace last row in ht.col 
        ht.col[ht.col.n, ] <-  ht0
  
        go.on <- FALSE
      }   # if (!last.row
    }     # while(go.on
  }       # if (P.hi.idx < ht.raw.n
                   

  #  Output as histogram
  breaks.col <- unname(c(ht.col[ , "x.lo"], tail(ht.col[ , "x.hi"], 1)))
  breaks.col.n <- length(breaks.col)
  counts.col <- unname(ht.col[ , "count"])

  hist.col <- list(breaks=breaks.col,
                    counts=counts.col,
                    density=counts.col/(diff(breaks.col)*x.n),
                    mids=(breaks.col[1:(breaks.col.n-1)]+
                          breaks.col[2:breaks.col.n])/2,
                    xname="x",
                    equidist=length(unique(diff(breaks.col)))==1)
  class(hist.col) <- "histogram"

  # cat("\n %%%% [CollapseBins V5]  End\n")

  return(hist.col)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (FALSE)
{
#stop(" ERROR: test sequence activated ")
  source("../func/F_Hist.table.R")

  # .........................................................................
  #  Simple test data
  counts <- c(0, 0, 10, 20, 100, 200, 100, 20, 10, 0, 0)
  breaks <- 1000 + seq(1:(length(counts)+1))

  ht.raw <- Hist.table(breaks, counts)

  ht.raw <- ht.raw[ ,c("x.lo", "x.hi", "count")]

  P.lo.idx <- 5           # idx refers to ht.raw
  P.hi.idx <- 7           # do not collapse inner range
  n.per.bin.min <- 135

hc <- CollapseBins(breaks, counts, P.lo.idx, P.hi.idx, n.per.bin.min)

Hist.table(hc$breaks, hc$counts)
}

if (FALSE)
{
stop(" ERROR: test sequence activated ")

  set.seed(73271)
  x      <- rnorm(100, mean=10.5, sd=1)
  #x.hist <- hist(x, breaks=seq(6, 14, by=1))
  x.hist <- hist(x, breaks=seq(4, 16, by=1))

  breaks <- x.hist$breaks
  counts <- x.hist$counts
  counts.cdf <- cumsum(counts/sum(counts))

  #   Corresponding indices in counts
  P.lo <- 0.20
  P.hi <- 0.80
  P.lo.idx <- min(which(counts.cdf > P.lo))
  P.hi.idx <- max(which(counts.cdf < P.hi))

  n.per.bin.min <- 25

  CB <- CollapseBins(breaks, counts, P.lo.idx, P.hi.idx, n.per.bin.min)

  cat("\n Original histogram\n")
  print(Hist.table(breaks, counts))
  cat("\n Collapsed histogram\n")
  print(Hist.table(CB$breaks, CB$counts))
}

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_color.palette.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Extract the component theta.est from theta
#  No checks for reasonable vector size and contents for time reasons!
 
#  CHANGE HISTORY
#  30.07.2021 Start
# ==============================================================================

color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

#color.palette(5)

#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ***************************************************************************
# F_CompileTabcn.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Compile output tables for TMC_seg100

#  CHANGE HISTORY

#  03.11.2021 Start 
# ===========================================================================
CompileTabcn <- function(meth.list, tabc.stra, tabn.stra,    
               tabc.npa, tabc.qqw, tabc.bha, tabc.tmc, tabc.tmu, tabc.tuk,
               tabn.npa, tabn.qqw, tabn.bha, tabn.tmc, tabn.tmu, tabn.tuk)
{
  #  Combine essential results from all methods for this stratum
  #  tab.stra has a dummy line from definition in the beginning here

  if ("npa" %in% meth.list)
  {
    tabc.stra <- rbind(tabc.stra, tabc.npa)
    tabn.stra <- rbind(tabn.stra, tabn.npa)
  }

  if ( "qqw" %in% meth.list )
  {
    tabc.stra <- rbind(tabc.stra, tabc.qqw)
    tabn.stra <- rbind(tabn.stra, tabn.qqw)
  }

  if ("bha" %in% meth.list)
  {
    #  @@@ subpopulation auswählen  - not yet updated
    tabc.stra["bha", "RL1"]     <- tabc.bha[["subpop"]]["comp01", "RL1.2"]
    tabc.stra["bha", "RL2"]     <- tabc.bha[["subpop"]]["comp01", "RL2.2"]
    tabn.stra["bha", "RL1"]     <- tabn.bha[["subpop"]]["comp01", "RL1.2"]
    tabn.stra["bha", "RL2"]     <- tabn.bha[["subpop"]]["comp01", "RL2.2"]
  }

  if ("tmc" %in% meth.list)
  {
    tabc.stra <- rbind(tabc.stra, tabc.tmc)
    tabn.stra <- rbind(tabn.stra, tabn.tmc)
  }

  if ("tmu" %in% meth.list)
  {
    tabc.stra <- rbind(tabc.stra, tabc.tmu)
    tabn.stra <- rbind(tabn.stra, tabn.tmu)
  }

  if ("tuk" %in% meth.list)
  {
    tabc.stra <- rbind(tabc.stra, tabc.tuk)
    tabn.stra <- rbind(tabn.stra, tabn.tuk)
  }            

  # tabc.stra, tabn.stra are complete, remove the dummy lines
  tabc.stra <- tabc.stra[-1, ] 
  tabn.stra <- tabn.stra[-1, ] 

  # Get rid of old row names
  row.names(tabc.stra) <- NULL
  row.names(tabn.stra) <- NULL

  return(list(tabc.stra=tabc.stra, tabn.stra=tabn.stra))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



# ****************************************************************************
#  F_Compose.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Composes the full parameter vector theta from components theta.est and
#  theta.fix
#  No checks for reasonable vector size and contents for time reasons!
 
#  CHANGE HISTORY
#  12.05.2020 Start
# ==============================================================================

Compose <- function(theta.fix, theta.est, idx.fix, idx.est)
{
  #  INPUT 
  #  theta.fix   the fixed component of the parameter vector. May be NULL.
  #  theta.est   the estimated component of the parameter vector
  #  idx.fix     the indices of the fix component(s) in theta. May be NULL.

  #  OUTPUT 
  #  theta.est   vector with length as idx.est
  # ---------------------------------------------------------------------------

  if (is.null(theta.fix))
  { # No fixed component
    theta <- theta.est
  }  else
  { # Fixed component(s) exist(s)
    theta.n <- length(theta.fix) + length(theta.est)
    theta   <- rep(NA, times=theta.n)
    theta[idx.fix] <- theta.fix
    theta[idx.est] <- theta.est
  } 
  
  return(theta)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#t.fix <- NULL
#t.est <- c(10, 20, 30)
#idx.fix <- NULL
#idx.est <- c(1,2,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))

#t.fix <- 10
#t.est <- c(20, 30)
#idx.fix <- 1
#idx.est <- c(2,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))


#t.fix <- 10
#t.est <- c(20, 30)
#idx.fix <- 2
#idx.est <- c(1,3)
#print(Compose(t.fix, t.est, idx.fix, idx.est))



# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Drift.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Generate a plot for drift identification for the monthly median
 
#  CHANGE HISTORY
#  13.04.2022 ef introduced
#  26.03.2021 Start
# ==============================================================================

Drift <- function(date, x, RL1, RL2, figA, xlabel, titletext, ef, 
                  alpha=0.05, gamcol="lightblue",
   subt="The blue confidence range should overlap the grey range at all times")
{
  #  INPUT
  #  date  date, format  "2021-03-26" generated by
  #        as.Date("26.03.2021", format = "%d.%m.%Y"),
  #  x     time series data
  #  RL1   Lower 95% reference limit ( 2.5%) 
  #  RL2   Upper 95% reference limit (97.5%) 
  #  figA  plot window
  #  xlabel label of x
  #  titletext text for plot title

  #  OUTPUT 
  # ---------------------------------------------------------------------------

  #  Median per unique time point

  x.median <- aggregate(x, by=list(date), median)
  colnames(x.median) <- c("date", "median")
 
  #  Fit a gam to the monthly medians, if enough points
  if (nrow(x.median) >= 12)
  {
    #cat("\n Drift estimation by GAM: median = time + s(time)\n")

    #  To simplify writing:
    time  <- as.numeric(x.median[ , "date"])
    x.med <- x.median[ , "median"]

    x.med.gam <- gam(x.med ~ time + s(time))
    #print(x.med.gam)
    #print(summary(x.med.gam))    
  
    x.med.gam.pred <- predict.gam(x.med.gam,se.fit=TRUE)
    # print(x.med.gam.pred)

    x.med.gam.pred.fit <- unlist(x.med.gam.pred$fit)
    x.med.gam.pred.se  <- unlist(x.med.gam.pred$se.fit)
    zalpha             <- qnorm(1-alpha/2)
    x.med.gam.pred.ciu <- x.med.gam.pred.fit - zalpha*x.med.gam.pred.se
    x.med.gam.pred.cio <- x.med.gam.pred.fit + zalpha*x.med.gam.pred.se

    # print(cbind(x.med,x.med.gam.pred.ciu,x.med.gam.pred.cio))

    #  Total median (of all data)
    #x.median.total <- median(x)

    x.median.total <- median(x.med)
    #  Permissible deviation around x.median.total
    x.median.total.pI <- PermissibleDiff(RL1, RL2, x.median.total, ef)

    #  Simplify writing
    xmt.lo <- x.median.total.pI[1, "xi.lo"]
    xmt.hi <- x.median.total.pI[1, "xi.hi"]

    yplotmin <- min(c(x.med, xmt.lo, x.med.gam.pred.ciu)) 
    yplotmax <- max(c(x.med, xmt.hi, x.med.gam.pred.cio)) 

    # Span the plot window
    dev.set(figA)
    bringToTop(figA)

    plot(c(x.median[1 ,"date"], tail(x.median[ ,"date"], 1)),
         c(yplotmin, yplotmax), 
         type="p", pch=3, col="transparent",
         xlab="Date", ylab=xlabel, main=titletext,
         sub=subt,
         cex.sub=0.8) 

    #  Permissible range of the total median
    polygon(c(x.median[1, "date"], tail(x.median[ , "date"], 1),
              tail(x.median[ , "date"], 1),  x.median[1, "date"] ), 
            c(xmt.lo, xmt.lo, xmt.hi, xmt.hi),
            col="gray95", border="gray95") 
    
    #  CI of drift
    polygon(c(time, rev(time)),
            c(x.med.gam.pred.ciu, rev(x.med.gam.pred.cio)), 
            col=gamcol, border="lightblue") 
    
    #  Monthly medians   
    points(x.median[ ,"date"], x.median[ ,"median"], type="p", col="black",
           pch=3)

    #   Total median
    abline(h=x.median.total, col="red")

    #   Total median
    abline(h=c(xmt.lo, xmt.hi), col="red", lty=2)

    #  GAM fit to monthly medians, + CI
    lines(time, x.med.gam.pred.fit, col="blue",lty=1)
    lines(time, x.med.gam.pred.ciu, col="blue",lty=2)
    lines(time, x.med.gam.pred.cio, col="blue",lty=2)

    #  Check for drift
    drift.likely <- any(x.med.gam.pred.ciu > xmt.hi) |
                    any(x.med.gam.pred.cio < xmt.lo)

  }   else
  {  #  No data for drift analysis
     drift.likely <- NA
  }

  return(drift.likely)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test: see Test_DriftAnalysis.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# F_DynAgeLimits.R
#  19.04.2022
#
#  (c) wwosniok@math.uni-bremen.de
# ============================================================================

DynAgeLimits <- function(age, x.n.min)
{ #  Find dynamic age limits auch that each class has at least size x.n.min
  #  Age classes do not overlap.
  #  
  #  INPUT 
  #  age         vector of ages 
  #
  #  OUTPUT
  #  age.class   Matrix of age class limits. All limits belong to the 
  #              interval
  # =========================================================================

  age.table <- table(age)
  age.table.n <- length(age.table)
  age.table.names <- as.numeric(names(age.table))

  age.class <- matrix(NA, nrow=age.table.n, ncol=3)  
  colnames(age.class) <- c("lo", "hi", "count")

  iage  <- 1   # index in age.class
  ilo   <- 1    # indices in age.table
  count <- 0

  while(ilo <= age.table.n)
  { 
    ihi <- ilo
    count <- sum(age.table[ilo:ihi])
    while ( (count < x.n.min) & (ihi < age.table.n) )
    {
      ihi <- ihi + 1
      count <- count + age.table[ihi]
    }
    age.class[iage, "lo"] <- age.table.names[ilo]
    age.class[iage, "hi"] <- age.table.names[ihi]
    age.class[iage, "count"] <- count

    iage  <- iage + 1
    ilo   <- ihi + 1
    count <- 0 
  }
  
  age.class <- age.class[!is.na(age.class[ ,"lo"]), ]
  age.class.n <- nrow(age.class)

  if (age.class.n > 1)
  {
    if (age.class[age.class.n, "count"] < x.n.min)
    { #  Last subset is too small
      age.class[age.class.n-1, "hi"] <- age.class[age.class.n, "hi"] 
      age.class[age.class.n-1, "count"] <- age.class[age.class.n-1, "count"] + 
                                         age.class[age.class.n, "count"] 
      age.class.n <- age.class.n - 1
      age.class   <- age.class[1:age.class.n, ]
    }
  }
  
  age.class <- RowToMatrix(age.class)
  
  return(age.class)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#source("../func/F_RowToMatrix.R")

#age <- c(
#       rep(17, times=4),
#       rep(18, times=5),
#       rep(19, times=6),
#       rep(20, times=7),
#       rep(21, times=8),
#       rep(22, times=9),
#       rep(23, times=10),
#       rep(24, times=11),
#       rep(25, times=12), 
#       rep(26, times=2), 
#       rep(27, times=1), 
#       rep(28, times=1)
#       )

#x.n.min <- 3
#ac <- DynAgeLimits(age, x.n.min)
#print(ac)


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_EstParentDist.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Estimate the parent normal distribution if only limits, mean and sd of a 
#  truncated subset are known. See @@@ for a method description. 

#  CHANGE HISTORY
#  Changes in functions are documented in the functions!

#  02.12.2020 RLs, CI(RL) addidionally calculated here, call changed
#  28.09.2020 ok range changed to { , )
#  08.09.2020 tmu installed, execution via do.tmu.  
#  30.03.2020 eps removed from parameter list
#  17.09.2019 p.fact as additional parameter
#  15.03.2019 Installation from directory ProgWW17_RH

# ============================================================================

EstParentDist <- function(x.tr.lo, x.tr.hi, x.tr, theta.parent.ini,
                          alpha, lambda, RL1.p, RL2.p)
{ 
  # x.tr.lo            lower truncation limit
  # x.tr.hi            upper truncation limit
  # x.tr               truncated data set
  # theta.parent.ini   initial guess for parent parameters

  # ==========================================================================

  # --------------------------------------------------------------------------
  # if (print.log.message) { cat("%%%   F_EstParentDist  Start\n") }
  # --------------------------------------------------------------------------

  #  Estimates of mue and sigma from the truncated, to ND transformed 
  #  interval
  x.tr.n   <- length(x.tr)
  x.tr.mea <- mean(x.tr)
  x.tr.sd  <- sd(x.tr)
 
  #  Improve parent parameter estimation iteratively  

  TMU.nlm <- nlm(CalcParDistance, theta.parent.ini, 
                 theta.trunc=c(x.tr.mea, x.tr.sd), 
                 a=x.tr.lo, b=x.tr.hi)
    
  #  Modified estimates for mue and sigma
  tab.tmu <- c(mue.tmu=TMU.nlm$estimate[1], sigma.tmu=TMU.nlm$estimate[2], 
              distance.tmu=TMU.nlm$minimum, rc.tmu=TMU.nlm$code, 
              iter.tmu=TMU.nlm$iterations)

  #  Meaning of "code" ("rc"):
  #  1: relative gradient is close to zero, current iterate is 
  #     probably solution.
  #  2: successive iterates within tolerance, current iterate is 
  #     probably solution.
  #  3: last global step failed to locate a point lower than estimate. 
  #     Either estimate is an approximate local minimum of the function 
  #     or steptol is too small.
  #  4: iteration limit exceeded.
  #  5: maximum step size stepmax exceeded five consecutive times. 
  #     Either the function is unbounded below, becomes asymptotic 
  #     to a finite value from above in some direction or stepmax is too small.

  # -------------------------------------------------------------------------
  #  Calculate RLs
  x.RL1.tmu <- unname(q.PN(RL1.p, lambda, tab.tmu["mue.tmu"], 
                          tab.tmu["sigma.tmu"]))
  x.RL2.tmu <- unname(q.PN(RL2.p, lambda, tab.tmu["mue.tmu"], 
                          tab.tmu["sigma.tmu"]))
  tab.tmu <- c(tab.tmu, x.RL1.tmu=x.RL1.tmu, x.RL2.tmu=x.RL2.tmu)

  # -------------------------------------------------------------------------
  #  Add asymptotic confidence intervals
  xc.RL1.tmu.CI <- CIQuant.PNV(x.tr.n, RL1.p, 
                               tab.tmu["x.RL1.tmu"],
                               lambda, 
                               tab.tmu["mue.tmu"], 
                               tab.tmu["sigma.tmu"],
                               alpha=alpha)
  tab.tmu <- c(tab.tmu, x.RL1.cilo=unname(xc.RL1.tmu.CI[1]),
                        x.RL1.cihi=unname(xc.RL1.tmu.CI[2]))

  xc.RL2.tmu.CI <- CIQuant.PNV(x.tr.n, RL2.p, 
                               tab.tmu["x.RL2.tmu"],
                               lambda, 
                               tab.tmu["mue.tmu"], 
                               tab.tmu["sigma.tmu"],
                               alpha=alpha)
  tab.tmu <- c(tab.tmu, x.RL2.cilo=unname(xc.RL2.tmu.CI[1]),
                        x.RL2.cihi=unname(xc.RL2.tmu.CI[2]))

  # --------------------------------------------------------------------------
  # if (print.log.message) { cat("%%%   F_EstParentDist  End\n") }
  # --------------------------------------------------------------------------

  return(tab.tmu)
}   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Expand.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Expand tied data (due to rounding) to fictive individual data. Counts per
#  bin are as in input, position in bin is log-equidistant, with half distance
#  to the margins
 
#  #  CHANGE HISTORY
#  16.02.2021 Special treatment for count = 1 added
#  06.01.2021 Step function as 2nd possibility of expansion added
#  30.12.2020 Start
#
# ==============================================================================

Expand <- function(breaks, counts, step=TRUE)
{
  #  INPUT 
  #  breaks          Jump positions of the step functions
  #                  Must be > 0 
  #                  Reasonable number of breaks must be checked outside 
  #                  this function
  #  counts          # of observations per bin
  #  step            TRUE: x.fict in bin uniformly distributed
  #                  FALSE: x.fict in bin logarithmically distributed
 
  #  OUTPUT 
  #  x.fict          fictive data points, counts per bin as in input,
  #                  position in bin depends on step
  # ---------------------------------------------------------------------------

  #  Allocate fictive cases in each (original) bin
  #  x positions for all data
  bins.n <- length(counts)

  x.fict.exists <- FALSE             # in place of exists()

  if (step)
  { 
    # Use a simple step function
    for (ibin in 1:bins.n)
    {
      if (counts[ibin] > 0)
      { 
        #  Avoid numerical trouble due to very small differences if count=1
         
        if (counts[ibin] == 1)
        {
          x.fict0 <-  (breaks[ibin+1] + breaks[ibin]) / 2
        } else
        {   
          dx <- (breaks[ibin+1] - breaks[ibin]) / counts[ibin]
          x.fict0 <- seq(breaks[ibin]+dx/2, breaks[ibin+1]-dx/2, by=dx)
        }

        if (!x.fict.exists)
        { 
          x.fict        <- x.fict0
          x.fict.exists <- TRUE 
        } else
        {
          x.fict <- c(x.fict, x.fict0)   
        }  
      }
    }
  }   else
  { # Locate points at logarithmic positions  

    for (ilo in 1:(bins.n))
    { ihi <- ilo + 1
      if (counts[ilo] > 0)
      { 
        delta.logx <- (log(breaks[ihi]) - log(breaks[ilo]))/ (counts[ilo]+1)  
        logx.fict  <- seq(log(breaks[ilo])+delta.logx/2,
                        log(breaks[ihi])-delta.logx/2, length.out=counts[ilo])

        if (!x.fict.exists)
        { 
          x.fict        <- exp(logx.fict)
          x.fict.exists <- TRUE 
        } else
        {
          x.fict <- c(x.fict, exp(logx.fict))   
        }
      }
    }
  } 

  return(x.fict) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
# see also Test.Expand.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_Extract.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Extract the component theta.est from theta
#  No checks for reasonable vector size and contents for time reasons!
 
#  CHANGE HISTORY
#  12.05.2020 Start
# ==============================================================================

Extract <- function(theta, idx.est)
{
  #  INPUT 
  #  theta   the full parameter vector
  #  idx.est the indices of components that are estimated (not fix)

  #  OUTPUT 
  #  theta.est   vector with length as idx.est
  # ---------------------------------------------------------------------------

  theta.est <- theta[idx.est]

  return(theta.est)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#t <- c(10, 20, 30)
#idx.est <- c(1,2,3)
#print(Extract(t, idx.est))
#idx.est <- c(1,3)
#print(Extract(t, idx.est))
#idx.est <- 2
#print(Extract(t, idx.est))
#idx.est <- NA
#print(Extract(t, idx.est))
#idx.est <- 4
#print(Extract(t, idx.est))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FC.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  Print values using format specifications as in formatC. 
#  formatC itself does not in all cases deal with NA correctly.

# x  <- c(123, 123.4, 123.45, NA)
# y  <- c(NA, NA, NA)
# formatC(x, width=8, digits=2)  
# " 1.2e+02" " 1.2e+02" " 1.2e+02" "      NA"
# formatC(x, format="f", width=8, digits=2)  
# "  123.00" "  123.40" "  123.45" "      NA"
# formatC(y, format="f", width=8, digits=2)
#Fehler in formatC(y, format = "f", width = 8, digits = 2) : 
#  unsupported type 

# formatC(as.numeric(y), format="f", width=8, digits=2)

#  CHANGE HISTORY
#  11.04.2021 Start
# ============================================================================ 

FC <- function(x, w, d)
{ 
  if (is.na(w)) 
  { Fx <- formatC(as.numeric(x), format="f", digits=d, flag="#") } else
  { Fx <- formatC(as.numeric(x), format="f", width=w, digits=d, flag="#") }

  return(Fx)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# x  <- c(123, 123.4, 123.45, NA)
# y  <- c(NA, NA, NA)
#FC(x,  8, 2)
#FC(x, NA, 2)
#FC(y, NA, 2)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv # *****************************************************************************
#  F_FilterWDay.R
#
#  Extract data from given weekdays
#
#  (c) wwosniok@math.uni-bremen.de
#
#  CHANGE HISTORY
#
#  29.10.2020 Start
#
# --------------------------------------------------------------------------
# INPUT
# wday      vector containing the weekdays of observations as numbers
#          (Monday=0)
# use.wday  vector containing the weekdays to select, 
#           either NA or "all" (==> use all data)
#           or a subset of c("mo", "tu", "we", "th", "fr", "sa", "su")

# OUTPUT
# ok        vector (length as wday) indicating the selected components of wday
# --------------------------------------------------------------------------

FilterWDay <- function(wday, use.wday)
{  
  #  If selection of weekdays is required: do it
  if (!is.na(use.wday) && use.wday!="all")
  { 
    ok <- rep(FALSE, times=length(wday))

    #  Check for valid entries

    for (day in use.wday)
    { 
      # print(day)

      idx <- which(tolower(day) == c("mo", "tu", "we", "th", "fr", "sa", "su"))
      # print(idx) 
      if (length(idx) == 1)
      { 
        ok  <- ok | (wday == (idx-1))
      }  else
      { cat("\n ++++++ Weekday selection ", day, " invalid  ++++++\n") } 
    }
  }  else
  {
    # No selection required
    ok <- rep(TRUE, times=length(wday))   
  }

  #  The calling programme has to handle the situation of no data left
  return(ok)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

#wday <- c(0, 1,1, 2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5,5, 6,6,6,6,6,6)

#use.wday <- c("mo", "su", "WE")

#select <- FilterWDay(wday, use.wday)
#print(data.frame(wday, select))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindBw_V4.R

FindBw <- function(x, kernel.n, round.unit, FigA, subtitle, kdecol,
                    min.bw=0.75, b.fact=1.05, inner.lo=0.025, inner.hi=0.975)
{  
  #  Adjust bandwidth of kernel density estimate such that the
  #  inner part of the density has only 1 maximum
  #
  #  INPUT
  #  x           data
  #  round.unit  rounding applied to the data
  #  kernel.n    number of support points for the kernel 
  #  FigA        window for plotting. Must exist. NA: no plot
  #  b.fact      factor by which bw.adjust is increased per iteration
  #  inner.lo    lower limit of the inner range (quantile of total data)
  #  inner.hi    upper limit of the inner range (quantile of total data)

  #  05.04.2021 temp[ , "lo"] <=  temp[ , "mid"]) allowed for maximum
  #  01.12.2020 Missing new calculation of density() added in loop 
  #  26.09.2020 Colours as function input
  #  15.09.2020 Bandwidth must be >= min.bw*round.unit
  #  13.09.2020 Additional condition: bandwidth must be >= round.unit
  #  22.05.2020 Target changed: smoothing is increased until all maxima
  #             are more than the rounding unit apart (or only 1 maximum left)
  #  21.05.2020 Target changed: up to max.n.acc maxima are accepted, the highest
  #             is taken as mode position
  #  20.05.2020 Start
  # ==========================================================================

  x.kde.mode <- NA
  bw.adjust  <- 1

  #  First try. Used to determine the support and to ensure that 
  #  bw > round.unit
  #  First guess. May have negative support.
  x.kde <- density(x, n=kernel.n, adjust=bw.adjust)

  #  Support points must be > 0 to avoid trouble with BC transformation
  x.from <- min(x.kde$x[x.kde$x>0])
  x.to   <- max(x.kde$x)

  x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)

  if (x.kde$bw < min.bw*round.unit)
  { #  Force bandwidth to be >= min.bw*round.unit
    bw.adjust <- min.bw*round.unit / x.kde$bw
    x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)
  }

  if (!is.na(FigA))
  { 
    dev.set(FigA)
    yplotmax <- max(x.kde$y)
    plot(x.kde, type="l", col="gray", ylim=c(0, yplotmax),
         main=paste("KDE with bandwidth adjustment, bw =", x.kde$bw, 
                    "  adj =", format(bw.adjust, digits=3)),          
         sub=subtitle, cex.sub=0.7)
    axis(1, at=seq(130, 160, by=round.unit), labels=FALSE)
  }

  while (is.na(x.kde.mode))
  {
    x.kde <- density(x, n=kernel.n, adjust=bw.adjust, from=x.from,to=x.to)

    if (!is.na(FigA))
    { 
      lines(x.kde, type="l", col="gray")
    }

    #  Check only the inner range of the density for maxima
    x.kde.cdf <- cumsum(x.kde$y)
    x.kde.cdf <- x.kde.cdf/tail(x.kde.cdf,1)
    Q1Q3      <- (inner.lo <= x.kde.cdf) & (x.kde.cdf < inner.hi)

    x.kde.Q1Q3.x <- x.kde$x[Q1Q3]
    x.kde.Q1Q3.y <- x.kde$y[Q1Q3]
    Q1Q3.n <- sum(Q1Q3)

    if (!is.na(FigA))
    { 
      abline(v=c(x.kde.Q1Q3.x[1], tail(x.kde.Q1Q3.x, 1)), col="gray")
    }

    #  Find maxima
    temp <- cbind(x=x.kde.Q1Q3.x[2:(Q1Q3.n-1)],
                  lo=x.kde.Q1Q3.y[1:(Q1Q3.n-2)],
                  mid=x.kde.Q1Q3.y[2:(Q1Q3.n-1)],
                  hi=x.kde.Q1Q3.y[3:(Q1Q3.n  )])
    is.max <- (temp[ , "lo"] <=  temp[ , "mid"]) & 
              (temp[ , "mid"]>   temp[ , "hi"])
    temp <- cbind(temp, is.max=is.max)

    max.n <- sum(is.max)
    # cat("\n bw.adjust: ", bw.adjust, " number of maxima: ", max.n)

    if (max.n == 0)
    { cat("\n [Find.bw] No maximum found\n")
      print(temp)
      stop("++++++ NO maximum found ++++++")
    }

    if (!is.na(FigA) & max.n > 0)
    { 
      abline(v=temp[is.max, "x"], col="orange", lty=2) 
    }

    #  If > 1 maximum: may be artificial due to strong rounding. If
    #  any distance between maximum positions ~~ round.unit, increase 
    #  bw.adjust, otherwise the solution is reached.
    #  Calculate pairwise distance between maximum positions. Slow method,
    #  but there should not be too may positions.

    if (max.n > 1)
    {
      temp.max <-  temp[is.max, ] 
      #cat("\n More than 1 maximum, maxima positions\n")
      #print(temp.max)

      x.diff <- abs(diff(temp.max[ , "x"]))

      #cat("\n Absolute maxima position differences \n")
      #print(x.diff)

      # Is there any distance ~~ round.unit outside the diagonal?
      d.eq.ru <- (0.55*round.unit <= x.diff) & (x.diff <= 1.45*round.unit)

      if (any(d.eq.ru))
      {  
        # at least 1 pair of maxima with distance ~~ r.u., increase bandwidth
        bw.adjust <- bw.adjust * b.fact
        #  do another round. x.kde.mode is still NA             
      } else   
      { #  Several maxima left, all distances > round.unit. 
        #  Take the one with max density
        #cat("\n\n")

        mode.idx <- which.max(temp[ ,"mid"])
        x.kde.mode <- temp[mode.idx, "x"]  #  this terminates while()

        #cat("\n Local maxima\n")
        #print(temp[is.max, ])
        #cat("\n Mode: ", x.kde.mode, "\n")
      } 
    }  else    
    { 
      # only 1 maximum left, calculate mode
      x.kde.mode <- temp[is.max, "x"]
      #cat("\n Mode: ", x.kde.mode, "\n")
    }
  }    #  while ...

  if (!is.na(FigA))  
  { 
    lines(x.kde, col=kdecol)
    abline(v=x.kde.mode, col="firebrick") 
    legend("topright", c("kde", "Q1, Q3", "mode"), 
           col=c(kdecol, "gray", "firebrick"), lty=1, lwd=2, cex=0.8)   
    abline(v=x.kde.mode, col="firebrick")
  }

  return(list(bw.adjust = bw.adjust,
              x.kde.mode=x.kde.mode,
              x.kde=x.kde))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  See Test_Bandwidth.R 

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_FindHisto_V13.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Find a histogram containing a specified number of bins with specified 
#  minimum size per bin
#  Generate also a reduced version and two smoothed versions
 
#  #  CHANGE HISTORY
#  08.09.2022 Equidistant breaks get non-eqidistant beginning, if DL exists
#  07.02.2022 If DL exists, first break always 2*fastnull
#  27.01.2022 V11 Start. First break always x.val[1]-round.unit/2, also for
#             h > 1
#  29.11.2021 V9: h fixed on basis of maximal # of bins while collapsing the 
#             extremes bins
#  28.11.2021 V8: h fixed depending on x.val.n, derived from V7
#  26.11.2021 V7: Search for h and collapsing in each cycle
#  17.11.2021 Strategy for collapsing changed
#             x.reso.red calculated as x.histo.eqd$mids
#  03.10.2021 Finer resolution in x.reso.red
#  09.03.2021 Treatment of values < DL extended to second break position
#  23.02.2021 Treatment of values < DL added
#  15.02.2021 bins.n2 corrected to integer
#  11.02.2021 Derivation of breaks corrected
#  07.02.2021 Start with new approach: eqd is determined by Freedman-Diaconis,
#             x.hist gets its breaks from quantiles of the kde-cdf
#
# ==============================================================================

FindHisto <- function(x, round.unit, detect.limits.max, x.val, x.kde,
                      n.per.bin.min, bins.n.min, bins.n.max,
                      smooth.hist1, smooth.hist2, red.limit, xlabel, subtitle, 
                      bordercol1, histcol1, bordercol2, histcol2, 
                      bordercol3, histcol3, kdecol, figB, 
                      fastnull)
{
  #  V13
  #  INPUT 
  #  x              data vector
  #  round.unit
  #  x.kde          kde of x
  #  x.kde.mode     mode of x.kde
  #  n.per.bin.min  min number of per bin
  #  red.limit      create a reduced data set if original size is > re.limit 
  #  xlabel         label for analyte 
  #  subtitle       subtitle in FigA 
  #  bordercol1     colour for histogram borders, initial
  #  histcol1       colour for bin area, initial
  #  bordercol2     colour for histogram borders, collapsed
  #  histcol2       colour for bin area, collapsed
  #  bordercol2     colour for histogram borders, collapsed and smoothed
  #  histcol2       colour for bin area, collapsed and smoothed
  #  kdecol         colour for kde
  #  figB           initial and collapsed histogram
  #  figC           collapsed and smoothed histogram

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------
  
  x.n     <- length(x)  
  x.val.n <- length(x.val)
  x.min   <- min(x)
  x.max   <- max(x)

  # ---------------------------------------------------------------------------
  #  Produce a set of reasonable histograms on the original scale
  # --------------------------------------------------------------------------
  #  Calculate equidistant histogram, standard approach using the 
  #  Freedman-Diaconis approach (which, like all others is not suitable for 
  #  rounded data)
  #  Parameter 'freq' not used here.

  #  All histograms have the same first and nearly the same last breaks. 
  #  Note that the data value zero is set to 0.1*round.unit in seg070. 
  #  first: max(x.val[1] - round.unit/2, 0.05*round.unit)
  #  last: for equidistant breaks: such that max(x) is just included
  #        for non-equidistant breaks: max(x) + round/2

  # ..........................................................................
  #  Equidistant histogram. Bin size depends on data between P20 and P80

  #  Maximal equidistant breaks. Make sure that max is included.

  #  Expansion factor
  h            <- 1

  if (is.na(detect.limits.max))
  { # No detection limit
    breaks.start <- max(x.val[1] - round.unit/2, 0.05*round.unit)
  }  else
  {
    # Detection limit exists. Fix breaks[1] at nearly 0 (not exactly 0,
    # makes trouble with some transformations).
    # More action below, after bin width has been determined
    breaks.start  <- 2*fastnull
  }

  breaks.end   <- x.max + 1.01*h*round.unit
  x.breaks0 <- c(breaks.start, 
                 seq(x.val[1]+round.unit/2, breaks.end, by=round.unit) )

  #  Maximal equidistant histogram, bin width > DL is round.unit

  x.hist0.eqd  <- hist(x, breaks=x.breaks0, right=FALSE, plot=FALSE)  

  x.hist0.eqd$cdf  <- cumsum(x.hist0.eqd$counts) / x.n

  #  Identify the range from P20 to P80
  P.lo <- 0.20
  P.hi <- 0.80

  #  Index of P.lo
  P.lo.idx <- min(which(x.hist0.eqd$cdf > P.lo))

  #  Index of P.hi
  P.hi.idx <- max(which(x.hist0.eqd$cdf < P.hi))

  #  As long as there are too many bins in the inner range 
  #  try to increase bin size (but stay above bins.n.min)

  inner.bins.n0 <- P.hi.idx - P.lo.idx + 1
  n.per.bin.ok0 <- all(x.hist0.eqd$counts[P.lo.idx:P.hi.idx] >= n.per.bin.min)

  #  If no improvement posssible (P.lo.idx and P.hi.idx are also needed)
  x.hist.eqd   <- x.hist0.eqd
  inner.bins.n <- inner.bins.n0

  #  Try improvement
  go.on         <- ((bins.n.min < inner.bins.n0) & (!n.per.bin.ok0)) | 
                   (inner.bins.n0 > bins.n.max)  
 
  #cat("\n (1) Bin range for h selection, inner.bins.n0,  n.per.bin.ok0:", 
  #      P.lo.idx, P.hi.idx, inner.bins.n0, n.per.bin.ok0, "\n")  
  #cat("\n (1.5) breaks, counts\n")  
  #print(x.hist0.eqd$breaks[P.lo.idx:P.hi.idx])
  #print(x.hist0.eqd$counts[P.lo.idx:P.hi.idx])

  while ( go.on )
  {
    #  Try reduction of bins.n
    h <- h + 1  

    breaks1.start <- breaks.start 
    breaks1.end   <- x.max + 1.01*h*round.unit
    x.breaks1     <- c(breaks1.start, 
                       seq(breaks1.start+h*round.unit/2, breaks1.end, 
                           by=h*round.unit) )

    #cat("\n (2) h, x.min, x.max, breaks1.start, breaks1.end, tail(x.breaks1,1)\n",
    #  h, x.min, x.max, breaks1.start, breaks1.end, tail(x.breaks1,1),
    #"\n") 

    x.hist1.eqd     <- hist(x, breaks=x.breaks1, right=FALSE, plot=FALSE)  
    x.hist1.eqd$cdf <- cumsum(x.hist1.eqd$counts) / x.n

    #  Identify the range from P20 to P80
    #  Index of P.lo
    P.lo.idx1 <- min(which(x.hist1.eqd$cdf > P.lo))

    #  Index of P.hi
    P.hi.idx1 <- max(which(x.hist1.eqd$cdf < P.hi))

    inner.bins.n1 <- P.hi.idx1 - P.lo.idx1 + 1
    n.per.bin.ok1 <- all(x.hist1.eqd$counts[P.lo.idx1:P.hi.idx1] >= n.per.bin.min)

    #cat("\n (3) Bin range for h selection, inner.bins.n1,  n.per.bin.ok1:", 
    #    P.lo.idx1, P.hi.idx1, inner.bins.n1, n.per.bin.ok1, "\n")  
    #cat("\n (3) # of inner bins, # total bins", P.hi.idx1-P.lo.idx1+1, 
    #    length(x.hist1.eqd$counts), "\n")  
    #cat("\n (3) x.hist1.eqd$counts[P.lo.idx1:P.hi.idx1] \n")  
    #print(x.hist1.eqd$counts[P.lo.idx1:P.hi.idx1])

    # if according to limits, x.hist1.eqd1 replaces x.hist1.eqd0 
    go.on        <- ((bins.n.min < inner.bins.n1) & (!n.per.bin.ok1)) | 
                     (inner.bins.n1 > bins.n.max)
    if (go.on)
    { 
      x.hist.eqd <- x.hist1.eqd
      inner.bins.n <- inner.bins.n1
      P.lo.idx <- P.lo.idx1
      P.hi.idx <- P.hi.idx1
    }
  
  }   #  while ( ...

  #  x.hist1.eqd is complete histogram with bin size factor h

  #  Original meaning of bins.n
  bins.n <- length(x.hist.eqd$counts)

  #  Derived and proposed bins.n
#  cat("\n [FindHisto] Proposals for the number of bins",
#      "\n             Sturges              ", nclass.Sturges(x),
#      "\n             Scott                ", nclass.scott(x),
#      "\n             Freedman-Diaconis    ", nclass.FD(x),
#      "\n             Target in FindHisto  ", bins.n.max,
#      "\n             Derived in FindHisto ", bins.n,
#      "\n             inner.bins.n         ", inner.bins.n,
#      "\n             Bin size factor h    ", h,
#      "\n             n.per.bin.min        ", n.per.bin.min,
#      "\n")   

  # -------------------------------------------------------------------------  
  #  Deal with detection limits

  if (!is.na(detect.limits.max))
  {
    #  Keep breaks[1], next break is the first of actual breaks which
    #  is >= detect.limits.max
    breaks.n <- length(x.hist.eqd$breaks)
    x.breaksDL <- c(x.hist.eqd$breaks[1], 
                    x.hist.eqd$breaks[x.hist.eqd$breaks >= detect.limits.max])

    # New histogram needed
    x.hist.eqd     <- hist(x, breaks=x.breaksDL, right=FALSE, plot=FALSE)  
  }

  # -------------------------------------------------------------------------  
  #  Smooth the histograms, if requested

  if ((!smooth.hist1) & (!smooth.hist2))
  { # No smoothing
    x.hist.eqd.smo <- x.hist.eqd
  }

  if (( smooth.hist1) & (!smooth.hist2))
  { # Smoothing by moving average
    x.hist.eqd.smo <- SmoothHist1(x.hist.eqd, figA=NA )
  }

  if ((!smooth.hist1) & ( smooth.hist2))
  { # Smoothing by kde fit
    x.hist.eqd.smo <- SmoothHist2(x.kde, x.hist.eqd, NA, 
                                  histcol2, bordercol2, 
                                  histcol3, bordercol3, kdecol, subtitle)
  } 

  # -------------------------------------------------------------------------  
  # x.hist.eqd.smo can contain a large number of outer bins with very few 
  # cases - unuseful for estimation and time consuming). Collapse these bins
  # if they exist
  #  To shorten the name: x.hist.sc = x.hist.smo.coll

  #  will be done by CollapseBins()
  x.hist.sc <- CollapseBins(x.hist.eqd.smo$breaks, x.hist.eqd.smo$counts, 
                             P.lo.idx, P.hi.idx, n.per.bin.min/2)

  # --------------------------------------------------------------------------
  #  Plot the initial histogram (h is already adjusted) and the collapsed
  #  histogram

  if (!is.na(figB))
  { dev.set(figB)
    bringToTop(figB)

    par(mfcol=c(3,1))

    x.clip <- c(0, 60)

    plot(x.hist.sc,
         border=bordercol1, col=histcol1, freq=FALSE, lty=1,
         xlim=x.clip, 
         main="Equidistant histogram, smoothed, collapsed",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=paste(subtitle, "  n =", x.n, 
                   "  bins.n =", bins.n), cex.sub=0.7)
    axis(1, at=seq(0, 50, by=2.5), labels=FALSE)

    plot(x.hist.eqd.smo,
         border=bordercol2, col=histcol2, freq=FALSE, lty=1,
         xlim=x.clip, 
         main="Equidistant histogram, smoothed",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=paste(subtitle, "  n =", x.n, 
                   "  bins.n =", bins.n), cex.sub=0.7)
    axis(1, at=seq(0, 50, by=2.5), labels=FALSE)

    plot(x.hist.eqd,
         border=bordercol3, col=histcol3, freq=FALSE, lty=1,
         xlim=x.clip, 
         main="Equidistant histogram, not smoothed",
         xlab=xlabel,
         ylab="pdf(y)",
         sub=paste(subtitle, "  n =", x.n, 
                   "  bins.n =", length(x.hist.eqd$counts)), cex.sub=0.7)
    axis(1, at=seq(0, 50, by=2.5), labels=FALSE)

  #  plot(x.hist0.eqd,
  #       border=bordercol3, col=histcol3, freq=FALSE, lty=1,
  #       xlim=x.clip, 
  #       main="Initial equidistant histogram, not smoothed",
  #       xlab=xlabel,
  #       ylab="pdf(y)",
  #       sub=paste(subtitle, "  n =", x.n, 
  #                 "  bins.n =", length(x.hist0.eqd$counts)), cex.sub=0.7)
  #  axis(1, at=seq(0, 50, by=2.5), labels=FALSE)

    legend("topright", c("x.hist.esc", "x.hist.eqd.smo", "x.hist.eqd"),
           col=c(bordercol1, bordercol2, bordercol3), lwd=3, cex=0.8)
           #  fill=c(histcol1, histcol2, histcol3),
  }

  # --------------------------------------------------------------------------

  # x.hist.eqd        equidistant, no smooth
  # x.hist.eqd.smo    equidistant, smoothed
  # x.hist.sc         (usually not equidistant), smoothed, collapsed
  # h                 bin size factor (1: no increase)
 
  return(list(x.hist.eqd=x.hist.eqd,
              x.hist.eqd.smo=x.hist.eqd.smo,
              x.hist.sc=x.hist.sc,
              inner.bins.n=inner.bins.n,
              h=h))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindKDE.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Find a reasonable kernel density estimate for rounded data. 
#  Only the inner interval (Q.lo, Q.hi) is approximated by kde, to avoid the 
#  need of extremely many support points due to extreme observations
#  Problematic: Q.lo if many values < DL. Choice of surrogate value is 
#  important here.
 
#  #  CHANGE HISTORY
#  08.10.2020 Start
#
# ============================================================================

FindKDE <- function(x, round.unit, kernel.n, 
                    figA, subtitle, kdecol,
                    Q.lo=0.025, Q.hi=0.975, b.fact=1.05)
{
  #  INPUT 
  #  x           data vector
  #  round.unit  unit used when rounding the data
  #  kernel.n    number of support points for kde
  #  figA        device for figure illustrating the steps of finding a 
  #              bandwidth. NA: no plot.
  #  subtitle    subtitle in FigA 
  #  kdecol      color for kde line in FigA 

  #  OUTPUT 
  #  
  # --------------------------------------------------------------------------

  KDE <- FindBw(x, kernel.n, round.unit, figA, subtitle, kdecol,
                    b.fact=1.05, inner.lo=0.025, inner.hi=0.975)

  x.kde      <- KDE$x.kde
  x.kde.mode <- KDE$x.kde.mode

  return(list(x.kde=x.kde, x.kde.mode=x.kde.mode)) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_FindModeIndex.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  In a histogram, find the index of the bin containing the mode 
 
#  #  CHANGE HISTORY
#  31.08.2020 Start
#
# ==============================================================================

FindModeIndex <- function(breaks, counts, mode)
{
  #  INPUT 
  #  breaks   breaks defining the bins of a histogram
  #           bins are right-open: [ , ) except the last bin, which is[ , ]
  #           hist() in R needs right=FALSE for this convention
  #  counts   number of values per bin
  #  mode     usually from a density estimation 

  #  OUTPUT 
  #  idx      the index of the bin (in count) containing the mode
  # ---------------------------------------------------------------------------

  counts.n <- length(counts)
  idx      <- NA

  #  Which bin contains the mode?
  if (mode < breaks[1])
  { # below the first break
    idx <- 1
  }
  if (mode > tail(breaks, 1))
  { # above the last break
    idx <- counts.n
  }
  if (is.na(idx))
  { # somewhere within breaks
    idx <- which( (breaks[1:counts.n] <= mode) & 
                  (mode < breaks[2:(counts.n+1)] ))
  }
  if (is.na(idx))
  {
    idx <- Round(counts.n/2, 1)
    cat("\ ++++++ [FindModeIndex]  Mode index not found ++++++",
        "\n       Set to ", idx,
        "\n")
  }
  return(idx) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_Floor.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Function for rounding down to the next smaller unit 

#  CHANGE HISTORY
#  08.01.2020  Last change
# ============================================================================

Floor <- function(x,unit)
{ 
  # INPUT
  # x    number to round down
  # unit rounding unit.
  x.floor <- trunc(x/unit) * unit

  return(x.floor)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_FprContour.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Create a contour plot background for the distribution of estimated RLs.
#  RL contours must be added externally.
#  A legend for contour colours is produced by FprLegend()
 
#  CHANGE HISTORY
#  01.08.2021 Device and file names changed
#  28.07.2021 Start
# ==============================================================================

FprContour <- function(theta, figA.dev, figB.dev, figA.file, figB.file, 
                       RL1.min, RL1.max, RL2.min, RL2.max,
                       show.extra.axes,
                       RL1.n=11, RL2.n=11, levels=5,
                       xlabel="xlabel", ylabel="ylabel", 
                       maintitle="Main title",
                       subtext="subtext",
                       legendtitle="Legend title",  
                       RL1.p=0.025, RL2.p=0.975,
                       xl=1, yl=10,
                       outfile=NA)
{
  #  INPUT 
  #  theta    the true parameter vector
  #  figA     plot device for contour plot, must exist
  #  figB     plot device for legend, must be NA (no legend) or exist,
  #           then defined preferrably by win.graph(3.5, 7)  
  #  RL1.max, RL1.min   RL1 interval to display
  #  RL2.max, RL1.min   RL2 interval to display
  #  RL1.n    number of points on RL1 axis
  #  RL2.n    number of points on RL2 axis
  #  levels   if scalar: number of levels, will be selected as equidistant
  #           if vector: countour levels (border of a colour)
  #  xlabel
  #  ylabel
  #  maintext
  #  legendtitle
  #  outfile  basic plotfile name, fully qualified, including file type as 
  #           suffix (last 3 characters).
  #           The file for the contour plot is produced if outfile != NA  
  #           The file for the legend is produced if  outfile != NA  
  #           The contour plot file name has  "cont" added at the end,
  #           the legend plot file has "legd". 

  #  OUTPUT 
  #  list with elements 
  #       $usr      user coordinates
  #       $borders  level of contour borders
  #       $col.list contour colours
  #  outfileCont as output file
  #  outfileLegd as output file
  # ---------------------------------------------------------------------------

  if (!is.na(figA.dev))  { dev.set(figA.dev) }

  #  Calculate the true RLs
  RL1.true <- q.PN(RL1.p, theta[1], theta[2], theta[3])
  RL2.true <- q.PN(RL2.p, theta[1], theta[2], theta[3])

  C <- matrix(NA, nrow=RL2.n, ncol=RL1.n)
  RL1.seq <- seq(RL1.min, RL1.max, length.out=RL1.n)
  RL2.seq <- seq(RL2.min, RL2.max, length.out=RL2.n)

  #  Difference between intended and actual false positive rate below
  delta.FPR1 <- cdf.PN(RL1.seq, theta[1], theta[2], theta[3]) - RL1.p
  
  #  Difference between intended and actual false positive rate below
  delta.FPR2 <- cdf.PN(RL2.seq, theta[1], theta[2], theta[3]) - RL2.p

  #  Sum of absolute differences, as %
  delta.FPR  <- matrix(NA, nrow=RL1.n, ncol=RL2.n)
  for (i in 1:RL1.n)
  {
    for (j in 1:RL2.n)
    {
      delta.FPR[i, j] <- 100*(abs(delta.FPR1[i]) + abs(delta.FPR2[j]))
    }
  }

  dfp <- sort(as.vector(delta.FPR))
  delta.FPR.qua <- quantile(dfp, probs=seq(0, 1, by=0.1)) 
 
  #  Next command produces plot with non-reported scales for the main plot.
  #  Not useful if additional points/lines are to be placed in the plot 
  #filled.contour(x=RL1.seq, y=RL2.seq, z=delta.FPR,
  #               xlim = range(RL1.seq, finite = TRUE),
  #               ylim = range(RL2.seq, finite = TRUE),
  #               zlim = range(delta.FPR, finite = TRUE),
  #               nlevels = 9)

  #  Do everything by hand

  #  Span the plot without real information

  # par(mgp=c(3, 1, 0))      # Default-Zeilen für axis title / label / line
  # par(mgp=c(2.5, 1, 0))    # Zeilen für axis title / label / line
  par(mgp=c(2.5, 1.25, 0))   # Zeilen für axis title / label / line
  # par(mgp=c(3, 1, -1))      # bringt die tickmarks direkt an die Konturfläche
  par(bty="n")             # kein Kasten um die Plot-Fläche

  if (show.extra.axes)
  {
    par(xaxt="n")
    par(yaxt="n")
  }

  par(xaxs="i")
  par(yaxs="i")

  plot(c(RL1.seq[1], RL1.seq[RL1.n]), 
       c(RL2.seq[1], RL2.seq[RL2.n]), type="p", col="transparent",
       xlim=c(RL1.seq[1], RL1.seq[RL1.n]), 
       ylim=c(RL2.seq[1], RL2.seq[RL2.n]),
       xlab=xlabel, ylab=ylabel, main=maintitle,
       sub=subtext, cex.sub=0.8)

  #  Plot contour levels

  if (length(levels) == 1)
  { #  Only number of levels is given. Mo precisely, this is the number of 
    #  borders between areas of the same colour. There is one area less
    #  than 'levels" .
    #  Partition the range of delta.FPR into 'levels' equidistant intervalls

    borders <- seq(min(delta.FPR), max(delta.FPR), length.out=levels+1) 
  } else
  { #  levels are explicitly given
    borders <- levels
  }      
 
  borders.n <- length(borders)
  color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

  col.list <- color.palette(borders.n-1)

  #  Do contour plot
  .filled.contour(RL1.seq, y = RL2.seq, z = delta.FPR, 
                  levels=borders, col=col.list)


  FprC <- list(usr=par("usr"), borders=borders, col.list=col.list) 

  if (!is.na(figA.file))
  {
    filetype <- substr(figA.file, nchar(figA.file)-2, nchar(figA.file))
    savePlot(file=figA.file, type=filetype)
  } 

  #  Produce legend

  if (!is.na(figB.dev)) 
  {
    FprLegend(FprC, figB.dev, legendtitle, xl=xl, yl=yl) 

    if (!is.na(figB.file))
    {
      filetype <- substr(figB.file, nchar(figB.file)-2, nchar(figB.file))
      savePlot(file=figB.file, type=filetype)
    } 
  }
  dev.set(figA.dev)

  return(FprC)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#  
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_FprLegend.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Produce a legend for contour colours. Function parameters are provided
#  as output from FprContour()
 
#  CHANGE HISTORY
#  29.07.2020 Start
# ==============================================================================

FprLegend <- function(FprC, figB, legendtitle, xl=1, yl=10) 
{
  #  INPUT 
  #  FprC     list from FprContour with elements 
  #           $usr     user coordinates
  #           borders  level of contour borders
  #           col.list contour colours
  #  figB     plot device, must exist
  #  xl       horizontal size of legend
  #  yl       vertical size of legend

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------

  dev.set(figB)

  #  No annotation (axes remain)
  par.old <- par()
 
  par(ann=FALSE)
  par(xaxt="n")
  par(yaxt="n")
  par(bty="n")

  hsize <- 2.5 * xl
  vsize <- 1.1 * yl

  plot(c(0, xl, xl, 0, 0), c(0, 0, yl, yl, 0), type="l", col="black",
       xlim=c(0, hsize), ylim=c(0, vsize))

  #  borders separate contour areas
  borders.n <- length(FprC$borders)
  cell.height <- yl/(borders.n-1)

  for (i in 1:(borders.n-1))
  {
    polygon(c(0, xl, xl, 0), 
            c((i-1)*cell.height, (i-1)*cell.height,
              i*cell.height, i*cell.height),
            col=FprC$col.list[i], border="black")

    text(1.8*xl, (i-0.5)*cell.height, 
         paste(formatC(FprC$borders[i], width=5, digits=2, format="f"), 
               "-", 
               formatC(FprC$borders[i+1], width=5, digits=2, format="f"))) 
  }
  text(hsize/2, 1.05*yl, legendtitle)
  
  par(par.old)

  # return(list(usr=par("usr"), borders=borders, col.list=col.list))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# *****************************************************************************
#  F_Hist.table.R
#  
#  Print histogram data in a readable form
#
#  22.11.2021 Name changed: hist.table -> Hist.table,
#             cum columns added 
#  09.09.2020 print() replaced by return()
#  04.03.2020
# =============================================================================

Hist.table <- function(breaks, counts)
{
  #  Print histogram data
  x.n      <- sum(counts)
  counts.n <- length(counts)
  tab <- data.frame(x.lo=breaks[1:counts.n],
                    x.hi=breaks[2:(counts.n+1)],
                    count=counts,
                    prop=counts/x.n,
                    cumcount=cumsum(counts),
                    cumprop=cumsum(counts/x.n) )
  return(tab) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
# F_Infoblock.R

#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Writes essential control parameters in results plot
#  Version for TMC with dynamic prop.tmc
#

#  TO DO       
#  -  

#  #  CHANGE HISTORY
#  09.04.2021 oh.sym, dev.sym added
#  23.01.2021 x.bins.min replaced by bins.n.min, bins.n.max added
#  15.01.2021 r.fact added to call
#  10.01.2021 pathol.position removed
#  11.12.2020 Doan.adjust removed from call
#  02.12.2020 second switch for smooth.hist added
#  11.05.2020 More variables added
#  27.02.2020 Changed to function
#  24.01.2017 No details into plot if no plot was generated
# =================================================================================
#
Infoblock <- function(RunId, yplotmin, yplotmax, 
                     datafile, x.n, subsample.n, oh.sym, dev.sym, 
                     round.unit, smooth.hist1, 
                     smooth.hist2, 
                     n.per.bin.min, bins.n.min, bins.n.max,
                     lambda.min, lambda.max,   
                     l.fact, p.fact, r.fact, s.fact, w.fact, 
                     x.tr.prop.min, x.tr.prop.max, p.fit.min)
{ 
  #  Find user coordinates of actual window
  usr.coord <- par("usr")
  xplotmin <- usr.coord[1]
  xplotmax <- usr.coord[2]

  #  Use them to position the text
  text(xplotmin+0.010*(xplotmax-xplotmin),
       yplotmin+0.30*(yplotmax-yplotmin),
       labels=paste("\nRun Id            ", RunId,
                    "\ndata              ", datafile,
                    "\ntotal size        ", x.n,"   subsample.n   ",subsample.n, 
                    "\nout-/inpatient    ", oh.sym,
                    "\ndevice            ", dev.sym,
                    "\nround.unit        ", round.unit,
                    "\nsmooth.hist1      ", smooth.hist1,
                    "\nsmooth.hist2      ", smooth.hist2,
                    "\nn.per.bin.min     ", n.per.bin.min,
                    "\nbins.n.min        ", bins.n.min,
                    "\nbins.n.max        ", bins.n.max,
                    "\nl.fact            ", l.fact,
                    "\np.fact            ", p.fact,
                    "\nr.fact            ", r.fact,
                    "\ns.fact            ", s.fact,
                    "\nw.fact            ", w.fact,
                    "\nx.tr.prop.min/max ", x.tr.prop.min, "/", x.tr.prop.max,
                    "\nlambda.min/max    ", lambda.min, "/", lambda.max,
                    "\np.fit.min         ", p.fit.min
                    ),
     cex=0.8,pos=4,family="mono",col="navy") 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#  F_IniMatrix.R
#
#  26.01.2021
# ==============================================================================

IniMatrix <- function(rnames, cnames, IniVal=NA)
{
  M <- matrix(IniVal, nrow=length(rnames), ncol=length(cnames))
  dimnames(M) <- list(rnames, cnames)
  return(M)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#rn <- c("r1", "r2", "r3")
#cn <- c("c1", "c2")

#A <- IniMatrix(rn, cn)
#A
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_IniTab.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Initialize output tables for TMC_seg100

#  CHANGE HISTORY

#  13.10.2022 Generation of2 vectors, separate for num and char values
#  07.10.2021 Start 
# ===========================================================================
IniTab <- function()
{
  tabc.names  <- c( "method", "Sex", "Age", "errcode")
  tabc        <- rep(NA, times=length(tabc.names))
  names(tabc) <- tabc.names

  tabn.names  <- c(     "irep"        , 
        "Age.mea"     , "subset.type"  , "n"         , "pct.lt.DL",
        "ilo"         , "ihi"          , "x.tr.lo"   , "x.tr.hi" ,
        "x.tr.n"      , "x.lt.tr.n"    , "x.ge.tr.n" ,   
        "x.tr.prop"   , "prop.lo"      , "prop.hi"   , "x.tr.bins",
        "TI.qual"     ,
        "lambda"      , "mue"          , "sigma"     ,
        "x.RL1"       , "x.RL2"        , 
        "x.RL1.cilo"  , "x.RL1.cihi"   , "x.RL2.cilo", "x.RL2.cihi", 
        "x.RL1.tilo"  , "x.RL1.tihi"   , "x.RL2.tilo", "x.RL2.tihi",
        "delta.RL1"   , "delta.RL2"    , 
        "prev.l"      , "prev.c"       , "prev.r"    , 
        "xl.n"        , "xc.n"         , "xr.n"      ,
        "r2"          , "rank.r2"      ,
        "chi2.total"  , "chi2.total.df", "chi2.trun" ,
        "chi2.trun.df", "chi2.path"    , "iter"      , "rc" ,  
        "crit1"       , "crit2"        , "crit3"     , "crit4" ,
        "crit5"       , "sol.score"    , "opt.crit"   , "p.fit" , 
        "p.rt"        , "rel.dist"  )

  tabn <- rep(NA, times=length(tabn.names))
  names(tabn) <- tabn.names

  return(list(tabc=tabc, tabn=tabn))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#tab <- IniTab()
#tab

#cc1 <- list(meth="??", n=11)
#cc1["meth"]
#cc1["n"]

#cc2 <- list(meth="!!", n=22)
#cc2["meth"]
#cc2["n"]

#cc3 <- list(NA, NA)
#names(cc3) <- c("meth", "n")
#cc3
#cc3["meth"] <- "xxx"
#cc3["n"]    <- 33

#cc4 <- list(rep(NA, times=2))  # rep erzeugt ein Object (den Vektor)
#names(cc4) <- c("meth", "n")
#cc4
#cc4["meth"] <- "yyy"
#cc4["n"]    <- 44

#cc5 <- list(NA, NA)
#names(cc5) <- c("meth", "n")
#cc5
#cc5["meth"] <- "zzz"
#cc5["n"]    <- 55
#cc5["n"]
#cc5["n"]    <- 66
#cc5["n"]
#unlist(cc5)  #  erzeugt char

#rbind(cc1, cc2, cc3)

#ttab           <- IniTab()
#ttab["n"]      <- 10000
#ttab["method"] <- "test"
#
#as.data.frame(ttab)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Interpolx.R

# Linear interpolation: x, y, y target are given, xtarget will be determined

# INPUT
# x        x values MUST BE sorted (ascending). Sorting is not checked!. 
# y        y values corresponding to x
# ytarget  y target value, the corresponding x will be calculated

# NAs are automatically removed from the input.
# y must be a strictly monotone sequence, otherwise there is no unique 
# solution. Monotonicity is not checked!

# OUTPUT 
# xtarget  required x, y(xtarget) = ytarget

# CHANGE HISTORY
# 19.12.2020 Superfluous checks removed
# 03.07.2015 Start
# =============================================================================

Interpolx <- function(x,y,ytarget)
{  
  # Remove (x,y) pairs containing NA
  ok <- (!is.na(x)) & (!is.na(y))
  xx <- x[ok]
  yy <- y[ok]

  xtarget <- NA
  arg.ok <- TRUE
 
  if (is.na(ytarget)) 
  { cat("\n++++++ [Interpolx]: ytarget is NA\n")
    arg.ok <- FALSE
  }
  n <- length(xx)

  if (n <= 1) 
  { cat("\n++++++ [Interpolx]: x too short\n")
    arg.ok <- FALSE
  }
  if (length(yy) != n) 
  { cat("\n++++++ [Interpolx]: length(x) != length(y)\n")
    arg.ok <- FALSE
  }
  #  Reasonable value for ytarget?
  if (arg.ok)
  {
    if ( (ytarget < min(yy)) | (ytarget > max(yy)) )
    { cat("\n++++++ [Interpolx]: ytarget not in [min(y), max(y)]\n")
      arg.ok <- FALSE
    }
  }
  # No result if
  # ytarget is NA
  # length(x) <= 1
  # length(x) != length(y)
  # ytarget not in [min(y), max(y)]

  if (arg.ok)
  {
    for (i in 2:n)
    {
      if ( ((yy[i-1] <= ytarget) & (ytarget < yy[i  ])) |
           ((yy[i  ] <= ytarget) & (ytarget < yy[i-1])) )
      {
        xtarget <- xx[i-1] + (xx[i]-xx[i-1])*(ytarget-yy[i-1])/(yy[i]-yy[i-1])
        break
      }
    }
    if (ytarget == yy[n]) xtarget <- xx[n]
  }
  return(xtarget)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#D <- matrix(c(
#000, 135.0, 0.000,
#135, 135.5, 0.005,
#136, 136.5, 0.045,
#137, 137.5, 0.075,
#138, 138.5, 0.160,
#139, 139.5, NA,
#140, 140.5, 0.470,
#141, 141.5, 0.645,
#142, 142.5, 0.785,
#143, 143.5, 0.890,
#144, 144.5, 0.940,
#145, 145.5, 0.970,
#146, 146.0, 1.000), byrow=TRUE, ncol=3, nrow=13)
#colnames(D) <- c("u", "x", "y")
#Interpolx(D[ ,"x"], D[ ,"y"], 0.00)
#Interpolx(D[ ,"x"], D[ ,"y"], 0.30)
#Interpolx(D[ ,"x"], D[ ,"y"], 1.00)
#Interpolx(D[ ,"x"], D[ ,"y"], 2.00)

# 0.295


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Interpoly.R

# Linear interpolation: x, y, x target are given, ytarget will be determined

# INPUT
# x        x values MUST BE sorted (ascending). Sorting is not checked!. 
# y        y values corresponding to x
# xtarget  x target value, the corresponding y will be calculated
# LOCF     "Last observation carried forward": use largest x as result if
#          xtarget > max(x) ?

# NAs are automatically removed from the input.
# (x, y) must be a strictly monotone sequence, otherwise there is no unique 
# solution. Monotonicity is not checked!

# OUTPUT 
# ytarget  required y, x(ytarget) = xtarget

# CHANGE HISTORY
# 18.01.2021 Superfluous checks removed
# 03.07.2015 Start
# =============================================================================

Interpoly <- function(x, y, xtarget, LOCF=FALSE)
{  
  # Remove (x,y) pairs containing NA
  ok <- (!is.na(x)) & (!is.na(y))
  xx <- x[ok]
  yy <- y[ok]

  ytarget <- NA

  arg.ok <- TRUE
  if (is.na(xtarget)) 
  { cat("++++++ [Interpoly]: xtarget is NA\n")
    arg.ok <- FALSE
  }
  n <- length(yy)
  if (n <= 1) 
  { cat("\n++++++ [Interpoly]: y too short\n")
    arg.ok <- FALSE
  }
  if (length(xx) != n) 
  { cat("\n++++++ [Interpoly]: length(x) != length(y)\n")
    arg.ok <- FALSE
  }
  #  Reasonable value for ytarget? 
  if (arg.ok & !LOCF)
  {
    if ( !(xx[1] <= xtarget & xtarget <= xx[n]) &
         !(xx[n] <= xtarget & xtarget <= xx[1]) )
    { cat("\n++++++ [Interpoly]: xtarget not in [min(x), max(x)]\n")
      arg.ok <- FALSE
    }
  }
  if (arg.ok & LOCF)
  {
    if ( !(xx[1] <= xtarget))
    { cat("\n++++++ [Interpoly]: xtarget < min(x)\n")
      arg.ok <- FALSE
    }
  }

  if (arg.ok)
  {
    for (i in 2:n)
    {
      if ( (xx[i-1] <= xtarget & xtarget < xx[i  ]) |
           (xx[i  ] <= xtarget & xtarget < xx[i-1]) )
      {
        ytarget <- yy[i-1] + (yy[i]-yy[i-1])*(xtarget-xx[i-1])/(xx[i]-xx[i-1])
        break
      }
    }
    if (xtarget == xx[n]) ytarget <- yy[n]

    if ((xtarget > xx[n]) & LOCF) {ytarget <- yy[n]}
  }
  return(ytarget)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test
#x0 <- c(1,2,3)
#y0 <- c(10,20,30)
#print(Interpoly(x0,y0,0.5))
#print(Interpoly(x0,y0,1.0))
#print(Interpoly(x0,y0,1.5))
#print(Interpoly(x0,y0,2.0))
#print(Interpoly(x0,y0,2.5))
#print(Interpoly(x0,y0,3.0))
#print(Interpoly(x0,y0,3.5))
#print(Interpoly(x0,y0,3.5,LOCF=TRUE))

#x0 <- c(1,2)
#y0 <- c(10,20,30)
#print(Interpoly(x0,y0,1.5))
#x0 <- c(1,2,3)
#y0 <- c(30,20,10)
#print(Interpoly(x0,y0,0.5))
#print(Interpoly(x0,y0,1.0))
#print(Interpoly(x0,y0,1.5))
#print(Interpoly(x0,y0,2.0))
#print(Interpoly(x0,y0,2.5))
#print(Interpoly(x0,y0,3.0))
#print(Interpoly(x0,y0,3.5))
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# F_modified.tukey.R


#  QQ-Plot zur Berechung von RL in versuchter Annäherung an G. Hoffmann
#  QQ-Plot modifizierte Form ähnlich G. Hoffmann, mail vom 3.5.2017,
#  Trillium Normalizer

#  Weitere Kommentare siehe F_QQPlotHW.R


#  Change History
#  22.08.2020 In dieser Form eingerichtet  
# ============================================================================


modified.tukey <- function(vNumeric=NULL, perc=2.5, log.mode=TRUE, n.cycles=0,
                           print.cycles=TRUE, print.boxplots=TRUE, ytext="")
{
  #returns a numeric vector without outliers
  #perc indicates, which percentage is removed from each side of a normal or lognormal distribution (see log.mode)
  #repeats the algorithm n.cycles times
  #if n.cycles==0: repeats algorithm until no further outliers are detected
  #option: plots boxplots

  if (is.null(vNumeric)){return(NULL)}

  #initializes variables

  x <- vNumeric

  if(log.mode){x <- log(x)}
  n0 <- length(x)
  n1 <- 0
  i <- 0
  if(print.cycles)
  {
    m <- c(min(x), max(x))
    if(log.mode){m <- exp(m)}
    print(paste("cycle",i,"n =", length(x), "min =", m[1], "max =", m[2])) 
  }

  #calculates a modified tukey factor from perc

  quartile.factor <- qnorm(perc/100)/qnorm(0.25)

  #eliminates outliers from vector x
  if (n.cycles > 0)
  { #performs predefined number of cycles

    for (i in 1:n.cycles) 
    {
      i <- i+1
      m <- median(x)
      var1 <- m - quantile(x, 0.25)
      var2 <- quantile(x, 0.75) - m
      var <- min(c(var1, var2))
      lim <- c(m - quartile.factor * var, m + quartile.factor * var)
      x <- subset(x, x >= lim[1] & x <= lim[2])
      if(print.cycles)
       {
        m <- c(min(x), max(x))
        if(log.mode){m <- exp(m)}
        print(paste("cycle",i,"n =", length(x), "min =", m[1], 
                    "max =", m[2]))
      }
    }
  } else 
  { #cycles until no further outliers are detected
    while (n0 > n1)
    {
      i <- i + 1
      n0 <- length(x)
      m <- median(x)
      var1 <- m - quantile(x, 0.25)
      var2 <- quantile(x, 0.75) - m
      var <- min(c(var1, var2))
      lim <- c(m - quartile.factor * var, m + quartile.factor * var)
      x <- subset(x, x >= lim[1] & x <= lim[2])
      n1 <- length(x)
      if(print.cycles)
      {
        m <- c(min(x), max(x))
        if(log.mode){m <- exp(m)}
        print(paste("cycle",i,"n =", length(x), "min =", m[1], "max =", m[2]))
      }
    }
  }
 
  if (log.mode){x <- exp(x)}

  #plots boxplots before and after outlier elimination

  if(print.boxplots)
  {
    par(mfrow=c(1,2))
    boxplot(vNumeric, xlab="before outlier elimination", ylab=ytext)
    boxplot(x, xlab="after outlier elimination")
    par(mfrow=c(1,1))
  }
  return(x)
}

#  F_npa.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Do nonparametric RL estimation. Only possible of group membership is known.

#  CHANGE HISTORY

#  18.10.2022 Start from TMC_seg102_Analysis_npa.R
# ===========================================================================

npa <- function(x, grp, round.unit, x.kde.mode, tabn.npa, 
                RL1.p, RL2.p, print.log.message=FALSE)
{
  # 
  #
  #
  # ==========================================================================

  if (print.log.message) { cat("%%%                    npa Start\n") }

  tabn.npa["xc.n"] <- sum(grp=="C") # this number would be known under perfect
                                 # identification of non-diseased subject
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

  if (print.log.message) { cat("%%%                         npa End\n") }

  return(tabn.npa)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_o2R.R

o2R <- function(theta,RB,fastnull=1.e-10)
{ #  Parameter der PNV auf Originalskala --> unbeschränkte Skala (ganz R)
  #  Transformation durch inverse logit-Transformation mit Grenzen in RB
  #  (offene Intervalle)
  #  16.12.2017 

  #  theta   - (lambda,mue,sigma) einer Power-Normal-Verteilung
  #            Werte müssen innerhalb der durch RB3 gegebenen Grenzen liegen 
  #  RB[1, ] - c(lambda.min,lambda.max)
  #  RB[2, ] - c(mue.min,mue.max)
  #  RB[3, ] - c(sig.min,sig.max) 
  #   

  theta.n <- length(theta)
  ttheta  <- rep(NA,times=theta.n) 

  #  Logistische Transformation nur möglich, wenn Parameter innerhalb der
  #  Randbedingungen liegen, sonst Meldung

  for (k in 1:theta.n)
  {
    if ( (RB[k,1] < theta[k]) & (theta[k] < RB[k,2]) ) 
    { ttheta[k] <- log((theta[k]-RB[k,1])/(RB[k,2]-theta[k])) } else
    { cat("\n++++++++  Invalid input to o2R:",
          "\n          Parameter no             ",k,           
          "\n          Parameter value          ",theta[k],
          "\n          Permissible lower limit >",RB[k,1],
          "\n          Permissible upper limit <",RB[k,2],
          "\n") 
    }
  }

  #cat("o2R: ttheta,theta\n")
  #print(cbind(ttheta,theta))

  return(ttheta)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# F_opt.crit2.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Calculates basic part of optimality criterion: 1 - goodness of fit 
#  - if the chi2 value is small enough (smaller than the quantile 
#    of 1-node.p): 1 - goodness of fit
#  - if not: find the chi2 quantile corresponding to   1-node.p,
#    continue the (1-p(chi2) curve from that quantile in with a 
#    sraight line with slope = slope of  1-p(chi2) at the node
#    This is differentiable, different from the 1st version from 06.01.2022
#
#  Note that opt.crit will be minimized
# 
#  Rechnung für Skalar ok,
#  für Vektor unnötig umständlich, weil nr wenige verschiedene df vorkommen,
#  aber für jede Komponente von chi2 wird node.chi neu berechnet
  
#  31.07.2022 Criterion made differentiable
#  06.01.2022 chi2 == NA, NaN treated correctly
#  18.12.2021 Start
# ------------------------------------------------------------------------------

oc2 <- function(chi2, df, node.p=0.1, na.rm=TRUE)
{ 
  #  INPUT
  #  chi2   term measuring distance between observed and predicted in 
  #         a chi square manner, may contain penalty components, 
  #         smaller is better
  #  df     corresponding degrees of freedom
  #  chi2 and df must have same length!
  #
  #  OUTPUT 
  #  oc2    optimality criterion (smaller is better)
  # -------------------------------------------------------------------------

  chi2.n <- length(chi2)
  df.n   <- length(df)

  if (chi2.n != df.n)
  {
    cat("\n +++ Error on call oc2(): arguments have different lengths",
        "\n     length(chi2) = ", chi2.n,       
        "\n     length(df)   = ", df.n,
        "\n +++ Execution stops",
        "\n")
     stop("+++ Execution stopped in oc2 +++")
  }  else
              
  { 
    #  If chi2 or df missing: remove if allowed
    if (na.rm)
    {    
      ok <- !is.na(chi2) & !is.na(df)
      if (sum(ok) > 0) 
      { 
        chi2 <- chi2[ok]
        df   <- df[ok]
      }  else
      { 
        #  Nothing left 
        cat("\n +++ Error on call oc2(): no non-missing arguments",
            "\n     chi2 ", chi2, "   df ", df, 
            "\n +++ Execution stops",
            "\n")
        stop("+++ Execution stopped in oc2 +++")
      }
    }

    #  Data is ok, execution
    #  Determine location for transition 1-P(chi2) -> straight line
    node.chi2 <- qchisq(1 - node.p, df)   # horizontal coord
                                          # vertical coord is node.p
    #  Calculate value of oc
    opt.crit <- rep(NA, times=chi2.n)

    chi2.le.node.chi2 <- as.numeric(chi2 <= node.chi2)

    if (sum(chi2.le.node.chi2) > 0)
    {
      #  There are chi2 values below the node

      opt.crit[chi2.le.node.chi2] <- pchisq(chi2[chi2.le.node.chi2], 
                                              df[chi2.le.node.chi2])
    }  

    chi2.gt.node.chi2 <- as.numeric(chi2 > node.chi2)  #  earlier: >=

    if (sum(chi2.gt.node.chi2) > 0)
    {
      #  There are chi2 values beyond the node

      # Determine slope of optimality criterion at node.chi2
      node.slope     <- dchisq(node.chi2[chi2.gt.node.chi2], 
                                      df[chi2.gt.node.chi2])
      node.intercept <- 1 - node.p - node.slope * 
                        node.chi2[chi2.gt.node.chi2]
      opt.crit[chi2.gt.node.chi2] <- node.intercept + 
                                     node.slope * chi2[chi2.gt.node.chi2]
    }
  }

  return(opt.crit)
}

# ------------------------------------------------------------------------------
#optcrit <- oc2(78.79624,56)

#chi2    <- seq(0, 80, by=10)
#df      <- rep(56, times=length(chi2))
#pfit    <- 1-pchisq(chi2, df)

#optcrit <- oc2(chi2, df)
#cbind(chi2, pfit, 1-pfit, optcrit)

#win.graph(7, 0.7*7)
#plot(chi2, optcrit, type="l", col="blue", lwd=2)
#lines(chi2, pfit, lwd=2)
#lines(chi2, 1-pfit, col="green4", lty=2, lwd=3)
#abline(h=c(0, 1), col="black")
#abline(h=0.9, col="red")

#legend("right", c("p.fit", "1-p.fit", "opt.crit", "node.p"),
#       col=c("black", "green4", "blue", "red"),
#       lty=c(1,2,1,1), lwd=2 )
#savePlot("chi2_vs_optcrit.bmp", type="bmp")
#savePlot("chi2_vs_optcrit.pdf", type="pdf")

#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_pdf.PN.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  CHANGE HISTORY
#  16.01.2021 Comments changed to English
# ============================================================================ 

pdf.PN <- function(x,lambda,mue,sigma,fastnull=1.e-10)
{
  #  Probability density function of the power normal distribution (PND) 
  #  Background: see 
  #  Freeman J, Modarres R. Inverse Box-Cox: the power normal 
  #  distribution. Stat Probabil Letters 2006;76; 764-72, particularly p 766

  #  INPUT 
  #  x       value at which to evaluate the pdf
  #  lambda  shape parameter. Real number. Special cases:
  #          1: PND is normal distribution with mean  and sd sigma
  #          0: PND is a logarithmic normal distribution with parameters
  #             mue and sigma
  #          If abs(lambda) < fastnull, lambda== 0 is assumed and the 
  #          intrinsic function plnorm is used
  #  mue     mean of the PND (of BoxCox(x, lambda))
  #  sigma   sd of the PND   (of BoxCox(x, lambda))
  #  fastnull (="nearly zero") absolute values < fastnull are treated as zero
   
  #  OUTPUT
  #  pdf     Probability density function  

  # ==========================================================================

  pdf  <- rep(fastnull, times=length(x))
  xgt0 <- (x > fastnull)
  T    <- 1/(lambda*sigma) + mue/sigma
 
  if (lambda > fastnull)
  { #  lambda > 0  in principle
    K <- pnorm(T,mean=0,sd=1)
    pdf[xgt0] <- (1/K) * (1/(sigma*sqrt(2*pi))) * (x[xgt0])^(lambda-1) *
                exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )
    pdf[!is.finite(pdf)] <- fastnull   
  } 
  if (abs(lambda) <=  fastnull)  
  { # lambda == 0 (in principle)
    #  K <- 1 
    #pdf[xgt0] <-         (1/(x[xgt0] * sigma*sqrt(2*pi))) * 
    #              exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )

    #  intrinsic function is faster
    pdf[xgt0] <- dlnorm(x[xgt0],meanlog=mue,sdlog=sigma)
  }
  if (lambda < -fastnull)
  { #  lambda < 0   in principle
    #  Potential problem: if T -> Inf , K -> 0, which generates an 
    #  overflow below. wenn T groß, dann aber auch R -> 0
    K <- pnorm(-T,mean=0,sd=1)
    K <- max(K,fastnull)

    pdf[xgt0] <- (1/K) * (1/(sigma*sqrt(2*pi))) * (x[xgt0])^(lambda-1) *
                  exp(-0.5*((BoxCox(x[xgt0],lambda)-mue)/sigma)^2 )
  }

  #  This is a definition. @@@ Check for appropriateness eg if lambda=1 
  pdf[!xgt0] <- 0

  if (any(!is.finite(pdf)))
  { cat("[pdf.PN] Infinite values:\n") 
    cat("lambda,mue,sigma: ",lambda,mue,sigma,"\n")
    cat("x, pdf(x) with non-finite pdf\n")
    control <- data.frame(x,pdf)
    print(control[!is.finite(control[ ,2]), ])   
  }
  return(pdf)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #  lambda = -0.09989172   mue = 11.25949   sigma = 0.01440878 
#  x 14 - 19
#lambda <- 0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#x.P50  <- q.PN(0.50, lambda, mue, sigma)
#x      <-  seq(1,20)
#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(x, x.pdf, x.cdf))

#  Comparison with R provided function
# install.packages("powdist", dependencies=TRUE)

#library(powdist)
#dpnorm(x, lambda = lambda, mu = mue, sigma = sigma, log = FALSE)

# Kundu, D. and Gupta, R. D. (2013) Power-normal distribution. 
#  Statistics 47, 110–125

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_PermissibleDiff.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate the permissible difference pD according to 
#  Rainer Haeckel*, Werner Wosniok and Eberhard Gurr
#  Diagnostic efficiency in models for permissible
#  measurement uncertainty
#  J Lab Med 2017; 41(6): 309–315
#  and 
#  2021_01_06_Zula__ssige_Messunsicherheit_DGKL.3.10.20_-1.xlsx
#  (from DGKL web site)

 
#  CHANGE HISTORY
#  13.04.2022 equivalence factor introduced
#  26.03.2021 Start
# ==============================================================================

PermissibleDiff <- function(RL1, RL2, xi, ef)
{
  #  INPUT 
  #  RL1   Lower 95% reference limit (2.5%)
  #  RL2   Upper 95% reference limit (97.5%) 
  #  xi    Point at which to calculate the pD
  #  ef    equivalence factor, default from publication = 1.28

  #  OUTPUT 
  #  xi.lo    lower limit of the permissible interval
  #  xi.hi    upper limit of the permissible interval
  # ---------------------------------------------------------------------------

  y.sd      <- (log(RL2) - log(RL1)) / 3.92  # SD of a log normal distribution
  x.med     <- sqrt(RL1 * RL2)               # Median of a log normal distribution
  CV.E.star <- 100*sqrt(exp(y.sd^2) - 1)     # CV of a log normal distribution

  pCV.A     <-  sqrt(CV.E.star - 0.25)       # see figure in paper
  ps.A.med  <- pCV.A * x.med * 0.01

  ps.A.xi <- 0.8 * ps.A.med * xi / x.med + 0.2 * ps.A.med
  pD.xi   <- ef * ps.A.xi

  xi.lo   <- xi - pD.xi
  xi.hi   <- xi + pD.xi

  #  Not used
  # pCV.A.RL2 <- 100 * ps.A.RL2 / RL2
  # pU.pct    <- 2.39 * pCV.A.RL2

  return(data.frame(xi.lo=xi.lo, xi.hi=xi.hi))
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test  (compare with  
#         2021_01_06_Zula__ssige_Messunsicherheit_DGKL.3.10.20_-1.xlsx
#         sodium)

#pI <- PermissibleDiff(135, 145, 135)
#pI

#pI <- PermissibleDiff(135, 145, 145)
#pI

#pI <- PermissibleDiff(135.9, 144.9, 144.9)
#pI

#pI <- PermissibleDiff(135, 145, c(135, 145))
#pI


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# F_plot.res.R

plot.res <- function(xsupp, comp.n, lambda.BC, lambda, mue, sig, col, lty)
{
  for (ic in 1:comp.n)
  { x.pdf <- lambda[ic] * pdf.PN(xsupp, lambda.BC, mue[ic], sig[ic])
    lines(xsupp, x.pdf, col=col, lty=lty)
  }
}

  # *****************************************************************************
#  F_PlotConfElli.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Plot a confidence ellipse
#
#  16.12.2018 Einrichtung
# =============================================================================

PlotConfElli <- function(mean,cov,alpha,figA)
{ #
  # mean   mean of the ellipse
  # cov    covariance matrix 
  # alpha  area outside the ellipse
  # figA   graphics window, NA: no plot
  # ---------------------------------------------------------------------------
  
  #  Start with ellipse for standard normal distribution
  r <- qnorm(1-alpha/2)
 
  phi <- seq(0,360,by=5) * pi/180
  x0  <- r * cos(phi)
  y0  <- r * sin(phi)

  #  Jordan decomposition of cov
  cov.decomp <- eigen(cov,TRUE)
  #print(cov.decomp)  

  #  Root
  cov.root <- cov.decomp$vectors %*% diag(sqrt(cov.decomp$values)) %*% 
              cov.decomp$vectors 
  #  Check root
  cov.root %*% cov.root

  cov.elli0 <- cbind(x0,y0) %*% cov.root 
  cov.elli  <- cov.elli0 + matrix(rep(1,times=length(phi)),ncol=1) %*%
                           matrix(mean,nrow=1)  

  if (!is.na(figA))
  { dev.set(figA)
    xplotmin <- min(c(x0,cov.elli))
    xplotmax <- max(c(x0,cov.elli))
    yplotmin <- min(c(y0,cov.elli))
    yplotmax <- max(c(y0,cov.elli))

    plot(x0,y0,type="l",col="gray",
         xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax))
    abline(h=0,col="gray")
    abline(v=0,col="gray")
    lines(cov.elli,col="green3")
  }

  return(cov.elli)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#alpha <- 0.05
#mean <- c(1,-2)
#cov <- matrix(c(10, 2,
#                 2, 5),byrow=TRUE,ncol=2) 

#c.e <- PlotConfElli(mean,cov,alpha,NA)

#xplotmin <- min(c(x0,cov.elli))
#xplotmax <- max(c(x0,cov.elli))
#yplotmin <- min(c(y0,cov.elli))
#yplotmax <- max(c(y0,cov.elli))

#plot(cov.elli,type="l",col="green3",
#     xlim=c(xplotmin,xplotmax),ylim=c(yplotmin,yplotmax))
#     abline(h=0,col="gray")
#     abline(v=0,col="gray")
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  PlotDS.R
#
#  Plot descriptive statistics into plot
#
#  28.02.2020
# ==============================================================================

PlotDS <- function(xvar)
{
  #  Plot must exist. Get extreme coordinates.
  usr.coord <- par("usr")
  xplotmin <- usr.coord[1]
  xplotmax <- usr.coord[2]
  yplotmin <- usr.coord[3]
  yplotmax <- usr.coord[4]

  #  Use them to position the text
  xvar.n    <- sum(!is.na(xvar))
  xvar.min  <- min(xvar, na.rm=TRUE)
  xvar.P025 <- Quantile(xvar, probs=0.025)
  xvar.med  <- median(xvar, na.rm=TRUE)
  xvar.mea  <- mean(xvar, na.rm=TRUE)
  xvar.P975 <- Quantile(xvar, probs=0.975)
  xvar.max  <- max(xvar, na.rm=TRUE)


  text(xplotmin+0.010*(xplotmax-xplotmin),
       yplotmin+0.72*(yplotmax-yplotmin),
       labels=paste("\n# values ", format(xvar.n,   digits=6),
                    "\nMinimum  ", format(xvar.min, digits=4),
                    "\nP025     ", format(xvar.P025,digits=4),
                    "\nMedian   ", format(xvar.med, digits=4),
                    "\nMean     ", format(xvar.mea, digits=4),
                    "\nP975     ", format(xvar.P975,igits=4),
                    "\nMaximum  ", format(xvar.max, digits=4)
                    ),
     cex=0.7,pos=4,family="mono",col="firebrick") 

  # Show data points
  points(xvar, rep(yplotmin+0.10*(yplotmax-yplotmin), times=length(xvar)),
         type="h", col="gray")  

  # Show quantiles etc.
  abline(v=xvar.mea, col="black")
  abline(v=c(xvar.P025, xvar.med,xvar.P975), col="black", lty=2)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_PlotHistFit.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Plot a histogram of the data plus additional components, see INPUT
 
#  #  CHANGE HISTORY
#  01.05.2021 kde plot, difference kde-predicted removed, residual histogram
#             added 
#  12.01.2021 x.tr.prop.range added to call
#  01.12.2020 Call changed (n.per.bin.mea removed, denlty added)
#  19.11.2020 Start
#
# ==============================================================================

PlotHistFit <- function(figA, figA.file, figA.type, 
                        x.hist, xlabel, title, subtitle, n.per.bin.min,
                        x.kde, x.val, counts.range, x.tr.prop.range,
                        x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                        x.clip.type, yplotmaxA, par.las, par.tcl,
                        bordercol, histcol, difcol, kdecol, denlty)
{
  #  =====================================================================
  #  Plot histogram for the subset of the x range defined by x.clip
  #  x plot limits can be forced to coincide with x.breaks 
  #  (x.clip.type=="coincide") or can be used as specified
  #  (x.clip.type=="as.specified")
  #  Additonal components may be added here, see INPUT

  #  INPUT
  #  figA
  #  figA.file 
  #  figA.type
  #
  #  OUTPUT  

  # ==========================================================================

  dev.set(figA)
  bringToTop(figA)
  
  #  x limits depend on x.clip settings

  if (is.na(x.clip.min))
  { #  Nothing defined
    xplotminA <- x.hist$breaks[1]
  } else
  { if (x.clip.type=="as.specified")
    { 
      xplotminA <- x.clip.min
    } else
    {
      #  Limits defined by x.clip are shifted such that the limits
      #  coincide with break points
      #  Find closest break <= x.clip.min
      x.breaks.cand   <- x.hist$breaks[x.hist$breaks <= x.clip.min]
      x.breaks.cand.n <- length(x.breaks.cand)
      if (x.breaks.cand.n > 0)
      { 
        xplotminA <- x.breaks.cand[x.breaks.cand.n] 
      } else
      {
        xplotminA <- x.clip.min
      }
    
    }  # x.clip.type=="as.specified" ...   
  }    # is.na(x.clip.min) ...

  if (is.na(x.clip.max))
  { #  Nothing defined
    xplotmaxA <- tail(x.hist$breaks, 1)
  } else
  { if (x.clip.type=="as.specified")
    { 
      xplotmaxA <- x.clip.max
    } else
    { 
      #  find closest break <= x.clip.max
      x.breaks.cand <- x.hist$breaks[x.hist$breaks >= x.clip.max]
      x.breaks.cand.n <- length(x.breaks.cand)
      if (x.breaks.cand.n > 0)
      { 
        xplotmaxA <- x.breaks.cand[1] 
      } else
      {
        xplotmaxA <- x.clip.max
      }
    }    # x.clip.type=="as.specified" ...
  }      # is.na(x.clip.max)

  # Use y maximum if specified, otherwise use histogram density maximum 
  if (is.na(yplotmaxA))
  {
    yplotmaxA <- 1.1 * max(x.hist$density)
  }

  #  Set yplotminA such that there is space for indicating observed values
  yplotminA <- -0.05 * yplotmaxA

  #  Start plot
  par(las=par.las)   
  par(tcl=par.tcl)
   
  if (!is.na(x.clip.by1))
  { par(xaxt="n") } else
  { par(xaxt="s") } 

  #  Histogram
  plot(x.hist, freq=FALSE, 
       xlim=c(xplotminA,xplotmaxA),
       ylim=c(yplotminA,yplotmaxA),
       col=histcol, border=bordercol,
       xlab=xlabel, ylab="Density",
       main=paste(title, 
                  " range", counts.range, 
                  " prop", x.tr.prop.range),
       sub=subtitle, cex.sub=0.7)

  par(xaxt="s")
  if (!is.na(x.clip.by1) )
  { 
    axis(1,at=seq(xplotminA, xplotmaxA, by=x.clip.by1), labels=TRUE)  
  }

  if (!is.na(x.clip.by2) )
  { 
    axis(1, at=seq(xplotminA, xplotmaxA,by= x.clip.by2), labels=FALSE) 
  }

  #  Indicate observed values
  x.val.n <- length(x.val)
  if (x.val.n > 1)
  {
    points(x.val, rep(0.5*yplotminA, times=x.val.n), type="h", col="gray10")
  }

  #  Density estimation
  if (length(x.kde) > 1)
  {
    lines(x.kde,col=kdecol,lty=denlty)  
  }

  #  Save plot, if requested
  if (!is.na(figA.file))
  {
    savePlot(file=figA.file,type=figA.type)
  }
  return(c(yplotmin=yplotminA,yplotmax=yplotmaxA))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#  F_PlotHistResid.R
#  former name: F_HistoResid.R
#  (c) wwosniok@math.uni-bremen.de

#  TMC: Truncated minimum chi-square estimation

#  Function calculating and plotting the difference
#  between 2 histograms (usually observed and expected)
 
#  #  CHANGE HISTORY
#  26.04.2021 Start
# ==============================================================================

PlotHistResid <- function(temp.met, figA, yplotminA, bordercol, fillcol)
{
  # INPUT
  #  temp.met  object containing observed and expected values of a histogram,
  #             produced by chi2trunc, see there for details
  #  xlabel     histogram label
  #  figA	    window to put the plot in (as addition)

  # OUTPUT
  # delta.hist    the difference histogram

  # ============================================================================

  var.list <- c("UG", "OG", "nobs", "nerw", "diff")

  histo.resid <- temp.met$tab[ , var.list] 
  if (length(temp.met$tab.lo) > 1)
  {
    histo.resid <- rbind(temp.met$tab.lo[ ,var.list],
                         histo.resid)
  }

  if (length(temp.met$tab.hi) > 1)
  {
    histo.resid <- rbind(histo.resid, 
                     temp.met$tab.hi[ ,var.list]  ) 
  }

  #  Create a histogram object
  #  Observed density
  x.n       <- sum(histo.resid[ , "nobs"])
  bin.width <- histo.resid[ , "OG"] - histo.resid[ , "UG"]
  obs.den   <- histo.resid[ , "nobs"] / ( x.n * bin.width)
  exp.den   <- histo.resid[ , "nerw"] / ( x.n * bin.width)
  delta.den <- histo.resid[ , "diff"] / ( x.n * bin.width)
  delta.de2 <- obs.den - exp.den

  histo.resid <- data.frame(histo.resid,
                            obs.den=obs.den,
                            exp.den=exp.den,
                            delta.den=delta.den)

  #  Truncate negative densities for aesthetic reasons 
  histo.resid[histo.resid[ , "delta.den"]<yplotminA, "delta.den"] <- yplotminA

  exp.hist <- list(breaks=c(histo.resid[ , "UG"],
                             tail(histo.resid[ , "OG"], 1)),
                   counts=histo.resid[ , "nerw"],
                   density=histo.resid[ , "exp.den"],
                   mids=(histo.resid[ , "UG"]+histo.resid[ , "OG"])/2,
                   xname="delta",
                   equidist=FALSE)
  class(exp.hist) <- "histogram" 

  delta.hist <- list(breaks=c(histo.resid[ , "UG"],
                             tail(histo.resid[ , "OG"], 1)),
                   counts=histo.resid[ , "diff"],
                   density=histo.resid[ , "delta.den"],
                   mids=(histo.resid[ , "UG"]+histo.resid[ , "OG"])/2,
                   xname="delta",
                   equidist=FALSE)
  class(delta.hist) <- "histogram"  

  #  Add the histogram to an existing plot
  if (!is.na(figA))
  {
    #plot(delta.hist, border=bordercol, col=fillcol,  density=20, angle=45,
    #     lwd.border=5, add=TRUE)
    plot(exp.hist, border=blue.64, col=green.16, add=TRUE)

    plot(delta.hist, border=bordercol, col=fillcol, add=TRUE)
  }

  #  Show estimated denity

  return(delta.hist)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test: see Test_HistoResid.R

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



  
  
# ****************************************************************************
#  F_PlotMetRes.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Plot results of truncated estimation by method met in an existing histogram plot
 
#  #  CHANGE HISTORY
#  01.05.2021 kde residuals replaced by histogram residuals
#  08.01.2021 mode added to plot
#  07.12.2020 xsupp added to call
#  19.11.2020 Start
#
# ==============================================================================

PlotMetRes <- function(figA, figA.file, figA.type, yplotminA,  
                       xsupp, prev.c.met, xc.pdf.met, metcol,
                       x.kde.mode, kdecol, x.kde.met, difcol,
                       difbordercol, difhistcol, 
                       x.RL1.met, x.RL2.met, y.RL, 
                       x.tr.lo, x.tr.hi, y.TI, trcol,
                       temp.met) 
{
  #  =====================================================================
  #  Plot histogram for the subset of the x range defined by x.clip
  #  x plot limits can be forced to coincide with x.breaks 
  #  (x.clip.type=="coincide") or can be used as specified
  #  (x.clip.type=="as.specified")
  #  Additonal components may be added here, see INPUT

  #  INPUT
  #  figA
  #  figA.file 
  #  figA.type
  #
  #  OUTPUT  

  # ==========================================================================

  dev.set(figA)
  bringToTop(figA)

  #  Density estimated by method met  
  lines(xsupp, prev.c.met * xc.pdf.met, col=metcol, lwd=2)

  if (length(x.kde.met) > 1)
  {
    #  Difference between overall kde and estimated density
    lines(xsupp, x.kde.met, col=difcol, lwd=2) 
  }

  #  Histogram of difference between observed and expected
  delta.hist <- PlotHistResid(temp.met, figA, yplotminA, difbordercol, 
                              difhistcol)

  #  Mode from kde
  abline(v=x.kde.mode, col=kdecol, lty=2)

  #  Estimated RLs  
  plotRL(x.RL1.met, y.RL, metcol, pos=2)
  plotRL(x.RL2.met, y.RL, metcol, pos=4)

  #  Truncation interval
  lines(c(x.tr.lo, x.tr.hi), c(y.TI, y.TI), col=trcol, lwd=2)  

  #  Save plot, if requested
  if (!is.na(figA.file))
  {
    savePlot(file=figA.file,type=figA.type)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_PlotRL.R
#
#  CHANGE HISTORY
#
#  07.01.2020 Start

#  -----------------------------------------------------------------------------
plotRL <- function(RL,y,farbe,pos=4)
{ 
  #  RL    - x-Position der RL
  #  y     - y-Position der RL
  #  farbe - Farbe für Linie und Text
  #  pos   - Position des Textes neben Linie 
  #          1,2,3,4 = unter, links, über, rechts
 
  lines(c(RL,RL),c(0,y),col=farbe)
  text(RL,y,pos=pos,labels=format(RL,digits=3),cex=0.8,col=farbe)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_ProcessDL.R
#
#  (c) wwosniok@math.uni-bremen.de

#  Process data entries of the form "< DL" (DL = detection/determination limit).
#  No action, if no such entries exist.
#  If >= 1 entries exist, several options available. The surrogate factor
#  0 <= s.fact < DL must be supplied. s.fact affects breaks, kde, 
#  QQ estimation.

#  If 1 unique DL exists:
#  option 1 (default):
#  - replace "< DL"                            by s.fact*DL
#    replace values < DL (should never exist)  by s.fact*DL
#    
#  option 2 (not recommended):
#  - replace "< DL"                            by s.fact*DL
#    leave values < DL (should never exist)    as they are

#  If > 1 unique DL exist:
#  option 1 (default):
#  - replace all "< DL"                        by s.fact*max(DL)
#    replace all values < max(DL)              by s.fact*max(DL)
#    
#  option 2 (not recommended):
#  - replace all "< DL"                        by s.fact*max(DL)
#    leave values < some DL (should never exist)    as they are

#  option 3 (not recommended):
#  - replace "< DL"                               by s.fact*DL
#    replace values < some DL (should never exist)  by s.fact*min_DL(DL>value)

#  option 4 (not recommended):
#  - replace "< DL"                               by s.fact*DL
#    leave values < some DL (should never exist)  as they are

#  Non-numeric entries other than "<DL" are left in the data unchanged.

#  Change History
#  27.04.2020 detect.limits.max added to return argument, is needed by
#             other functions
#             Recommendation: use option 1 generally 
#  14.04.2020 Start from 070_ReadData.R
# =============================================================================

ProcessDL <- function(v, s.fact, option=1)
{
  #
  # INPUT
  #
  # v - data vector
  # s.fact - surrogate factor, \in [0, DL]
  # option - controls how values < DL are processed, see above

  # OUTPUT
  # w <- vector with DLs processed.
  # ===========================================================================

  #  Check for permissible option
  if (option %in% c(1,2,3,4))
  {
    detect.limits.max <- NA
    w <- v

    # Are there non-numerical entries?
    options(warn=-1)           #  suppress warning regarding NAs resulting from
                               #  transformation 
    char.val.table <- table(v[is.na(as.numeric(v))]) 
    options(warn=0)

    if (length(char.val.table) == 0)
    { 
      cat("\n No non-numerical values found\n")
    } else
    {
      cat("\n Table of non-numerical entries in the data\n")
      print(char.val.table)
    
      #  Is "<" in the string?
      char.val <- names(char.val.table)
      has.lt <- str_detect(char.val,"<")

      #  Answer may be NA, means here "no"
      has.lt[is.na(has.lt)] <- FALSE
      has.lt.n <- sum(has.lt)        # # of values in char.val with form "< DL"

      if (has.lt.n > 0)
      {
        #  Remove non-numerical values other than "< ..." from char.val.table
        char.val.table <- char.val.table[has.lt]
        char.val       <- char.val[has.lt]

        cat("\n Table of '< DL' entries in the data\n")
        print(char.val.table)

        #  Determine the value(s) of DL. May look different, but be the same 
        #  numerically ("<5" "<5.0", " 5", ...),
        #  so change to numerical before counting  
        DL <- rep(NA, times=has.lt.n)

        for (i.lt in 1:has.lt.n)
        {
          lt.idx <- str_locate(char.val[i.lt], "<")
          DL[i.lt] <- char.val[i.lt]
          str_sub(DL[i.lt], start = lt.idx[ ,1], end = lt.idx[ ,2]) <- " "
        }
        DL <- as.numeric(DL)
        detect.limits.max <- max(DL)

        char.val.table <- rbind(char.val.table, DL)
        char.val.table <- matrix(as.numeric(char.val.table), nrow=2)
        colnames(char.val.table) <- char.val
        rownames(char.val.table) <- c("count", "DL")

        cat("\n Table of '< DL' entries in the data and their numerical values\n")
        print(char.val.table)

        #  Make v numeric, where possible, otherwise (non-numeric gets NA)
        options(warn=-1)      #  suppress warning regarding NAs resulting from
                              #  transformation
        options(warn=-1)
        v.num <- as.numeric(v)
        options(warn=0)       #  return to usual warning style

        #  In the data, replace "< DL" by their surrogate values according to
        #  the chosen option

        if (option == 1)
        { #  Global surrogate for DL and values < max. DL
          for (i.lt in 1:has.lt.n)
          {
            w[v==char.val[i.lt]] <-  s.fact*detect.limits.max 
          }

          #  Real values < DL are replaced by s.fact*detect.limits.max
          subset <- (!is.na(v.num)) & (v.num < detect.limits.max)
          #print(data.frame(v, subset))

          subset[is.na(subset)] <- FALSE         
          w[subset]   <-  s.fact*detect.limits.max
        }

        if (option == 2)
        { #  Global surrogate only for DL
          for (i.lt in 1:has.lt.n)
          {
            w[v==char.val[i.lt]] <-  s.fact*detect.limits.max 
          }
          #  Real values < DL are unchanged
        }

        if (option == 3)
        { #  Individual surrogate for DL and values < indiv. DL
          for (i.lt in seq(has.lt.n, 1, by=-1))
          {
            w[v==char.val[i.lt]] <-  s.fact*char.val.table[2, i.lt]
            #cat("\n option=3 (1), i.lt=", i.lt) 
            #print(data.frame(v, v.num, w))   
 
            #  Real values < actual DL are replaced by s.fact*actual DL
            # look for values < the actual DL
            subset <-  (!is.na(v.num)) & (v.num < char.val.table[2, i.lt])
            subset[is.na(subset)] <- FALSE
            if (sum(subset) > 0)          
            {
              w[subset] <-  s.fact*char.val.table[2, i.lt]
              #cat("\n option=3 (2)") 
              #print(data.frame(v, v.num, w))   
            }
          }
        }

        if (option == 4)
          { #  Individual surrogate for DL only
          for (i.lt in seq(has.lt.n, 1, by=-1))
          {
            w[v==char.val[i.lt]] <-  s.fact*char.val.table[2, i.lt] 
          }
          #  Real values < DL are unchanged
        }

        options(warn=-1)      #  suppress warning regarding NAs resulting from
                              #  transformation
        w0 <- as.numeric(w)   #  as.numeric generates NA from non-numeric
                              #  data that may still be in w. Will be deleted
                              #  later, but shall stay here unchanged.
        is.na.w0 <- is.na(w0)
        w0[is.na.w0] <- w[is.na.w0]
        w <- w0 

        options(warn=0)       #  return to usual warning

      }   # has.lt.n > 0 ...

    }     # length(char.val.table) ...
  
    return(list(data=w, detect.limits.max=detect.limits.max))

  } else       # wrong option
  { cat("\n ++++++ [ProcessDL] Wrong option: ", option,
        "\n ++++++ Values below DL are not processed",
        "\n")
    return(list(data=v, detect.limits.max=NA))
  }         
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_ProcessPropLimits_V1.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Process the proportion contained in the truncation intervals.
#  Default depends on settings in 030, in turn dependent on whether a
#  detection limit exists. User input overrides.

#  CHANGE HISTORY

#  29.09.2022 Start
# ============================================================================

ProcessPropLimits <- function(x.tr.prop.min, x.tr.prop.max,
                              x.tr.prop.min.DL, x.tr.prop.max.DL,
                              x.tr.prop.min.noDL, x.tr.prop.max.noDL,
                              detect.limits.max, 
                              logfile=NA, print.log.message=FALSE)
{ 
  # ProcessPropLimits_V1
  # INPUT
  # x.tr.prop.min       lower limit for proportion, default or user input
  # x.tr.prop.max       upper limit for proportion, default or user input
  # x.tr.prop.min.DL    lower limit for proportion, DL present
  # x.tr.prop.max.DL    upper limit for proportion, DL present
  # x.tr.prop.min.noDL  lower limit for proportion, DL not present
  # x.tr.prop.max.noDL  upper limit for proportion, DL not present
  # logfile             logfile

  # OUTPUT
  # x.tr.prop.limits    effective values for x.tr.prop.min, x.tr.prop.max
  # ==========================================================================     

  # ==========================================================================
  if (print.log.message) { cat("%%%   ProcessPropLimits_V9   Start\n") }
  # ==========================================================================

  #  Process user input for x.tr.prop.min, max. If the parameter is not given
  #  by the user, its value is NA
  if ( is.na(x.tr.prop.min) )
  { # No user input
    if (is.na(detect.limits.max) ) 
    { x.tr.prop.min <- x.tr.prop.min.noDL } else
    { x.tr.prop.min <- x.tr.prop.min.DL }
  }

  if ( is.na(x.tr.prop.max) )
  { # No user input
    if (is.na(detect.limits.max) )
    { x.tr.prop.max <- x.tr.prop.max.noDL } else
    { x.tr.prop.max <- x.tr.prop.max.DL }
  }

  # --------------------------------------------------------------------------
  #  Check for reasonable truncation limits

  if ( !( (0 <= x.tr.prop.min) &
        (x.tr.prop.min < x.tr.prop.max) & 
        (x.tr.prop.max <= 1) ))
  { if (!is.na(logfile)) { sink(file=logfile, split=TRUE, append=TRUE) }
    cat("\n ++++++ [ProcessPropLimits] Invalid trunation limits:",
        "\n        x.tr.prop.min ", x.tr.prop.min,   
        "\n        x.tr.prop.max ", x.tr.prop.max,
        "\n ++++++ Execution stops +++++++",
        "\n")
    if (!is.na(logfile)) { sink() }
    stop("++++++ Invalid truncation limits ++++++")
  }    

  return(c(x.tr.prop.min, x.tr.prop.max))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ProcessPropLimits(NA, NA,
#                  0.40, 1.00,
#                  0.55, 1.00,
#                  NA)

#ProcessPropLimits(NA, NA,
#                  0.40, 1.00,
#                  0.55, 1.00,
#                  3)

#ProcessPropLimits(0.35, 0.65,
#                  0.40, 1.00,
#                  0.55, 1.00,
#                  NA)

#ProcessPropLimits(0.35, 0.65,
#                  0.40, 1.00,
#                  0.55, 1.00,
#                  3)

#ProcessPropLimits(0.65, 0.35,
#                  0.40, 1.00,
#                  0.55, 1.00,
#                  3)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# rm(list=ls())

# ****************************************************************************
#  F_Proportions_V9.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Function determining the proportions of observations contained in an 
#  interval with limits ilo, ihi.
#  ilo and ihi are indices in counts (e.g. from a histogram)


#  CHANGE HISTORY

#  29.09.2022 Parameter TI.left.only introduced
#  10.08.2022 TI definition changed. More TIs as a second line
#  28.01.2022 Additional requirement: TI must contain (deltap.hi-deltap.lo)
#             of data around mode
#             p.seq is used to exclude extreme bins
#             Ranges for QQ calculation developed here
#  18.01.2022
#  07.01.2022 Exclusion of TI candidates with extreme values
#  17.11.2021 Start from proportions()

# ============================================================================

Proportions <- function(counts, xbreaks.mode.idx, bins.n.min, prop.min, 
                        prop.max, TI.left.only, deltap.lo, deltap.hi, 
                        print.log.message=FALSE)
{ 
  # Proportions_V9
  # counts          vector of counts (typically from a histogram)
  # xbreaks.mode.idx index of the interval containing the mode of the 
  #                  distribution = index of the left breakpoint
  # bins.n.min       minimum number of bins in the interval
  # prop.min         minimal proportion included in the interval
  # prop.max         proportion in the interval must be < prop.max
  #                  Exception: prop.max == 1
  # deltap.lo        Proportion of data left of mode that must lie in TI
  # deltap.hi        Proportion of data right of mode that must lie in TI
  # ==========================================================================     

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   Start\n") }
  # ==========================================================================

  counts.n  <- length(counts)
  x.n       <- sum(counts)
  x.cdf     <- cumsum(counts)/x.n

  prop.max.w <- prop.max  #  working value
  if (prop.max == 1) { prop.max.w <- 1.1 }

  #  Try to find breaks with probabilities in the empirical cdf corresponding
  #  to p \in p.seq

  # p.seq <- seq(0.03, 0.97, by=0.02)    # bis 02.02.2022 48
  # p.seq <- seq(0.04, 0.96, by=0.04)    # 03.02.2022  24
  # p.seq <- seq(0.04, 0.96, by=0.03)    # 03.02.2022  31
  p.seq <- c(0.0, seq(0.04, 0.96, by=0.03), 1.00)    # 13.08.2022  33

  p.n   <- length(p.seq)  

  idx <- rep(NA, times=p.n)
  i <- 0
  for (p in p.seq)
  { i <- i + 1
    idx[i] <- which.min(abs(p - x.cdf))[1]  # [1] for the case of multiple
                                            #     minima
  }
  D <- data.frame(p.seq, idx, cdf=x.cdf[idx])

  #cat("\n [Proportions_V9] All indices of p.seq in counts\n")
  #print(D)

  #  Reduce to unique quantiles
  D <- D[ , -1]
  D <- unique(D)
  D <- D[order(D[ ,"cdf"]), ]
  p.n <- nrow(D)

  #cat("\n [Proportions_V9] Unique indices of p.seq in counts\n")  
  #print(D)
  
  xbreaks.mode.p <- x.cdf[xbreaks.mode.idx]

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   100\n") }
  # ==========================================================================
  # -------------------------------------------------------------------------
  #  New from here 10.08.22
  #  Generate all intervals with bin limits that have indices in D

  #  List all intervals of sufficient size with indices from D[ ,"idx"]
  #  ilo and ihi are indices in breaks
  #  If TI.left.only then consider only ilo == 1

  tab.names      <- c("ilo", "ihi", "prop.lo", "x.tr.prop", "prop.hi", 
                      "x.lt.tr.n", "x.tr.n", "x.ge.tr.n")
  tab0           <- matrix(NA, nrow=1, ncol=length(tab.names))
  colnames(tab0) <- tab.names

  #  ilo, ihi are indices in counts
  #  Dlo, Dhi are indices in D
  #  itab is index in tab
  #  tab lists interval indices in counts

  if (TI.left.only) 
  { D.lo.seq <- 1 } else
  { D.lo.seq <- 1:(p.n-1)}

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   200\n") }
  # ==========================================================================

  tab.n <- 0

  for (Dlo in D.lo.seq)
  {
    ilo <- D[Dlo, "idx"] 
    x.lt.tr.n <- 0
    if (ilo > 1) { x.lt.tr.n <- sum(counts[1:(ilo-1)]) }

    for (Dhi in (Dlo+1):p.n)
    { 
      ihi    <- D[Dhi, "idx"]

      x.ge.tr.n <- 0
      if (ihi < counts.n) { x.ge.tr.n <- sum(counts[(ihi+1):counts.n]) }
      x.tr.n <- sum(counts[ilo:ihi])

      if ( abs(x.lt.tr.n + x.tr.n + x.ge.tr.n - x.n) > 1.e-6 )
      {
        cat("\n [Proportions_V9] Error in count decomposition",
            "\n x.lt.tr.n + x.tr.n + x.ge.tr.n =", 
                x.lt.tr.n + x.tr.n + x.ge.tr.n,
            "\n                               x.n =", x.n,
            "\n")

        stop(" +++ [Proportions_V9] Error in count decomposition +++")
      }

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   300\n") }
  # ==========================================================================

      prop      <- x.tr.n / x.n
      bins.n    <- ihi - ilo + 1

      #  TIs must contain the required number of bins and proportion of cases
      if ( (bins.n >= bins.n.min) & 
           (prop.min <= prop) & (prop < prop.max.w)
         ) 
      { 
        # Put information into tab
     
        tab0[1, "ilo"]       <- ilo
        tab0[1, "ihi"]       <- ihi
        tab0[1, "prop.lo"]   <- x.lt.tr.n/x.n
        tab0[1, "x.tr.prop"] <- prop
        tab0[1, "prop.hi"]   <- x.ge.tr.n/x.n
        tab0[1, "x.lt.tr.n"] <- x.lt.tr.n
        tab0[1, "x.tr.n"]    <- x.tr.n
        tab0[1, "x.ge.tr.n"] <- x.ge.tr.n

        tab.n <- tab.n + 1
  
        if (tab.n > 1) { tab <- rbind(tab, tab0) } else
                       { tab <- tab0 }
      }
    }
  }

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   400\n") }
  # ==========================================================================

  # Mode contained in TI?
  contains.mode <- (tab[ , "ilo"] <= xbreaks.mode.idx) &
                   (xbreaks.mode.idx <= tab[ , "ihi"])

  # ==========================================================================
  # if (print.log.message) { cat("%%%   Proportions_V9   410\n") }
  # ==========================================================================

  #  Identify bin indices that contain deltap.lo below and deltap.hi
  #  above the mode
  #  Find probability (cdf) of mode bin (in counts)
  xbreaks.mode.p <- x.cdf[xbreaks.mode.idx]

  # ==========================================================================
  # if (print.log.message) { cat("%%%   Proportions_V9   420\n") }
  # ==========================================================================

  #  Find bin in counts which has cdf = xbreaks.mode.p - deltap.lo
  idx.lo <- which(x.cdf <= xbreaks.mode.p - deltap.lo)

  if (length(idx.lo) == 0)
  { # No such index
    idx.lo <- NA
  } else
  {
    idx.lo <- tail(idx.lo, 1)
  }

  # ==========================================================================
  # if (print.log.message) { cat("%%%   Proportions_V9   430\n") }
  # ==========================================================================

  #  Find bin in counts which has cdf = xbreaks.mode.p + deltap.hi
  idx.hi <- which(x.cdf >= xbreaks.mode.p + deltap.hi)

  # ==========================================================================
  # if (print.log.message) { cat("%%%   Proportions_V9   440\n") }
  # ==========================================================================

  if (length(idx.hi) == 0)
  { # No such index
    idx.hi <- NA
  } else
  {
    idx.hi <- idx.hi[1]
  }

  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   500\n") }
  # ==========================================================================

  # Add TI quality descriptor to tab
  # 3: interval (deltap.lo, deltap.hi) in TI (==> mode in TI), 
  # 2: mode in TI, but not the interval (deltap.lo, deltap.hi)
  # 1: mode not in TI

  TI.quality <- rep(NA, times=tab.n)
  TI.quality[!contains.mode] <- 1
  TI.quality[contains.mode]  <- 2
  TI.quality[(tab[ ,"ilo"] <= idx.lo) & (tab[ ,"ihi"] >= idx.hi)]  <- 3

  tab <- data.frame(tab, TI.qual=TI.quality)
  tab <- tab[order(-tab[ , "TI.qual"], -tab[ , "x.tr.prop"] ), ]

  #cat("\n [Proportions_V9] xbreaks.mode.idx", xbreaks.mode.idx,
  #    "\n                  xbreaks.mode.p  ", xbreaks.mode.p,
  #    "\n                  idx.lo, idx.hi  ", idx.lo, idx.hi,
  #    "\n")


  # ==========================================================================
  if (print.log.message) { cat("%%%   Proportions_V9   End\n") }
  # ==========================================================================

  return(tab)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#counts <- x.hist$counts
#breaks <- x.hist$breaks
#save(breaks, counts, file="Proportions.Testdata.RData")
#source("../func/F_RowToMatrix.R")

#load(file="Proportions.Testdata.RData")

#temp <- Proportions(counts, which.max(counts), 8, 0.6, 
#                        1.00, 0.15, 0.15)

#print(temp)

#cbind(data.frame(counts=counts, x.cdf=cumsum(counts)/sum(counts)))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# ****************************************************************************
#  F_q.PN.R

#  Truncated minimum chi square estimation

# (c) wwosniok@math.uni-bremen.de

#  CHANGE HISTORY
#  29.01.2021 Treatment of p=1 for lambda < 0 addded
#  28.01.2021 Error in negative branch corrected (Definition of T erroneously
#             deleted)
#  23.01.2021 Error in equation (1) of Freeman, Modarres corrected, see
#             PND_properties.tex / pdf
#             Note that for lambda < 0 the quantile function can attain 
#             indefinite values, because a power of a negative quantity arises. 
#  18.01.2021 Comments changed to English
# ============================================================================ 

q.PN <- function(p,lambda,mue,sigma,fastnull=1.e-10)
{
  #  Quantile of the PND, original scale (x)
  #  mue, sigma gelten innerhalb der PNV
  # siehe Freeman, S. 767
  #  Probability density function of the power normal distribution (PND) 
  #  Background: see 
  #  Freeman J, Modarres R. Inverse Box-Cox: the power normal 
  #  distribution. Stat Probabil Letters 2006;76; 764-72, particularly p 766

  #  INPUT   
  #  p       probability for which to determine the quantile, 0 < p < 1
  #  lambda  shape parameter. Real number. Special cases:
  #          1: PND is normal distribution with mean  and sd sigma
  #          0: PND is a logarithmic normal distribution with parameters
  #             mue and sigma
  #          If abs(lambda) < fastnull, lambda== 0 is assumed and the 
  #          intrinsic function plnorm is used
  #  mue     mean of the PND (of BoxCox(x, lambda))
  #  sigma   sd of the PND   (of BoxCox(x, lambda))
  #  fastnull (="nearly zero") absolute values < fastnull are treated as zero
   
  #  OUTPUT
  #  q.PN    Quantiles  

  # ==========================================================================
 
  q.PN <- rep(NA, times=length(p))
 
  if (lambda > fastnull)
  { 
    # cat("\n Positive lambda\n")
    T <- 1/(lambda*sigma) + mue/sigma
    V <- 1 - (1-p)*pnorm(T)    
    qnormV                  <- qnorm(V)
    
    #  qnormV  might be -Inf or Inf, set replacements to avoid overflow or NaN 
    qnormV[qnormV == -Inf]  <-  -10     #  -40 previously - why?
    qnormV[qnormV ==  Inf]  <-   10
    
    W <- lambda*(sigma*qnormV+mue) + 1     # kann durch numerische Fehler
                                               # negativ werden
    #  W might be < 0, avoid trouble

    W.gt.0        <- W > 0
    q.PN[W.gt.0]  <- (W[W.gt.0])^(1/lambda)
    q.PN[!W.gt.0] <- 0     
  }

  if (abs(lambda) < fastnull)
  { # LNV
    # cat("\n lambda = 0\n")
    q.PN  <- qlnorm(p,meanlog=mue,sdlog=sigma)
  } 

  if (lambda < -fastnull)
  { 
    # cat("\n Negative lambda\n")
    #  Calculation according to  F&M  (is wrong)
    # W <- lambda*(sigma*qnorm(p)+mue) + 1   # can get < 0
    # W.gt.0        <- W > 0
    # q.PN[W.gt.0]  <- (W[W.gt.0])^(1/lambda)
    # q.PN[!W.gt.0] <- 0     

    #  Corrected calculation, see PND_properties.tex / pdf  23.01.2021
    #  Holds for p in [0, 1). p = 1 and negative C need special treatment.
    q.PN <- rep(NA, times=length(p))

    T <- 1/(lambda*sigma) + mue/sigma
    A <- p * pnorm(-T)
    B <- lambda * sigma * qnorm(A) + lambda * mue
    C <- (1 + B)

    p.eq.1 <- (p > (1-fastnull)) 
    C.lt.0 <- (C < 0)    
    #  
    ok <- (!p.eq.1) & (!C.lt.0)
 
    q.PN[ok]  <- C[ok]^(1/lambda)    
    
    #  The final expression tends to zero^(1/lambda), which is 0, if p -> 1. 
    #  However, the sequence of q.PN for p -> 1 tends to Inf.
    if (any(p.eq.1)) { q.PN[p.eq.1] <- Inf } 

    #  Next condition seems to coincide with first condition
    #if (any(C.lt.0))
    #{ # This is an error condition
    #  cat("\n [q.PN] ++++++ No q.PN solution for  following p values",
    #      "\n        ++++++ lambda, mue, sigma",
    #      "\n        ++++++", lambda, mue, sigma,
    #      "\n")
    #  print(p[C.lt.0])
    #}  
  }

  return(unname(q.PN))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test

#source("F_BoxCox.R")
#source("F_BoxCoxInv.R")
#source("F_cdf.PN.R")
#source("F_pdf.PN.R")

#lambda <- -0.09989172
#mue    <- 11.25949
#sigma  <-  0.01440878
#pseq   <- seq(0, 1, by=0.1)
#x      <- q.PN(pseq, lambda, mue, sigma)

#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(pseq, x, x.pdf, x.cdf))


#lambda <- -0.1
#mue    <- 2.184875
#sigma  <-  0.3562688
#pseq   <- c(seq(0, 1, by=0.1), seq(0.91, 1.00, by=0.01), 
#            0.999, 0.9999, 0.99999, 0.999999)
#x      <- q.PN(pseq, lambda, mue, sigma)

#x.pdf  <- pdf.PN(x, lambda, mue, sigma)
#x.cdf  <- cdf.PN(x, lambda, mue, sigma)
#cat("\n x, x.pdf")
#print(cbind(pseq, x, x.pdf, x.cdf))

#q.PN(0.975, -1.000000e-01, 2.521592, 0.2969054)  # 12.30498  41.08872
#cdf.PN(41.08872, -1.000000e-01, 2.521592, 0.2969054)  # 

#Q <- c()
#lambda.seq <- seq(0.99, 1.01, by=0.001)
#for (lambda in lambda.seq)
#{ 
#  Q <- c(Q, q.PN(0.5, lambda, 140, 2.551027))
#}
#print(cbind(lambda.seq, Q))
#plot(lambda.seq, Q, type="o")
#abline(h=140)
#abline(v=1)

#for (lambda in seq(0.90, 1.10, by=0.01))
#{ cat("\n", lambda, cdf.PN(140, lambda, 140, 2.551027) ) }

# ************************************************************************
#  F_QQEstimateIni_V1.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Calculate initial values for TMC-theta by a modified QQ approach for
#  fixed lambda

#  CHANGE HISTORY

#  30.09.2022 Start
# =========================================================================

QQEstimateIni <- function(x, lambda, prop.lo, prop.hi, subset,
                          figA, print.log.message=FALSE)
{ 
  # ProcessPropLimits_V1
  # INPUT
  # x        data
  # lambda   parameter for BoxCox transformation of x, for title text only
  # prop.lo  proportion of data left of the data subset to use (the TI)
  # prop.hi  proportion of data right of the data subset to use (the TI)
  # subset   index vector for the part of x to use in estimation
  # figA     if not NA, show QQ plot. Graph window must exist.

  # OUTPUT
  # mue      ) estimated PND parameters
  # sigma    )
  # rsq      adjusted residual sum of squares for this solution

  # =======================================================================     

  # =======================================================================
  if (print.log.message) { cat("%%%   QQEstimateIni_V1   Start\n") }
  # =======================================================================

  #  QQnorm for whole data
  x.qq   <- qqnorm(x, plot.it=!is.na(figA), 
                   main=paste("lambda =", lambda ))
 
  #  Do regression in subset
  x.qq.lm      <- lm(x.qq$y[subset] ~ x.qq$x[subset])
  x.qq.lm.sum  <- summary(x.qq.lm) 

  if (!is.na(figA)) 
  { 
    points(x.qq$x[subset], x.qq$y[subset], col="blue")   
    abline(y.qq2.qq.lm.sum$coefficients[1:2], col="blue") 
  }  

  return(c(x.qq.lm.sum$adj.r.squared,
           x.qq.lm.sum$coefficients[1], x.qq.lm.sum$coefficients[2]) )
     
  # =======================================================================
  if (print.log.message) { cat("%%%   QQEstimateIni_V1   End\n") }
  # =======================================================================
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#n <- 100
## x <- qnorm((0:(n-1) + 0.5)/n, 5, 2)
#x <- rnorm(n, 5, 2)
#X <- cbind(BoxCox(x, 0), BoxCox(x, 0.25), BoxCox(x, 1))

#prop.lo <- 0
#prop.hi <- 0
#subset  <- rep(TRUE, times=n)
#figA    <- NA

#il <- 0
#for (lambda in c(0, 0.5, 1))
#{
#  il <- il+1
#  QQE <- QQEstimateIni(X[ ,il], lambda, prop.lo, prop.hi, subset,
#                          figA, print.log.message=FALSE)
#  print(QQE)
#}

#  Should be faster
#QQE.tab <- apply(X, MARGIN=2, FUN=QQEstimateIni,
#                 prop.lo=prop.lo, prop.hi=prop.hi, subset=subset, 
#                 figA=figA)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv#  F_QQPlotIMSS.R 
#  (c) wwosniok@math.uni-bremen.de
#
#  Standort: /user/wwosniok/TMC11/func
#  Testprogramm siehe C:\user\wwosniok\RLE2\archive\QQplot

#  QQ-Plot zur Berechung von RL in versuchter Annäherung an G. Hoffmann
#  QQ-Plot modifizierte Form ähnlich G. Hoffmann, mail vom 3.5.2017,
#  Trillium Normalizer

#  Verkürzte Fassung von F_QQPlotHW: nur modified Tukey und QQ-Plot, der
#  auf dessen Ergebnisse angewendet wird.. QQ-Parameterschätzung aus
#  allen Punkten, die aus modified Tukey hervorgehen. 
#  Keine Datentransformation

#  Kommentar zu F_QQPlotHW:
#  Estimate RL and parameters from the QQ plot following G. Hoffmann,
#  using the modified Tukey procedure to identify the unaffected part of the 
#  data
#  Background: 
#  Hoffmann et al 2015 points to 
#  Horn (CC 2001) (describes Tukey method, values outside (Q1- 1.5 IQR,
#  Q1 + 1.5 IQR) are outliers and removed,  no iteration, RLs are  calculated 
#  from the rest, RIs from the regression parameters plus compensation for 
#  using a reduced data set)
#  and Shaw (CB 2014) (flow chart, identification of linear portion by
#  piecewise regression, regression parameters provide distribution parameters
#  without any compensation)
#  In all cases, outliers are extreme points - no mixing of normal and outlier
#  populations.
#   

#  (c) of function modified.tukey: G. Hoffmann (?)
#  (c) the rest:  wwosniok@math.uni-bremen.de

#  CHANGE HISTORY

#  30.05.2022 Reduktion auf einfache Form, siehe oben
#  25.08.2020 Beschränkte Optimierung durch nlm
#  24.08.2020 Erste Schätzung via modified.tukey ergänzt um Kalkulation der
#             erzeugenden Verteilung unter Ausnutzung der Beziehung zwischen
#             Parametern der vollen und der gestutzten Verteilung 
#  06.07.2020 Treatment of tied values changed.
#  Verbesserungsvorschlag WW zur Schätzung: nach Tukey-Reduktion die 
#  erwartete Anzahl xc ausrechnen. Neuen QQ-Plot machen, dabei nur die 
#  reduzierten Daten plus fiktive Werte oben und unten entsprechend der 
#  geschätzten Anzahl von xc und den gekappten Anteilen oben und unten.
#  Dann Regression nur mit den gekappten Daten im QQ-Plot. 

#  05.07.2020 Umbau als Bestandteil für Manuskript IMSS
#  13.11.2018 Einrichtung aus F_QQPlotHoffmann.R in progexp

# ===========================================================================

# ===========================================================================

QQPlotIMSS <- function(x, RL1.p, RL2.p, figA, figB, fastnull=1.e-10)
{
  #  x      sorted data (no transformation done here)
  #  figA   window for Box plots, NA: no Box plot
  #  figB   window for QQ plot
  # ------------------------------------------------------------------------

  if (is.na(figA)) { print.boxplots <- FALSE } else
  { 
    print.boxplots <- TRUE
    dev.set(figA)
  }

  #  Reduce y. To deal appropriately with ties, add a name to each value
  y     <- x            # work on a copy of x that can be modified
  y.n   <- length(y)
  names(y) <- 1:y.n

  y.red <- modified.tukey(vNumeric=y, perc=2.5,log.mode=FALSE, n.cycles=0,
                          print.cycles=TRUE, print.boxplots=print.boxplots,
                          ytext="")

  #print(y.red)
  
  y.tr.lo <- min(y.red)
  y.tr.hi <- max(y.red)

  #  Use the names to identify the reduced dataset
  y.red.names <- as.numeric(names(y.red))
  y.red.lo    <- y.red.names[1] - 1           # number removed lo
  y.red.hi    <- y.n - tail(y.red.names, 1)   # number removed hi 
 
  y.red.lo.pct <- 100 * y.red.lo / y.n
  y.red.hi.pct <- 100 * y.red.hi / y.n
  abline(h=c(y.tr.lo,y.tr.hi),col="gray",lty=2)

  cat("Extreme values after Tukey reduction: ",y.tr.lo,y.tr.hi,"\n")

  #  Index vector
  red.idx <- rep(FALSE, times=y.n)
  red.idx[y.red.names[1] : tail(y.red.names, 1)] <- TRUE 

  dev.set(figB)
  if (is.na(figB)) { plot.it <- FALSE } else
  { 
    plot.it <- TRUE
    dev.set(figB)
  }

  #  Schätzung der Verteilungsparameter
  #  Variante 0: 
  #  MW und SD sind MW und SD des Tukey-reduzierten Datensatzes

  y.mue.qq0 <- mean(y.red)
  y.sig.qq0 <- sd(y.red)

  #  Schätzung der Verteilungsparameter
  #  Variante 1:
  #  x-Positionen aus dem Tukey-reduzierten Datensatz berechnet
  #  y-Positionen aus dem Tukey-reduzierten Datensatz berechnet
  #  Regression benutzt den gesamten Tukey-reduzierten Datensatz

  y.qq <- qqnorm(y.red,col="black", 
                 plot.it=plot.it,
                 main="Tukey-reduced data",                 
                 sub=paste("Originaldaten: n =",y.n,
                             "  reduziert: n =",length(y.red),
                             " li:",y.red.lo.pct,
                             " re:",y.red.hi.pct),
                   cex.sub=0.8)
  abline(h=c(y.tr.lo,y.tr.hi),col="gray",lty=2)

  #  Regression im reduzierten Bereich,
  y.qq.lm1 <- lm(y.qq$y ~ y.qq$x)
  y.qq.lm1.sum <- summary(y.qq.lm1)
  cat("Variante 1\n")
  print(names(y.qq.lm1.sum))
  print(y.qq.lm1.sum)

  abline(c(y.mue.qq0, y.mue.qq0),lwd=2, col="orange")
  abline(y.qq.lm1$coefficients,lwd=1, col="green4")
                    
  return(c(y.mue.qq0=y.mue.qq0,
           y.sig.qq0=y.sig.qq0,
           y.mue.qq1=unname(y.qq.lm1$coefficients[1]),
           y.sig.qq1=unname(y.qq.lm1$coefficients[2]),
           rsq=unname(y.qq.lm1.sum$adj.r.squared),
           x.tr.n=length(y.red)) )
}

# ===========================================================================
# ****************************************************************************
#  F_QQW.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Calculate a QQ-plot based RL estimate
#  Method similar to the modified Tukey based method by G. Hoffmann (?),
#  but modified by WW in several respects.
 
#  #  CHANGE HISTORY
#  17.12.2020 Emergency action if only nonacceptable solutions changed
#  13.12.2020 prev.acc.lo/hi added to call
#  07.12.2020 bins.n.min replaced by bins.n.min.act
#  05.12.2020 Result plot (J) done by standard function
#  01.12.2020 x.red.hist removed from call
#  28.11.2020 Calculation of estimated prevalences moved to CalcPre()
#  16.11.2020 Call of TIQQ changed
#             Calculation of x.poly changed 
#  14.11.2020 bins.n.min added to call and to TIQQ( )
#  08.10.2020 Start
#
# ==============================================================================

QQW <- function(x, x.red, round.unit, lambda.seq, 
                x.hist, x.kde, x.kde.mode, x.kde.mode.idx, x.val, 
                x.unaffected.lo, x.unaffected.hi,   
                n.per.bin.min, bins.n.min, x.tr.bins.min, 
                counts.range, x.tr.prop.range, 
                x.tr.prop.min, x.tr.prop.max, prev.acc.lo, prev.acc.hi, 
                xsupp,
                x.Q1, x.Q2, RL1.p, RL2.p,
                lambda.gen, mue.gen, sig.gen,
                df.est,df.con,
                l.fact, p.fact, r.fact, w.fact, 
                xlabel, subtitle,
                figD, figE, figF, figG, figH, 
                figI, figI.file, 
                figJ, figJ.file, figtype,
                x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                x.clip.type, par.las, par.tcl,
                bordercol2, histcol2, kdecol, qqrawcol, polycol, gencol, 
                difcol, denlty,
                print.log.message)
{
  #  INPUT 
  #  x              data vector
  #  round.unit
  #  lambda.seq     candidate lambda values 
  #  x.hist         histogram that will later be used by tmc, breaks are used
  #                 here as candidate limits of the truncation interval
  #  x.tr.bins.min  min number of bins in the truncation interval (TI), 
  #  x.tr.prop.min  min proportion of values in the TI
  #  x.tr.prop.max  max proportion of values in the TI
  #  lambda.gen     if simulated data: lambda for generation 
  #  mue.gen        if simulated data: mue for generation 
  #  sig.gen        if simulated data: sig for generation
  #  xlabel         label for analyte 
  #  subtitle       subtitle in FigA 
  #  figD           in TIQQ: BC transformed collapsed histogram
  #  figE           in TIQQ: QQ plot for BC transformed data, polygon 
  #                       approximation of the QQ plot
  #  figF           in TIQQ: All QQ plot based parameter estimates
  #  figG           in TIQQ: All QQ plot based RL estimates
  #  figH           in TIQQ: QQ plot based relative errors in RL estimates
  #  figI           Collapsed histogram and QQ plot based fits
  #  bordercol2     colour for histogram borders, collapsed
  #  histcol2       colour for bin area, initial, collapsed
  #  kdecol         colour for kde
  #  qqrawcol       symbol color in raw QQ plot 
  #  polycol        symbol color in untied QQ plot
  #  gencol         colour for plot objects showing values from generation

  #  OUTPUT 
  #  
  # ---------------------------------------------------------------------------

  #cat("\n *** Start QQW() \n")

  x.n <- length(x)

  #  The original QQ plot method as in Hoffmann 1963 and others is done
  #  by @@@@ elsewhere (execute this programm on the raw data)   
 
  #  Loop over lambda sequence
  #
  ilambda <- 0
  for (lambda in lambda.seq)
  { ilambda <- ilambda + 1 
    
    #  Transform data, histogram breaks and kde
    y        <- BoxCox(x, lambda)
    y.red    <- BoxCox(x.red, lambda)
    y.breaks <- BoxCox(x.hist$breaks, lambda)

    y.red.hist <- hist(y.red, breaks=y.breaks, right=FALSE, plot=FALSE)
    y.kde      <- TransformKDE(x.kde, BoxCox(x.kde$x, lambda))

    # TIs must contain the transformed value of x.kde.mode, not the mode
    # of y.red. Plots contain this transformed mode.
    y.kde.mode <- BoxCox(x.kde.mode, lambda)
                    
    #  Changed: Find the bin (in y.red.hist$counts) containing the 
    #  the transformation of x.kde. This concentrates TI intervals around
    #  the mode of the x distribution on the original scale. 

    # y.red.hist.mode.idx is the same as x.kde.mode.idx
    
    # Get a proposal for reasonable truncation intervals, given lambda

    tab1.qqw <- TIQQ(x.hist, x.kde.mode, 
                     y.red, y.red.hist, y.kde, y.kde.mode, 
                     x.kde.mode.idx,                       
                     lambda, round.unit, kernel.n, bins.n.min, 
                     x.tr.bins.min, x.tr.prop.min, x.tr.prop.max,
                     df.est,df.con,
                     x.Q1, x.Q2, RL1.p, RL2.p,
                     l.fact, p.fact, r.fact, w.fact, 
                     xlabel, subtitle,
                     figD, figE, figF, figG, figH,
                     lambda.gen, mue.gen, sig.gen,
                     histcol2, bordercol2, kdecol, qqrawcol, polycol, gencol,
                     print.log.message)

    #  There may be no proposal
    if (length(tab1.qqw) == 1)
    {
      #  No result from TIQQ

      cat("\n [QQW] TIQQ without result\n")
      tab.qqw <- NA
    } else
    {
      #  There is a result
    
      tab1.qqw    <- RowToMatrix(tab1.qqw)
      tab1.qqw.n <- nrow(tab1.qqw)

      #  CalcPrev provides general names, no reference to calculating method
      tab2.qqw.names <- c("n.l", "prev.l", "n.c", "prev.c", 
                        "n.r", "prev.r", "neg.prev.sum")
      tab2.qqw       <- matrix(NA, nrow=tab1.qqw.n, ncol=length(tab2.qqw.names))
      colnames(tab2.qqw) <- tab2.qqw.names

      #  Next calculation by row 
      for (i in 1:tab1.qqw.n)
      {
        x.tr.n <- sum( (tab1.qqw[i, "x.tr.lo"] <= x) & 
                     (x < tab1.qqw[i, "x.tr.hi"])  )
        x.lt.tr.n <- sum(x <   tab1.qqw[i, "x.tr.lo"])
        x.ge.tr.n <- sum(x >=  tab1.qqw[i, "x.tr.hi"])

        tab2.qqw[i, ] <- CalcPrev(x.n, 
                                x.tr.n, 
                                x.lt.tr.n,
                                x.ge.tr.n, 
                                tab1.qqw[i, "x.tr.lo"], tab1.qqw[i, "x.tr.hi"], 
                                lambda, 
                                tab1.qqw[i, "mue.qqw"], 
                                tab1.qqw[i, "sigma.qqw"])
      }

      #  Make names in tab2.qqw method-specific. Check with first definitions a 
      #  few  lines higher.
      colnames(tab2.qqw) <- c("n.l.qqw", "prev.l.qqw", 
                            "n.c.qqw", "prev.c.qqw",
                            "n.r.qqw", "prev.r.qqw", "neg.prev.sum.qqw")
   
      if (ilambda==1)
      {
        tab.qqw <- data.frame(tab1.qqw, tab2.qqw, stringsAsFactors=FALSE)
      }  else
      {
        tab.qqw <- rbind(tab.qqw,
                       data.frame(tab1.qqw, tab2.qqw, stringsAsFactors=FALSE))
      }
    }
  }

  if (length(tab.qqw) > 1)
  {

    # --------------------------------------------------------------------------
    #  Decide for the best (some of the best) TIQQ result

    #cat("\n [QQW] Initial result from TIQQ, sorted by lambda, start, end\n")
    #print(tab.qqw[order(tab.qqw[ ,"lambda.qqw"],
    #                    tab.qqw[ ,"x.tr.lo"],
    #                    tab.qqw[ ,"x.tr.hi"]),  ])

    #cat("\n [QQW] Initial result from TIQQ, sorted by RL1, RL2 \n")
    #print(tab.qqw[order(tab.qqw[ ,"x.RL1.qqw"],
    #                    tab.qqw[ ,"x.RL2.qqw"]),  ])
  
    #  If necessary, turn vector into matrix
    tab.qqw    <- RowToMatrix(tab.qqw)

    #  Needed if no acceptable solutions
    tab.qqw.backup <- tab.qqw

    #  Provide overview over estimates
    #temp <- matrix(c(min(tab.qqw.backup[ , "prev.l.qqw"]),
    #               max(tab.qqw.backup[ , "prev.l.qqw"]),
    #               min(tab.qqw.backup[ , "prev.c.qqw"]),
    #               max(tab.qqw.backup[ , "prev.c.qqw"]),
    #               min(tab.qqw.backup[ , "prev.r.qqw"]),
    #               max(tab.qqw.backup[ , "prev.r.qqw"]) ), byrow=TRUE,
    #               nrow=3, ncol=2)
    #dimnames(temp) <- list(c("l", "c", "r"), c("min", "max"))
    #cat("\n Characteristics of QQW prevalence estimates\n")
    #print(temp)

    #  Remove estimates with implausible prevalence estimates
    #  Check prev.l estimates
    tab.qqw.n0  <- nrow(tab.qqw)
    ok <- (prev.acc.lo <= tab.qqw[ ,"prev.l.qqw"]) & 
          (tab.qqw[ ,"prev.l.qqw"] <= prev.acc.hi)
    tab.qqw.n <- sum(ok)
    tab.qqw    <- tab.qqw[ok, ] 
    tab.qqw    <- RowToMatrix(tab.qqw)

    #  Check for nothing more left  
  
    if (tab.qqw.n > 0)
    { #  Check prev.c estimates
      tab.qqw.n0  <- nrow(tab.qqw)
      ok <- (prev.acc.lo <= tab.qqw[ ,"prev.c.qqw"]) & 
            (tab.qqw[ ,"prev.c.qqw"] <= prev.acc.hi)
      tab.qqw.n <- sum(ok)
      tab.qqw   <- tab.qqw[ok, ] 
      tab.qqw   <- RowToMatrix(tab.qqw)
    }
  
    if (tab.qqw.n > 0)
    { #  Check prev.r estimates
      tab.qqw.n0  <- nrow(tab.qqw)
      ok <- (prev.acc.lo <= tab.qqw[ ,"prev.r.qqw"]) & 
          (tab.qqw[ ,"prev.r.qqw"] <= prev.acc.hi)
      tab.qqw.n <- sum(ok)
      tab.qqw    <- tab.qqw[ok, ] 
      tab.qqw    <- RowToMatrix(tab.qqw)
    }
  
    #  Emergency action, if no results are left
  
    if (tab.qqw.n == 0)
    { cat("\n ++++++  Method QQW found no acceptable solution",
          "\n ++++++  Modify prev.acc.lo and / or prev.acc.hi in",
          "\n ++++++  TMC_seg015_DefaultSettings.R",
          "\n ++++++  Below: the first up to 10 prevalence estimates by QQW,",
          "\n ++++++  sorted by the sum of negative prevalences and rsq",
          "\n\n")
      tab.qqw.backup <- tab.qqw.backup[order(-tab.qqw.backup[ ,"neg.prev.sum.qqw"],
                                           -tab.qqw.backup[ ,"r2"]), ] 
      imax <- min(10, nrow(tab.qqw.backup))
      print(tab.qqw.backup[1:imax, c("x.RL1.qqw", "x.RL2.qqw",
                                 "prev.l.qqw", "prev.c.qqw", 
                                 "prev.r.qqw", "neg.prev.sum.qqw", "r2")])

      cat("\n ++++++  Below: the first up to 10 prevalence estimates by QQW,",
          "\n ++++++  sorted by rsq",
          "\n ++++++  First solution willl be taken",
          "\n\n")
      tab.qqw.backup <- tab.qqw.backup[order(-tab.qqw.backup[ ,"r2"]), ]

      print(tab.qqw.backup[1:imax, c("x.RL1.qqw", "x.RL2.qqw",
                                 "prev.l.qqw", "prev.c.qqw", 
                                 "prev.r.qqw", "neg.prev.sum.qqw", "r2")])

      tab.qqw <- tab.qqw.backup
      tab.qqw.n <- nrow(tab.qqw)

    } else
    { 
      #  No problems so far
      rm(tab.qqw.backup)
    }  

    #  Sort by decreasing r2
    tab.qqw    <- tab.qqw[order(-tab.qqw[ ,"r2"]), ]

    #  Sort by increasing opt.crit
    #tab.qqw    <- tab.qqw[order(tab.qqw[ ,"opt.crit"]), ]
  
    #cat("\n [QQW] tab.qqw before reduction\n")
    #print(tab.qqw[1:30, ])

    #  Reduce to the first at most ... candidates.
    tab.qqw.n <- min(5, tab.qqw.n)
    tab.qqw   <- tab.qqw[1:tab.qqw.n, ]
    tab.qqw   <- RowToMatrix(tab.qqw)

    # -------------------------------------------------------------------------
    #  Add asymptotic confidence intervals. Calculation uses the actually
    #  used numbers, which may come from a reduced histogram
    for (i in 1:tab.qqw.n)
    {
      xc.RL1.qqw.CI <- CIQuant.PNV(tab.qqw[i, "x.tr.n"], RL1.p, 
                                 tab.qqw[i, "x.RL1.qqw"],
                                 tab.qqw[i, "lambda.qqw"], 
                                 tab.qqw[i, "mue.qqw"], 
                                 tab.qqw[i, "sigma.qqw"],
                                 alpha=alpha)
      tab.qqw[i, "x.RL1.cilo"] <- xc.RL1.qqw.CI[1]
      tab.qqw[i, "x.RL1.cihi"] <- xc.RL1.qqw.CI[2]

      xc.RL2.qqw.CI <- CIQuant.PNV(tab.qqw[i, "x.tr.n"], RL1.p, 
                                 tab.qqw[i, "x.RL2.qqw"],
                                 tab.qqw[i, "lambda.qqw"], 
                                 tab.qqw[i, "mue.qqw"], 
                                 tab.qqw[i, "sigma.qqw"],
                                 alpha=alpha)
      tab.qqw[i, "x.RL2.cilo"] <- xc.RL2.qqw.CI[1]
      tab.qqw[i, "x.RL2.cihi"] <- xc.RL2.qqw.CI[2]
    }

    # -------------------------------------------------------------------------
    if (!is.na(figI))
    { dev.set(figI)
      bringToTop(figI)

      #  Show QQ plot for optimal truncation interval
      lambda.qqw.opt <- tab.qqw[1, "lambda.qqw"]

      qqnorm(BoxCox(x.red, lambda.qqw.opt),
           main=paste("Truncation interval according to QQ plot, lambda =",
                      lambda.qqw.opt),
           sub=subtitle, cex.sub=0.7)
      abline(h=BoxCox(tab.qqw[1, "x.tr.lo"], lambda.qqw.opt), col="blue")   
      abline(h=BoxCox(tab.qqw[1, "x.tr.hi"], lambda.qqw.opt), col="blue")

      abline(c(tab.qqw[1, "mue.qqw"], tab.qqw[1, "sigma.qqw"]), col="red")   
      abline(h=y.kde.mode, col=kdecol, lty=2)   

      #  savePlot(file=figI.file, type=figtype)  #  01.05.2021
    }

    # -------------------------------------------------------------------------
    #  Plot figJ moved to calling programme  01.05.2021
    #  @@@ move to calling program?

  }  

  #cat("\n *** End QQW() \n")
  return(tab.qqw=tab.qqw) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_Quantile.R
#  
#  (c) wwosniok@math.uni-bremen.de


#  19.05.2021 WW  Check for empty input vector
#  16.03.2021 WW  Determination of round.unit is always necessary
#  14.03.2021 WW  Linear interpolation instead of parabola interpolation
#                 Continuity connection added
#  19.12.2020 WW  New interpretation of input data: x is considered
#                 as a rounded value ==> for a given x, the true value
#                 lies in the interval [x-ru/2, x+ru/2), where ru is the
#                 rounding unit. Due to rounding, ties are likely, therefore
#                 quantiles are calculated from the polygon connecting
#                 (x[1]-ru/2, 0), (x[1], cdf[1]),  (x[2], cdf[2]), ...
#                 (x[n], cdf[n]), where cdf[i] is the # of x values <= x[i].
#                 x[1]-ru/2 may be negative - no precaution for using negative
#                 results so far.  
#  09.07.2019 WW  Behandlung des ersten und letzten Intervalls geändert.
#                 Interpolation verwendet jetzt die korrekte Definition
#                 der empirischen Verteilungsfunktion. Kein Problem mehr,
#                 wenn p-Werte unterhalb / oberhalb der empirischen ps liegen.
#  17.08.2018 WW  Stetigkeitskorrektur eingeführt. Die interpolierende
#                 Linie verbindet die Mittelpunkte der Treppenstufen,
#                 nicht die Eckpunkte.
#  22.11.2017 WW  Einrichtung

# ===========================================================================
Quantile <- function(x, probs=c(0.25,0.50,0.75), na.rm=FALSE, type=0)
{ #  Berechnet Quantile durch Interpolation der empirischen 
  #  Verteilungsfunktion, auch bei Bindungen (das macht quantile() nicht!)
  #  Mit Stetigkeitskorrektur.
  #  EINGABE
  #  x     - Daten
  #  probs - Wahrscheinlichkeiten \in [0,1] zu den gewünschten Quantilen
  #  na.rm - falls TRUE, werden fehlenden Werte aus x entfernt
  #  type  - ohne Bedeutung, nur zur Kompatitibilität mit Aufruf von quantile()
  #  AUSGABE
  #  Quantile mit den Wahrscheinlichkeiten in % als Namen
  # ==========================================================================

  probs.n <- length(probs)

  if (length(x) > 1)
  { # At least 2 values, try execution 
    
    #  Remove missing vales, if required
    if (na.rm) { x <- x[!is.na(x)] }


    #  If there still are NAs, stop, otherwise execute
    if (any(is.na(x)) & !na.rm ) 
    { 
      cat("In Quantile(): x has NA and na.rm == FALSE \n")
      q <- rep(NA, times=probs.n) 
    } else
    { #  No (more) NAs
      x.table <- table(x)
      x.val   <- as.numeric(names(x.table))
      x.val.n <- length(x.val)
      x.val.cdf <- cumsum(x.table)/length(x)

      q <- rep(NA,times=probs.n)
      if (x.val.n == 1)
      {
        # Only 1 unique value  
        q <- x.val 
      } else
      { #  x.val.n > 1
        #  More than 1 unique value
        #  Calculate coordinates for interpolation

        intx <- x.val
        inty <- x.val.cdf

        #  Continuity correction
        round.unit <- min(diff(intx))
        intx       <- intx + round.unit/2 
      
        # If a quantile for a probability < x.val.cdf[1] is requested,
        # we need to construct x[1] - ru
        if (any(probs < x.val.cdf[1]))
        { # yes, find the rounding unit and  add left limit
          intx <- c(intx[1] - round.unit/2, intx)            
          inty <- c(0, inty)
        }   

        for (i in 1:probs.n)
        { 
          #  Interpolate between the shifted jump positions
         
          q[i] <- Interpolx(intx, inty, probs[i]) 
        }
      }     #  if (x.val.n == 1)
    }       #  if (any(is.na(x)) & !na.rm ) 
  
  } else
  {
    #  Too few data
    q <- rep(NA, times=probs.n) 

  }         # if (length(x) > 1)

  names(q) <- paste(100*probs,"%",sep="")
  return(q)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test see TMC6/dev/Test_Quantile.R 
#
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_R2o.R

R2o <- function(ttheta,RB)
{ #  Auf ganz R  transformierte Parameter --> beschränkte Originalskala 
  #  16.12.2017 

  theta.n <- length(ttheta)
  theta   <- rep(NA,times=theta.n)
  
  for (k in 1:theta.n) 
  { 
    #  Schutz vor Überlauf
    ttheta.eff <- min(ttheta[k],50)
    ttheta.eff <- max(-50,ttheta.eff)
    theta[k] <- RB[k,1] + (RB[k,2]-RB[k,1])/(1+exp(-ttheta.eff))
  }

  #cat("R2o: ttheta,theta\n")
  #print(cbind(ttheta,theta))

  return(theta)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test
#RB <- matrix(c(-1.e-10,1+1.e10),nrow=1,ncol=2)
#theta     <- 1-1.e-5
#ttheta    <- o2R(theta,RB,fastnull=1.e-10)
#theta.rec <- R2o(ttheta,RB)
#cat(theta,ttheta,theta.rec,"\n")

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_reldist.R
#
#  (c) wwosniok@math.uni-bremen.de  
#
#  Calculate relative distance between xtest and xref
#
#  04.01.2020

# =============================================================================
reldist <- function(xtest,xref)
{ 
  temp <- (xtest-xref)/xref
  reldist <- sqrt(temp %*% temp)
  return(reldist)
} 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_RL.age.spline.R

#
#  Show course of RL1, RL2 estimated by method RL.age.meth, vs 
#  mean of age group
#  Smooth by spline function
#  Separation by sex
#  Show tolerance intervals, if present, otherwise asymptotic CIs
#
#  22.03.2022  age.limits removed from call
#  01.03.2022  Call changed (RL.age.meth moved)
#  28.02.2022  y label added as input
#  08.05.2021  Estimation method added as call parameter
#  30.04.2021  x.clip.@@ correctly processed
#              par.tcl and par.las as parameters   
#  02.12.2020  Names x.RL, ..  changed to RL, marker for unsafe results
#              changed to "!!"
#  21.04.2020  ylab changed back to 'Reference limits'
#  24.01.2020  Indicator for successful plotting addded to output
#  10.01.2020  External parameters for the RL axis
#  08.01.2020  Var names changed for TMC4
#  12.12.2019  Mark plots in gray, if solution type != 2
#  26.07.2019  Plotten auch, wenn errcode nicht " "
#  08.05.2019  Plotten nur, wenn errcode nicht " "
#              Definition des errcode in analysis geändert
#  17.04.2019  Plotten auch dann, wenn errrcode nicht " "
#              (Überlegen, ob wackelige Schätzungen irgendwie markiert
#               werden können)  
#  02.03.2019  ylab geändert
#  02.05.2018
#  
#  --------------------------------------------------------------------------

RL.age.spline <- function(fig,gtab,sexcode,Verlauf.Meth,shift,
                           xplotmin,xplotmax,yplotmin,yplotmax,
                           yplotby1,yplotby2,  
                           xlabel, ylabel, RL.age.meth,
                           x.clip.min, x.clip.max, x.clip.by1, x.clip.by2, 
                           par.tcl, par.las,
                           outname.stra,methcol,methpch,
                           arrow.length=0.03,newplot=TRUE, 
                           unsafe.col="dimgray")
{ #  
  #  fig
  #  figfile
  #  gtab
  #  sexcode       
  #  Verlauf.Meth
  #  xplotmin
  #  xplotmax
  #  yplotmin
  #  yplotmax
  #  outname.stra
  #  methcol
  #  methpch
  #  age.limits
  #  arrow.length  arrow length in inches! Default figure width is 7 in

  #  Additional style parameters
  arrlwd <- 1    #  line width for arrows (confidence bars) 
  pchlwd <- 1    #  line width for point symbols 
  spllwd <- 1    #  line width for spline lines 
  pchcex <- 1    #  point symbol size
     
  dev.set(fig)
  bringToTop(fig)

  #  Select data. Note that some values are formatted (=> of character type).
  errcode  <- gtab[ ,"errcode"]
  age.mea  <- as.numeric(gtab[ ,"Age.mea"])
  RL1.hat  <- as.numeric(gtab[ ,"x.RL1"])
  RL2.hat  <- as.numeric(gtab[ ,"x.RL2"])
  RL1.cilo <- as.numeric(gtab[ ,"x.RL1.cilo"])
  RL1.cihi <- as.numeric(gtab[ ,"x.RL1.cihi"])
  RL1.tilo <- as.numeric(gtab[ ,"x.RL1.tilo"])
  RL1.tihi <- as.numeric(gtab[ ,"x.RL1.tihi"])
  RL2.cilo <- as.numeric(gtab[ ,"x.RL2.cilo"])
  RL2.cihi <- as.numeric(gtab[ ,"x.RL2.cihi"])
  RL2.tilo <- as.numeric(gtab[ ,"x.RL2.tilo"])
  RL2.tihi <- as.numeric(gtab[ ,"x.RL2.tihi"])

  use     <- (!is.na(age.mea)) & 
             (!is.na(RL1.hat)) & (!is.na(RL2.hat))
  n       <- sum(use)

  plot.ok   <- (n >= 2)
  spline.ok <- (n >= 4 )

  if (plot.ok)
  { #  Plotting is worthwhile
    #  Select data

    errcode  <- errcode[use]
    age.mea  <- age.mea[use]
    RL1.hat  <- RL1.hat[use]
    RL2.hat  <- RL2.hat[use]
    RL1.cilo <- RL1.cilo[use]
    RL1.cihi <- RL1.cihi[use]
    RL1.tilo <- RL1.tilo[use]
    RL1.tihi <- RL1.tihi[use]
    RL2.cilo <- RL2.cilo[use]
    RL2.cihi <- RL2.cihi[use]
    RL2.tilo <- RL2.tilo[use]
    RL2.tihi <- RL2.tihi[use]

    if (newplot)
    { 
      #  Punkte eintragen
      #  RL1
      #  No automatic x legend

      par(las=par.las)
      par(tcl=par.tcl)
      par(xaxt="n") 
      if (!is.na(x.clip.by1))  
      {  # User values are given 
        par(xaxt="n") 
        xplotmin.eff <- x.clip.min
        xplotmax.eff <- x.clip.max
      } else 
      {  # User values are not given 
        xplotmin.eff <- xplotmin
        xplotmax.eff <- xplotmax
      }
  
      if (!is.na(yplotby1))    { par(yaxt="n") }  # User values are given 
 
      #  Standard y label
      ylab <- paste("Reference limits, method =", RL.age.meth)

      #  User defined label
      if (!is.na(ylabel)) ylab <- ylabel

      plot(age.mea+shift,RL1.hat,type="p",
          col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd,
          xlim=c(xplotmin.eff,xplotmax.eff),
          ylim=c(yplotmin,yplotmax),
          xlab=xlabel,
          ylab=ylab )

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
      
      #  User-defined tickmarks for x
      if (!is.na(x.clip.by1))  
      {
        par(xaxt="s")
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by1),labels=TRUE)
        axis(1,at=seq(x.clip.min,x.clip.max,by=x.clip.by2),labels=FALSE)
      }

      #  User-defined tickmarks for y
      if (!is.na(yplotby1))  
      { 
        par(yaxt="s") 
        axis(2,at=seq(yplotmin,yplotmax,by=yplotby1),labels=TRUE)
        if (!is.na(yplotby2))
        {
          axis(2,at=seq(yplotmin,yplotmax,by=yplotby2),labels=FALSE)
        }
      }

    }  else
    {
      #  Plot existiert schon
      #  Punkte eintragen
      #  RL1 

      points(age.mea+shift,RL1.hat,type="p",
             col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

      # Mark unsafe estimates with errcode == "!!" in extra color (gray) 
      points(age.mea[errcode=="!!"]+shift,RL1.hat[errcode=="!!"],
             type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

      #  Konfidenz- / Toleranzintervslle
      if (all(is.na(RL1.tilo)))
      { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
        #cat("[SplineVerlauf3] (1) Asymptotisches CI für RL1 verwendet\n")
        ci.type <- "CI"
        arrows(age.mea+shift,RL1.cilo,
               age.mea+shift,RL1.cihi,
               code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe  
        arrows(age.mea[errcode=="!!"]+shift,RL1.cilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.cihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }  else
      { #  Bootstrap-Ergebnis vorhanden - nehmen
        # cat("[SplineVerlauf3] (2) Bootstrap-TI für RL1 verwendet\n")
        ci.type <- "TI"
        arrows(age.mea+shift,RL1.tilo,
               age.mea+shift,RL1.tihi,
               code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

        #  Mark CI in gray, if solution unsafe
        arrows(age.mea[errcode=="!!"]+shift,RL1.tilo[errcode=="!!"],
               age.mea[errcode=="!!"]+shift,RL1.tihi[errcode=="!!"],
               code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
      }
    }   # newplot

    #  Punkte eintragen
    #  RL2 

    points(age.mea+shift,RL2.hat,type="p",
           col=methcol,pch=methpch,cex=pchcex,lwd=pchlwd)

    # Mark unsafe estimates  in extra color (gray) 
    points(age.mea[errcode=="!!"]+shift,RL2.hat[errcode=="!!"],
           type="p",col=unsafe.col,pch=methpch,cex=pchcex,lwd=pchlwd) 

    #  Konfidenz- / Toleranzintervslle
    if (all(is.na(RL2.tilo)))
    { #  Kein Bootstrap-Ergebnis - asymptotische Intervalle nehmen
      #cat("[SplineVerlauf3] (3) Asymptotisches CI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.cilo,
             age.mea+shift,RL2.cihi,
             code=3,angle=90,length=arrow.length,col=methcol,lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.cilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.cihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)

    }  else
    { #  Bootstrap-Ergebnis vorhanden - nehmen
      #cat("[SplineVerlauf3] (4) Bootstrap-TI für RL2 verwendet\n")
      arrows(age.mea+shift,RL2.tilo,
             age.mea+shift,RL2.tihi,
             code=3,angle=90,length=arrow.length,col="magenta",lwd=arrlwd)

      #  Mark CI in gray, if solution type != 2  
      arrows(age.mea[errcode=="!!"]+shift,RL2.tilo[errcode=="!!"],
             age.mea[errcode=="!!"]+shift,RL2.tihi[errcode=="!!"],
             code=3,angle=90,length=arrow.length,col=unsafe.col,lwd=arrlwd)
    }

    # Legende am Ende, wenn CI-Typ bekannt ist

  } # if (plot.ok)

  # .........................................................................
  #  Spline-Berechnung

  if (spline.ok)
  { #  Spline ist möglich  

    age.out <- seq(Floor(age.mea[1],1),Ceiling(age.mea[length(age.mea)],1),
                         by=1)

    # RL1.sspline.gcv <- smooth.spline(age.mea, RL1.hat)

    if (( 4 <= n) & (n <= 10)) { df.user <- Round(0.9 * n, 1) }
    if ((11 <= n) & (n <= 12)) { df.user <- Round(0.8 * n, 1) }
    if ((13 <= n) & (n <= 14)) { df.user <- Round(0.6 * n, 1) }
    # if ((15 <= n)            ) { df.user <- RL1.sspline.gcv$df }
    if ((15 <= n)            ) { df.user <- Round(0.5 * n, 1) }

    RL1.sspline.df      <- smooth.spline(age.mea, RL1.hat,df=df.user)
    RL1.sspline.df.pred <- predict(RL1.sspline.df,age.out)

    #  Spline eintragen
    lines(RL1.sspline.df.pred$x+shift,RL1.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

 
    # .........................................................................
    #  RL2 

    RL2.sspline.df      <- smooth.spline(age.mea, RL2.hat,df=df.user)
    RL2.sspline.df.pred <- predict(RL2.sspline.df,age.out)

    #  Spline eintragen
    lines(RL2.sspline.df.pred$x+shift,RL2.sspline.df.pred$y,
          col=methcol,lwd=spllwd)

    return(list(RL.emp=data.frame(age=age.mea,RL1=RL1.hat,RL2=RL2.hat,
                                  RL1.cilo=RL1.cilo,RL2.cilo=RL2.cilo,
                                  RL1.tilo=RL1.tilo,RL2.tilo=RL2.tilo),
                RL.smo=data.frame(age=age.out,
                                  RL1.hut=RL1.sspline.df.pred$y,
                                  RL2.hut=RL2.sspline.df.pred$y),
                plot.ok=plot.ok,
                spline.ok=spline.ok))

  } else  #  if spline.ok ..
  { 
    #  Spline nicht möglich, da zu wenig Stützpunkte
    cat(  "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++",
        "\n  No spline smoothing for sex =",sexcode,
        ", because < 4 support points",
        "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

    return(list(RL.emp=data.frame(age=NA, RL1=NA,RL2=NA,
                                  RL1.cilo=NA, RL2.cilo=NA,
                                  RL1.tilo=NA, RL2.tilo=NA),
                RL.smo=data.frame(age=NA,
                                  RL1.hut=NA,
                                  RL2.hut=NA),
                plot.ok=plot.ok,
                spline.ok=spline.ok))
  }
} 

# ============================================================================
# *****************************************************************************
#  F_Round.R

Round <- function(x,unit)
{ # Rounds in standard manner (round() proceeds differently!)
  
  x.sign  <- sign(x)
  x.round <- unit * trunc((abs(x)+0.5*unit)/unit + 1.e-15)
  x.round <- x.round * x.sign
  return(x.round)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#xxx <- seq(-0.50,0.50,by=0.01)
#xxxR <- Round(xxx,0.10)
#cbind(xxx,xxxR)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# ****************************************************************************
#  F_RowToMatrix.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Make input object a matrix / data frame, if possible. Needed after 
#  extracting a single row from a matrix / data frame, which produces a vector.
#  If input object is NULL, return NULL
#  If input object is a vector, return a 1-row matrix with column names taken from 
#  the input
#  If input object is a matrix / data.frame, return it as it is
 
#  #  CHANGE HISTORY
#  09.09.2020 Treatmetn of data frames added 
#  09.03.2020 Start
#
# ****************************************************************************

RowToMatrix <- function(obj)
{
  #  INPUT 
  #  obj   a scalar, a vector, or a matrix

  #  OUTPUT 
  #  obj.out   a 1-row matrix, if possible, otherwise NA
  # ---------------------------------------------------------------------------

  #cat("\n [RowToMatrix] mode(obj)      ", mode(obj) ) 
  #cat("\n [RowToMatrix] class(obj)     ", class(obj) ) 
  #cat("\n [RowToMatrix] typeof(obj)    ", typeof(obj) ) 
  #cat("\n [RowToMatrix] length(obj)    ", length(obj) ) 
  #cat("\n [RowToMatrix] is.matrix(obj) ", is.matrix(obj), "\n" ) 

  if ( (length(obj) > 0) && (is.matrix(obj) & !is.data.frame(obj)) )  
  { 
    obj.out <- obj
  }

  if ( (length(obj) > 0) && (!is.matrix(obj) & is.data.frame(obj)) )  
  { 
    obj.out <- obj
  }

  if ( (length(obj) > 0) && (!is.matrix(obj) & !is.data.frame(obj)) )
  { 
    cnames  <- names(obj)
    obj.out <- matrix(obj, nrow=1)
    colnames(obj.out) <- cnames
  }

  if ( is.null(obj) )
  { 
    obj.out <- obj
  }

  # cat("\n [RowToMatrix] class(obj.out)    ", class(obj.out), "\n") 

  return(obj.out)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test
#A <- matrix(c(11, 12, 13,
#              21, 22, 23), byrow=TRUE, nrow=2)
#dimnames(A) <- list(c("r1", "r2"), c("c1", "c2", "c3"))
## A <- data.frame(A)
#print(A)

#length(A)
#is.matrix(A)
#is.data.frame(A)

#B <- RowToMatrix(A)
#print(B)

#C <- A[1, ]
#print(C)

#D <- RowToMatrix(C)
#print(D)

#E <- NA
#is.null(E)
#print(E)
#F <- RowToMatrix(E)
#print(F)

#G <- NULL
#print(G)
#H <- RowToMatrix(G)
#print(H)
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#  F_ShowTiProps.R
#
ShowTiProps <- function(tabn.tmc.all, p.fit.min, p.rt.min, text)
{
  #  F_ShowTiProps
  #
  #  (c) wwosniok@matrh.uni-bremen.de
  #
  #  Show truncation interval properties
  # 
  #  CHANGE HISTORY
  #  
  #  14.11.2022 Start
  # =========================================================================

  tabn.tmc.all.n <- nrow(tabn.tmc.all)

  # col.iti <- rainbow(tabn.tmc.all.n)
  # col.iti <- heat.colors(tabn.tmc.all.n)   #  red -> yellow
  # col.iti <- topo.colors(tabn.tmc.all.n)   #  blue -> yellow
  # col.iti <- terrain.colors(tabn.tmc.all.n)   #  green -> very light brown
  # col.iti <- cm.colors(tabn.tmc.all.n)   #  cyan -> pink, middle invisible
  
  #  hcl palettes, hcl.pals()
  hcl.palettes <- hcl.pals()
  hcl.palettes.n <- length(hcl.palettes)
#  [1] "Pastel 1"      "Dark 2"        "Dark 3"        "Set 2"        
#  [5] "Set 3"         "Warm"          "Cold"          "Harmonic"     
#  [9] "Dynamic"       "Grays"         "Light Grays"   "Blues 2"      
# [13] "Blues 3"       "Purples 2"     "Purples 3"     "Reds 2"       
# [17] "Reds 3"        "Greens 2"      "Greens 3"      "Oslo"         
# [21] "Purple-Blue"   "Red-Purple"    "Red-Blue"      "Purple-Orange"
# [25] "Purple-Yellow" "Blue-Yellow"   "Green-Yellow"  "Red-Yellow"   
# [29] "Heat"          "Heat 2"        "Terrain"       "Terrain 2"    
# [33] "Viridis"       "Plasma"        "Inferno"       "Rocket"       
# [37] "Mako"          "Dark Mint"     "Mint"          "BluGrn"       
# [41] "Teal"          "TealGrn"       "Emrld"         "BluYl"        
# [45] "ag_GrnYl"      "Peach"         "PinkYl"        "Burg"         
# [49] "BurgYl"        "RedOr"         "OrYel"         "Purp"         
# [53] "PurpOr"        "Sunset"        "Magenta"       "SunsetDark"   
# [57] "ag_Sunset"     "BrwnYl"        "YlOrRd"        "YlOrBr"       
# [61] "OrRd"          "Oranges"       "YlGn"          "YlGnBu"       
# [65] "Reds"          "RdPu"          "PuRd"          "Purples"      
# [69] "PuBuGn"        "PuBu"          "Greens"        "BuGn"         
# [73] "GnBu"          "BuPu"          "Blues"         "Lajolla"      
# [77] "Turku"         "Hawaii"        "Batlow"        "Blue-Red"     
# [81] "Blue-Red 2"    "Blue-Red 3"    "Red-Green"     "Purple-Green" 
# [85] "Purple-Brown"  "Green-Brown"   "Blue-Yellow 2" "Blue-Yellow 3"
# [89] "Green-Orange"  "Cyan-Magenta"  "Tropic"        "Broc"         
# [93] "Cork"          "Vik"           "Berlin"        "Lisbon"       
# [97] "Tofino"        "ArmyRose"      "Earth"         "Fall"         
#[101] "Geyser"        "TealRose"      "Temps"         "PuOr"         
#[105] "RdBu"          "RdGy"          "PiYG"          "PRGn"         
#[109] "BrBG"          "RdYlBu"        "RdYlGn"        "Spectral"     
#[113] "Zissou 1"      "Cividis"       "Roma" 

  ipal <- 113
  col.iti <- hcl.colors(tabn.tmc.all.n, palette = hcl.palettes[ipal])

  col.p.fit.min <- "firebrick"
  col.p.rt.min  <- "steelblue"

  lty.p.fit.min <- 2
  lty.p.rt.min  <- 3

  pch.RL1 <- 1
  pch.RL2 <- 4
  pch.p.fit <- 1
  pch.p.rt  <- 4
  pch.ok     <- rep(4, times=tabn.tmc.all.n)
  pch.ok[ (tabn.tmc.all[ , "p.fit"]>=p.fit.min) &
          (tabn.tmc.all[ , "p.rt"] >=p.rt.min) ] <- 1
  cex.ok  <- rep(1, times=tabn.tmc.all.n)
  cex.ok[pch.ok==1] <- 1.5

  par(mfrow=c(2,2))

  # -------------------------------------------------------------------------- 
  #  Plot 1: TI limits
  iti <- 1
  plot(c(tabn.tmc.all[iti, "x.tr.lo"], 
         tabn.tmc.all[iti, "x.tr.hi"]), 
         c(1,1) * iti,
       type="l", col=col.iti[iti], lwd=2,
       xlim=c(min(tabn.tmc.all[ , "x.tr.lo"]), max(tabn.tmc.all[ , "x.tr.hi"])),
       ylim=c(0, tabn.tmc.all.n+1),
       xlab="Truncation interval (TI)",
       ylab="TI index",
       main=text,
       sub=hcl.palettes[ipal])
  points(c(tabn.tmc.all[iti, "x.tr.lo"], 
         tabn.tmc.all[iti, "x.tr.hi"]), 
         c(1,1) * iti, 
         col=col.iti[iti], pch=pch.ok[iti], cex=cex.ok[iti],lwd=2)

  if (tabn.tmc.all.n > 1)
  {
    for (iti in 2:tabn.tmc.all.n)
    {
      lines(c(tabn.tmc.all[iti, "x.tr.lo"], 
              tabn.tmc.all[iti, "x.tr.hi"]), 
            c(1,1) * iti,
            type="l", col=col.iti[iti], lwd=2)
      points(c(tabn.tmc.all[iti, "x.tr.lo"], 
               tabn.tmc.all[iti, "x.tr.hi"]), 
             c(1,1) * iti,
             type="p", col=col.iti[iti], pch=pch.ok[iti], cex=cex.ok[iti], lwd=2)
    }
  }

  # -------------------------------------------------------------------------- 
  #  Plot 2: Estimated RLs vs x.tr.prop
  plot(tabn.tmc.all[ , "x.tr.prop"], tabn.tmc.all[ , "x.RL2"], 
       type="p", pch=pch.RL2, col=col.iti, lwd=2,
       ylim=c(0, max(tabn.tmc.all[ , "x.RL2"])),
       xlab="Proportion of values in TI",
       ylab="RLx")
  points(tabn.tmc.all[ , "x.tr.prop"], tabn.tmc.all[ , "x.RL1"], 
         type="p", pch=pch.RL1, col=col.iti, lwd=2)
  legend("right", c("RL1", "RL2"), pch=c(pch.RL1, pch.RL2), lwd=2, cex=0.8)
  abline(h=c(10,35), col="magenta")

  # -------------------------------------------------------------------------- 
  #  Plot 3: p.fit, p.rt vs x.tr.prop
  plot(tabn.tmc.all[ , "x.tr.prop"], tabn.tmc.all[ , "p.fit"], 
       type="p", pch=pch.p.fit, col=col.iti, lwd=2,
       ylim=c(0, 1),
       xlab="Proportion of values in TI",
       ylab="p.fit, p.rt")
  points(tabn.tmc.all[ , "x.tr.prop"], tabn.tmc.all[ , "p.rt"], 
         type="p", pch=pch.p.rt, col=col.iti, lwd=2)
  abline(h=c(p.fit.min, p.rt.min), col=c(col.p.fit.min, col.p.rt.min),
             lty=c(lty.p.fit.min, lty.p.rt.min), lwd=c(1, 2))
  legend("topright", c("p.fit.min", "p.rt.min"), 
         col=c(col.p.fit.min, col.p.rt.min),
         pch=c(pch.p.fit, pch.p.rt),  
         lty=c(lty.p.fit.min, lty.p.rt.min), lwd=2, cex=0.8)

  # -------------------------------------------------------------------------- 
  #  Plot 4: p.fit, vs p.rt
  plot(tabn.tmc.all[ , "p.fit"], tabn.tmc.all[ , "p.rt"], type="p", pch=1,
       xlim=c(0, 1), ylim=c(0, 1), col=col.iti, lwd=2,
       xlab="p.fit",
       ylab="p.rt")
  abline(h=c(p.fit.min), col=c(col.p.fit.min),
             lty=c(lty.p.fit.min, lwd=1 ))
  abline(v=c(p.rt.min), col=c(col.p.rt.min),
             lty=c(lty.p.rt.min, lwd=1 ))

}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load(file="tabn.tmc.all.RData")
#ls()

# "tabn.tmc.all"
# View(tabn.tmc.all)

#p.fit.min <- 0.05
#p.rt.min  <- 0.05

#ShowTiProps(tabn.tmc.all, p.fit.min, p.rt.min, "Test data tabn.tmc.all")



# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_SmoothHist1_V7.R
#
#  Smooth a histogram by 5 point moving average
#
#  29.01.2022 Smoothing range extended
#  23.01.2021 Negative predictions replaced by zero
#  23.02.2021 Five point weighted parabola smoothing 
# ============================================================================

SmoothHist1 <- function(x.hist,
                        weights5=c(1, 1, 1, 1, 1),
                        weights7=c(1, 1, 1, 1, 1, 1, 1),
                        figA=NA )
{
  density1   <- x.hist$density
  density1.n <- length(x.hist$density)

  #  No smoothing if less than 5 bins

  if (density1.n < 5)
  { 
    x.hist2 <- x.hist
  } else
  { 
    x.n       <- sum(x.hist$counts) 
    breaks1   <- x.hist$breaks
    breaks1.n <- length(x.hist$breaks)
 
    if (5 <= density1.n  & density1.n < 20 )
    {
      #   Normalize weights
      weights5 <- weights5/sum(weights5)

      #  Leave first and last 2 points unchanged
      #  Smooth interior counts
      density2 <- rep(NA, times=density1.n)
      density2[1:2] <- x.hist$density[1:2]
      density2[(density1.n-1):density1.n] <- tail(x.hist$density, 2)

      for (i in 3:(density1.n-2))
      {
        subset <- (i-2):(i+2)
        x.sub  <- x.hist$mids[subset]
        x.sub2 <- x.sub^2
        y.sub  <- x.hist$density[subset]
        xy.sub.lm <- lm(y.sub ~ x.sub + x.sub2, weight=weights5)
        xy.sub.lm.pred <- predict(xy.sub.lm)
        # density2[i] <- xy.sub.lm.pred[3] 
        density2[i] <- max(0, xy.sub.lm.pred[3])   # 23.01.2021     
      }
    }  else
    {
      #   Normalize weights
      weights7 <- weights7/sum(weights7)

      #  Leave first and last 3 points unchanged
      #  Smooth interior counts
      density2 <- rep(NA, times=density1.n)
      density2[1:3] <- x.hist$density[1:3]
      density2[(density1.n-2):density1.n] <- tail(x.hist$density, 3)

      for (i in 4:(density1.n-3))
      {
        subset <- (i-3):(i+3)
        x.sub  <- x.hist$mids[subset]
        x.sub2 <- x.sub^2
        y.sub  <- x.hist$density[subset]
        xy.sub.lm <- lm(y.sub ~ x.sub + x.sub2, weight=weights7)
        xy.sub.lm.pred <- predict(xy.sub.lm)
        density2[i] <- max(0, xy.sub.lm.pred[4])
      }
    }

    #  Normalize smoothed density
    bin.width    <- breaks1[2:breaks1.n] - breaks1[1:(breaks1.n-1)]
    density2.sum <- sum(density2 * bin.width)
    density2     <- density2 / density2.sum

    density2.sum <- sum(density2 * bin.width)

    #  Calculate smoothed counts
    counts2   <- x.n * density2 * bin.width 

    x.hist2 <- list(breaks=breaks1, 
                    counts=counts2, 
                    density=density2, 
                    mids=(breaks1[2:breaks1.n] + 
                          breaks1[1:(breaks1.n-1)])/2 )
    class(x.hist2) <- "histogram"
  }

  if (!is.na(figA))
  { 
    dev.set(figA)
    plot(x.hist, col="cornsilk")
    plot(x.hist2, border="green3", add=TRUE)
  }

  return(x.hist2)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#text <- "a"
#if (text == "a") { cat( "Text is a") } else 
#if (text == "b") { cat( "Text is b") } else 
#if (text == "c") { cat( "Text is c") } else 
#{ cat( "Text nicht a, b, c") }
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_SmoothHist2.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Smooth a histogram by using the kernel density estimate as
#  smoothing tool

#  CHANGE HISTORY
#  04.01.2021 Concept is problematic, if many bins / small differences
#             between density values. Stop using it for the moment.
#  20.11.2020 Name changed to SmoothHist2 (TMC has already SmoothHist2) 
#  10.10.2020 Calculation of x.n added
#  26.09.2020 Colours as function input
#  25.09.2020 Start
# ==========================================================================

SmoothHist2 <- function(x.kde, x.hist, FigA, histcol2, bordercol2, 
                       histcol3, bordercol3, kdecol, subtitle)
{  
  #
  #  INPUT
  #  x.kde       kernel density estimate of the raw data
  #  x.hist      the histogram to smooth
  #  FigA        window for plotting. Must exist. NA: no plot
  #  subtitle

  # ==========================================================================
  #

  #  Produce a smoothed version of the input histogram
  x.n <- sum(x.hist$counts)

  x.kde.cdf <- cumsum(x.kde$y)
  x.kde.cdf <- x.kde.cdf / tail(x.kde.cdf, 1)

  x.hist.tab   <- hist.table(x.hist$breaks, x.hist$counts)
  x.hist.tab.n <- nrow(x.hist.tab)
  prop.kde     <- rep(NA, times=x.hist.tab.n)
  
  for (i in 1:x.hist.tab.n)
  {
    prop.kde[i] <- max(x.kde.cdf[x.kde$x < x.hist.tab[i, "x.hi"]])
  }
  prop.kde[2:x.hist.tab.n] <- prop.kde[2:x.hist.tab.n] - 
                              prop.kde[1:(x.hist.tab.n-1)] 

  prop.kde[x.hist.tab.n] <- 1 - sum(prop.kde[1:(x.hist.tab.n-1)])

  x.hist.tab <- cbind(x.hist.tab, prop.kde=prop.kde) 
  x.hist.tab <- cbind(x.hist.tab, count.kde=prop.kde*x.n) 
  x.hist.tab <- cbind(x.hist.tab, density.kde=prop.kde/
                             (x.hist.tab[ , "x.hi"] - x.hist.tab[ , "x.lo"]) )

  x.hist.smo <- x.hist
  x.hist.smo$counts <- x.hist.tab[ , "count.kde"]
  x.hist.smo$density <- x.hist.tab[ , "density.kde"]

  if (!is.na(FigA))
  {
    #  Show original (unsmoothed) histogram
    dev.set(FigA)
    plot(x.hist, freq=FALSE, xlim=c(0, 100), col=histcol1, border=bordercol1,
         main="Original and smoothed histogram", 
         sub=subtitle, cex=0.7)
    lines(x.kde,col=kdecol) 

    plot(x.hist.smo, freq=FALSE, add=TRUE, border=bordercol3, 
         col=histcol3)

   legend("topright", c("kde", "Original histogram", "Smoothed histogram"), 
          col=c("black", bordercol2, bordercol3), 
          lwd=3, lty=1, cex=0.7)
  }

  return(x.hist.smo)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test   

#hs <- SmoothHist2(x.kde, x.hist.eqd, 5, histcol3, bordercol3, 
#                  histcol1, bordercol1, "chartreuse", subtitle)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_tmc.gridsearch.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Grid seach to improve the initial value for tmc/nlm estimation
#
#  CHANGE HISTORY
#
#  06.10.2022 Start
# ==============================================================================

tmc.gridsearch <- function(idx, t23.ini, 
                           x.hist, ilo, ihi,
                           x.lt.tr.n, x.ge.tr.n, 
                           x.Q1, x.Q2, RL1.p, RL2.p,
                           t1.ini,
                           df.est,df.con,
                           l.fact, p.fact, r.fact, 
                           rt.fact, w.fact, 
                           opt.crit.only, fastnull)
{                   

  crit.ini <- chi2trunc( x.hist, ilo, ihi,
                         x.lt.tr.n, x.ge.tr.n, 
                         x.Q1, x.Q2, RL1.p, RL2.p,
                         t1.ini,t23.ini[idx, 1], t23.ini[idx,2],
                         df.est,df.con,
                         l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                         rt.fact=rt.fact, w.fact=w.fact, 
                         opt.crit.only=TRUE, fastnull=fastnull)

  return(crit.ini)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_tmc.master_V8.R
#
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Master function for TMC estimation. This function does the TMC estimation 
#  for a single truncation interval (iti). Initial values for tmc are 
#  produced here by a modified qq plot approach.  
#  
#  Action here, for the truncation interval iti:
#
#  - find initial values via a QQ plot approach in a loop over lambda values
#  - use the QQ plot estimate as inital value for TMC, thereby estimating a 
#    possibly new lambda, but keeping the truncation interval fixed  
#  - Do bias reduction by tmu, switched on/ off by do.tmu
#  - Reported results: from QQ plot (tab.qqw) and TMC (tab.tmc). Note that
#    qqw here does not mean the same as in earlier versions. 
#  

#  To do      This version produces 1 tmc result for 1 TI, for which 1
#             vector of initial valules (from Tukey/QQW) is provided. This 
#             makes the iqq loop superfluous as well as athe search for the
#             optimal tmc solution. There is only 1. ==> Simplify!

#  CHANGE HISTORY
#
#  07.10.2022 qqw component moved from seg 100 to here
#             Grid search reorganized, use of apply, optimum determined
#             at loop end, summary informaiton only upon request

# ============================================================================

tmc.master <- function(x, Y.red, subset, xsupp,
                       round.unit, detect.limits.max, 
                       x.Q1, x.Q2, x.RL1.npa, x.RL2.npa, 
                       x.hist,
                       idx.fix, idx.est,
                       lambda.min, lambda.max, lambda.seq, 
                       l.fact, p.fact, r.fact, rt.fact, w.fact,
                       df.est, df.con, 
                       RL1.p, RL2.p,
                       prev.acc.lo, prev.acc.hi                      ,
                       gencol, kdecol,tmccol, figA,
                       tabn.qqw, tabn.tmc,
                       fastnull, fastnull.chi,
                       print.summary=FALSE, print.log.message=FALSE)
{ #
  #  tmc.master_V8 

  #  INPUT

  #  x
  #  Y.red   
  #  subset
  #  xsupp
  #  x.Q1
  #  x.Q2
  #  x.RL1.npa
  #  x.RL2.npa
  #  x.hist
  #  idx.fix
  #  idx.est
  #  lambda.min 
  #  lambda.max
  #  lambda.seq
  #  l.fact
  #  p.fact
  #  r.fact
  #  rt.fact
  #  w.fact
  #  df.est
  #  df.con 
  #  RL1.p
  #  RL2.p                      
  #  gencol
  #  kdecol
  #  tmccol
  #  figA
  #  tabn.qqw
  #  tabn.tmc
  #  fastnull
  #  fastnull.chi
  #  print.summary
  #  print.log.message

  #  OUTPUT
  #  a list containing
  #  tabn.tmc  vector containing the TMC result. For components see
  #           tabn.tmc.names
  #  tabn.qqw  vector containing the modified qq result. For components see
  #           tabn.qqw.names
 
  # ==========================================================================

  # --------------------------------------------------------------------------
  if (print.log.message) { cat("%%%   tmc.master   Start\n") }
  # --------------------------------------------------------------------------

  #  Find and keep optimal rsq from QQ plot estimate
  rsq.opt.qqw <- 0

  #  qq loop over lambda candidates 

  r2.tab <- apply(Y.red, MARGIN=2, FUN=QQEstimateIni,
                  prop.lo=prop.lo, prop.hi=prop.hi, subset=subset, 
                  figA=figA)
  r2.tab <- t(r2.tab)

  #  Add lambda
  r2.tab <- cbind(lambda.seq, r2.tab)
  colnames(r2.tab) <- c("lambda", "r2", "mue", "sigma")

  #  Find optimum
  idx <- which.max(r2.tab[ ,"r2"])

  tmc.master.t1 <- Sys.time()

  # --------------------------------------------------------------------  
  #  Calculate TMC estimate for this truncation interval

  if (print.log.message) { cat("%%%   tmc.master   100\n") }

  #  Put qqw result in tabn.qqw, is itself a result and will be used 
  #  as initial value for tmc
  tabn.qqw["lambda"] <- r2.tab[idx ,"lambda"]
  tabn.qqw["mue"]    <- r2.tab[idx ,"mue"]
  tabn.qqw["sigma"]  <- r2.tab[idx ,"sigma"]
  tabn.qqw["r2"]     <- r2.tab[idx ,"r2"]

  tabn.qqw["x.RL1"]  <- q.PN(RL1.p, tabn.qqw["lambda"], 
                                    tabn.qqw["mue"], 
                                    tabn.qqw["sigma"])
  tabn.qqw["x.RL2"]  <- q.PN(RL2.p, tabn.qqw["lambda"], 
                                    tabn.qqw["mue"], 
                                    tabn.qqw["sigma"] )

  if (print.log.message) { cat("%%%   tmc.master   200\n") }

  tmc.master.t2 <- Sys.time()
  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n tmc.master; t2-t1;", 
          format(difftime(tmc.master.t2, tmc.master.t1)))
  sink()


  # -----------------------------------------------------------------------
  #  From here: essential beginning of earlier tmc.master_V7

  #  Set constraints for TMC estimation. Take initial values from qqw. 
  #  Transform all components, including those that are fixed. 

  if (print.log.message) { cat("%%%   tmc.master   300\n") }
 
  y.ini.first   <- BoxCox(x[1],             tabn.qqw["lambda"])
  y.ini.last    <- BoxCox(x[tabn.tmc["n"]], tabn.qqw["lambda"])

  RB3           <- matrix(NA,nrow=3,ncol=2)
  dimnames(RB3) <- list(c("lambda","mue","sigma"),c("min","max"))

  RB3[1, ]  <- c( lambda.min-fastnull,lambda.max+fastnull)

  if (y.ini.first <  0) { RB3[2, 1] <- 2*y.ini.first } 
  if (y.ini.first == 0) { RB3[2, 1] <- 2*fastnull }   #  25.05.2020

  if (y.ini.first >  0) { RB3[2, 1] <- 0.5*y.ini.first } 

  RB3[2, 2]  <- 2.0*y.ini.last                # (y is sorted ascending)

  RB3[3, ]  <- c( 0.001 * tabn.qqw["sigma"], 
                      2 * tabn.qqw["sigma"] )

  if (print.log.message) { cat("%%%   tmc.master   400\n") }

  # -----------------------------------------------------------------------
  #  Prepare grid search to improve initial values from qqw: oc for qqw

  opt.crit.qqw <-  chi2trunc(x.hist, tabn.qqw["ilo"], tabn.qqw["ihi"],
                          tabn.qqw["x.lt.tr.n"], tabn.qqw["x.ge.tr.n"], 
                          x.Q1, x.Q2, RL1.p, RL2.p,
                          tabn.qqw["lambda"], tabn.qqw["mue"], 
                          tabn.qqw["sigma"],
                          df.est,df.con,
                          l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                          rt.fact=rt.fact, w.fact=w.fact, 
                          opt.crit.only=TRUE, fastnull=fastnull)

  #  Create summary of qqw, grid, tmc, if requested
  if (print.summary)
  { comp <- matrix(NA, nrow=4, ncol=7)
    dimnames(comp) <- list(c("ini", "grid", "tmc", "diff"), 
             c("p.fit", "p.rt", "opt.crit", "mue", "sigma", 
               "iter", "rc")) 

    comp["ini", "opt.crit"]     <- opt.crit.qqw
    comp["ini", "mue"]          <- tabn.qqw["mue"]
    comp["ini", "sigma"]        <- tabn.qqw["sigma"]
  }

  #  Keep present estimate as initial value
  opt.crit.ini  <- opt.crit.qqw
  theta.ini     <- c(tabn.qqw["lambda"], tabn.qqw["mue"], tabn.qqw["sigma"])

  tmc.master.t3 <- Sys.time()
  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n tmc.master; t3-t2;", 
          format(difftime(tmc.master.t3, tmc.master.t2)))
  sink()


  # -----------------------------------------------------------------------
  #  Grid search to improve initial values from qqw
  #  Only if initial opt.crit bad 

  if (print.log.message) { cat("%%%   tmc.master   500\n") }

  if (opt.crit.qqw > 2.0)  
  {   
    if (print.log.message) { cat("%%%   tmc.master   600\n") }

    #  Additional grid search for mue, sigma, keeping lambda fixed
    #  activated 07.10.2022
    fseq <- c(0.950, 0.970, 0.980, 0.990, 0.995,
              1.005, 1.010, 1.020, 1.030, 1.050) 

    t2.ini.seq <- fseq * tabn.qqw["mue"]
    t3.ini.seq <- fseq * tabn.qqw["sigma"]

    t2.ini.seq.n <- length(t2.ini.seq)
    t3.ini.seq.n <- length(t3.ini.seq)

    #  Put parameters in index vector for apply
    t23.ini <- cbind(rep(t2.ini.seq, each=t3.ini.seq.n),
                     rep(t3.ini.seq, times=t2.ini.seq.n))

    idx <- matrix(rep(1:nrow(t23.ini)), ncol=1)

    opt.crit.grid <- apply(idx, 1, tmc.gridsearch, t23.ini=t23.ini,
                           x.hist=x.hist,
                           ilo=tabn.tmc["ilo"], ihi=tabn.tmc["ihi"],
                           x.lt.tr.n=tabn.tmc["x.lt.tr.n"], 
                           x.ge.tr.n=tabn.tmc["x.ge.tr.n"], 
                           x.Q1=x.Q1, x.Q2=x.Q2, RL1.p=RL1.p, RL2.p=RL2.p,
                           t1.ini=tabn.qqw["lambda"],
                           df.est=df.est,df.con=df.con,
                           l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                           rt.fact=rt.fact, w.fact=w.fact, 
                           opt.crit.only=TRUE, fastnull=fastnull)

    if (print.log.message) { cat("%%%   tmc.master   700\n") }

    #  Find optimum position in grid
    k <- which.min(opt.crit.grid)

    ##  which() gives index in column-wise vector representaton  
    #iy <- floor(k/t2.ini.seq.n) + 1    
    #ix <- k - (iy-1) * t2.ini.seq.n

    #  Reduce grid to optimal point
    opt.crit.grid <- opt.crit.grid[k]

    if (opt.crit.grid < opt.crit.qqw)
    { 
      # grid search found better optimum than qqw
      opt.crit.ini  <- opt.crit.grid
      theta.ini     <- c(tabn.qqw["lambda"], t23.ini[k, 1], theta.ini[3])
    }
  }   #    if (opt.crit.qqw > 2.0 ...

  theta.ini <- unlist(theta.ini)

  if (print.log.message) { cat("%%%   tmc.master   800\n") }

  tmc.master.t4 <- Sys.time()
  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n tmc.master; t4-t3;", 
          format(difftime(tmc.master.t4, tmc.master.t3)))
  sink()


  # -----------------------------------------------------------------------
  #  Start TMC estimation

  #  Transform initial parameter to real scale
  ttheta.ini <- o2R(theta.ini, RB3, fastnull=fastnull)

  #  Estimate better theta by chi^2 minimization
  #  Calculate directly, using parts of tmc0, tmc1

  # ttheta contains all 3 components of theta (lambda, mue, sigma),
  # but possibly some of them are fixed by the user. Pass only components
  # to nlm that are to be estimated. 

  ttheta.ini.est <- Extract(ttheta.ini, idx.est)

  if (is.null(idx.fix[1]))
  { # No fixed component
    ttheta.ini.fix <- NULL
  }  else
  { # >= 1 components in idx.fix
    ttheta.ini.fix <- ttheta.ini[idx.fix]
  }

  if (print.log.message) { cat("%%%   tmc.master  1000\n") }

  tmc.nlm <- nlm(chi2.PNV.tr.lms, ttheta.ini.est, iterlim=500, print.level=0,
                 gradtol=1.e-8, steptol=1.e-8, 
                 RB=RB3, 
                 ttheta.ini.fix=ttheta.ini.fix, 
                 idx.fix=idx.fix, idx.est=idx.est,
                 x.hist=x.hist, 
                 ilo=tabn.tmc["ilo"], ihi=tabn.tmc["ihi"],
                 x.lt.tr.n=tabn.tmc["x.lt.tr.n"], 
                 x.ge.tr.n=tabn.tmc["x.ge.tr.n"],
                 l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, rt.fact=rt.fact,
                 w.fact=w.fact,
                 x.Q1=x.Q1, x.Q2=x.Q2, RL1.p=RL1.p, RL2.p=RL2.p,
                 df.est=df.est, df.con=df.con, 
                 opt.crit.only=TRUE, fastnull=fastnull.chi)

  if (print.log.message) { cat("%%%   tmc.master  1100\n") }

  #  Transform estimated and fixed parameters to original scale
  #  First build complete parameter vector
  ttheta.hut <- Compose(ttheta.ini.fix, tmc.nlm$estimate, idx.fix, idx.est)

  theta.tmc.hut <- R2o(ttheta.hut,RB3)

  #  Put nlm result in tabn.tmc
  
  tabn.tmc["lambda"]   <- theta.tmc.hut[1]
  tabn.tmc["mue"]      <- theta.tmc.hut[2]
  tabn.tmc["sigma"]    <- theta.tmc.hut[3]
  tabn.tmc["opt.crit"] <- tmc.nlm$minimum
  tabn.tmc["iter"]     <- tmc.nlm$iterations
  tabn.tmc["rc"]       <- tmc.nlm$code

  #  Meaning of "code" ("rc"):
  #  1: relative gradient is close to zero, current iterate is 
  #     probably solution.
  #  2: successive iterates within tolerance, current iterate is 
  #     probably solution.
  #  3: last global step failed to locate a point lower than estimate. 
  #     Either estimate is an approximate local minimum of the function 
  #     or steptol is too small.
  #  4: iteration limit exceeded.
  #  5: maximum step size stepmax exceeded five consecutive times. 
  #     Either the function is unbounded below, becomes asymptotic 
  #     to a finite value from above in some direction or stepmax is too small.

  tmc.master.t5 <- Sys.time()
  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n tmc.master; t5-t4;", 
          format(difftime(tmc.master.t5, tmc.master.t4)))
  sink()

  # --------------------------------------------------------------------------
  #  RLs and other derived quantities

  if (print.log.message) { cat("%%%   tmc.master  1200\n") }

  tabn.tmc["x.RL1"] <- q.PN(RL1.p, tabn.tmc["lambda"], tabn.tmc["mue"], 
                         tabn.tmc["sigma"])
  tabn.tmc["x.RL2"] <- q.PN(RL2.p, tabn.tmc["lambda"], tabn.tmc["mue"],
                         tabn.tmc["sigma"])

  #  Relative distance between (P025.c, P975.c) and (x.RL1.tmc, x.RL2.tmc)
  #  This is available only for test data!

  if (is.na(x.RL1.npa))
  { # No test data
    tabn.tmc["rel.dist"] <- NA
  } else
  { #  Test data  
    tabn.tmc["rel.dist"] <- reldist(c(tabn.tmc["x.RL1"], tabn.tmc["x.RL2"]), 
                                   c(x.RL1.npa, x.RL2.npa))
  }

  if (print.log.message) { cat("%%%   tmc.master  1300\n") }

  # Calculate more properties of the actual estimate
  TMC.prop <-  chi2trunc(x.hist, 
                         tabn.tmc["ilo"], tabn.tmc["ihi"],
                         tabn.tmc["x.lt.tr.n"], tabn.tmc["x.ge.tr.n"],
                         x.Q1, x.Q2, RL1.p, RL2.p,
                         tabn.tmc["lambda"], tabn.tmc["mue"],tabn.tmc["sigma"],
                         df.est,df.con,
                         l.fact=l.fact, p.fact=p.fact, r.fact=r.fact, 
                         rt.fact=rt.fact, w.fact=w.fact, 
                         opt.crit.only=FALSE, fastnull=fastnull)

  if (print.log.message) { cat("%%%   tmc.master   1400\n") }

  tabn.tmc["chi2.total"]     <- TMC.prop$res["chi2.total"]
  tabn.tmc["chi2.total.df"]  <- TMC.prop$res["chi2.total.df"]
  tabn.tmc["chi2.trun"]      <- TMC.prop$res["chi2.trun"]
  tabn.tmc["chi2.trun.df"]   <- TMC.prop$res["chi2.trun.df"]
  tabn.tmc["chi2.path"]      <- TMC.prop$res["chi2.path.w"]  

  tabn.tmc["p.fit"]          <- TMC.prop$res["chi2.trun.p"]
  tabn.tmc["p.rt"]           <- TMC.prop$res["p.rt"]
  tabn.tmc["prev.l"]         <- TMC.prop$res["prev.l.tmc"]
  tabn.tmc["prev.c"]         <- TMC.prop$res["prev.c.tmc"]
  tabn.tmc["prev.r"]         <- TMC.prop$res["prev.r.tmc"]
  tabn.tmc["xl.n"]           <- TMC.prop$res["xl.n.tmc"]
  tabn.tmc["xc.n"]           <- TMC.prop$res["xc.n.tmc"]
  tabn.tmc["xr.n"]           <- TMC.prop$res["xr.n.tmc"]

  if (print.log.message) { cat("%%%   tmc.master  1500\n") }

  #  Calculate qualitative assessment criteria
  #  Estimated RLs in [x.Q1, x.Q2]?
  #  Consider crit 1 only if no values < DL, because x.Q1 is
  #  very unreliable if values < DL exist
  if (is.na(detect.limits.max))
  { 
    tabn.tmc["crit1"] <- unname(x.Q1-round.unit <= tabn.tmc0["x.RL1"]) # 01.03.2022
  } else
  { tabn.tmc["crit1"] <- TRUE } 

  #  Test 28.02.2022: crit1 always TRUE
  tabn.tmc["crit1"] <- TRUE

  tabn.tmc["crit2"] <- unname(x.Q2+round.unit >= tabn.tmc["x.RL2"])  # 01.03.2022

  #  Estimated prevalences in [prev.acc.lo, prev.acc.hi]?

  tabn.tmc["crit3"] <- unname((prev.acc.lo <= tabn.tmc["prev.l"])) &
                      unname((tabn.tmc["prev.l"] <= prev.acc.hi))  
  tabn.tmc["crit4"] <- unname((prev.acc.lo <= tabn.tmc["prev.c"])) &
                      unname((tabn.tmc["prev.c"] <= prev.acc.hi))  
  tabn.tmc["crit5"] <- unname((prev.acc.lo <= tabn.tmc["prev.r"])) &
                      unname((tabn.tmc["prev.r"] <= prev.acc.hi))

  #  Basic solution score
  tabn.tmc["sol.score"] <- sum(tabn.tmc["crit1"], tabn.tmc["crit2"],
                               tabn.tmc["crit3"], tabn.tmc["crit4"],
                               tabn.tmc["crit5"], na.rm=FALSE)
 
  if (print.log.message) { cat("%%%   tmc.master   1600\n") }

  #  If requested, calculate summary
  if (print.summary)
  { comp["tmc", "p.fit"]        <- tabn.tmc["p.fit"]
    comp["tmc", "p.rt"]         <- tabn.tmc["p.rt"]
    comp["tmc", "opt.crit"]     <- tabn.tmc["opt.crit"]
    comp["tmc", "mue"]          <- tabn.tmc["mue"]
    comp["tmc", "sigma"]        <- tabn.tmc["sigma"]
    comp["tmc", "iter"]         <- tabn.tmc["iter"]
    comp["tmc", "rc"]           <- tabn.tmc["rc"]

    #  Was there an improvement by TMC?
    comp["diff", "p.fit"]    <- comp["tmc", "p.fit"] - comp["grid", "p.fit"]
    comp["diff", "p.rt"]     <- comp["tmc", "p.rt"]  - comp["grid", "p.rt"]
    comp["diff", "opt.crit"] <- 
                             comp["tmc", "opt.crit"] - comp["grid", "opt.crit"]
    comp["diff", "mue"]   <- comp["tmc", "mue"] - comp["grid", "mue"]
    comp["diff", "sigma"] <- comp["tmc", "sigma"] - comp["grid", "sigma"]

    cat("\n [tmc.master] Change due to TMC (tmc - ini) \n")
    print(comp)   
  }

  tmc.master.t6 <- Sys.time()
  sink(file="../temp/UsedTime.txt", append=TRUE)
  cat("\n tmc.master; t6-t5;", 
          format(difftime(tmc.master.t6, tmc.master.t5)))
  sink()

  #  qqw and tmc result for this TI finalized

  if (print.log.message) { cat("%%%   tmc.master   End\n") }

  return(list(tabn.qqw=tabn.qqw, tabn.tmc=tabn.tmc))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_TransformKDE.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Transform the kde of x to the kde of y <- f(x)
 
#  #  CHANGE HISTORY
#
# 18.11.2020 Start
# ============================================================================

TransformKDE <- function(x.kde, y)
{
  #  INPUT 
  #  x.kde          kernel density estimate of x, has components $x, $y
  #  y              transformed x.kde$x

  #  OUTPUT 
  #  y.kde          kernel density estimate of y, has components $x, $y
  # =========================================================================

  kernel.n <- length(x.kde$x)

  #  Calculate area of each trapezoid in x.kde
  x.area <- diff(x.kde$x) * (x.kde$y[1:(kernel.n-1)] + x.kde$y[2:kernel.n])/2
  # cat("\n sum(x.area): ", sum(x.area), "\n")

  y.diff <- diff(y)
  y.kde.y <- rep(NA, times=kernel.n)
  y.kde.y[1] <- x.kde$y[1]     #  should be zero in both cases

  for (i in 2:kernel.n)
  { 
    y.kde.y[i] <- 2*x.area[i-1]/y.diff[i-1] -  y.kde.y[i-1]
  }
  y.area <- y.diff * (y.kde.y[1:(kernel.n-1)] + y.kde.y[2:kernel.n])/2
  # cat("\n sum(y.area): ", sum(y.area), "\n")
  
  return(y.kde=list(x=y, y=y.kde.y))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#  Test
#x.kde <- list(x=seq(6, 14, by=0.1), y=dnorm(seq(6, 14, by=0.1), mean=10, sd=1))
#kern.n <- length(x.kde$x)

#y <- log(x.kde$x)
##y <- x.kde$x

#y.kde <- TransformKDE(x.kde, y) 

#cbind(x.kde$x, y.kde$x)
#cbind(x.kde$y, y.kde$y)

#z <- exp(y.kde$x)
#z.kde <- TransformKDE(y.kde, z)

#cbind(x.kde$y, z.kde$y)
#  Hin- und Rücktransformation ist ok
# ............................................................................

#x <- rnorm(1000, mean=10, sd=2)
#x.kde <- density(x)

#y <- log(x)
#y.kde1 <- density(y)

#y.kde2 <- TransformKDE(x.kde, log(x.kde$x))

#plot(y.kde1, type="l", col="red")
#lines(y.kde2$x, y.kde2$y, type="l", col="blue")

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# ****************************************************************************
#  F_TukeyWW.R 
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  'Tukey' method of removing 'outliers'
#  Find a subset of the input data which has data only in the range
#  Q1 - 1.5*IQR, Q3 + 1.5*IQR, calculate mean and SD, adjust SD estimate
#  because of using a truncated sample

#  Source for 'Tukey' method: Horn, Pesce: Reference intervals: an update.
#  Clinica Chimica Acta 334(2003) 5-23, p.13
#  Expansion by WW: iterated solution
#  Source for truncation adjustment: WW

#  'Fences' in Tukey terminology correspond to truncation interval limits  

#  NEEDS
#  Quantile()   function to calculate quantiles
 
#  #  CHANGE HISTORY
#  22.05.2021 Iteration stop criterion x.tr.prop.min added
#  16.05.2021 Iteration made optional, default 1 step
#  25.04.2021 Iteration added
#  22.04.2021 Start  
# ============================================================================

TukeyWW <- function(x, xlabel, IQR.fact=1.5, iter.max=1, x.tr.prop.min=0.50, 
                    RL1.p=0.025, RL2.p=0.975, 
                    print.details=FALSE, figA=NA)
{
  #  INPUT 
  #  x              data vector. NAs are not checked, not allowed
  #  IQR.fact       IQR range factor 
  #  iter.max       controls iteration: 1 : basic (original) Tukey approach
  #                                    >1 : extended version
  #  x.tr.prop.min  iteration emergency break, prevents generation of too small
  #                 Tukey data sets  
  #  RL1.p, RL2.p   probabiilities for RLs       
  #  print.details  print details of the calculation and iteration process

  #  OUTPUT 
  #  mue.hat        estimate for mean
  #  sigma.hat      estimate for standard deviation 
  #  --> see more terms in the return statement
  # ---------------------------------------------------------------------------
  
  #  For a standard normal distribution, the fences below correspond to a 
  #  probability of ...
  Q1.z    <- qnorm(0.25)
  Q3.z    <- -Q1.z
  IQR.z   <- Q3.z - Q1.z
  z.lo    <- Q1.z - IQR.fact*IQR.z     #  lower fence, theoretical
  z.hi    <- Q3.z + IQR.fact*IQR.z     #  upper fence, theoretical
  p.trunc <- pnorm(z.hi) - pnorm(z.lo) #  probability between fences
  range.z <- z.hi - z.lo

  #  Calculate fences for input data
  x.n <- length(x)
  ok  <- rep(TRUE, times=x.n)

  #  Make the loop start at all
  out.lo    <- 1
  out.hi    <- 1
  x.tr.prop <- 1
  iter      <- 0

  #  Iterate until no more cases lie outside the fences or the data 
  #  set is reduced to less than x.tr.prop.min

  #  Problem: with rounded data, the decrease in IQR can be smaller than
  #  the rounding unit, causing that no cases are removed and IQR does no more
  #  change. 
  while ((out.lo > 0 | out.hi > 0) & (iter < iter.max) & 
          (x.tr.prop > x.tr.prop.min))
  {
    iter <- iter + 1
    Q1 <- Quantile(x[ok], probs=0.25)
    Q3 <- Quantile(x[ok], probs=0.75)
    IQR <- Q3 - Q1

    #  The fences x.lo and x.hi are part of the accepted interval
    #  Note that these limits may coincide with  existing values if the data
    #  is rounded.
    x.lo <- Q1 - IQR.fact*IQR
    x.hi <- Q3 + IQR.fact*IQR

    #  Shift the theoretical limits to really existing values
    # x.lo <- min(x[x >= x.lo])
    # x.hi <- max(x[x <= x.hi])
    #  No good idea. Changes the proportion of data between fences (p.trunc)
    #  to an unknown level

    out.lo <- sum(x[ok] < x.lo)
    out.hi <- sum(x[ok] > x.hi)
    ok <- (x.lo <= x) & (x <= x.hi)
    x.tr.prop <- sum(ok) / x.n

    if (print.details)
    {
      cat("\n Iteration ", iter, " x.lo", x.lo, " x.hi", x.hi,
                                 " out.lo", out.lo, " out.hi", out.hi, 
                                 " x.tr.prop", x.tr.prop)
    }
  }
  cat("\n")

  #  No more values outside the fences
  mean.hat  <- mean(x[ok])
  sigma.raw <- sd(x[ok])

  # Adjusted estimate of sigma: leave the mean as it is and look for the sd
  # that assigns a probability p.trunc to (x.lo, x.hi)
  
  range.x   <- x.hi - x.lo
  sigma.adj <- range.x/range.z  

  # Calculate xc.n. Needs total numbers of removed.
  out.lo <- sum(x < x.lo)
  out.hi <- sum(x > x.hi)
  xc.n   <- (x.n - out.lo - out.hi)/p.trunc

  #  Distribution based  estimate for RLx
  x.RL1 <- qnorm(RL1.p, mean=mean.hat, sd=sigma.adj)
  x.RL2 <- qnorm(RL2.p, mean=mean.hat, sd=sigma.adj)

  if (print.details)
  {
    cat("\n Total # of values        :", x.n,
        "\n Range of input data      :", min(x), max(x),  
        "\n Range of truncated data  :", min(x[ok]), max(x[ok]),  
        "\n Tukey fences             :", x.lo, x.hi,  
        "\n # of values removed below:", out.lo, 
        "\n # of values removed above:", out.hi,
        "\n Tukey estimate for mean  :", mean.hat,
        "\n Raw Tukey estimate for sd:", sigma.raw,
        "\n Q1.z, Q3.z               :", Q1.z, Q3.z,
        "\n p.trunc                  :", p.trunc,
        "\n range.z                  :", range.z,
        "\n range.x                  :", range.x,
        "\n Adjusted estimate for sd :", sigma.adj,
        "\n Estimated xc.n           :", xc.n,
        "\n Allowed # of iterations  :", iter.max,
        "\n")
  }

  return(c(mean=mean.hat, 
           sd.raw=sigma.raw, sd.adj=unname(sigma.adj), 
           x.tr.lo=unname(x.lo), x.tr.hi=unname(x.hi), 
           z.tr.lo=unname(z.lo), z.tr.hi=unname(z.hi), 
           out.lo=out.lo, out.hi=out.hi,            
           x.RL1=x.RL1, x.RL2=x.RL2, 
           xc.n=xc.n)) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Test: see ../dev/Test_Tukey.R 


# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ****************************************************************************
#  F_WaWoTest.R
#
#  Wald-Wolfowitz test (runs test)
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  CHANGE HISTORY
#  30.08.2022  Start
#  ============================================================================

WaWo <-  function(sequence, exact=FALSE, alternative="two.sided", 
                  fastnull=1.e-10)
{
  # sequence    input sequence, only allowed members are 0 and 1 
  # exact       do exact test? (if FALSE: do asymptotic test)
  # alternative type of alternative (two.sided, less, greater)
  #             "less": less than mean
  # =========================================================================

  sequence.n <- length(sequence)

  #  Initialize output vector
  WaWo <- list(n=sequence.n, exact=exact, alternative=alternative, 
            n0=NA, n1=NA, runs=NA, mean=NA, var=NA, z=NA, p.value=NA)

  ok        <- TRUE

  # More than 1 value?

  if (sequence.n < 2)
  {
    cat("\n [WaWo] < 2 observations \n")
    print(sequence)
    ok <- FALSE
  }

  if (ok)
  {
    s.table       <- table(sequence)
    s.table.n     <- length(s.table)
    s.table.names <- names(s.table)
  }

  #  Less than 2 different values in the sequence?
  if (s.table.n < 2)
  {
    #  No. Put available information into WaWo
    WaWo["runs"]    <- 1
    WaWo["z"]       <- qnorm(fastnull)
    WaWo["p.value"] <- fastnull
    ok <- FALSE
  }

  #  More than 2 different values in the sequence?
  if (s.table.n > 2)
  {
    cat("\n [WaWo] Unpermissible value in input \n")
    print(s.table)
    ok <- FALSE

    #  No further entries in WaWo
  }

  if (ok)
  {
    idx0 <- which(s.table.names == "0")
    WaWo["n0"] <- s.table[idx0]

    idx1 <- which(s.table.names == "1")
    WaWo["n1"] <- s.table[idx1]

    # Find number of runs
    r <- 1
    for (i in 2:sequence.n)
    { 
      r <- r + (sequence[i-1] != sequence[i])
    }
    WaWo["runs"] <- r

    if (!exact)
    {
      r.mea <- 2 * as.numeric(WaWo["n0"]) * as.numeric(WaWo["n1"]) / 
               sequence.n + 1
      r.var <- (r.mea-1) * (r.mea-2) / (sequence.n-1)
      z     <- (r - r.mea) / sqrt(r.var)

      if (alternative == "two.sided") {p.value <- 2 * (1-pnorm(abs(z), 0, 1)) }
      if (alternative == "less")      {p.value <-        pnorm(z, 0, 1) }
      if (alternative == "greater")   {p.value <-      1-pnorm(z, 0, 1) }

      WaWo["mean"]     <- r.mea
      WaWo["var"]      <- r.var
      WaWo["z"]        <- z
      WaWo["p.value"]  <- p.value
    }   
  }

  return(WaWo)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Test

if (FALSE)
{
  sign0 <- c(0,0,0)  # only 1 run, too few
  rt0t <- WaWo(sign0, exact=FALSE, alternative="two.sided")
  rt0l <- WaWo(sign0, exact=FALSE, alternative="less")
  rt0g <- WaWo(sign0, exact=FALSE, alternative="greater")
  print(rt0t)
  print(rt0l)
  print(rt0g)

  sign1 <- c(1,1,1)   # only 1 run, too few
  rt1t <- WaWo(sign1, exact=FALSE, alternative="two.sided")
  rt1l <- WaWo(sign1, exact=FALSE, alternative="less")
  rt1g <- WaWo(sign1, exact=FALSE, alternative="greater")
  print(rt1t)
  print(rt1l)
  print(rt1g)

  sign1 <- c(0,0,0,0,1,1,1,1)  #  only 2 runs, too few
  rt1t <- WaWo(sign1, exact=FALSE, alternative="two.sided")
  rt1l <- WaWo(sign1, exact=FALSE, alternative="less")
  rt1g <- WaWo(sign1, exact=FALSE, alternative="greater")
  print(rt1t)
  print(rt1l)
  print(rt1g)


  sign2 <- c(0,1,0,1,0,1,0,1)  #  too many runs
  rt2t <- WaWo(sign2, exact=FALSE, alternative="two.sided")
  rt2l <- WaWo(sign2, exact=FALSE, alternative="less")
  rt2g <- WaWo(sign2, exact=FALSE, alternative="greater")
  print(rt2t)
  print(rt2l)
  print(rt2g)
  rt2g$z
  rt2g["z"]
  rt2g[["z"]]
}

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



xxx <- list(a="a", b = 3)
xxx["a"] <- "c"
# *****************************************************************************
# F_x.mean.PN.R

# ----------------------------------------------------------------------       

x.mean.PN <- function(lambda,mue,sigma,fastnull=1.e-10)
{
  #  Mean of the power normal distribution, numerically
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarres, SPL, 2006, p. 766
  #  Dichte ist für x>0 positiv, Null für x = 0, sonst Null

  #  #  CHANGE HISTORY
  #  03.03.2021 q.X replaced by q.PN
  #             q.PN replaced by 0, Inf
  #  31.07.     Integration limits fixed from quantiles
  #  27.07.2019 Installation
  # ===========================================================================

  #  Calculate integral numerically

  #  Using 0, Inf as integration limits fails sometimes
  #  Use quantiles instead
  #x.lo <- q.PN(0.0001,lambda,mue,sigma)
  #x.hi <- q.PN(1-0.0001,lambda,mue,sigma)

  EX <- integrate(x.pdf.PN, 0, Inf, lambda=lambda, mue=mue, sigma=sigma)

  return(EX[[1]]) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_x.mode.PN.R

# --------------------------------------------------------------------------
x.mode.PN <- function(lambda,mue,sigma,fastnull=1.e-10)
{
  #  Mode of x with x ~ PND
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarres, SPL, 2006, p. 766
  #  Numerical calculation (no case distinction necessary

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================

  #  Determine integration interval

  x.min <- q.X(  0.001, lambda, mue, sigma)
  x.max <- q.X(1-0.001, lambda, mue, sigma)

  x <- seq(x.min, x.max, length.out=1000)
  f <- pdf.PN(x, lambda, mue, sigma)
  x.mode <- x[which.max(f)]
  if (length(x.mode) > 1)
  { 
    cat("\n x has multiple modes - first is used \n")
    print(x.mode)
    x.mode <- x.mode[1] 
  } 
  return(x.mode) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  F_x.pdf.PN.R

# -----------------------------------------------------------------------------

x.pdf.PN <- function(x, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Auxiliary function for mean.PN
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================

  return(x * pdf.PN(x,lambda,mue,sigma))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv# *****************************************************************************
#  F_x.Var.PN.R

# ----------------------------------------------------------------------       

x.Var.PN <- function(EX, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Variance of the power normal distribution, numerically
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766
  #  Dichte ist für x>0 positiv, Null für x = 0, sonst Null

  #  #  CHANGE HISTORY
  #  03.03.2021 q.X replaced by q.PN
  #             Quantiles replaced by 0, Inf
  #             EX required as argument 
  #  31.07.     Integration limits fixed from quantiles
  #  27.07.2019 Installation
  # ===========================================================================

  #  Calculate integral numerically

  #  Using 0, Inf as integration limits fails sometimes
  #  Use quantiles instead
  # x.lo <- q.PN(0.0001,lambda,mue,sigma)
  # x.hi <- q.PN(1-0.0001,lambda,mue,sigma)

  VarX <- integrate(xx.pdf.PN, 0, Inf, EX=EX, lambda=lambda, mue=mue, 
                    sigma=sigma)
  return(VarX[[1]]) 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# *****************************************************************************
#  xx.pdf.PN.R

xx.pdf.PN <- function(x, EX, lambda,mue,sigma,fastnull=1.e-10)
{
  #  Auxiliary function for Var.PN
  #  Parameters of the PND: lambda, mue, sigma
  #  siehe Freeman & Modarra, SPL, 2006, p. 766

  #  #  CHANGE HISTORY
  #  27.07.2019 Installation
  # ===========================================================================
  
  return( (x-EX)^2 * pdf.PN(x,lambda,mue,sigma))
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
