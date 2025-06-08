#  TMC_seg090_Plot_Age_Val.R

#  (c) wwosniok@math.uni-bremen.de

#  Contour plot of the data distribution, by M / F
#  Countours obtained from density estimate per age

#  TODO       
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# =====================================================================
if (print.log.message) { cat("%%%   TMC_seg090_Plot_Age_Val.R  Start\n") }
# ===========================================================================

#  Execution only if age is given in the data

if (!is.na(spalte.a))
{ 
# ------------------------------------------------------------------
  #  Graph windows
  source("TMC_seg089_OpenWindows.R")

  # ------------------------------------------------------------------

  # Determine x range for common support of pdf
  x.min <- min(dataset[ ,spalte.w])
  x.max <- max(dataset[ ,spalte.w])

  for (isex in 1:sex.val.n)
  {
    # Wenn Trennung nach Sex:
    if (sex.val[isex] != "All")
    {
      # Quantile der Verteilung aller Werte für sex = sex.val[isex]
      x.sub   <- dataset[dataset[ ,spalte.s] == sex.val[isex], ]
      x.sub.n <- nrow(x.sub)
    } else
    {
      # Quantile der Verteilung aller Werte
      x.sub   <- dataset
      x.sub.n <- nrow(x.sub)
    }

    # KDE per age
    # Density estimation only if more than 25 cases

    age.list   <- sort(unique(x.sub[ ,spalte.a]))
    age.tab    <- table(x.sub[ ,spalte.a])
    age.tab    <- age.tab[age.tab > 25]

    #  The previous command may have left no data!
    
    age.list   <- as.numeric(names(age.tab))
    age.list.n <- length(age.list)

    if (age.list.n > 0)
    {
      # # of support points
      kde.n      <- 256

      x.kde.mat <- matrix(NA,nrow=kde.n,ncol=age.list.n)

      if (isex == 1) 
      { 
        if (plot.fig090.010)
        {
          dev.set(fig090.010)
        }       
      }

      if (isex == 2) 
      { 
        if (plot.fig090.020)
        { 
          dev.set(fig090.020)
        }       
      }
 
      for (iage in 1:age.list.n)
      { 
        x.sub.age <- x.sub[x.sub[ ,spalte.a]==age.list[iage], spalte.w]
        x.kde <- density(x.sub.age, from=x.min, to=x.max, n=kde.n)
 
        #  Replace very small values > 0 by 0 (numerical problems)
        x.kde$y[x.kde$y < fastnull] <- 0 
        x.kde.mat[ ,iage] <- x.kde$y
      }
      #  Transform to a vector (columnwise)  
      x.kde.vec <- as.vector(x.kde.mat)

      #  Quantiles of observed values
      x.kde.quant   <- quantile(x.kde.vec,probs=seq(0.00,1.00,by=0.05)) 

      #  Collapse identical quantiles
      x.kde.quant   <- sort(unique(x.kde.quant)) 
      x.kde.quant.n <- length(x.kde.quant)

      #  Set up plot coordinates
      xplot  <- rep(age.list,each=kde.n)
      yplot  <- rep(x.kde$x,times=age.list.n)

    # xx.mat <- matrix(c(11,12,21,22),byrow=TRUE,ncol=2)
    # xx.vec <- as.vector(xx.mat) # wird spaltenweise zum Vektor
    # xx.mat2 <- matrix(xx.vec,ncol=2)   # Originalmatrix

      #  Assign quantile group to each component in x.kde.mat
      x.qgrp.mat <- x.kde.mat

      for (iqgrp in 1:(x.kde.quant.n-1))
      {  
        x.qgrp.mat[(x.kde.quant[iqgrp] <= x.kde.mat) & 
                   (x.kde.mat <  x.kde.quant[iqgrp+1])] <- iqgrp
      }
      x.qgrp.mat[(x.kde.quant[x.kde.quant.n] == x.kde.mat)] <- 
                 x.kde.quant.n-1

      #col.list <- heat.colors(x.kde.quant.n)
      col.list <- topo.colors(x.kde.quant.n)
      #col.list <- cm.colors(x.kde.quant.n)
      #col.list <- terrain.colors(x.kde.quant.n)

      xplotmin <- age.list[1]
      xplotmax <- tail(age.list, 1) + 10
      yplotmin <- min(x.kde$x)
      yplotmax <- max(x.kde$x)

      par(las=par.las)
      par(tcl=par.tcl)

      plot(xplot, yplot, type="p", pch=15, 
         col=col.list[as.vector(x.qgrp.mat)], cex=0.8,
         xlim=c(xplotmin,xplotmax), ylim=c(yplotmin,yplotmax),
         xlab="Age", ylab="x", main=paste("Sex =",sex.val[isex]),
         sub=paste(date()," / ",outname," / n =",x.sub.n),cex.sub=0.7)

      abline(v=age.limits,col="red")

      # Legende
      cont.levels <- paste(1:x.kde.quant.n,"\n", sep="") 
      points(rep(xplotmax,times=x.kde.quant.n),
           seq(x.kde$x[1],tail(x.kde$x,1),length.out=x.kde.quant.n),
           col=col.list[1:x.kde.quant.n], pch=15, cex=1.5)
      text(rep(xplotmax*1.00,times=x.kde.quant.n), 
         seq(x.kde$x[1],tail(x.kde$x,1),length.out=x.kde.quant.n)-
             0.025*(yplotmax-yplotmin),
         labels=cont.levels,
         pos=2, cex=0.90,col="firebrick")

      if (print.age.contours)  
      { 
        cat("\nContour levels for sex = ",sex.val[isex],"\n")
        print(data.frame(Level.num=1:x.kde.quant.n, Value=x.kde.quant))
      }

      if (isex == 1) { savePlot(file=fig090.010.file,type=figtype) }
      if (isex == 2) { savePlot(file=fig090.020.file,type=figtype) }    
    }  else
    { 
      cat("\n+++  Insufficient data for a Age * Value contour plot for sex =",
           sex.val[isex], "\n" ) 
    }  
  }
}  

# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg090_Plot_Age_Val.R  Ende\n") }
# ===========================================================================

