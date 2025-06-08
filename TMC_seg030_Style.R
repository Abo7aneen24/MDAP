#  TMC_seg030_Style.R
#  (c) wwosniok@math.uni-bremen.de
#
#  Truncated minimum chi-square estimation
#
#  Define style (colors, lines, ...)

#  CHANGE HISTORY

##  08.11.2022 Last change (clean up)
# ============================================================================

if (print.log.message) { cat("%%%   Start TMC_seg30_Style\n") }

#  Define transparent colors
rgb.black <- col2rgb("black")
black.64  <- rgb(rgb.black["red",1], rgb.black["green",1], rgb.black["blue",1], 
                 64, names = NULL, maxColorValue = 255)

black.32  <- rgb(rgb.black["red",1], rgb.black["green",1], rgb.black["blue",1], 
                 32, names = NULL, maxColorValue = 255)

black.16  <- rgb(rgb.black["red",1], rgb.black["green",1], rgb.black["blue",1], 
                 16, names = NULL, maxColorValue = 255)


rgb.blue <- col2rgb("blue")
blue.64  <- rgb(rgb.blue["red",1], rgb.blue["green",1], rgb.blue["blue",1], 
                 64, names = NULL, maxColorValue = 255)

blue.32  <- rgb(rgb.blue["red",1], rgb.blue["green",1], rgb.blue["blue",1], 
                 32, names = NULL, maxColorValue = 255)

blue.16  <- rgb(rgb.blue["red",1], rgb.blue["green",1], rgb.blue["blue",1], 
                 16, names = NULL, maxColorValue = 255)

rgb.green4 <- col2rgb("green4")
green.64   <- rgb(rgb.green4["red",1], rgb.green4["green",1], rgb.green4["blue",1], 
                 64, names = NULL, maxColorValue = 255)
green.32   <- rgb(rgb.green4["red",1], rgb.green4["green",1], rgb.green4["blue",1], 
                 32, names = NULL, maxColorValue = 255)
green.16   <- rgb(rgb.green4["red",1], rgb.green4["green",1], rgb.green4["blue",1], 
                 16, names = NULL, maxColorValue = 255)

rgb.red <- col2rgb("red")
red.64  <- rgb(rgb.red["red",1], rgb.red["green",1], rgb.red["blue",1], 
                 64, names = NULL, maxColorValue = 255)

red.32  <- rgb(rgb.red["red",1], rgb.red["green",1], rgb.red["blue",1], 
                 32, names = NULL, maxColorValue = 255)

red.16  <- rgb(rgb.red["red",1], rgb.red["green",1], rgb.red["blue",1], 
                 16, names = NULL, maxColorValue = 255)

rgb.yellow <- col2rgb("yellow")
yellow.64  <- rgb(rgb.yellow["red",1], rgb.yellow["green",1], rgb.yellow["blue",1], 
                 64, names = NULL, maxColorValue = 255)

yellow.32  <- rgb(rgb.yellow["red",1], rgb.yellow["green",1], rgb.yellow["blue",1], 
                 32, names = NULL, maxColorValue = 255)

yellow.16  <- rgb(rgb.yellow["red",1], rgb.yellow["green",1], rgb.yellow["blue",1], 
                 16, names = NULL, maxColorValue = 255)

rgb.violet <- col2rgb("violet")
violet.64  <- rgb(rgb.violet["red",1], rgb.violet["green",1], rgb.violet["blue",1], 
                 64, names = NULL, maxColorValue = 255)

violet.32  <- rgb(rgb.violet["red",1], rgb.violet["green",1], rgb.violet["blue",1], 
                 32, names = NULL, maxColorValue = 255)

violet.16  <- rgb(rgb.violet["red",1], rgb.violet["green",1], rgb.violet["blue",1], 
                 16, names = NULL, maxColorValue = 255)

#  Colors for data points 
xlcol <- "blue"               # left distribution
xccol <- "green3"             # central distribution 
xrcol <- "red"                # right distribution 

#  Colors for histograms
#  Test colors
bordercol  <- "chocolate3"   # "gray"         # "steelblue"
histcol    <- "linen"        # "cornsilk"     # "goldenrod"
bordercol1 <- "black"         #  FD 
histcol1   <- black.16      
bordercol2 <- "blue"          #  x.hist, collapsed
histcol2   <- blue.16  
bordercol3 <- "green4"        # x.hist.s1
histcol3   <- green.16  
bordercol4 <- "firebrick"     # x.hist.s2 
histcol4   <- red.16  

#  Regular colors
bordercol  <- "gray40"   # "chocolate3"   # "gray"         # "steelblue"
histcol    <- "gray97"  # linen        # "cornsilk"     # "goldenrod"
bordercol1 <- "black"         #  FD 
histcol1   <- black.16      
bordercol2 <- "blue"          #  x.hist, collapsed
# histcol2   <- blue.16  
histcol2   <- yellow.32   
bordercol3 <- "green4"        # x.hist.s1
histcol3   <- green.16  
bordercol4 <- "firebrick"     # x.hist.s2 
histcol4   <- red.16  
difbordercol <- "violet"        # residuals histogram 
difhistcol   <- violet.64 
polycol      <- "black" 


#  See ColorPalettes.R for a selection of colors
col.list <- brewer.pal(9, "Set1")  # 9 different colors from palette Set1

#  Colours for densities and other curves
bhacol      <- col.list[8]
boscol      <- "magenta"       # bootstrap
difcol      <- "violet"
gencol      <- "magenta"
inicol      <- "chartreuse"
kdecol      <- "black"
nlmcol      <- "firebrick"
npacol      <- col.list[9]
patcol      <- "firebrick"
qqcol       <- col.list[7]
qqrawcol    <- "gray"
qqtcol      <- "pink"
qqwcol      <- col.list[5]
tmccol      <- col.list[2]
tmlcol      <- col.list[3]
tmucol      <- col.list[1]
tukcol      <- col.list[4]

#  Colours and symbols for RLs and densities, 
#  depending on the estimation method

# meth.list  : "bha" "npa" "qqw" "tmc"  "tml" "tmu" "tuk"
  
meth.dsg <- matrix(c("bha",bhacol,2,
                     "npa",npacol,0,
                     "qqw",qqwcol,5,
                     "tmc",tmccol,4,
                     "tml",tmlcol,1,
                     "tmu",tmucol,3,
                     "tuk",tukcol,3
                     ),byrow=TRUE,ncol=3) 
meth.dsg <- data.frame(meth.dsg,stringsAsFactors=FALSE)

#  Line types
denlty      <- 1
RLxlty      <- 2
modlty      <- 3

#  Colours and symbols for plotting RLs vs age
sexcol     <- c("red","blue")   # F, M
sexcolfill <- c(red.32, blue.32)
sexpch     <- c(1,4)

sexcolFM <- "green4"
sexpchFM <- 0

if (print.log.message) { cat("%%%   End TMC_seg30_Style\n") }
