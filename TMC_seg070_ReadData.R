#  TMC_seg070_ReadData.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation
#
#  Read the input data and produce the "first value" data set,
#  if requested

#  TODO      
#  -  

#  CHANGE HISTORY
#  08.11.2022 Last change (clean up)
# ===========================================================================

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg70_ReadData  Start\n") }
# ===========================================================================

source("TMC_seg080_NamesPerFile.R")  #  

# ==========================================================================
#  Datensatz lesen
if (file.exists(infile)) 
{ 
  cat("\n***  Reading input data ...")

  #  Eingabe-Daten existieren, lesen
  if (nchar(dec)==0)  { dec <- "."  }
  if (is.null(dec))   { dec <- "."  }
 
  dataset <- read.table(infile, sep=";", dec=dec, header=TRUE,
                      comment.char="", as.is=TRUE, stringsAsFactors=FALSE)

  cat("   input data read \n")
  cat("\n [seg070] dataset has ", nrow(dataset), " rows and  ", 
           ncol(dataset), " columns\n")

  # dimnames(dataset)    # is a list, first component=vector of row names,
  #                      # second component = vector of x values if x has 
  #                      # only 1 column

  #  Workaround, if only 1 column: then dataset is changed to a vector later,
  #  produces trouble with 2dimensional index [ ,1]. Therefore add row 
  #  number as vector thus producing a matrix data frame. Add at end to 
  #  maintain column numbers of existing variables.

  if (ncol(dataset) == 1)
  { 
    dataset <- cbind(dataset, recno=rownames(dataset))
  }

  # -------------------------------------------------------------------------

  x.n0 <- nrow(dataset)

  cat("\n", rep("-", times=76), "\n", sep="")

  #  erste 10 Zeilen ausgeben
  cat("\nOriginal (unscaled) lines 1 - 10 from ",infile,"\n")
  print(dataset[1:10, ])

  #  letzte 10 Zeilen ausgeben
  cat("\nOriginal (unscaled) lines",x.n0-9,"-",x.n0," from ",infile,"\n")
  print(dataset[(x.n0-9):x.n0, ])
  cat("\n", rep("-", times=76), "\n", sep="")

} else
{ 
  #  Eingabe-Datei existiert nicht, Abbruch
  cat("\n", rep("=", times=76), "\n", sep="")
  cat("\n\n ++++++  Input data file",
        "\n",infile,
        "\n ++++++  does not exist - programme stops \n\n")
  cat("\n", rep("=", times=76), "\n", sep="")
  stop("++++++  Programme stop because of missing input file  ++++")
}

# =================================================================
#  Initial filter: use only values below / above the given split date
#  Split date belongs to first split group.

cat("\n***  Filter input data ...\n")

x.n1 <- nrow(dataset)

#  resolve date-time information  ........................................

x.n2 <- x.n1
x.n3 <- x.n1
x.n4 <- x.n1

#  Information is given in column spalte.dz (SpalteDatumZeit in 
#  Dateiinformation.csv), if at all
if ((length(spalte.dz) > 0) &&  !is.na(spalte.dz))
{ 
  #  Standard or alternative date time format?
  #  Standard      dd.mm.yyyy HH.MM
  #                13.04.2022 16:35
  #                1234567890123456 
  #  Alternative 1 yyyy.mm.dd HH.MM
  #                2022.04.13 16:35
  #                1234567890123456 
  #  as.numeric("13.0") > 1000  FALSE
  #  as.numeric("2022") > 1000  TRUE
  #  Other alternatives not covered, e.g American format mm.dd.yyyy

  if (as.numeric(substr(dataset[1 ,spalte.dz], 1, 4)) < 1000)
  { #  Standard
    ana.day  <- as.numeric(substr(dataset[ ,spalte.dz], 1, 2))
    ana.mon  <- as.numeric(substr(dataset[ ,spalte.dz], 4, 5))
    ana.yea  <- as.numeric(substr(dataset[ ,spalte.dz], 7, 10))
    ana.hour <- as.numeric(substr(dataset[ ,spalte.dz], 12, 13))
    ana.minu <- as.numeric(substr(dataset[ ,spalte.dz], 15, 16))
  } else
  { #  Alternative 1
    ana.day  <- as.numeric(substr(dataset[ ,spalte.dz], 9,10))
    ana.mon  <- as.numeric(substr(dataset[ ,spalte.dz], 6, 7))
    ana.yea  <- as.numeric(substr(dataset[ ,spalte.dz], 1, 4))
    ana.hour <- as.numeric(substr(dataset[ ,spalte.dz], 12, 13))
    ana.minu <- as.numeric(substr(dataset[ ,spalte.dz], 15, 16))
  }

  #  Put all variables in 'dataset' to ensure proper selection of records  
  #  Put date components in dataset
  dataset <- data.frame(dataset, ana.day, ana.mon, ana.yea, ana.hour,
                        ana.minu, stringsAsFactors=FALSE) 

  #  Delete single date / time variables 
  rm(ana.day, ana.mon, ana.yea, ana.hour, ana.minu)

  #  If partitioning the data before/after a splitting date 
  #  or a time range selection is requested 
  #  Either a split date or strat.date + end.date must be given

  if ( ( is.na(start.date) & !is.na(end.date)) |
       (!is.na(start.date) &  is.na(end.date)) )
  { #  Incomplete request
    cat("\n +++ Incomplete range for analysis time range - is ignored +++\n",
        "\n     Input was start.date =", start.date, " end.date =", end.date,
        "\n")
  }

  if (!is.na(split.date) | (!is.na(start.date) & !is.na(end.date)))
  {
    #  Recode sampling date as # of days since 01Jan1960. Needs library 'date'.
    ana.date.n  <- as.numeric(mdy.date(dataset[ ,"ana.mon"],
                                       dataset[ ,"ana.day"], 
                                       dataset[ ,"ana.yea"]))

    ana.time.n  <- 60 * as.numeric(dataset[ ,"ana.hour"]) +
                        as.numeric(dataset[ ,"ana.minu"])
 
    ana.dt <- 1440 * ana.date.n + ana.time.n

    dataset <- data.frame(dataset, ana.date.n)
 
    if (!is.na(split.date))
    { #  Split was requested
      # Numerical version of the split date, input format  "dd.mm.yyyy". 
      # Needs library 'date'.
      split.date.n <- as.numeric(
                      mdy.date(as.numeric(substr(split.date, 4, 5)),
                               as.numeric(substr(split.date, 1, 2)),
                               as.numeric(substr(split.date, 7,10))))
 
      # Numerical version of the split time, format  "hh.mm"
      split.time.n <- 60 * as.numeric(substr(split.time, 1, 2)) +
                           as.numeric(substr(split.time, 4, 5))

      split.dt     <- 1440 * split.date.n + split.time.n
      split.group  <- 1 + (ana.dt > split.dt) 

      #  Do the split
      ok <- split.group == use.split
      x.n2 <- sum(ok)

      if (x.n2 == 0)
      {
        cat("\n", rep("=", times=76), "\n", sep="")
        cat("\n ++++++ No more data after filtering split group ++++++\n\n")
        cat("\n", rep("=", times=76), "\n", sep="")
        stop("++++++ Filtering left no data ++++++") 
      }
      dataset  <- dataset[ok, ]

      #  clean up 
      rm(ana.date.n, ana.time.n, ana.dt) 

    }   # split requested?

    if (!is.na(start.date) & !is.na(start.date))
    {
      # Numerical version of the start date, input format  "dd.mm.yyyy". 
      start.date.n <- as.numeric(
                      mdy.date(as.numeric(substr(start.date, 4, 5)),
                               as.numeric(substr(start.date, 1, 2)),
                               as.numeric(substr(start.date, 7,10))))
      end.date.n <- as.numeric(
                      mdy.date(as.numeric(substr(end.date, 4, 5)),
                               as.numeric(substr(end.date, 1, 2)),
                               as.numeric(substr(end.date, 7,10))))

      #  Do the selection. Limits are included.
      ok <- (start.date.n <= dataset[ , "ana.date.n"]) &
            (dataset[ , "ana.date.n"] <= end.date.n)
      x.n2    <- sum(ok)

      if (x.n2 == 0)
      {
        cat("\n", rep("=", times=76), "\n", sep="")
        cat("\n ++++++ No more data after filtering time interval ++++++\n\n")
        cat("\n", rep("=", times=76), "\n", sep="")
        stop("++++++ Filtering left no data ++++++") 
      }
      dataset  <- dataset[ok, ]

    }   # time interval requested

  }     # split or time interval requested
 
  # --------------------------------------------------------------------------
  #  Special need 20.04.2020
  #write.table(dataset, file="../temp/tnth_UKSH_Kiel_split1.csv",
  #            row.names=FALSE, col.names=TRUE, sep=";", quote=FALSE)

  # --------------------------------------------------------------------------
  #  Filter for weekday or daytime 
  ana.wday <- as.numeric(format(as.Date(paste(dataset[ ,"ana.yea"],"-", 
                                              dataset[ ,"ana.mon"],"-", 
                                              dataset[ ,"ana.day"], 
                                              sep="")
                                       ), "%w" )
                                ) #  weekday, number, 0 = Sunday

  # Recode ana.wday: 0 = Monday
  ana.wday <- ana.wday - 1
  ana.wday[ana.wday==-1] <- 6
  
  dataset <- data.frame(dataset, ana.wday)
  rm(ana.wday)
 
  if (print.daytime.distr)
  {
    cat("\n", rep("-", times=76), "\n", sep="")
    cat("Distribution of 'analysis day' - all data\n")
    print(table(dataset[ ,"ana.day"]))
    cat("\n")

    cat("Distribution of 'analysis weekday' - all data\n")
    print(table(dataset[ ,"ana.wday"]))
    cat("\n")

    cat("Distribution of 'analysis month' - all data\n")
    print(table(dataset[ ,"ana.mon"]))
    cat("\n")  

    cat("Distribution of 'analysis year' - all data\n")
    print(table(dataset[ ,"ana.yea"]))
    cat("\n")  

    cat("Distribution of 'analysis hour' - all data\n")
    print(table(dataset[ ,"ana.hour"]))
    cat("\n")  

    cat("Distribution of 'analysis minu' - all data\n")
    print(table(dataset[ ,"ana.minu"]))
    cat("\n")  
  }

  #  If selection of daytime is required: do it

  if (!is.na(ana.hour.min))
  {
    # Analysis daytime in minutes
    dt.minutes <- dataset[ ,"ana.hour"]*60 + dataset[ ,"ana.minu"]
    
    dt.min     <- ana.hour.min*60 + ana.minu.min
    dt.max     <- ana.hour.max*60 + ana.minu.max

    ok <- (dt.min <= dt.minutes) & (dt.minutes < dt.max)
    ok[is.na(ok)] <- FALSE

    x.n3 <- sum(ok)

    rm(dt.minutes)  

    if (x.n3 == 0)
    {
      cat("\n", rep("=", times=76), "\n", sep="")
      cat("\n ++++++ No more data after filtering daytime ++++++\n\n")
      cat("\n", rep("=", times=76), "\n", sep="")
      stop("++++++ Filtering left no data ++++++") 
    }
 
    dataset <- dataset[ok, ]

    cat("\n", rep("-", times=76), "\n", sep="")
    cat("Distribution of 'analysis day' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.day"]))
    cat("\n")

    cat("Distribution of 'analysis weekday' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.wday"]))
    cat("\n")

    cat("Distribution of 'analysis month' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.mon"]))
    cat("\n")  

    cat("Distribution of 'analysis year' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.yea"]))
    cat("\n")  

    cat("Distribution of 'analysis hour' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.hour"]))
    cat("\n")  

    cat("Distribution of 'analysis minu' - data in selected day-time range\n")
    print(table(dataset[ ,"ana.minu"]))

    cat("\n", rep("-", times=76), "\n", sep="")
  }

  # --------------------------------------------------------------------------
  #  If a selection of weekdays is required: do it
  ok <- FilterWDay(dataset[ ,"ana.wday"], use.wday)
  x.n4 <- sum(ok)

  if (x.n4 == 0)
  {
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after filtering weekdays ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }
  
  dataset <- dataset[ok, ]
}

# =================================================================
#  Filter invalid records: missing or invalid values, missing sex, age

#  valid values ..................................................

#  Messwert muss numerisch und > 0 sein. Problem: Einträge "< 0.03"

v <- dataset[ ,spalte.w]

#  values < DL ..................................................

# Process non-numerical entries. See ProcessDL().
vDL <- ProcessDL(v, s.fact, option=DL.option)

#  v may still contain invalid entries

v <- as.numeric(vDL[["data"]])

detect.limits.max <- vDL[["detect.limits.max"]]    

#  invalid values ..................................................

#  Remove records with missing analysis value from the dataset
ok                  <- !is.na(v) 
x.n5                <- sum(ok)

if (x.n5 == 0)
{  
  cat("\n", rep("=", times=76), "\n", sep="")
  cat("\n ++++++ No more data after removing NAs ++++++\n\n")
  cat("\n", rep("=", times=76), "\n", sep="")
  stop("++++++ Filtering left no data ++++++") 
}

dataset[ ,spalte.w] <- v
dataset             <- dataset[ok, ]

#  negative  values ..................................................
#  Remove records with negative analysis value from the dataset
ok                <- dataset[ ,spalte.w] >= 0    
x.n6              <- sum(ok)

if (x.n6 == 0)
{  
  cat("\n", rep("=", times=76), "\n", sep="")
  cat("\n ++++++ No more data after removing records with negative value ++++++\n\n")
  cat("\n", rep("=", times=76), "\n", sep="")
  stop("++++++ Filtering left no data ++++++") 
}

dataset             <- dataset[ok, ]
 
#  Scale values  ...................................................
#  (Zeroes remain zeroes)
#  done before rounding

v <- dataset[ ,spalte.w]
cat("\nValue range before scaling             ", min(v,na.rm=TRUE), " - ", 
    max(v,na.rm=TRUE))
v <- v * scale.fact

cat("\nValue range after  scaling             ", min(v,na.rm=TRUE), " - ", 
    max(v,na.rm=TRUE))

dataset[ ,spalte.w] <- v

#  Round values  ...................................................

v <- Round(dataset[ ,spalte.w], round.unit)
dataset[1:10 ,spalte.w]

cat("\n", rep("-", times=76), "\n", sep="")

cat("\nValue range after  rounding            ", min(v,na.rm=TRUE), " - ", 
    max(v,na.rm=TRUE))

dataset[ ,spalte.w] <- v

#  Check for inconsistent rounding .................................

ldt <- CheckRounding(dataset[ ,spalte.w])

cat("\n", rep("-", times=76), "\n", sep="")
cat("\n Frequencies of last digit \n")
print(ldt)

ldt.p <- 1-pchisq(sum(ldt[ ,"chi2"]), df=9) 
cat("\n p value for H0 'Last digits are equally likely' ", ldt.p,"\n")
if (ldt.p < 0.05)
{
  cat("\n *** Uniform rounding for the whole data set is questionable ***\n")
}      

#  Restrict to allowed value range  .........................................
#  Limits apply to scaled values

ok <- rep(TRUE,times=nrow(dataset))  
if (!is.na(x.lo.limit))
{ 
  ok <- ok & (x.lo.limit <= dataset[ ,spalte.w] ) 
}

if (!is.na(x.hi.limit))
{ 
  ok <- ok & (dataset[ ,spalte.w] <= x.hi.limit) 
}
x.n7  <- sum(ok)

if (x.n7 == 0)
{  
  cat("\n", rep("=", times=76), "\n", sep="")
  cat("\n ++++++ No more data after restriction to allowed range ++++++\n\n")
  cat("\n", rep("=", times=76), "\n", sep="")
  stop("++++++ Filtering left no data ++++++") 
}

#  Delete records outside the specified range
dataset <- dataset[ok, ]

#  Treatment of zeros  ...................................................
#  (if any are left)

v <- dataset[ ,spalte.w]

#  Replace zeroes by fastnull
ist.null <- (v==0)

cat("\n", rep("-", times=76), "\n", sep="")

cat("\n# of values == 0 before zero removal   ",sum(ist.null,na.rm=TRUE))
cat("\nValue range before zero removal        ", min(v,na.rm=TRUE), " - ", 
    max(v,na.rm=TRUE))

# v[ist.null] <- fastnull    # changed 20.11.2020
v[ist.null] <- 0.1 * round.unit

ist.null <- (v==0)
cat("\n# of values == 0 after zero removal    ",sum(ist.null,na.rm=TRUE))
cat("\nValue range after zero removal         ", min(v,na.rm=TRUE), " - ", 
    max(v,na.rm=TRUE),"\n\n")

dataset[ ,spalte.w] <- v

rm(ist.null)

# ...................................................................

x.n8 <- nrow(dataset)

#  Wenn Sex angegeben, Sätze mit fehlenden Werten von Sex löschen ...
if (!is.na(spalte.s))
{ #  remove records with NA
  ok <- !is.na(dataset[ ,spalte.s]) 
  x.n8 <- sum(ok) 
  if (x.n8 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after removal of sex == NA ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }

  dataset <- dataset[ok, ]
  
  #  Eventuelle Eigenmächtigkeit von Excel beheben
  if (!is.na(spalte.s))
  { 
    dataset[dataset[ ,spalte.s] == "FALSE", spalte.s] <- "F" 
  }
}

#  sex range  .............................................................

x.n9 <- nrow(dataset)

# Remove unwanted sex categories
if (!is.na(spalte.s))
{ 
  #  Provide overview
  sex.table <- table(dataset[ ,spalte.s])
  cat("\n", rep("-", times=76), "\n", sep="")
  cat("Sex distribution, sex codes as in data, prior to filter 3\n")
  print(sex.table)
  cat("\n")

  #  Filtern - wenn gültige Werte zum Filtern angegeben
  #  Wenn nicht, ist use.sex.n == 0
  #  Filtering uses the original sex codes in the data
  if (use.sex.n > 0)
  { ok   <- dataset[ ,spalte.s] == use.sex[1]
    isex <- 1
    while(isex < use.sex.n)
    { isex <- isex + 1
      ok   <- ok | (dataset[ ,spalte.s] == use.sex[isex])
    }
  }
  x.n9  <- sum(ok)
  if (x.n9 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after filtering desired sex codes ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }
  dataset <- dataset[ok, ]
}

#  recode sex codes  ...................................................

if (!is.na(spalte.s))

{ #  Transform sex code for females to 'F', if the coding in the data 
  #  is not already 'F', and the analogue for males 
  #  Applies only if use.sex is not NA (which is the case in simulations)

  if ( !((length(use.sex) == 1) && is.na(use.sex)) )
  { 
    dataset[dataset[ ,spalte.s] == info[idx,"SexCodeF"], spalte.s] <- "F"
    dataset[dataset[ ,spalte.s] == info[idx,"SexCodeM"], spalte.s] <- "M"

    # Update sex codes in use.sex 
    use.sex <- c("F", "M")
    sex.val <- use.sex

    if (eval.rep) { sex.val <- "All" }
  }

  sex.table <- table(dataset[ ,spalte.s])

  cat("\n", rep("-", times=76), "\n", sep="")
  cat("Sex distribution, sex codes transformed to standard, after filter 3\n")
  print(sex.table)
  cat("\n")
}

#  valid age  .........................................................

x.n10 <- nrow(dataset)

#  Wenn Alter angegeben, Sätze mit fehlenden Werten von Alter löschen
if (!is.na(spalte.a))
{ #  Falsche Einträge bei Alter könnten dazu führen, dass Alter als 
  #  character-Variable interpretiert wird. In numerisch umwandeln,
  #  dann auf NA prüfen.
  v                   <- as.numeric(as.character(dataset[ ,spalte.a]))
  dataset[ ,spalte.a] <- v

  ok <- (!is.na(dataset[ ,spalte.a])) 
  x.n10 <- sum(ok)

  if (x.n10 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after removal of age == NA ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }

  #  Sätze mit fehlenden Werten tatsächlich rauswerfen
  dataset <- dataset[ok, ]
}

#  age range  ...................................................

x.n11 <- nrow(dataset)

# Vierte Filterung: nicht erwünschte Werte von Alter entfernen 
if (!is.na(spalte.a))
{ 
  #  Erstmal einen Überblick
  age.table <- table(dataset[ ,spalte.a])

  if (print.age.distr)
  {
    cat("\n", rep("-", times=76), "\n", sep="")
    cat("Age distribution prior to age filter\n")
    print(age.table)
    cat("\n", rep("-", times=76), "\n", sep="")

    cat("\nCumulative age distribution prior to age filter\n")
    print(data.frame(age=as.numeric(names(age.table)),
                   count=as.numeric(age.table), 
                   cumcount=cumsum(age.table),
                   cdf=format(cumsum(age.table)/x.n11, digits=4) ))
    cat("\n")
  }

  ok <- rep(TRUE,times=x.n11)

  #  First replace surrogate values in age.class, if present
  if (is.na(age.class[1, "lo"]) |  age.class[1, "lo"] == -99999) 
  { 
    age.class[1, "lo"] <- as.numeric(names(age.table)[1])
  }
  if (is.na(age.class[age.class.n, "hi"]) | 
            age.class[age.class.n, "hi"] == 99999) 
  { 
    age.class[age.class.n, "hi"] <- as.numeric(tail(names(age.table),1))
  }

  #  Select data
  ok <- ok & (age.class[1, "lo"] <= dataset[ ,spalte.a]) 
  ok <- ok & (dataset[ ,spalte.a] <= age.class[age.class.n, "hi"])

  x.n11 <- sum(ok)
  if (x.n11 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after filtering desired age range ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }
  #  Remove unwanted data
  dataset <- dataset[ok, ]
}

#  outpatient / hospitalized  .............................................

x.n12 <- nrow(dataset)

# Filter outpatient / hospitalized filtern, if requested
if ( (length(spalte.oh) > 0) &&  !is.na(spalte.oh))
{ 
  #  Erstmal einen Überblick
  oh.table <- table(dataset[ ,spalte.oh])

  cat("\n", rep("-", times=76), "\n", sep="")
  cat("Distribution of outpatient / hospitalized, prior to filter 5\n")
  print(oh.table)
  cat("\n")

  ok <- rep(TRUE,times=x.n12)

  #  Filtern - wenn gültige Werte zum Filtern angegeben
  #  Wenn nicht, ist use.oh.n == 0
  if (use.oh.n > 0)
  { ok   <- dataset[ ,spalte.oh] == use.oh[1]
    ioh <- 1
    while(ioh < use.oh.n)
    { ioh <- ioh + 1
      ok   <- ok | (dataset[ ,spalte.oh] == use.oh[ioh])
    }
  }
  x.n12  <- sum(ok)
  if (x.n12 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after filtering outpatient/ hospitalized ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }
  dataset <- dataset[ok, ]
}

#  device  ...................................................

x.n14 <- nrow(dataset)

# Filter device filtern, if requested
if ((length(spalte.dev) > 0) &&  !is.na(spalte.dev))
{ 
  #  Erstmal einen Überblick
  dev.table <- table(dataset[ ,spalte.dev], useNA="ifany")

  cat("\n", rep("-", times=76), "\n", sep="")
  cat("Distribution of 'device', prior to filter 6\n")
  print(dev.table)
  cat("\n")

  ok <- rep(TRUE,times=x.n14)

  #  Filtern - wenn gültige Werte zum Filtern angegeben
  #  Wenn nicht, ist use.dev.n == 0
  if (use.dev.n > 0)
  { ok   <- dataset[ ,spalte.dev] == use.dev[1]
    idev <- 1
    while(idev < use.dev.n)
    { idev <- idev + 1
      ok   <- ok | (dataset[ ,spalte.dev] == use.dev[idev])
    }
  }
  dataset <- dataset[ok, ]
  x.n14  <-  nrow(dataset)

  if (x.n14 == 0)
  {  
    cat("\n", rep("=", times=76), "\n", sep="")
    cat("\n ++++++ No more data after filtering device ++++++\n\n")
    cat("\n", rep("=", times=76), "\n", sep="")
    stop("++++++ Filtering left no data ++++++") 
  }
}

# ---------------------------------------------------------------------------
#  Generate subsample, if requested. 
#  Default setting in  seg015_DefaultSettings.R
#  Actual settng should be done in the start file
if (!is.na(subsample.n))
{ #  Random sample. Note that dataset might be sorted by spalte.w.
  #  Make random ordering reproducible by setting a seed for runif()
  #  set.seed(73271)   #  no, gives misleading selection    02.08.2021

  ok <- data.frame(1:x.n14, runif(x.n14))
  ok <- ok[order(ok[ ,2]), ]

  #  Data subset might be smaller than requested subsample size
  subsample.n.eff <- min(x.n14, subsample.n)
  if (subsample.n.eff < subsample.n)
  {
    cat("\n", rep("-", times=76), "\n", sep="")
    cat("\n",
        "\n *** Requested subsample of size ", subsample.n,
        "\n     is larger than this stratum (size =", x.n14, ").",
        "\n     The complete stratum will be used.",
        "\n")
    cat("\n", rep("-", times=76), "\n", sep="")
  }

  dataset <- dataset[ok[1:subsample.n.eff, 1], ] 
}
x.n   <- nrow(dataset)

cat("\n", rep("-", times=76), "\n", sep="")
cat("\n***  ... Global filters applied to input data\n")

# ---------------------------------------------------------------------------

#  Wenn verlangt, den Datensatz auf First Value reduzieren.
#  Braucht die Angabe der Spalte für Pat-Idendifikation (spalte.i)
#  und den Namen der Ausgabedatei (DataFileFirstValue). 
#  Beides zunächst noch in speziellen Optionen untergebracht.
#  If requested: reduce data to first values. Only possible, if
#  spalte.i  column for patient identification
#  spalte.d  column for analysis day
#  DataFileFirstValue  name of data file for first values 
#  are all !NA.
#  Default setting in  seg015_DefaultSettings.R
#  Actual setting should be done in the start file

if (!is.na(DataFileFirstValue))
{ 

  #  Nach Patid und Alter / Datum sortieren
  dataset <- dataset[order(dataset[ ,spalte.i], dataset[ ,spalte.d]), ]

  cat("Rohdaten, 1:100, sortiert nach PatId und Alter/Datum \n")
  print(dataset[1:100, ])

  #  Nur den ersten Wert von jeder PatId behalten
  pid1 <- dataset[1:(x.n-1),2]
  pid2 <- dataset[2:x.n,2]

  ok <- (pid1 != pid2)
  sum(ok)

  pid1.table <- table(pid1)
  sum(pid1.table == 1)
  length(pid1.table)

  dataset <- dataset[ok, ]

  write.table(dataset,file=DataFileFirstValue,sep=";",
            col.names=TRUE,row.names=FALSE)
}

# ----------------------------------------------------------------------------
#  Final parameter settings 

source("TMC_seg075_FinalParSettings.R")

# ----------------------------------------------------------------------------
#  Wertebereich 

x.min <- min(dataset[ ,spalte.w])
x.max <- max(dataset[ ,spalte.w])

#  Wie differenziert sind die Daten (wieviel verschiedene Werte) ?
x.table <- table(dataset[ ,spalte.w], useNA="ifany")
x.val.n <- length(x.table)

cat("\n", rep("-", times=76), "\n", sep="")
cat("\n  Data input summary\n")
cat("Data file ",infile,"\n")
cat("Data lines read:                                      ",x.n1,"\n")
cat("# in split range:                                     ",x.n2,"\n")
cat("# in daytime range:                                   ",x.n3,"\n")
cat("# in weekday range:                                   ",x.n4,"\n")
cat("# with non-missing measurements:                      ",x.n5,"\n")
cat("# with nonnegative measurements:                      ",x.n6,"\n")
cat("# in accepted value range:                            ",x.n7,"\n")
cat("# with non-missing values for 'sex':                  ",x.n8,"\n")
cat("# with accepted values for 'sex':                     ",x.n9,"\n")
cat("# with non-missing values for 'age':                  ",x.n10,"\n")
cat("# with accepted values for 'age':                     ",x.n11,"\n")
cat("# with accepted values for 'outpatient/hospitalized': ",x.n12,"\n")
cat("# with accepted values for 'device':                  ",x.n14,"\n")
cat("# left after subsampling:                             ",x.n,"\n")
cat("Final range of measurement values:                    ",x.min,x.max,"\n")
cat("# of distinct measurement values:                     ",x.val.n,"\n")
cat("\n", rep("-", times=76), "\n", sep="")

#  Die Häufigkeiten der 10 kleinsten Werte angucken: Hinweis auf
#  nicht deklarierte Werte unter NWG?

cat("\nFrequencies of the ten smallest values\n")
print(x.table[1:10])
cat("\n")

#  Look at largest 10 values
cat("\nFrequencies of the ten largest values\n")
print(tail(x.table, 10))
cat("\n")

#  Most frequent value
cat("\nMost frequent 10 values\n")
print(-sort(-x.table)[1:10])
cat("\n")

cat("\nQuantiles of the dataset\n")
print(Quantile(dataset[ ,spalte.w], probs=c(0.01, 0.025, 0.050, 0.25, 0.50, 
                                            0.75, 0.95, 0.975, 0.990)))
cat("\n")
cat("\n", rep("-", times=76), "\n", sep="")

#  For test data with known group membership
if (!is.na(spalte.g))
{ 
  grp.table <- table(dataset[ ,spalte.g])
  cat("\n", rep("-", times=76), "\n", sep="")
  cat("Distribution of 'Group' after global filtering\n")
  print(grp.table)
  cat("\n")
  cat("\n", rep("-", times=76), "\n", sep="")
  
  #  If spalte.g is given, more variables describing the data generation 
  #  should be given - see seg050_master
}

# .................................................................
#  Clean up. 'dataset' is the essential remaining data frame

rm(ok,v)

# =================================================================
#  Sort data. Is required for further processing-

dataset <- dataset[order(dataset[ ,spalte.w]), ]

# =================================================================
####  Globale Bedingungen für die Auswertung sind festgelegt

# ===========================================================================
if (print.log.message) { cat("%%%   TMC_seg70_ReadData  End\n") }
# ===========================================================================

