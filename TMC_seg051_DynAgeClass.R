# TMC_seg051_DynAgeClass.R
#  (c) wwosniok@math.uni-bremen.de

#  Truncated minimum chi-square estimation

#  Find sex-age strata based on minimum stratum size = x.n.min, set in 030

#  #  CHANGE HISTORY

#  22.03.2022
# ==============================================================================

# names(dataset)
# [1] "PatId"    "Value"    "Age"      "Sex"      "DateTime" "OH"      
# [7] "Device"   "Grp"      "ValueRaw" "for.npa"  "ana.day"  "ana.mon" 
#[13] "ana.yea"  "ana.hour" "ana.minu" "ana.wday" "DM"  

#dataset[1, ]
#    PatId Value Age Sex            DateTime OH Device Grp ValueRaw for.npa
#687   687   133  53   F 22.01.2017 15:31:00  s    D02   C 132.9755   FALSE
#    ana.day ana.mon ana.yea ana.hour ana.minu ana.wday         DM
#687      22       1    2017       15       31        6 2017-01-15

#gtab[1, ]
#    irep Sex sexlabel     Age agelabel Age.num Age.mea method
#npa    1 F+M      F+M  18-100   18-100      83 58.9393    npa
#    subset.type     n pct.lt.DL x.tr.n prop n.per.bin.min n.per.bin.min.eff
#npa           1 10000         0     NA    1            NA                NA
# ....
