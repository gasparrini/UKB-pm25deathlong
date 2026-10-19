################################################################################
# R code for reproducing the analysis in:
#
# Vanoli J, et al. Long-term associations between time-varying exposure to 
#   ambient PM2.5 and mortality: an analysis of the UK Biobank. Epidemiology. 
#   2024;36(1):1-10. DOI: 10.1097/EDE.0000000000001796
# http://www.ag-myresearch.com/2025_vanoli_epidemiol.html
#
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/UKB-pm25deathlong
################################################################################

################################################################################
# DLM ANALYSIS
################################################################################

# DERIVE THE CROSS-BASES FOR PM2.5
# LINEAR EXP-RESP, SPLINE AND STRATA LAG-RESP
argvar <- list(fun="lin") 
arglag1 <- list(fun="ns", df=4, int=T)
arglag2 <- list(fun="strata", breaks=c(1,3,5), int=T)
cbpm1 <- crossbasis(subset(data, select=paste0("pm25_",0:lag)),
  lag=lag, argvar=argvar, arglag=arglag1)
cbpm2 <- crossbasis(subset(data, select=paste0("pm25_",0:lag)),
  lag=lag, argvar=argvar, arglag=arglag2)

# FIT THE MODELS UPDATING THE MAIN MODEL
dlnmmod1 <- update(mainmod, . ~ . - pm25_ma + cbpm1)
dlnmmod2 <- update(mainmod, . ~ . - pm25_ma + cbpm2)

# OBTAIN THE PREDICTIONS FOR A 10-UNIT INCREASE IN PM2.5
cppm1 <- crosspred(cbpm1, dlnmmod1, cen=0, at=pminc, bylag=0.1)
cppm2 <- crosspred(cbpm2, dlnmmod2, cen=0, at=pminc)

# EXTRAC THE ESTIMATES AND 95%CI FOR THE OVERALL CUMULATIVE HR
cppm1$allRRfit ; cppm1$allRRlow ; cppm1$allRRhigh
cppm2$allRRfit ; cppm2$allRRlow ; cppm2$allRRhigh
