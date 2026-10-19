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
# PREPARE THE DATA
################################################################################

# DOWNLOAD SYNTHETIC DATA FROM THE ZENODO REPOSITORY (IF NEEDED)
files <- c("synthbdcohortinfo", "synthbdbasevar", "synthpmdata", "synthoutdeath")
for(x in files) if(! paste0(x,".RDS") %in% list.files("data"))
  download_zenodo("10.5281/zenodo.13983169", path="data", files=paste0(x,".RDS"))

# LOAD THE DATA
bdcohortinfo <- readRDS("data/synthbdcohortinfo.RDS") |> as.data.table()
bdbasevar <- readRDS("data/synthbdbasevar.RDS") |> as.data.table()
pmdata <- readRDS("data/synthpmdata.RDS") |> as.data.table()
outdeath <- readRDS("data/synthoutdeath.RDS") |> as.data.table()

# DEFINE THE ANNUAL LAGS FOR PM2.5
pmdata[, paste0("pm25_",0:lag):=shift(pm25, 0:lag), by=eid]

# TRANSFORM BASELINE VARIABLES IN UNORDERED FACTORS (FOR REGRESSION MODEL)
ordvar <- names(bdbasevar)[sapply(bdbasevar, is.ordered)]
bdbasevar[, (ordvar):=lapply(.SD, factor, ordered=F), .SDcols=ordvar]

# MERGE COHORT DATA AND BASELINE VARIABLES
maindata <- merge(bdcohortinfo, bdbasevar[, c("eid","asscentre","sex")])

# CREATE YEAR AND MONTH OF BIRTH
maindata[, birthyear:=year(dob)]
maindata[, birthmonth:=paste(year(dob), month(dob), sep="-")]
