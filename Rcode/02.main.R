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
# PERFORM THE MAIN ANALYSIS
################################################################################

# SELECT NON-ACCIDENTAL MORTALITY AS AN OUTCOME (CHANGE FOR OTHER OUTCOMES)
icd <- icdcode[[2]]
len <- icdlen[2]
outdata <- outdeath[substr(icd10,1,len) %in% icd]
setkey(outdata, eid, devent)

# MERGE WITH MAIN DATA
data <- merge(maindata, outdata[, c("eid","devent")], all.x=T)

# DEFINE THE EVENT AND EXIT TIME
data[, event:=(!is.na(devent) & devent<=dendfu) + 0]
data[, dexit:=fifelse(event==1, devent, dendfu)]

# EXCLUDE SUBJECTS WITH EVENT BEFORE THE START OF THE FOLLOW-UP
data <- data[dstartfu<dexit]

# SPLIT THE DATA BY CALENDAR YEAR USING FIRST DAY OF THE YEAR
cut <- year(range(data$dstartfu)[1]):year(range(data$dendfu)[2]) |>
  paste0("-01-01") |> as.Date() |> as.numeric()
data[, `:=`(dstartfu=as.numeric(dstartfu), dexit=as.numeric(dexit))]
data <- survSplit(Surv(dstartfu, dexit, event) ~., data, cut=cut) |> 
  as.data.table()

# ASSIGN THE YEAR AND THE YEAR OF LAG-0 EXPOSURE (MINUS ONE)
data[, year:= year(as.Date(dstartfu, origin=as.Date("1970-01-01")))]
data[, yearexp:= year-1]

# MERGE WITH PM DATA USING LAG-0 YEAR DEFINITION
# NB: OMITTING MISSING TO KEEP FULL LAG HISTORIES
setkey(data, eid, year)
setkey(pmdata, eid, year)
data <- merge(data, na.omit(pmdata), by.x=c("eid","yearexp"), 
  by.y=c("eid","year"))

# COMPUTE MOVING AVERAGE OF TIME-VARYING PM DATA
data[, pm25_ma:=rowMeans(.SD), .SDcols=paste0("pm25_",0:lag)]

# MERGE WITH BASELINE VARS
data <- subset(data, select=-c(sex,asscentre)) |> merge(bdbasevar, by="eid")
setkey(data, eid, year)

# CREATE THE MODEL FORMULA
fmod <- paste("Surv(dstartfu,dexit,event) ~ strata(asscentre,sex,birthyear) +", 
  paste(conf, collapse="+"), "+ pm25_ma") |> as.formula()

# FIT THE MODEL
# NB: THIS TAKES SLIGHTLY LONGER THAN 3 MIN WITH 1.30GHz/16GbRAM LAPTOP
mainmod <- coxph(fmod, data=data, ties="efron")
