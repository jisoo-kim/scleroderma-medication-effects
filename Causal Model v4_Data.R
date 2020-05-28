##### Packages ##### 
library("easypackages")
my_packages = c("plyr","tidyr","ggplot2","gridExtra","openxlsx","lubridate","data.table",
                "VIM","poLCA","lme4","splines","MCMCglmm","nnet","zoo",
                "plot.matrix")
#packages(my_packages)
libraries(my_packages)

##### Functions ##### 

datf <- function(x){ x <- ifelse(x == 0, NA, x)
x <- as.Date(x, origin = "1899-12-30")
return(x)}

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}
#scalevec <- function(mvec){mvec = mvec[!is.na(mvec)];qnorm(order(order(mvec))/length(mvec))}
#scalevec <- function(mvec){qnorm(ppoints(length(mvec)))[order(order(mvec))]}

backscale <- function(mvec, obs){
  obs = obs[!is.na(obs)]; obs = sort(obs); 
  mvec1 = mvec[!is.na(mvec)]; ind = which(!is.na(mvec))
  mvec[ind] = obs[ceiling(pnorm(mvec1)*length(obs))]
  return(mvec)}

##### Data Read-in ##### 
visdat <- read.xlsx(filename, sheet = 1)
meds <- read.xlsx(filename, sheet = 3)

mrss <- read.xlsx(filename, sheet = 4)
echo <- read.xlsx(filename, sheet = 11)
pft <- read.xlsx(filename, sheet = 14)


##### Medication data ##### 

# Medication Date

meds$Meds.Date= datf(meds$Meds.Date)

# Extract the Seven Medications

## data frame -> data table
setDT(meds) 
## Missingness set to zero
mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
tmp = meds[, mv, with=F]; allNAind = which(apply(tmp, 1, function(x) all(is.na(x))))
## Extract med = (ID, Date, meds)
med = meds[-allNAind,]
for(m in mv){
  med[is.na(get(m)), (m) := 0]
}
med = med[, c("Patient.ID","Meds.Date",mv), with=F] # 11859 13


########################################################### 

##### Combine Data (Med, Lung) ##### 

# Process pft date, ID; turn pft into data table
colnames(pft)[2] = "pftDate"; pft$pftDate = datf(pft$pftDate)
pft$Pt.ID = as.numeric(pft$Pt.ID)
setDT(pft)

# Remove outliers in pft
pft[is.na(stppFVC) & stppFVC >= 200, stppFVC := NA]
pft[is.na(stppDLCO) & stppDLCO >= 200, stppDLCO := NA]

# Remove missing entries
pft = pft[is.na(stppFVC) + is.na(stppDLCO) < 2, ] # 8076 -> 7920

# Define bmi
pft[, bmi := Weight/(Height)^2 ]

# Combine med and pft (unique patients: med 1201, pft 1041)
## Create temporary IDs
med[, tmpid := paste0(Patient.ID, "_", Meds.Date)]
pft[, tmpid := paste0(Pt.ID, "_", pftDate)]
## Merge med and pft into dat
setkey(med, tmpid); setkey(pft, tmpid)
dat = merge(med, pft, all = TRUE)
## Define common Date
dat$Date = dat$pftDate; dat[is.na(pftDate), Date := Meds.Date]
## Process ID and sort
dat[is.na(Patient.ID), Patient.ID := Pt.ID]; dat[, c("Pt.ID") := NULL]
dat = dat[order(Patient.ID, Date),]
## Define indicator of data source
dat[, pft := !is.na(pftDate)]; dat[, med := !is.na(Meds.Date)]
dat[, (c("pftDate", "Meds.Date")) := NULL]
# dat[,.N, by=list(pft, med)]
## 2061 perfect merge, 9798 with only med, 5859 with only pft

# dat: 17718 records, 1207 patients


##### Combine Data (Med, Lung, Heart) ##### 

# Process echo date, ID; turn echo into data table
setDT(echo)
colnames(echo)[1] = "Pt.ID"
echo$echoDate = datf(echo$Date.of.ECHO); echo[, "Date.of.ECHO":=NULL]

# Remove missing entries
echo = echo[is.na(RVSP) + is.na(Ejection.Fraction) < 2, ] # 5599-> 5462

# Combine dat and echo (unique patients: dat 1207, echo 983)
## Create temporary IDs
echo[, tmpid := paste0(Pt.ID, "_", echoDate)]
## Merge (med, lung) and heart into dat
setkey(dat, tmpid)
setkey(echo, tmpid)
dat = merge(dat, echo, all=TRUE)
## Process Date and ID
dat[is.na(Date), Date := echoDate]
dat[is.na(Patient.ID), Patient.ID := Pt.ID]; dat[, c("Pt.ID") := NULL]
## Define indicator of data source
dat[, echo := !is.na(echoDate)]; dat[, ("echoDate") := NULL]

# dat: 21567 records, 1207 patients

##### Combine Data (Med, Lung, Heart, Skin) ##### 

# Process skin date, ID; turn mrss into data table
setDT(mrss)
colnames(mrss)[1] = "Pt.ID"
colnames(mrss)[2] = "skinDate"; mrss$skinDate = datf(mrss$skinDate)

# Remove missing entries
mrss = mrss[!is.na(Total.Skin.Score), ] # 11851 -> 11818

# Combine dat and mrss (unique patients: dat 1207, mrss 1205)
## Create temporary IDs
mrss[, tmpid := paste0(Pt.ID, "_", skinDate)]
## Merge (med, lung) and heart into dat
setkey(dat, tmpid)
setkey(mrss, tmpid)
dat = merge(dat, mrss, all=TRUE)
## Process Date and ID
dat[is.na(Date), Date := skinDate]
dat[is.na(Patient.ID), Patient.ID := Pt.ID]; dat[, c("Pt.ID") := NULL]
## Define indicator of data source
dat[, mrss := !is.na(skinDate)]; dat[, ("skinDate") := NULL]

# dat: 21870 records, 1211 patients

##### Combine Data (Med, Lung, Heart, Skin, Visits) ##### 

setDT(visdat)

# Extract variables from visit data ("DateOnset", "DOB","Sex","ethnic","Race1")

visdat[, DateOnset := pmin(datf(DateOnRP), datf(Date1stSymptom), na.rm = T)]
visdat[, DateOnset1 := datf(DateOnRP)]
visdat[, DateOnset0 := datf(Date1stSymptom)]

# Further define baseline variables
visdat[,.N,by = 'SCL-70']; visdat[, SCL70 := get('SCL-70')]; visdat[SCL70 == 3, SCL70 := 1]; visdat[,.N,by = SCL70]
visdat[,.N,by = ACA]; visdat[ACA == 3, ACA := 1]; visdat[,.N,by = ACA]
visdat[,.N,by = RNAPol]; visdat[RNAPol == 3, RNAPol := 1]; visdat[,.N,by = RNAPol]

##### Select relevant covariates #####

tmp = visdat[, c("Patient.ID", "DateOnset", "DateOnset0", "DateOnset1", "DOB","Sex","ethnic","Race1", "ACA", "SCL70", "RNAPol")]; setDT(tmp)

# Combine the extracted visit data into med
setkey(tmp, Patient.ID); setkey(dat, Patient.ID)
dat = merge(dat, tmp, all.x=TRUE)
## Define Time since Onset
#dat[, YTime1 := round(as.numeric(Date - DateOnset1)/365,2)]
#dat[, YTime0 := round(as.numeric(Date - DateOnset0)/365,2)]
#dat[, YTime := as.numeric(Date - DateOnset)]
dat[, Time := round(as.numeric(Date - DateOnset0)/365,2)]# for skin

## Define Age based on Date of Birth; age: 16 ~ 96, median 55
dat$DOB = datf(dat$DOB)
#dat[, age := floor(as.numeric(Date - DOB)/365) ]
dat[, age := floor(as.numeric(DateOnset0 - DOB)/365) ]
dat[, DOB := NULL]

dat = dat[order(Patient.ID, Date),]

##### Clean Data for modeling ##### 
d = dat[, c("Patient.ID", "Date", "Time", "age", "Sex", "ethnic", "Race1", "bmi", 
            "ACA", "SCL70", "RNAPol",
            "stppFVC", "stppDLCO", "RVSP", "Ejection.Fraction", "Total.Skin.Score",
            "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
            "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
            "pft", "mrss", "echo", "med"), with=F]
colnames(d) =  c("Patient.ID", "Date", "Time",  "age", "Sex", "ethnic", "Race", "bmi", 
                 "ACA", "SCL70", "RNAPol",
                 "FVC", "DLCO", "RVSP", "EF", "mRSS",
                 "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
                 "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
                 "pft", "mrss", "echo", "med")

d[, pft := 1*(!is.na(pft) & pft == TRUE)]
d[, mrss := 1*(!is.na(mrss) & mrss == TRUE)]
d[, echo := 1*(!is.na(echo) & echo == TRUE)]
d[, med := 1*(!is.na(med) & med == TRUE)]

# 1211 patients from 21,870 rows

# which baseline characteristic(s) contain missingness?
PDlist = c()
for(v in c("age", "Sex", "ethnic", "Race", "ACA", "SCL70", "RNAPol")){
  numNA = sum(d[,  is.na(get(v))])
  print(paste0(v, " ", numNA) )
  if(numNA > 0) PDlist = c(PDlist, unique(d[is.na(get(v)), Patient.ID])  )
}
PDlist = unique(PDlist)
`%notin%` <- Negate(`%in%`)
d = d[Patient.ID %notin% PDlist,] # 1202 patients from 21,749 rows

# only keep lung, skin, med
d = d[pft + mrss + med > 0,] ## 1157 patients from 17,694 rows

d = d[order(Patient.ID, YTime)]

##### Define First Medication ##### 

library(BBmisc)
d[, tmp := which.first(med==1) , by = Patient.ID]
cvec = function(x,y) {x[y] = 1; x}
d[,firstmed := 0]; d[, firstmed := cvec(firstmed, tmp[1]), by = Patient.ID]
# table(d[, sum(med), by=Patient.ID]$V1) # 10 people do not have medication, delete them
d[, tmp := 1*(sum(med)==0),by = Patient.ID]
d = d[tmp == 0,]; d[, tmp := NULL]
#length(unique(d$Patient.ID)) # 1202 -> 1192 (17,887 rows)

##### Set Time = 0 to the 1st Medication Record #####
d[, MTime := as.numeric(Date - Date[firstmed == 1]), by = Patient.ID]


##### Define Outcome as cb (FVC, mRSS) #####
# Carry backward to the nearest medication record within 6 months

# Index medication records within patient
d[, ind := cumsum(med), by = Patient.ID]
# Define tmpid
d[, tmpid := paste0(Patient.ID, "_", ind)]

# Define cb outcome within 6 months of med
for(var in c("FVC", "mRSS")){
  v1 = paste0("cb", var)
  d[, (v1) := na.locf(get(var),na.rm = F, fromLast = TRUE), by = tmpid]
  d[, vart := Time]; d[is.na(get(var)), vart := NA]
  d[, tmp := na.locf(vart,na.rm = F, fromLast = TRUE), by = tmpid]
  d[, tmp := tmp - Time]
  # only cb for those <= 6 months
  d[!is.na(tmp) & tmp > 180 & is.na(FVC) & !is.na(cbFVC), cbFVC := NA] 
}
d[, c("vart", "tmp") := NULL]
#reduce missingness at medication records:
#FVC: 9750 -> 7027 ##sum(is.na(d$FVC[d$med==1])); sum(is.na(d$cbFVC[d$med==1]))
#mRSS: 413 -> 384 ##sum(is.na(d$mRSS[d$med==1])); sum(is.na(d$cbmRSS[d$med==1]))

##### Define Baseline Covariates dosage and cf (FVC, mRSS) #####

# Define dosage as the number of times a subject went on trt from enrollment to a med record
mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
d$nonMMF = apply(d[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)
d[, type := 1] # no trt
d[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
d[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
d[MMF == 1 & nonMMF == 0, type := 2] 
d[MMF == 1 & HCQ == 1, type := 2] # 2 MMF only OR MMF + HCQ
d[med==0, type := NA]

d[, tmp := 1*(type == 2)]; d[is.na(tmp), tmp:=0]
d[, MMFdos0 := cumsum(tmp), by = Patient.ID]
d[MMFdos0 > 5, MMFdos0 := 5]

# Define cf prev (FVC, mRSS)
## cfFVC0t: how many days ago was the cfFVC0 measured?
## cfmRSS0t: how many days ago was the cfmRSS0 measured?
for(var in c("cbFVC", "cbmRSS")){
  v1 = paste0("cf", var,"0"); v1t = paste0(v1,"t")
  d[, (v1) := na.locf(get(var),na.rm = F), by = Patient.ID]
  d[, tmp := Time]; d[is.na(get(var)), tmp:=NA]
  d[, (v1t) := na.locf(tmp,na.rm = F), by = Patient.ID]
  d[, (v1t) :=  Time - get(v1t)]
}

##### Extract Data and Variables #####
# Keep only those with 
# 1. Time >= 0 (post On Set)
# 2. Have at least 3 medication records before data closure

d1 = d[MTime >= 0,] #17664 -> 15985 rows
d1 = d1[med == 1,] #11643 rows

d1 = d1[, c("Patient.ID", "Time", "MTime", "cbFVC", "cbmRSS", "type", "MMFdos0",
            "age", "Sex", "ethnic", "Race", "YTime",
            "ACA", "SCL70", "RNAPol","cfcbFVC0", "cfcbmRSS0","cfcbFVC0t", "cfcbmRSS0t"), with=F]

d1[, N := .N, by=Patient.ID]

### Only look at PAIRS!!
d1 = d1[N>=2,] # 11472 rows from 977 people; keep those having at least 3 (med) records
d1 = d1[order(Patient.ID, Time)]

# Define FVCj, mRSSj, typej for j = 1,2
vn = c("FVC", "mRSS", "type", "Time")
vn1 = c("cbFVC", "cbmRSS", "type", "Time")
for(j in 1:4){
  var = vn1[j]
  d1[, (paste0(vn[j], 1)) := c(get(var)[-1], NA), by = Patient.ID] # next 1
  d1[, (paste0(vn[j], 2)) := c(get(paste0(vn[j], 1))[-1], NA), by = Patient.ID] # next 1
}
colnames(d1)[ which(colnames(d1)%in% c("Time","cbFVC", "cbmRSS","type"))] = c("Time0","FVC0", "mRSS0", "type0")

### Only take PAIRS!!!
d1 = d1[,c("Patient.ID", paste0("Time",0:1),
           paste0("FVC",0:1),paste0("mRSS",0:1),paste0("type",0:1),
           "MMFdos0",
           "age", "Sex", "ethnic", "Race", "MTime",
           "ACA", "SCL70", "RNAPol",
           "cfcbFVC0", "cfcbmRSS0","cfcbFVC0t", "cfcbmRSS0t"), with=F]


##### Identify Qualifed Baselines for PAIRS #####
# Criteria: subsequent visits each have gap times < 1yr
d1 = d1[!is.na(Time1),]# 11472 -> 10495
d1[, gap := Time1 - Time0]
d1[, qlf := 1*(gap<=1 & gap>90/365)] # on year scale
d1 = d1[qlf == 1,] # 9044 rows from 944 ppl

##### Variable Processing Before Modeling #####
d2 = copy(d1)
d2[, qlf := NULL]

# baseline health outcomes on year scale
for(v in c("cfcbFVC0t", "cfcbmRSS0t")){
  d2[, (v) := round(get(v)/365,2)]
}
# Sex
d2$Sex = as.numeric(d2$Sex)

# Race (1 white 2 black, use Race_1 and Race_2)
d2[is.na(Race), Race := 9]; d[Race == 9, Race := 0]
tmp = class.ind(d2[, Race]); colnames(tmp) = paste0("Race_",colnames(tmp))
d2 = cbind(d2, tmp[, -ncol(tmp)]); d2[, Race := NULL]

# Ethnic (use ethnic_0 and ethnic_1; indicator of latino)
d2[is.na(ethnic), ethnic := 9]
d2[, ethnic_0 := 1*(ethnic == 0)]; d2[, ethnic_1 := 1*(ethnic == 1)]; d2[, ethnic := NULL] 

# Define ILD
#d2[, ILD := 1*(cfcbFVC0 < 70)]

#V = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q","cfcbFVC0t", "cfcbmRSS0t", "Race_1", "Race_2", "ethnic_0", "ethnic_1")
V = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0", "cfcbmRSS0", "Race_1", "Race_2", "ethnic_0", "ethnic_1", "Time0")

misind = c()
for(v in V){
  print(v)
  print(sum(is.na(d2[,get(v)])))
  print(d2[,.N,by = get(v)])
  misind = c(misind, which(is.na(d2[,get(v)])))
}
misind = unique(misind)

# Only consider 3-visit windows that have no missingness in V
d3 = d2[-misind, ] # 8840 pairs from 874 

# Only look at binary trt for now
for(j in 0:1)
  d3[, (paste0("A",j)) := 1*(get(paste0("type",j)) ==2)]

###############################################
###############################################
###############################################
### Restrict to only those taking either no trt or MMF
d3 = d3[type0 %in% c(1,2) & type1 %in% c(1,2),]
#5577 pairs from 675 ppl

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

backscale <- function(mvec, obs){
  obs = obs[!is.na(obs)]; obs = sort(obs); 
  mvec1 = mvec[!is.na(mvec)]; ind = which(!is.na(mvec))
  mvec[ind] = obs[ceiling(pnorm(mvec1)*length(obs))]
  return(mvec)}

# transform health outcomes; DO THIS AT EXACTLY THE DATA FOR MODELING!!!!!!!!!!!
for(v in c(paste0("FVC",0:1),paste0("mRSS",0:1),"cfcbFVC0", "cfcbmRSS0")){
  d3[, (paste0(v,"q")) := scalevec(get(v))]
}

dat = copy(d3)
save(dat, file = "tmp_JK1.RData")


