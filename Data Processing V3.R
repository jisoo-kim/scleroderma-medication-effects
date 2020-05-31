##### Packages ##### 
library("easypackages")
my_packages = c("plyr","tidyr","ggplot2","gridExtra","openxlsx","lubridate","data.table",
                "VIM","poLCA","lme4","splines","MCMCglmm","nnet","zoo",
                "plot.matrix","zoo")
#packages(my_packages)
libraries(my_packages)

##### Functions ##### 

datf <- function(x){ x <- ifelse(x == 0, NA, x)
x <- as.Date(x, origin = "1899-12-30")
return(x)}

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

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

d = d[order(Patient.ID, Time)]

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

# Sex
d$Sex = as.numeric(d$Sex)

# Race (1 white 2 black, use Race_1 and Race_2)
d[is.na(Race), Race := 9]; d[Race == 9, Race := 0]
tmp = class.ind(d[, Race]); colnames(tmp) = paste0("Race_",colnames(tmp))
d = cbind(d, tmp[, -ncol(tmp)]); d[, Race := NULL]

# Ethnic (use ethnic_0 and ethnic_1; indicator of latino)
d[is.na(ethnic), ethnic := 9]
d[, ethnic_0 := 1*(ethnic == 0)]; d[, ethnic_1 := 1*(ethnic == 1)]; d[, ethnic := NULL] 

save(d, file = "tmp_plot.RData")

