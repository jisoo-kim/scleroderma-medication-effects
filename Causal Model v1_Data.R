
##### Functions ##### 

datf <- function(x){ x <- ifelse(x == 0, NA, x)
x <- as.Date(x, origin = "1899-12-30")
return(x)}

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

backscale <- function(mvec, obs){obs = sort(obs); obs[ceiling(pnorm(mvec)*length(mvec))]}

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
for(m in mv){
  meds[is.na(get(m)), (m) := 0]
}
## Extract med = (ID, Date, meds)
med = meds[, c("Patient.ID","Meds.Date",mv), with=F]


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

# Combine med and pft (unique patients: med 1201, pft 1051)
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
#table(dat$pft, dat$med) 
## 2062 perfect merge, 9800 with only med, 5858 with only pft

# dat: 17720 records, 1207 patients


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

# dat: 21569 records, 1207 patients

##### Combine Data (Med, Lung, Heart, Skin) ##### 

# Process echo date, ID; turn echo into data table
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

# dat: 21872 records, 1211 patients

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
dat[, YTime1 := round(as.numeric(Date - DateOnset1)/365,2)]
dat[, YTime0 := round(as.numeric(Date - DateOnset0)/365,2)]
dat[, YTime := round(as.numeric(Date - DateOnset)/365,2)]

## Define Age based on Date of Birth; age: 16 ~ 96, median 55
dat$DOB = datf(dat$DOB)
dat[, age := round(as.numeric(Date - DOB)/365) ]
dat[, DOB := NULL]

dat = dat[order(Patient.ID, Date),]

##### Clean Data for modeling ##### 
d = dat[, c("Patient.ID", "Date", "YTime", "YTime0", "YTime1", "age", "Sex", "ethnic", "Race1", "Height", "Weight", 
            "ACA", "SCL70", "RNAPol",
            "stppFVC", "stppDLCO", "RVSP", "Ejection.Fraction", "Total.Skin.Score",
            "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
            "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
            "pft", "mrss", "echo", "med"), with=F]
colnames(d) =  c("Patient.ID", "Date", "YTime", "YTime0", "YTime1", "age", "Sex", "ethnic", "Race", "Height", "Weight", 
                 "ACA", "SCL70", "RNAPol",
                 "FVC", "DLCO", "RVSP", "EF", "mRSS",
                 "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
                 "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
                 "pft", "mrss", "echo", "med")

d[, pft := 1*(!is.na(pft) & pft == TRUE)]
d[, mrss := 1*(!is.na(mrss) & mrss == TRUE)]
d[, echo := 1*(!is.na(echo) & echo == TRUE)]
d[, med := 1*(!is.na(med) & med == TRUE)]

# Delete missing covariates in antibodies (ACA, SCL70, RNAPol) (6ppl) and Sex/ethnic/race (3ppl)

d = d[!is.na(RNAPol) & !is.na(ethnic),] # 21872 -> 21751 rows

# Number of outcome observation per person
vlist = c("FVC", "DLCO", "RVSP", "EF", "mRSS")
for(v in  vlist){
  d[, (paste0("n",v)) := sum(!is.na(get(v))), by = Patient.ID]
}
d[, nmed := sum(med), by = Patient.ID]
d[, nmin := pmin(nFVC, nDLCO, nRVSP, nEF, nmRSS, nmed)]
#d[, keep := 1*(nmin >=4)]; d[, .N, by = keep] # keep 13193 discard 8679, 
#length(unique(d[keep==1, Patient.ID])) # keep 371 patients out of 1211 

# scale outcomes
for(v in c("FVC", "DLCO", "RVSP", "EF", "mRSS")){
  if(v %in% c("RVSP", "mRSS")){
    d[, (paste0(v,"t")) := scalevec(-get(v))]
  } else {
    d[, (paste0(v,"t")) := scalevec(get(v))]
  }
  
}

###################################################
##### Define Baseline (time = 0) and Lay Grid #####

# Order by ID and YTime
d = d[order(Patient.ID, YTime),]
#tmpd = copy(d) # make a copy

# Identify first medical record
library(BBmisc)
d[, tmp := which.first(med==1) , by = Patient.ID]
cvec = function(x,y) {x[y] = 1; x}
d[,firstmed := 0]; d[, firstmed := cvec(firstmed, tmp[1]), by = Patient.ID]
# table(d[, sum(med), by=Patient.ID]$V1) # 10 people do not have medication, delete them
d[, tmp := 1*(sum(med)==0),by = Patient.ID]
d = d[tmp == 0,]; d[, tmp := NULL]
#length(unique(d$Patient.ID)) # 1202 -> 1192

# Define time based on t = 0 being the first medication
d[, tmp := Date[firstmed == 1], by = Patient.ID]
d[, Time := round(as.numeric(Date - tmp)/365,3)]
# sum(d$Time == 0) == length(unique(d$Patient.ID)) # Need this to be TRUE
# range(d$Time) # -27.592 ~ 17.570

# Mark last medication record
d[, tmp := which.last(med==1) , by = Patient.ID]
cvec = function(x,y) {x[y] = 1; x}
d[,lastmed := 0]; d[, lastmed := cvec(lastmed, tmp[1]), by = Patient.ID]

# Length of intervals
int_days = 0.5

# Lay grid between the first and the last medication record
## Discard patients with less than half a year (0.5) of span between the first and the last medication records
d[, span := as.numeric(Time[lastmed == 1] - Time[firstmed == 1]), by = Patient.ID]
d = d[span > int_days,]
# dim(d) # 21710 -> 21134 
# length(unique(d$Patient.ID))# 1192 -> 984

## Define grid data and merge
### HOW MANY GRID POINT FOR EACH PERSON?
length_in_care = d[, span[1], by = Patient.ID]
num_grid = floor(length_in_care$V1/int_days) # range(num_grid): 1-35

### DEFINE GRID DATA
grid_id = rep(length_in_care$Patient.ID, num_grid+1) # plus 1 for the zeros (baseline) 
grid_times = unlist(lapply(num_grid,function(x) c(0,(1:x)*int_days)))
grid_ind = unlist(lapply(num_grid,function(x) c(0,1:x))) 
grid_data = as.data.frame(cbind(  "Patient.ID" = grid_id, Time = grid_times, grid_ind = grid_ind))
grid_data = setDT(grid_data)
remove(grid_id,grid_times,grid_ind,num_grid,length_in_care)

### SET KEY AND MERGE 
d[, tmpid := paste0(Patient.ID,"_",Time)]
grid_data[, tmpid := paste0(Patient.ID,"_",Time)]
setkey(d,tmpid)
setkey(grid_data, tmpid)
d1 = merge(d, grid_data, all = TRUE)
#dim(d1)# 34812    56

remove(grid_data);gc()

### CLEAN SOME VARIABLES FROM MERGING, AND SORT
d1[, Patient.ID  := Patient.ID.x]
d1$Patient.ID[is.na(d1$Patient.ID)] = d1$Patient.ID.y[is.na(d1$Patient.ID)]

d1[,Time := Time.x]
d1[is.na(Time), Time := Time.y]

d1 = d1[order(Patient.ID, Time)]

d1[, c("Patient.ID.x", "Patient.ID.y", "Time.x", "Time.y", "tmpid", "tmp") := NULL]

rm(list = ls.str(mode = 'numeric'))

### DEFINE TMPID = PATIENT_ID X GRID_ID 
d1[, grid_ind_cb := na.locf(grid_ind,na.rm=F,fromLast = TRUE), by = Patient.ID] # grid intervals (,]
d1[, tmpgridid := paste0(Patient.ID,"_",grid_ind_cb)]
d1[, totrec := sum(!is.na(span)), by = tmpgridid]
d1[, totmed := sum(med,na.rm=T), by = tmpgridid]

### Covariates
#### age at baseline
d1[, age0 := age[Time==0], by = Patient.ID]
d1[, span := span[Time==0], by = Patient.ID]
d1[, bmi := 703*Weight/(Height)^2 ]
d1[, bmit := scalevec(bmi)]

### CARRY FORWARD FEATURES -> EXTRACT GRIDS -> DEFINE OUTCOME -> DEFINE COVARIATES FOR ANALYSIS
PatientX = c("Sex", "ethnic", "Race", "ACA", "SCL70", "RNAPol")
GridX = c("bmit", "bmi",
          "FVCt", "DLCOt", "RVSPt", "EFt", "mRSSt",
          "FVC", "DLCO", "RVSP", "EF", "mRSS",
          "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
          "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF")

#Xs = c("Patient.ID","Time","age0", PatientX, GridX,
#       "firstmed", "premed", "lastmed", "span", "grid_ind", "grid_ind_cb","tmpgridid", "totrec", "totmed" )
Xs = c("Patient.ID","Time","age0", PatientX, GridX,
       "span", "grid_ind", "tmpgridid", "totrec", "totmed" )

for(vn in PatientX){
  d1[, (vn) := na.locf(get(vn),na.rm = F), by = Patient.ID]
}

for(vn in GridX){
  d1[, (vn) := na.locf(get(vn),na.rm = F), by = tmpgridid]
}

d2 = d1[, Xs, with=F]

### Extract grids
d2 = d2[order(Patient.ID, Time)]
d2 = d2[ Time >=0, ]
d2 = d2[!is.na(grid_ind),]

### Take only the first two intervals
d3 = d2[grid_ind <=2,] #984 Patients
#### Define treatment type
mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
nl = 1:4; levlist = c("no trt","MMF","MMF+","others")

d3$nonMMF = apply(d3[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)
d3[, type := 1] # no trt
d3[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
d3[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
d3[MMF == 1 & nonMMF == 0, type := 2] 
d3[MMF == 1 & HCQ == 1, type := 2] # 2 MMF only OR MMF + HCQ
d3[is.na(nonMMF), type := NA]

#### Extract baseline characteristics, treatment type, outcomes, and time-varying confounders
PatientX = c("age0","Sex", "ethnic", "Race", "ACA", "SCL70", "RNAPol")
GridX = c("bmit", "bmi")
GridY = c("FVCt", "mRSSt", "FVC", "mRSS")
GridA = "type"
Xs = c("Patient.ID","Time",PatientX, GridX, GridY, GridA, "totmed")
d3 = d3[, Xs, with=F]
d3 = d3[order(Patient.ID, Time)]

#### Summary of visit (med) counts in each of the two intervals
tmp = d3[,.N, by =  list(totmed, Time)]
tmp[order(totmed, Time),]

##### PREPARE FOR CAUSAL MODELING AND SUMMARIES #####

# Discard those with missing baseline outcome measurement (keep only those w. no missingness in Y0)
d3[, ind := 1*(!is.na(FVCt[Time == 0]))*(!is.na(mRSSt[Time == 0])) , by = Patient.ID]
d3 = d3[ind == 1,] # 626 ppl, 2872 rows

# Delete the 33 people whose (Weight, Height, FVC) or (mRSS, medication) are not updated simultaneously
d3[, .N, list(is.na(mRSSt),is.na(type))]
d3[, .N, list(is.na(bmit), is.na(FVCt))]

tmp = d3[(is.na(bmit) != is.na(FVCt)) | (is.na(mRSSt) != is.na(type)),] #33 ppl
`%notin%` <- Negate(`%in%`)
d3 = d3[Patient.ID %notin% unique(tmp$Patient.ID),] # 593 ppl, 1751 rows

d3 = d3[, c("Patient.ID", "Time","bmit", "bmi","FVCt", "mRSSt", "FVC", "mRSS","type",
            "age0", "Sex", "ethnic", "Race", "ACA", "SCL70", "RNAPol"), with=F]

# Keep only those with at least one year of followup (two intervals)
tmp = d3[,.N,Patient.ID]
d4 = d3[Patient.ID %in% tmp[N==3, Patient.ID],] # 565 patients 1695 rows

# Summary of (A0, A1, A2)
tmp = d4[, list(type[Time==0],type[Time==0.5],type[Time==1]), by = Patient.ID]
tmp = tmp[,-1]
tmp[tmp==1] = "no trt"; tmp[tmp==2] = "MMF"; tmp[tmp==3] = "MMF+"; tmp[tmp==4] = "others"
see = tmp[,.N,by = list(V1,V2,V3)]
colnames(see) = c(paste0("med_type", c(0,0.5,1)),"N")

see[order(get(paste0("med_type", c(0,0.5,1))))]
 
# Plot of (Y0,Y1,Y2)|(A0,A1)
tmp = d4[, list(type[Time==0],type[Time==0.5],type[Time==1]), by = Patient.ID]
g1 = tmp[V1 ==1 & V2 == 1 , Patient.ID]; g2 = tmp[V1 ==2 & V2 == 2 , Patient.ID]
#g1 = tmp[V1 ==1 & V2 == 1 & V3 == 1, Patient.ID]; g2 = tmp[V1 ==2 & V2 == 2 & V3 == 2, Patient.ID]
pd = d4[Patient.ID %in% c(g1,g2),]
pd$col = 0; pd[Patient.ID %in% g2,col :=1] # col being 0s for (A0,A1) = (0,0), 1s for (A0,A1) = (1,1)

pcurve = function(pd, g1, g2, var){
  pd = d4[Patient.ID %in% g1,]; dat1 = matrix(pd[,get(var)],nrow = 3)
  pd = d4[Patient.ID %in% g2,]; dat2 = matrix(pd[,get(var)],nrow = 3)
  dat = cbind(dat1, dat2); cv = c(rep(1, ncol(dat1)), rep(2, ncol(dat2)))
  show(matplot(dat,type = c("b"),pch=1,col=cv,main=var))
  legend("topright", legend = c("(0,0)", "(1,1)"), col=1:2, pch=1)
}

par(mfrow = c(1,3))
pcurve(pd, g1, g2, "FVC")
pcurve(pd, g1, g2, "mRSS")
pcurve(pd, g1, g2, "bmi")


save(d4, file = "GridData.RData")

remove(dat, tmp)

remove(echo, mrss, meds, med, pft, visdat)

rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))