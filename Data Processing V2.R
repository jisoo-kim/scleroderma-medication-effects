
##### Functions ##### 

datf <- function(x){ x <- ifelse(x == 0, NA, x)
x <- as.Date(x, origin = "1899-12-30")
return(x)}

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

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

#visdat[, DateOnset1 := pmin(datf(DateOnRP), datf(Date1stSymptom), na.rm = T)]
visdat[, DateOnset1 := datf(DateOnRP)]
visdat[, DateOnset0 := datf(Date1stSymptom)]

# Further define baseline variables
visdat[,.N,by = 'SCL-70']; visdat[, SCL70 := get('SCL-70')]; visdat[SCL70 == 3, SCL70 := 1]; visdat[,.N,by = SCL70]
visdat[,.N,by = ACA]; visdat[ACA == 3, ACA := 1]; visdat[,.N,by = ACA]
visdat[,.N,by = RNAPol]; visdat[RNAPol == 3, RNAPol := 1]; visdat[,.N,by = RNAPol]

##### Select relevant covariates #####

tmp = visdat[, c("Patient.ID", "DateOnset0", "DateOnset1", "DOB","Sex","ethnic","Race1", "ACA", "SCL70", "RNAPol")]; setDT(tmp)

# Combine the extracted visit data into med
setkey(tmp, Patient.ID); setkey(dat, Patient.ID)
dat = merge(dat, tmp, all.x=TRUE)
## Define Time since Onset
dat[, YTime1 := round(as.numeric(Date - DateOnset1)/365,2)]
dat[, YTime0 := round(as.numeric(Date - DateOnset0)/365,2)]

## Define Age based on Date of Birth; age: 16 ~ 96, median 55
dat$DOB = datf(dat$DOB)
dat[, age := round(as.numeric(Date - DOB)/365) ]
dat[, DOB := NULL]

##### Clean Data for modeling ##### 
d = dat[, c("Patient.ID", "Date", "YTime0", "YTime1", "age", "Sex", "ethnic", "Race1", "Height", "Weight", 
            "ACA", "SCL70", "RNAPol",
            "stppFVC", "stppDLCO", "RVSP", "Ejection.Fraction", "Total.Skin.Score",
            "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
            "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
            "pft", "mrss", "echo", "med"), with=F]
colnames(d) =  c("Patient.ID", "Date", "YTime0", "YTime1", "age", "Sex", "ethnic", "Race", "Height", "Weight", 
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

#x = d$FVC; tmp = x[!is.na(x)]
#b = qqnorm(tmp)$x
#bb = qnorm(order(order(tmp))/length(tmp)) # this is good
#bb  = qnorm((order(order(tmp)) - 0.375)/ (length(tmp)+0.25)) # textbook
#bb = qnorm(ppoints(length(tmp)))[order(order(tmp))] # original qqnorm code
#par(mfrow = c(1,2));hist(b, breaks = 30);hist(bb, breaks = 30)

save(d, file = "Data.RData")

remove(dat, tmp)

remove(echo, mrss, meds, med, pft, visdat)

rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))