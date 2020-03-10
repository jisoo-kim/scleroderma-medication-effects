
##### Functions ##### 

datf <- function(x){ x <- ifelse(x == 0, NA, x)
x <- as.Date(x, origin = "1899-12-30")
return(x)}


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
## 2068 perfect merge, 9794 with only med, 6008 with only pft

# dat: 17870 records, 1207 patients


##### Combine Data (Med, Lung, Heart) ##### 

# Process echo date, ID; turn echo into data table
setDT(echo)
colnames(echo)[1] = "Pt.ID"
echo$echoDate = datf(echo$Date.of.ECHO); echo[, "Date.of.ECHO":=NULL]

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


##### Combine Data (Med, Lung, Heart, Skin) ##### 

# Process echo date, ID; turn echo into data table
setDT(mrss)
colnames(mrss)[1] = "Pt.ID"
colnames(mrss)[2] = "skinDate"; mrss$skinDate = datf(mrss$skinDate)

# Combine dat and mrss (unique patients: dat 1207, echo 1205)
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

# dat: 22130 records, 1211 patients


##### Combine Data (Med, Lung, Heart, Skin, Visits) ##### 

# Extract variables from visit data ("DateOnRP", "DOB","Sex","ethnic","Race1")
tmp = visdat[, c("Patient.ID", "DateOnRP", "DOB","Sex","ethnic","Race1")]; setDT(tmp)

# Combine the extracted visit data into med
setkey(tmp, Patient.ID); setkey(dat, Patient.ID)
dat = merge(dat, tmp, all.x=TRUE)
## Process onset date 
dat$DateOnRP = datf(dat$DateOnRP)
## Define Time since Onset
dat[, YTime := round(as.numeric(Date - DateOnRP)/365,2)]
## Define Age based on Date of Birth; age: 16 ~ 96, median 55
dat$DOB = datf(dat$DOB)
dat[, age := round(as.numeric(Date - DOB)/365) ]
dat[, DOB := NULL]

##### Clean Data for modeling ##### 
d = dat[, c("Patient.ID", "Date", "YTime", "age", "Sex", "ethnic", "Race1", "Height", "Weight", 
            "stppFVC", "stppDLCO", "RVSP", "Ejection.Fraction", "Total.Skin.Score",
            "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
            "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
            "pft", "mrss", "echo", "med"), with=F]
colnames(d) =  c("Patient.ID", "Date", "YTime", "age", "Sex", "ethnic", "Race", "Height", "Weight", 
                 "FVC", "DLCO", "RVSP", "EF", "mRSS",
                 "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
                 "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
                 "pft", "mrss", "echo", "med")

d[, pft := 1*(!is.na(pft) & pft == TRUE)]
d[, mrss := 1*(!is.na(mrss) & mrss == TRUE)]
d[, echo := 1*(!is.na(echo) & echo == TRUE)]
d[, med := 1*(!is.na(med) & med == TRUE)]

save(d, file = "Data Processing V1.RData")

remove(dat, tmp)

remove(echo, mrss, meds, med, pft, visdat)

rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))