load("Data Processing V1.RData") #load d

# Delete records from people without date of onset or before onset
d = d[!is.na(YTime),] #21841 records from 1202 patients
d = d[YTime >= 0,] # reduces to 21802 records

pd = d[, c("Patient.ID","YTime", 
            "FVC", "DLCO", "RVSP", "EF", "mRSS",
            "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
            "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
            "pft", "mrss", "echo", "med"), with=F]

mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")

nl = 1:4; levlist = c("no trt","MMF","MMF+","others")

pd$nonMMF = apply(pd[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)
pd[, type := 1] # no trt
pd[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
pd[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
pd[MMF == 1 & nonMMF == 0, type := 2] 
pd[MMF == 1 & HCQ == 1, type := 2] # 2 MMF only OR MMF + HCQ

pd[med == 0, type :=NA] # type 9956 NA's

# Define indicator of records BEFORE the 1st medication, medtype_n1 

library(BBmisc)
pd[, tmp := which.first(!is.na(type)) -1, by = Patient.ID]
pd[, medtype_n1 := 0]
pd[, num := 1:(.N), by = Patient.ID]; pd[num <= tmp[1], medtype_n1 := 1, by = Patient.ID]
pd[, num := NULL]

# Carry forward medication type; NA if the observed is NA and the gap > 180 days
pd[, cbtype :=na.locf(type,na.rm = F,fromLast = TRUE), by = Patient.ID]
pd[, mTime := YTime]; pd[med == 0, mTime := NA]
pd[, cbTime := na.locf(mTime,na.rm = F,fromLast = TRUE), by = Patient.ID]; pd[, mTime := NULL]

# medtype: 5 missing before 1st med, 0 missing after 1st med, 1 no trt, 2 MMF only or (MMF, HCQ), 3 MMF+, 4 others
pd[, medtype := cbtype]
pd[, gap := cbTime -YTime]; pd[gap >= (180/365), medtype := NA]
pd[medtype_n1 == 1, medtype := 5]
pd[is.na(medtype) , medtype :=0] # all of these has medtype_n1 == 0, i.e. post 1st med

#pd[, c("gap", "cftype", "cfTime") := NULL]

pd = pd[, c("Patient.ID", "YTime", "FVC", "DLCO", "RVSP", "EF", "mRSS", "medtype", "medtype_n1")]

########################################################################

# create medtype dummy variables

tmp = class.ind(pd[, medtype]); colnames(tmp) = paste0("medtype_", colnames(tmp))
pd = cbind(pd, tmp[,c(1,3,4,5)]); pd[, medtype :=NULL]

save(pd, file = "Data Processing for Outcome Treatment Models V1.RData")
