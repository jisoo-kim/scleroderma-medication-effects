load("Data Processing V1.RData") #load d

# Delete records from people without date of onset or before onset
d = d[!is.na(YTime),] #11661 records from 1155 patients
d = d[YTime >= 0,] # reduces to 21576 records

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

pd[med == 0, type :=NA] # type 10268 NA's

# Carry forward medication type; NA if the observed is NA and the gap > 180 days
pd[, cftype :=na.locf(type,na.rm = F), by = Patient.ID]
pd[, mTime := YTime]; pd[med == 0, mTime := NA]
pd[, cfTime := na.locf(mTime,na.rm = F), by = Patient.ID]; pd[, mTime := NULL]

pd[, medtype := cftype]
pd[, gap := YTime - cfTime]; pd[gap >= (180/365), medtype := NA]
#pd[, c("gap", "cftype", "cfTime") := NULL]

pd = pd[, c("Patient.ID", "YTime", "FVC", "DLCO", "RVSP", "EF", "mRSS", "medtype")]

########################################################################

# medtype: 0 missing, 1 no trt, 2 MMF only or (MMF, HCQ), 3 MMF+, 4 others
# create medtype dummy variables

pd[is.na(medtype), medtype := 0]
tmp = class.ind(pd[, medtype]); colnames(tmp) = paste0("medtype_", colnames(tmp))
pd = cbind(pd, tmp[,-5]); pd[, medtype :=NULL]

# quantitilize outcomes
if(0){
  qfun = function(x){
    
    ind = which(!is.na(x)); x1 = x[ind]
    sx1 = sort(x1); a = min(sx1, na.rm=T); b = max(sx1, na.rm=T)
    q = 1; l = length(sx1)-1
    for(i in 1:l){
      sx1[i] = q
      if(i < l & sx1[i] != sx1[i+1]) q = q + 1
      if(i == l){
        if(sx1[i] != sx1[i+1]) sx1[l+1] = q + 1
        if(sx1[i] == sx1[i+1]) sx1[l+1] = q 
      } 
    } # P(sx1 <= sx1[i])
    sx1 = sx1/max(sx1)
    x2 = rep(NA, length(x1)); x2[order(x1)] = sx1
    res = x; res[ind] = x2; return(res)
    #t # original vector
    #t[order(t)] # sorted
    #tt = rep(NA, length(t)); tt[order(t)] = sort(t) # goes from order back to t, all(tt==t) is TRUE
  }
  
  normalizefun = function(x){
    a = min(x,na.rm=T); b = max(x,na.rm=T)
    return((x-a)/b)
  }
  tmp =qfun(pd$DLCO)
  hist(normalizefun(pd$FVC))
  hist(normalizefun(pd$RVSP))
  hist(normalizefun(log(pd$mRSS + 1)))
  
  par(mfrow = c(1,3))
  hist(qfun(pd$FVC))
  hist(qfun(pd$RVSP))
  hist(qfun(pd$mRSS))
  
}

save(pd, file = "Data Processing for Outcome Treatment Models V1.RData")
