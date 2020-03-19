load("Data.RData") #load d from Data Processing V2.R

# Use non-RP as onset
d$YTime = d$YTime0

## Delete records from people without date of onset or before onset
pd = d[!is.na(YTime) & YTime >= 0 & nmed > 0,] #21751 records from 1202 patients reduces to 21373 records from 1148 patients 

#pd = d[, c("Patient.ID","YTime", 
#           "FVC", "DLCO", "RVSP", "EF", "mRSS",
#           "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA",
#           "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
#           "pft", "mrss", "echo", "med"), with=F]

mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")

nl = 1:4; levlist = c("no trt","MMF","MMF+","others")

pd$nonMMF = apply(pd[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)
pd[, type := 1] # no trt
pd[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
pd[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
pd[MMF == 1 & nonMMF == 0, type := 2] 
pd[MMF == 1 & HCQ == 1, type := 2] # 2 MMF only OR MMF + HCQ

pd[med == 0, type :=NA] # type 9765 NA's

# Define indicator of records BEFORE the 1st medication, medtype_n1 

library(BBmisc)
pd[, tmp := which.first(!is.na(type)) -1, by = Patient.ID]
pd[, medtype_n1 := 0]

cvec = function(x,y) {x[1:y] = 1; x}
pd[,medtype_n1 := 0]; pd[, medtype_n1 := cvec(medtype_n1, tmp[1]), by = Patient.ID]

# Define indicator of the 1st medication record, T0

cvec= function(x,y) {x[y+1] = 1; x}
pd[, T0 := 0]; pd[, T0 := cvec(T0, tmp[1]), by = Patient.ID]

pd[, tmp := NULL]

#View(pd[Patient.ID == 78, c("type","T0", "medtype_n1")]) #check

# Percentage of mRSS missingness at the 1st medication record; 49/1148 = 4% of missingness

pd[T0 == 1, .N, by = is.na(mRSS)]

# Percentage of mRSS missingness BEFORE or at the 1st medication record; 

# Most recent outcomes' measurements (strictly) becfore the current record time, with the gap time in between 
for(v in c("FVC", "DLCO","RVSP","EF","mRSS","type")){
  pd[, (paste0(v,"cf")) := c(NA, get(v)[-.N]), by = Patient.ID]
  pd[, (paste0(v,"cf")) := na.locf(get(paste0(v,"cf")),na.rm = F),  by = Patient.ID]
  
  pd[, (paste0(v,"gap")) := YTime]; pd[is.na(get(v)), (paste0(v,"gap")) := NA]
  pd[, (paste0(v,"gap")) := c(NA, get(paste0(v,"gap"))[-.N]), by = Patient.ID]
  pd[, (paste0(v,"gap")) := na.locf(get(paste0(v,"gap")),na.rm = F),  by = Patient.ID]
  pd[, (paste0(v,"gap")) := YTime - get((paste0(v,"gap")))]
}

## check FVC
# tmp = pd[FVCgap ==0,]; tmp1 = pd[Patient.ID %in% tmp$Patient.ID,]
# View(tmp1[, c("Patient.ID", "Date","YTime", "FVC","med", "FVCgap", "FVCcf"), with=F])

## check trt type
# View(pd[, c("Patient.ID", "Date", "YTime", "type", "typecf", "typegap", "medswon", "medswoff"), with=F])

# Define indicator of the current med being an observed switch-on

pd[, medswon := 1*(!is.na(type) & type > 1 & typecf == 1)]; pd[is.na(medswon), medswon := 0]

# Define indicator of the current med being an observed switch-off

pd[, medswoff := 1*(!is.na(type) & type == 1 & typecf > 1)]; pd[is.na(medswoff), medswoff := 0]

# Extract variables and subdata from meds

var = c("Patient.ID", "YTime", "type", "T0",
        "age", "Sex", "ethnic", "Race", "Height","Weight",
        "ACA","SCL70","RNAPol",
        "FVCcf","FVCgap","DLCOcf","DLCOgap","RVSPcf","RVSPgap","EFcf",
        "EFgap", "mRSScf","mRSSgap","typecf", "typegap", "medswon", "medswoff")

pd = pd[med == 1 , var, with=F]

# Indicators of the PREVIOUS med being an observed switch-on/off

pd[, medswon := c(0, medswon[-.N]), by = Patient.ID]
pd[, medswoff := c(0, medswoff[-.N]), by = Patient.ID]

# View(pd[, c("Patient.ID","YTime", "type", "typecf","medswon","medswoff")])


# Carry forward and also backward for height and weight
for(var in c("Height", "Weight")){
  pd[, (paste0(var,"_M")) := 1*(is.na(get(var)))]
  pd[, var1 := na.locf(get(var),na.rm = F), by = Patient.ID] #cf
  pd[, var1 := na.locf(get(var),na.rm = F, fromLast = TRUE), by = Patient.ID] #cb
  pd[is.na(get(var)), (var) := var1]
  pd[is.na(get(var)), (var) := -99] # missing after carrying f/b coded as -99
}
pd[, var1 :=NULL]

# Missingness
vv = c("FVCcf","FVCgap","DLCOcf","DLCOgap","RVSPcf","RVSPgap","EFcf",
       "EFgap", "mRSScf","mRSSgap")
for(var in vv){
  if(var %in% grep("cf",vv,value = T)){#cf values
    pd[, (paste0(var,"_M")) := 1*(is.na(get(var)))]
    pd[is.na(get(var)), (var) := -99]
  } else {#gap values
    pd[is.na(get(var)), (var) := -99]
  }
}

# create indicators
## Sex
pd$Sex = as.numeric(pd$Sex)
## Race
pd[is.na(Race), Race := 9]; d[Race == 9, Race := 0]
tmp = class.ind(pd[, Race]); colnames(tmp) = paste0("Race_",colnames(tmp))
pd = cbind(pd, tmp[, -ncol(tmp)]); pd[, Race := NULL]
## ethnic
pd[is.na(ethnic), ethnic := 9]
pd[, ethnic_0 := 1*(ethnic == 0)]; pd[, ethnic_1 := 1*(ethnic == 1)]; pd[, ethnic := NULL] 
## previous outcome
tmp = class.ind(pd[, typecf]); colnames(tmp) = paste0("pretype_", colnames(tmp))
pd = cbind(pd, tmp[,-1]); pd[,typecf :=NULL]

pd = pd[T0 == 0, ] # discard the initial medical record for modeling; it will be given as the starting info for predictions
pd[, T0 :=NULL]

# NO MISSINGNESS IN pd!!!

save(pd, file = "pd.RData")

##############################################
load("pd.RData")


library(data.table)
library(nnet)
library(GcompBART)

# MPBART Alg3
KD = F; nuC = 1; Vm = diag(3); by = "1"; da = FALSE
nd = 9000; nb = 5000; nt = 100; p=length(unique(pd$type)) # number of outcome levels
seed = 1

nsub = nrow(pd); ntr = round(nsub*0.7)
set.seed(seed)
samp = sample(1:nsub, nsub, replace = F)

trid = samp[1:ntr]; trid = sort(trid)
trd = copy(pd[trid,])

teid = samp[(ntr+1):nsub]; teid = sort(teid)
ted = copy(pd[teid,])

RHS = colnames(pd)[c(2,4:42)]
fml = paste0("type ~ ",paste(RHS,collapse = " + "))

source("BartGcomp_19-11-22_Functions.R")
mod = model_bart(as.formula(fml), data = trd, type = "multinomial",
                 base = by,
                 Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                 Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                 correction = FALSE, Kindo = KD, do_alpha2_prior = da)
save(mod, file = "mod.RData")

load("mod.RData")
ted = ted[, c("type",RHS), with=F]; trd = trd[, c("type",RHS),with=F]
setDF(ted); setDF(trd)
colnames(ted) = colnames(trd) = c("type",mod$xcolnames)
acc = accuracyfun(bmpfit = mod,testd = ted,trdy = trd$type,tedy = ted$type)

save(acc, file = "acc.RData")

acc

BTr = mod$samp_y # train prediction
BTe = predict_bart(obj = mod, newdata = ted)$samp_y # test prediction
PmodeTr = apply(BTr, 1, ymode)
PmodeTe = apply(BTe, 1, ymode)

trd = copy(d[trid,]); ted = copy(d[teid,])
trd = trd[, c("Patient.ID", "YTime", "type"), with=F]
ted = ted[, c("Patient.ID", "YTime", "type"), with=F]
trd$pred = as.numeric(PmodeTr)
ted$pred = as.numeric(PmodeTe)

res = rbind(trd, ted); pd = pd[order(Patient.ID, YTime)]
save(res, file = "res.RData")