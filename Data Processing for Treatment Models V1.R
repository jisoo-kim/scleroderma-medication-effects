#library(MNP)
#https://github.com/yizhenxu/GcompBART

# Delete records from people without date of onset
d = dat[!is.na(YTime),] #11661 records from 1155 patients

# data: 21753 records, 1163 patients (48 people missing DateOnRP)

# Need to delete negative time since onset
#summary(d$YTime)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-10.55    4.41    9.71   12.96   18.05   71.46 

d = d[YTime >= 0,] # reduces to 21576 records

# carry forward (and also backward for height and weight)
for(var in c("Height", "Weight", "FVC", "DLCO", "RVSP", "EF", "mRSS")){
  d[, (paste0(var,"_M")) := 1*(is.na(get(var)))]
  d[, var1 := na.locf(get(var),na.rm = F), by = Patient.ID]
  if(var %in% c("Height", "Weight")){
    d[, var1 := na.locf(get(var),na.rm = F, fromLast = TRUE), by = Patient.ID]
  }
  d[is.na(get(var)), (var) := var1]
  d[is.na(get(var)), (var) := -99] # missing after carrying f/b coded as -99
}
d[, var1 :=NULL]

# previous outcome
d[, pretype := c(NA, type[-.N]), by=Patient.ID]
d[is.na(pretype), pretype := 0]
tmp = class.ind(d[, pretype]); colnames(tmp) = paste0("pretype_", colnames(tmp))
d = cbind(d, tmp[,-1]); d[,pretype :=NULL]

# create indicators
# Race
d[is.na(Race), Race := 9]; d[Race == 9, Race := 0]
tmp = class.ind(d[, Race]); colnames(tmp) = paste0("Race_",colnames(tmp))
d = cbind(d, tmp[, -ncol(tmp)]); d[, Race := NULL]
# ethnic
d[is.na(ethnic), ethnic := 9]
d[, ethnic_0 := 1*(ethnic == 0)]; d[, ethnic_1 := 1*(ethnic == 1)]; d[, ethnic := NULL] 
# Sex
d$Sex = as.numeric(d$Sex)

# keep those with non-negative time since onset (with onset date info non-missing)
d = d[!is.na(YTime) & YTime >= 0,]
# keep only medication records and delete 1 person with no gender info
d = d[!is.na(type) & !is.na(Sex),]

save(d, file = "d.RData")

# MPBART Alg3
KD = F; nuC = 1; Vm = diag(3); by = "1"; da = FALSE
nd = 3000; nb = 5000; nt = 100; p=length(unique(d$type)) # number of outcome levels
seed = 1

nsub = nrow(d); ntr = round(nsub*0.7)
set.seed(seed)
samp = sample(1:nsub, nsub, replace = F)

trid = samp[1:ntr]; trid = sort(trid)
trd = copy(d[trid,])

teid = samp[(ntr+1):nsub]; teid = sort(teid)
ted = copy(d[teid,])

RHS = colnames(d)[c(3,5:32)]
fml = paste0("type ~ ",paste(RHS,collapse = " + "))

source("BartGcomp_19-11-22_Functions.R")
mod = model_bart(as.formula(fml), data = trd, type = "multinomial",
                 base = by,
                 Prior = Prior_mult(p = p, nu_C = nuC, Vmat = Vm, ntree = nt),
                 Mcmc = Mcmc_mult(p = p,nb = nb, nd = nd),
                 correction = FALSE, Kindo = KD, do_alpha2_prior = da)
save(mod, file = "mod.RData")

ted = ted[, c("type",RHS), with=F]; trd = trd[, c("type",RHS),with=F]
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

pd = rbind(trd, ted); pd = pd[order(Patient.ID, YTime)]
save(pd, file = "pd.RData")