# Fit PENCOMP
library(nlme)
### number of imputed datasets (#bootstrap)
numRun = 3 
minK = 15 # at least minK number of knots
V = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1")

load("tmp.RData")
### Restrict to only those taking either no trt or MMF
d3 = d3[type0 %in% c(1,2) & type1 %in% c(1,2) & type2 %in% c(1,2),]
#4608 windows from 572 ppl

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

backscale <- function(mvec, obs){
  obs = obs[!is.na(obs)]; obs = sort(obs); 
  mvec1 = mvec[!is.na(mvec)]; ind = which(!is.na(mvec))
  mvec[ind] = obs[ceiling(pnorm(mvec1)*length(obs))]
  return(mvec)}

# transform health outcomes; DO THIS AT EXACTLY THE DATA FOR MODELING!!!!!!!!!!!
for(v in c(paste0("FVC",0:2),paste0("mRSS",0:2),"cfcbFVC0", "cfcbmRSS0")){
  d3[, (paste0(v,"q")) := scalevec(get(v))]
}

# checks
#all(backscale(d3$mRSS1q, d3$mRSS1) == d3$mRSS1, na.rm=T)
#all(backscale(d3$mRSS0q, d3$mRSS0) == d3$mRSS0, na.rm=T)

############## without excluding nonoverlaps
estFinal_ps11_all=numeric(numRun)
sdFinal_ps11_all=numeric(numRun)
CIFinal_ps11_all=matrix(0, nrow=numRun,ncol=2)

estFinal_ps10_all=numeric(numRun)
sdFinal_ps10_all=numeric(numRun)
CIFinal_ps10_all=matrix(0, nrow=numRun,ncol=2)

estFinal_ps01_all=numeric(numRun)
sdFinal_ps01_all=numeric(numRun)
CIFinal_ps01_all=matrix(0, nrow=numRun,ncol=2)

############## Estimate Propensity Model and Predict Propensities (to be used as covariates)

### P(A1 | V, A0, Y0)
RHS_A1 = c(V, "A0") #Y0 coded as "cfcbFVC0"   "cfcbmRSS0" in V
fml_A1 = paste0("A1 ~ ", paste0(RHS_A1,collapse = "+"), collapse = "")
model2a=glm(as.formula(fml_A1), data=d3, family="binomial") 
d3[,temp2a :=(A1)*model2a$fitted.values + (1-A1)*(1-model2a$fitted.values)]
# Summary
summary(model2a)
# Accuracy
mean((model2a$fitted.values<0.5 & d3$A1==0) | (model2a$fitted.values>=0.5 & d3$A1==1))#0.95

if(0){ # 2 outcomes
  ### P(A2 | V, A0, A1, Y0, Y1) # ONLY USING 3153 SAMPLES DUE TO MISSINGNESS IN (FVC1, mRSS1)
  RHS_A2 = c(V, "Time0", "A0", "A1", "FVC1", "mRSS1")
  fml_A2 = paste0("A2 ~ ", paste0(RHS_A2,collapse = "+"), collapse = "")
  model2b=glm(as.formula(fml_A2), data=d3, family="binomial") 
  A2 = d3[!is.na(FVC1) & !is.na(mRSS1), A2]
  d3$temp2b = NA
  d3[!is.na(FVC1) & !is.na(mRSS1), temp2b:=(A2)*model2b$fitted.values + (1-A2)*(1-model2b$fitted.values)]
  # Summary
  summary(model2b)
  # Accuracy 
  mean((model2b$fitted.values<0.5 & A2==0) | (model2b$fitted.values>=0.5 & A2==1))#0.933
  if(0){# the high accuracy is pretty much predicting A2 to be A1...
    mean(d3[!is.na(FVC1) & !is.na(mRSS1), A2] == d3[!is.na(FVC1) & !is.na(mRSS1), A1]) # 0.934
  }
}

# 1 outcome 
Y1q = "mRSS1q"; Y1 = "mRSS1"
### P(A2 | V, A0, A1, Y0, Y1) # ONLY USING 3153 SAMPLES DUE TO MISSINGNESS IN (mRSS1)
RHS_A2 = c(V,  "A0", "A1", Y1)
fml_A2 = paste0("A2 ~ ", paste0(RHS_A2,collapse = "+"), collapse = "")
model2b=glm(as.formula(fml_A2), data=d3, family="binomial") 
A2 = d3[!is.na(mRSS1), A2]
d3$temp2b = NA; d3[, temp2b := as.numeric(temp2b)]
d3[!is.na(mRSS1), temp2b:= (A2)*model2b$fitted.values + (1-A2)*(1-model2b$fitted.values)]
# Summary
summary(model2b)
# Accuracy 
mean((model2b$fitted.values<0.5 & A2==0) | (model2b$fitted.values>=0.5 & A2==1))#0.956
if(0){# the high accuracy is pretty much predicting A2 to be A1...
  mean(d3[!is.na(mRSS1), A2] == d3[!is.na(mRSS1), A1]) # 0.956
}

### probability of observed treatment at first time point
d3[,pslogit := log(temp2a/(1-temp2a)) ]
### probability of observed treatment at second time point
d3$sw = NA;d3[, sw := as.numeric(sw)]
#d3[!is.na(FVC1) & !is.na(mRSS1), sw:=log((temp2a*temp2b)/(1-temp2a*temp2b))] 
d3[!is.na(mRSS1), sw:=log((temp2a*temp2b)/(1-temp2a*temp2b))] 

############## For each regime, Define Propensity Knots -> Fit Outcome Model -> Impute Outcome

### column for random intercept in outcome models
d3$id = 1

### P(Y1 | V, A0, Y0, pA1)  # avoid singularity by leave-1-out variable selection
RHS_Y1 = c(V, "A0","pslogit")
fml_Y1 = paste0(Y1q," ~ ", paste0(RHS_Y1,collapse = "+"), collapse = "")

############## A1 = 0

a1 = 0

###use equally spaced fixed knots assuming K knots
#d0=d3[A0 == 0 & A1==0 & !is.na(get(Y1)),] # 6038 rows
d0=d3[A1==a1 & !is.na(get(Y1)),] # 6241 rows
d0 = as.data.frame(d0)
K1=min(c(floor(0.25*dim(d0)[1]), minK))
space0=(max(d0$pslogit)-min(d0$pslogit))/(K1+1)
knots0=(min(d0$pslogit)+space0*(1:K1))
if(0){
  hist(d0$pslogit, breaks = 60); abline(v = knots0, col="red")
}

###assume a truncated linear basis
linear0=NULL
for (var in 1:K1) {
  temp=(d0$pslogit-knots0[var])
  temp2=(as.numeric(d0$pslogit<knots0[var]))*0 + (as.numeric(d0$pslogit>=knots0[var]))*temp
  linear0=cbind(linear0, temp2)
}

colnames(linear0)=paste("basis", 1:K1, sep="")

### Do leave-one-out model selection
RHS_Y1 = V
ee = vector('list', (length(RHS_Y1)+1));
fmll = vector('list', (length(RHS_Y1)+1));
RHSl = vector('list', (length(RHS_Y1)+1)); 
ll = rep(NA, (length(RHS_Y1)+1))
for(i in 1:(length(RHS_Y1)+1)){
  if(i == (length(RHS_Y1)+1)){
    RHSl[[i]] = c("pslogit", "A0",RHS_Y1)
  } else {
    RHSl[[i]] = c("pslogit","A0",RHS_Y1[-i])
  }
  
  fmll[[i]] = paste0(Y1q," ~ ", paste0(RHSl[[i]],collapse = "+"), collapse = "")
  ee[[i]] =  tryCatch(lme(as.formula(fmll[[i]]), random=list(id=pdIdent(~0+linear0)), data=d0),
                     error=function(e) "error")
  if(ee[[i]]!="error") ll[i] = ee[[i]]$logLik
}
RHS_Y1 = RHSl[[which.max(ll)]]; fml_Y1 = fmll[[which.max(ll)]]
pspp0 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear0)), data=d0)

#pspp0 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear0)), data=d0)
summary(pspp0)

### Impute Y1 | A1 = 0
newData0 = d3[, c(unique(c(RHS_Y1,RHS_A1)) ), with=F]; newData0$id = 1
predict2a=1-predict(model2a, newData0, type="response")
newData0$pslogit=log(predict2a/(1-predict2a))

linear0=NULL
for (var in 1:K1) {
  temp=(newData0$pslogit-knots0[var])
  temp2=(as.numeric(newData0$pslogit<knots0[var]))*0 + (as.numeric(newData0$pslogit>=knots0[var]))*temp
  linear0=cbind(linear0, temp2)
}

colnames(linear0)=paste("basis", 1:K1, sep="")

newData0t = newData0[, RHS_Y1, with=F]; newData0t$id = 1
newData0t = as.data.frame(newData0t)
### Predict without grouping
#predictedM0=predict(pspp0, newData0t, level = 0) + rnorm(dim(newData0)[1], 0, summary(pspp0)$sigma)
### Predict with grouping
newData0t = cbind(linear0, newData0t)
predictedM0=predict(pspp0, newData0t, level = 1) + rnorm(dim(newData0)[1], 0, summary(pspp0)$sigma)

newData0$mRSS1q=predictedM0
newData0[!is.na(d3$mRSS1q) & d3$A1 == 0, mRSS1q := d3[!is.na(mRSS1q) & A1 == 0, mRSS1q]]

############## A1 = 1

a1 = 1

###use equally spaced fixed knots assuming K knots

d1=d3[A1==a1 & !is.na(get(Y1)),] # 1257 rows
d1 = as.data.frame(d1)
K1=min(c(floor(0.25*dim(d1)[1]), minK))
space1=(max(d1$pslogit)-min(d1$pslogit))/(K1+1)
knots1=(min(d1$pslogit)+space1*(1:K1))
if(0){
  hist(d1$pslogit, breaks = 60); abline(v = knots1, col="red")
}

###assume a truncated linear basis
linear1=NULL
for (var in 1:K1) {
  temp=(d1$pslogit-knots1[var])
  temp2=(as.numeric(d1$pslogit<knots1[var]))*0 + (as.numeric(d1$pslogit>=knots1[var]))*temp
  linear1=cbind(linear1, temp2)
}

colnames(linear1)=paste("basis", 1:K1, sep="")

### Do leave-one-out model selection
RHS_Y1 = V
ee = vector('list', (length(RHS_Y1)+1));
fmll = vector('list', (length(RHS_Y1)+1));
RHSl = vector('list', (length(RHS_Y1)+1)); 
ll = rep(NA, (length(RHS_Y1)+1))
for(i in 1:(length(RHS_Y1)+1)){
  if(i == (length(RHS_Y1)+1)){
    RHSl[[i]] = c("pslogit", "A0",RHS_Y1)
  } else {
    RHSl[[i]] = c("pslogit","A0",RHS_Y1[-i])
  }
  
  fmll[[i]] = paste0(Y1q," ~ ", paste0(RHSl[[i]],collapse = "+"), collapse = "")
  ee[[i]] =  tryCatch(lme(as.formula(fmll[[i]]), random=list(id=pdIdent(~0+linear1)), data=d1),
                      error=function(e) "error")
  if(ee[[i]]!="error") ll[i] = ee[[i]]$logLik
}
RHS_Y1 = RHSl[[which.max(ll)]]; fml_Y1 = fmll[[which.max(ll)]]
pspp1 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear1)), data=d1)

summary(pspp1)

### Impute Y1 | A0 = 0, A1 = 0
newData1 = d3[, c(unique(c(RHS_Y1,RHS_A1)) ), with=F]; newData1$id = 1
predict2a=predict(model2a, newData0, type="response")
newData1$pslogit=log(predict2a/(1-predict2a))

linear1=NULL
for (var in 1:K1) {
  temp=(newData1$pslogit-knots1[var])
  temp2=(as.numeric(newData1$pslogit<knots1[var]))*0 + (as.numeric(newData1$pslogit>=knots1[var]))*temp
  linear1=cbind(linear1, temp2)
}

colnames(linear1)=paste("basis", 1:K1, sep="")

newData1t = newData1[, RHS_Y1, with=F]; newData1t$id = 1
newData1t = as.data.frame(newData1t)
### Predict without grouping
#predictedM1=predict(pspp1, newData1t, level = 0) + rnorm(dim(newData1)[1], 0, summary(pspp1)$sigma)
### Predict with grouping
newData1t = cbind(linear1, newData1t)
predictedM1=predict(pspp1, newData1t, level = 1) + rnorm(dim(newData1)[1], 0, summary(pspp1)$sigma)

newData1$mRSS1q=predictedM1
newData1[!is.na(d3$mRSS1q) & d3$A1 == 1, mRSS1q := d3[!is.na(mRSS1q) & A1 == 1, mRSS1q]]

######################
if(0){
  a = backscale(newData1$mRSS1q, d3$mRSS1)
  b = backscale(newData0$mRSS1q, d3$mRSS1)
  mean(a - b) #1.272135
  mean(a[d3$A0 ==1 & d3$A1==1]) - mean(b[d3$A0 ==1 & d3$A1==1])#-0.5869797
  mean(a[d3$A0 ==0 & d3$A1==1]) - mean(b[d3$A0 ==0 & d3$A1==1])#0.6870229
  mean(a[d3$A0 ==1 & d3$A1==0]) - mean(b[d3$A0 ==1 & d3$A1==0])#0.1020408
  mean(a[d3$A0 ==0 & d3$A1==0]) - mean(b[d3$A0 ==0 & d3$A1==0])#1.833818
  
  mean(a[d3$A0 ==1]) - mean(b[d3$A0 ==1])#-0.5217391
  mean(a[d3$A0 ==0]) - mean(b[d3$A0 ==0])#1.791772
  
  # compare ppl w (1,1)  to ppl w (0,0)
  mean(a[d3$A0 ==1 & d3$A1==1] - d3[A0 ==1 & A1==1,mRSS0], na.rm=T) - mean(b[d3$A0 ==0 & d3$A1==0]-d3[A0==0 & A1==0,mRSS0], na.rm=T)
  #-1.024538
  mean(d3[A0 ==1 & A1==1,mRSS1-mRSS0], na.rm=T) - mean(d3[A0==0 & A1==0,mRSS1-mRSS0], na.rm=T)
  #-1.066997
}
######################

Y2 = "mRSS2"; Y2q = "mRSS2q"

### P(Y2 | V, A0, Y0, Y1, pA2)  # avoid singularity by leave-1-out variable selection
### A0 is part of the baseline distribution, pA2 = P(A1,A2)
RHS_Y2 = c(V, "A0", Y1,"sw")
fml_Y2 = paste0(Y2q," ~ ", paste0(RHS_Y2,collapse = "+"), collapse = "")

############## A1 = 0, A2 = 0

a1 = 0; a2 = 0
a1 = 1; a2 = 1

###use equally spaced fixed knots assuming K knots

d00=d3[A1==a1 & A2==a2 & !is.na(get(Y1)) & !is.na(get(Y2)),] # 3245 rows
d00 = as.data.frame(d00)
K1=min(c(floor(0.25*dim(d00)[1]), minK))
space00=(max(d00$sw)-min(d00$sw))/(K1+1)
knots00=(min(d00$sw)+space00*(1:K1))
if(0){
  hist(d00$sw, breaks = 60); abline(v = knots00, col="red")
}

###assume a truncated linear basis
linear00=NULL
for (var in 1:K1) {
  temp=(d00$sw-knots00[var])
  temp2=(as.numeric(d00$sw<knots00[var]))*0 + (as.numeric(d00$sw>=knots00[var]))*temp
  linear00=cbind(linear00, temp2)
}

colnames(linear00)=paste("basis", 1:K1, sep="")

pspp_t2 = lme(as.formula(fml_Y2), random=list(id=pdIdent(~0+linear00)), data=d00)

### Impute Y2 | A1 = 0, A2 = 0, V, A0, Y1(A1=0)
newData00 = d3[, RHS_Y2, with=F]; newData00$id = 1; newData00$A1 = a1
newData00$mRSS1 = backscale(newData0$mRSS1q, d3$mRSS1) # check passed: sum(a[newData0$A1 ==0]!=b1[newData0$A1 ==0],na.rm=T)
### sw = P(A1 = a1, A2 = a2) accounts for A2 = a2
ta = predict(model2a, newData00, type="response")
predict2a= a1*ta + (1-a1)*(1-ta)
tb = predict(model2b, newData00, type="response")
predict2b= a2*tb + (1-a2)*(1-tb)
newData00$sw=log((predict2a * predict2b)/(1-predict2a * predict2b))

linear00=NULL
for (var in 1:K1) {
  temp=(newData00$sw-knots00[var])
  temp2=(as.numeric(newData00$sw<knots00[var]))*0 + (as.numeric(newData00$sw>=knots00[var]))*temp
  linear00=cbind(linear00, temp2)
}

colnames(linear00)=paste("basis", 1:K1, sep="")

newData00t = newData00[, c(RHS_Y2,"id"), with=F]
newData00t = as.data.frame(newData00t)
### Predict without grouping
#predictedM1=predict(pspp1, newData00t, level = 0) + rnorm(dim(newData00)[1], 0, summary(pspp1)$sigma)
### Predict with grouping
newData00t = cbind(linear00, newData00t)
impute_t2=predict(pspp_t2, newData00t,level=1) + rnorm(dim(newData00)[1], 0, summary(pspp_t2)$sigma)

##########################
#summary(pspp_t2)
if(a1 == 0 & a2 == 0){
  pspp00 = pspp_t2; imputed00 = impute_t2
}
if(a1 == 1 & a2 == 1){
  pspp11 = pspp_t2; imputed31 = impute_t2
}
summary(pspp00); summary(pspp11)
# chekc p values
a = summary(pspp0)$tTable; a[,5] = round(a[,5],3); a
b = summary(pspp1)$tTable; b[,5] = round(b[,5],3); b

a = summary(pspp00)$tTable; a[,5] = round(a[,5],3); a
b = summary(pspp11)$tTable; b[,5] = round(b[,5],3); b

y00 = imputed00
y00[which(d3$A1==0 & d3$A2==0 & !is.na(d3$mRSS2q))] = d3$mRSS2q[which(d3$A1==0 & d3$A2==0 & !is.na(d3$mRSS2q))]
y00 = backscale(y00,d3$mRSS2)
y11 = imputed31
y11[which(d3$A1==1 & d3$A2==1 & !is.na(d3$mRSS2q))] = d3$mRSS2q[which(d3$A1==1 & d3$A2==1 & !is.na(d3$mRSS2q))]

mean(y11-y00) #-5.115173
mean(y11[d3$A0 == 1]) - mean(y00[d3$A0 == 1]) #-9.62013
mean(y11[d3$A0 == 0]) - mean(y00[d3$A0 == 0]) #-3.81021

d3[, dlt := y11 - y00]
library(ggplot2)

#tmp1 = data.frame(Time0 = d3$Time0, dlt = dlt)
tmpACA1 = d3[ACA ==1, list(Time0, dlt) ]; tmpACA1$note = "ACA = 1"
tmpACA0 = d3[ACA ==0, list(Time0, dlt) ]; tmpACA0$note = "ACA = 0"

tmpRC1 = d3[Race_2 ==1, list(Time0, dlt) ]; tmpRC1$note = "Race_2 = 1"
tmpRC0 = d3[Race_2 ==0, list(Time0, dlt) ]; tmpRC0$note = "Race_2 = 0"

tmp000 = d3[A0==1 & A1 ==0 & A2 ==0, list(Time0, dlt) ]; tmp000$note = "A0 = 0, A1=0, A2=0"
tmp111 = d3[A0==1 & A1 ==0 & A2 ==0, list(Time0, dlt) ]; tmp111$note = "A1 = 1, A1=1, A2=1"
tmp00 = d3[A1 == 0 & A2 == 0, list(Time0, dlt) ]; tmp00$note = "A1=0, A2=0"
tmp11 = d3[A1 == 1 & A2 == 1, list(Time0, dlt) ]; tmp11$note = "A1=1, A2=1"
tmp0 = d3[A0 == 0 , list(Time0, dlt) ]; tmp0$note = "A0=0"
tmp1 = d3[A0 == 1, list(Time0, dlt) ]; tmp1$note = "A0=1"
tmpall = d3[, list(Time0, dlt) ]; tmpall$note = "all"
#tmp = rbind(tmpall, tmp0,tmp00,tmp1,tmp11,tmp000,tmp111)

#tmp = rbind(tmpall, tmp0,tmp00,tmp1,tmp11,tmpACA1,tmpACA0)

tmp = rbind(tmp0,tmp1,tmpRC0, tmpRC1,tmpACA1,tmpACA0)

tmp[, Time0 := floor(Time0/1)*1]
tmp$Time0 = as.factor(tmp$Time0)
ggplot(tmp, aes(x=Time0, y=dlt, fill=note)) + 
  geom_boxplot() +
  facet_wrap(~note, scale="free")

if(0){
  ##### Summary of treatment pattern #####
  # Freq table of (A0, A1, A2)
  tmp0 = tmp = d3[, .N, list(type0, type1,type2)]
  tmp[tmp==1] = "no trt"; tmp[tmp==2] = "MMF"; tmp[tmp==3] = "MMF+"; tmp[tmp==4] = "others"
  tmp$N = tmp0$N
  tmp[order(get(paste0("type",0:2)))]
  ## 852 (1,1,1), 3441 (0,0,0)
  
  # Plot Delta/Gap
  
  par(mfrow=c(1,2))
  # plot FVC
  tmp1 = d3[type1 == 2 & type2 == 2, list(Time0, (FVC2 - FVC0)/(Time2 - Time0)) ]
  tmp0 = d3[type1 == 1 & type2 == 1, list(Time0, (FVC2 - FVC0)/(Time2 - Time0)) ]
  # MMF red (2), no trt black (1)
  tmp1$col = 2; tmp0$col=4; tmp = rbind(tmp1, tmp0); tmp = tmp[!is.na(V2)]
  
  plot(tmp[,V2], x=tmp[,Time0], col = alpha(tmp$col, 0.3), pch = tmp$col+18,cex=0.2, ylab =expression(Delta~'FVC'[2]),xlab = "Start Day of 3-Visit Window", main = paste0("[1,1] (N=",sum(tmp$col==2),") vs [0,0] (N=",sum(tmp$col==4),")"))
  pd = tmp[col == 4]; pd = pd[order(Time0)]
  lo <- loess(V2 ~ Time0, data=pd, span=0.30)
  lines(pd$Time0,predict(lo), col='blue', lwd=3)
  pd = tmp[col == 2]; pd = pd[order(Time0)]
  lo <- loess(V2 ~ Time0, data=pd, span=0.30)
  lines(pd$Time,predict(lo), col='red', lwd=3)
  
  # plot mRSS
  tmp1 = d3[type1 == 2 & type2 == 2, list(Time0, (mRSS2 - mRSS0)/(Time2 - Time0)) ]
  tmp0 = d3[type1 == 1 & type2 == 1, list(Time0, (mRSS2 - mRSS0)/(Time2 - Time0)) ]
  
  # MMF red, no trt black
  tmp1$col = 2; tmp0$col=4; tmp = rbind(tmp1, tmp0); tmp = tmp[!is.na(V2)]
  
  plot(tmp[,V2], x=tmp[,Time0], col = alpha(tmp$col, 0.3), pch = tmp$col+18,cex=0.2, ylab =expression(Delta~'mRSS'[2]),xlab = "Start Day of 3-Visit Window", main = paste0("[1,1] (N=",sum(tmp$col==2),") vs [0,0] (N=",sum(tmp$col==4),")"))
  pd = tmp[col == 4]; pd = pd[order(Time0)]
  lo <- loess(V2 ~ Time0, data=pd, span=0.30)
  lines(pd$Time0,predict(lo), col='blue', lwd=3)
  pd = tmp[col == 2]; pd = pd[order(Time0)]
  lo <- loess(V2 ~ Time0, data=pd, span=0.30)
  lines(pd$Time,predict(lo), col='red', lwd=3)
  
}
