if(0){
#!/bin/bash
#SBATCH -J btstrp
#SBATCH --time 12:00:00
#SBATCH --array=[1-100]%10
#SBATCH --output=job_%A_%a.out
#SBATCH --error=job_%A_%a.err
#SBATCH --mail-type=all
#SBATCH --mail-user=yizhen_xu@brown.edu

#Rscript Causal_Model_v2_Bootstrap.R
}



rm(list=ls()) 

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid = Sys.getenv('SLURM_ARRAY_TASK_ID')

h = as.numeric(slurm_arrayid) # index of simulation, h in 1:numLoop

# Fit PENCOMP
library(nlme)
library(data.table)

### number of imputed datasets (#bootstrap)

minK = 15 # at least minK number of knots

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

# define interactions with years since onset YTime
intvlist = c("A0","MMFdos0", "Sex", "ACA", "SCL70", "RNAPol", "Race_1", "Race_2", "ethnic_0", "ethnic_1")
for( v in intvlist){
  intv = paste0(v, "_YTime")
  d3[, (intv) := get(v)*YTime]
}

V = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "Race_1", "Race_2", "ethnic_0", "ethnic_1", "YTime", paste0(intvlist, "_YTime"))
Y0  = c("cfcbFVC0q", "cfcbmRSS0q")

# checks
#all(backscale(d3$mRSS1q, d3$mRSS1) == d3$mRSS1, na.rm=T)
#all(backscale(d3$mRSS0q, d3$mRSS0) == d3$mRSS0, na.rm=T)


d3$id2=1:dim(d3)[1]

num00=d3$id2[d3$A1==0 & d3$A2==0]
num01=d3$id2[d3$A1==0 & d3$A2==1]
num10=d3$id2[d3$A1==1 & d3$A2==0]
num11=d3$id2[d3$A1==1 & d3$A2==1]

#numLoop = 3

#for(h in 1:numLoop){
  set.seed(h)
  
  simdatFit=d3[c(sample(x=num00, size=length(num00),replace=T),sample(x=num01, size=length(num01),replace=T),
                 sample(x=num10, size=length(num10),replace=T),sample(x=num11, size=length(num11),replace=T)),]
  simdatFit = simdatFit[order(Patient.ID, Time0),]
  
  ############## Estimate Propensity Model and Predict Propensities (to be used as covariates)
  
  ### P(A1 | V, A0, Y0)
  RHS_A1 = c(V, "A0", Y0) 
  fml_A1 = paste0("A1 ~ ", paste0(RHS_A1,collapse = "+"), collapse = "")
  model2a=glm(as.formula(fml_A1), data=simdatFit, family="binomial") 
  simdatFit[,temp2a :=(A1)*model2a$fitted.values + (1-A1)*(1-model2a$fitted.values)]
  # Summary
  #summary(model2a)
  # Accuracy
  #mean((model2a$fitted.values<0.5 & simdatFit$A1==0) | (model2a$fitted.values>=0.5 & simdatFit$A1==1))#0.95
  
  # 1 outcome 
  Y1q = "mRSS1q"; Y1 = "mRSS1"
  
  ### P(A2 | V, A0, A1, Y0, Y1) # ONLY USING PART OF SAMPLES DUE TO MISSINGNESS IN (mRSS1)
  RHS_A2 = c(V, "A0", "A1", Y0, Y1)
  fml_A2 = paste0("A2 ~ ", paste0(RHS_A2,collapse = "+"), collapse = "")
  model2b=glm(as.formula(fml_A2), data=simdatFit, family="binomial") 
  A2 = simdatFit[!is.na(mRSS1), A2]
  simdatFit$temp2b = NA; simdatFit[, temp2b := as.numeric(temp2b)]
  simdatFit[!is.na(mRSS1), temp2b:= (A2)*model2b$fitted.values + (1-A2)*(1-model2b$fitted.values)]
  # Summary
  #summary(model2b)
  # Accuracy 
  #mean((model2b$fitted.values<0.5 & A2==0) | (model2b$fitted.values>=0.5 & A2==1))#0.956
  
  ### probability of observed treatment at first time point
  simdatFit[,pslogit := log(temp2a/(1-temp2a)) ]
  ### probability of observed treatment at second time point
  simdatFit$sw = NA;simdatFit[, sw := as.numeric(sw)]
  #simdatFit[!is.na(FVC1) & !is.na(mRSS1), sw:=log((temp2a*temp2b)/(1-temp2a*temp2b))] 
  simdatFit[!is.na(mRSS1), sw:=log((temp2a*temp2b)/(1-temp2a*temp2b))] 
  
  ############## For each regime, Define Propensity Knots -> Fit Outcome Model -> Impute Outcome
  
  ### column for random intercept in outcome models
  simdatFit$id = 1
  
  ### P(Y1 | V, A0, Y0, pA1)  # avoid singularity by leave-1-out variable selection
  RHS_Y1 = c(V, "A0", Y0, "pslogit")
  fml_Y1 = paste0(Y1q," ~ ", paste0(RHS_Y1,collapse = "+"), collapse = "")
  
  ############## A1 = 0
  
  a1 = 0
  
  ###use equally spaced fixed knots assuming K knots
  d0=simdatFit[A1==a1 & !is.na(get(Y1)),] # 6241 rows
  d0 = as.data.frame(d0)
  K1=min(c(floor(0.25*dim(d0)[1]), minK))
  space0=(max(d0$pslogit)-min(d0$pslogit))/(K1+1)
  knots0=(min(d0$pslogit)+space0*(1:K1))
  
  ###assume a truncated linear basis
  linear0=NULL
  for (var in 1:K1) {
    temp=(d0$pslogit-knots0[var])
    temp2=(as.numeric(d0$pslogit<knots0[var]))*0 + (as.numeric(d0$pslogit>=knots0[var]))*temp
    linear0=cbind(linear0, temp2)
  }
  
  colnames(linear0)=paste("basis", 1:K1, sep="")
  
  ### Do leave-one-out model selection
  ee = vector('list', (length(V)+1));
  fmll = vector('list', (length(V)+1));
  RHSl = vector('list', (length(V)+1)); 
  ll = rep(NA, (length(V)+1))
  #for(i in 1:(length(V)+1)){
  for(i in 1:(length(V))){
    if(i == (length(V)+1)){
      RHSl[[i]] = c("pslogit", "A0", Y0, V)
    } else {
      RHSl[[i]] = c("pslogit","A0",Y0, V[-i])
    }
    
    fmll[[i]] = paste0(Y1q," ~ ", paste0(RHSl[[i]],collapse = "+"), collapse = "")
    ee[[i]] =  tryCatch(lme(as.formula(fmll[[i]]), random=list(id=pdIdent(~0+linear0)), data=d0),
                        error=function(e) "error")
    if(ee[[i]]!="error") ll[i] = ee[[i]]$logLik
  }
  RHS_Y1 = RHSl[[which.max(ll)]]; fml_Y1 = fmll[[which.max(ll)]]
  pspp0 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear0)), data=d0)
  
  #pspp0 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear0)), data=d0)
  #summary(pspp0)
  
  ### Impute Y1 | A1 = 0
  newData0 = simdatFit[, c(unique(c(RHS_Y1,RHS_A1)) ), with=F]; newData0$id = 1
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
  newData0[!is.na(simdatFit$mRSS1q) & simdatFit$A1 == 0, mRSS1q := simdatFit[!is.na(mRSS1q) & A1 == 0, mRSS1q]]
  
  ############## A1 = 1
  
  a1 = 1
  
  ###use equally spaced fixed knots assuming K knots
  
  d1=simdatFit[A1==a1 & !is.na(get(Y1)),] # 1257 rows
  d1 = as.data.frame(d1)
  K1=min(c(floor(0.25*dim(d1)[1]), minK))
  space1=(max(d1$pslogit)-min(d1$pslogit))/(K1+1)
  knots1=(min(d1$pslogit)+space1*(1:K1))
  
  ###assume a truncated linear basis
  linear1=NULL
  for (var in 1:K1) {
    temp=(d1$pslogit-knots1[var])
    temp2=(as.numeric(d1$pslogit<knots1[var]))*0 + (as.numeric(d1$pslogit>=knots1[var]))*temp
    linear1=cbind(linear1, temp2)
  }
  
  colnames(linear1)=paste("basis", 1:K1, sep="")
  
  ### Do leave-one-out model selection
  ee = vector('list', (length(V)+1));
  fmll = vector('list', (length(V)+1));
  RHSl = vector('list', (length(V)+1)); 
  ll = rep(NA, (length(V)+1))
  #for(i in 1:(length(V)+1)){
  for(i in 1:(length(V))){
    if(i == (length(V)+1)){
      RHSl[[i]] = c("pslogit", "A0", Y0, V)
    } else {
      RHSl[[i]] = c("pslogit","A0",Y0, V[-i])
    }
    
    fmll[[i]] = paste0(Y1q," ~ ", paste0(RHSl[[i]],collapse = "+"), collapse = "")
    ee[[i]] =  tryCatch(lme(as.formula(fmll[[i]]), random=list(id=pdIdent(~0+linear1)), data=d1),
                        error=function(e) "error")
    if(ee[[i]]!="error") ll[i] = ee[[i]]$logLik
  }
  RHS_Y1 = RHSl[[which.max(ll)]]; fml_Y1 = fmll[[which.max(ll)]]
  pspp1 = lme(as.formula(fml_Y1), random=list(id=pdIdent(~0+linear1)), data=d1)
  #library(lme4)
  #paste0(colnames(linear1), collapse = ",")
  #d1tmp = cbind(d1,linear1)
  #tmp = lmer(mRSS1q ~ pslogit+A0+cfcbFVC0q+cfcbmRSS0q+Time0+MMFdos0+age+Sex+ACA+SCL70+RNAPol+Race_1+Race_2+ethnic_0+ethnic_1+YTime+A0_YTime+MMFdos0_YTime+Sex_YTime+ACA_YTime+RNAPol_YTime+Race_1_YTime+Race_2_YTime+ethnic_0_YTime+ethnic_1_YTime + 
  #             (1|basis1,basis2,basis3,basis4,basis5,basis6,basis7,basis8,basis9,basis10,basis11,basis12,basis13,basis14,basis15), data = d1tmp)
  summary(pspp1)
  
  ### Impute Y1 | A0 = 0, A1 = 0
  newData1 = simdatFit[, c(unique(c(RHS_Y1,RHS_A1)) ), with=F]; newData1$id = 1
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
  ### Predict with grouping
  newData1t = cbind(linear1, newData1t)
  predictedM1=predict(pspp1, newData1t, level = 1) + rnorm(dim(newData1)[1], 0, summary(pspp1)$sigma)
  
  newData1$mRSS1q=predictedM1
  newData1[!is.na(simdatFit$mRSS1q) & simdatFit$A1 == 1, mRSS1q := simdatFit[!is.na(mRSS1q) & A1 == 1, mRSS1q]]
  
  ######################
  
  Y2 = "mRSS2"; Y2q = "mRSS2q"
  
  ### P(Y2 | V, A0, Y0, Y1, pA2)  # avoid singularity by leave-1-out variable selection
  ### A0 is part of the baseline distribution, pA2 = P(A1,A2)
  RHS_Y2 = c(V, "A0", Y0, Y1,"sw")
  fml_Y2 = paste0(Y2q," ~ ", paste0(RHS_Y2,collapse = "+"), collapse = "")
  
  ############## 
  
  #A1 = 0; A2 = 0
  #A1 = 1; A2 = 1
  
  a1a2 = list(c(0,0), c(1,1))
  for(k in 1:2){
    a1 = a1a2[[k]][1]; a2 = a1a2[[k]][2]
    
    
    ###use equally spaced fixed knots assuming K knots
    
    d00=simdatFit[A1==a1 & A2==a2 & !is.na(get(Y1)) & !is.na(get(Y2)),] # 3245 rows
    d00 = as.data.frame(d00)
    K1=min(c(floor(0.25*dim(d00)[1]), minK))
    space00=(max(d00$sw)-min(d00$sw))/(K1+1)
    knots00=(min(d00$sw)+space00*(1:K1))
    
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
    newData00 = simdatFit[, RHS_Y2, with=F]; newData00$id = 1; newData00$A1 = a1
    newData00$mRSS1 = backscale(newData0$mRSS1q, simdatFit$mRSS1) # check passed: sum(a[newData0$A1 ==0]!=b1[newData0$A1 ==0],na.rm=T)
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
      pspp11 = pspp_t2; imputed11 = impute_t2
    }
  }
  
  y00 = imputed00
  y00[which(simdatFit$A1==0 & simdatFit$A2==0 & !is.na(simdatFit$mRSS2q))] = simdatFit$mRSS2q[which(simdatFit$A1==0 & simdatFit$A2==0 & !is.na(simdatFit$mRSS2q))]
  y00 = backscale(y00,d3$mRSS2)
  y11 = imputed11
  y11[which(simdatFit$A1==1 & simdatFit$A2==1 & !is.na(simdatFit$mRSS2q))] = simdatFit$mRSS2q[which(simdatFit$A1==1 & simdatFit$A2==1 & !is.na(simdatFit$mRSS2q))]
  y11 = backscale(y11,d3$mRSS2)
  simdatFit[, dlt := y11 - y00]
  
  res = simdatFit[, c("id2", "dlt"),with=F]
  
  write.table(res, file =paste0("Bootstrap_",h,".txt"),
              append = FALSE, sep = " ", dec = ".", row.names = F, col.names = T)
  
#}

