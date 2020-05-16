numLoop = 100

library(data.table)
load("tmp.RData")

### Restrict to only those taking either no trt or MMF
d3 = d3[type0 %in% c(1,2) & type1 %in% c(1,2) & type2 %in% c(1,2),]
#4608 windows from 572 ppl

scalevec <- function(mvec){as.data.frame(qqnorm(mvec))[, 1]}

# transform health outcomes; DO THIS AT EXACTLY THE DATA FOR MODELING!!!!!!!!!!!
for(v in c(paste0("FVC",0:2),paste0("mRSS",0:2),"cfcbFVC0", "cfcbmRSS0")){
  d3[, (paste0(v,"q")) := scalevec(get(v))]
}

d3$id2=1:dim(d3)[1]

idl = dtl = matrix(NA, nrow(d3), numLoop)

for(h in 1:numLoop){
  destfile = paste0("Bootstrap_",h,".txt")
  if(file.exists(destfile)){
    res = read.table(destfile, header = TRUE)
    idl[,h] = res[,1]
    dtl[,h] = res[,2]
  }
}

save(idl, dtl, file = "btstrp_res.RData")
load("btstrp_res.RData")


##########################
### K : Number of bootstrap
### Wd: Variance within, Wd = sum(SE^2)/K, where SE[k] is the standard error of mu[k] from the kth bootstrap
### SE[k]: because mu[k] = sum(mu[k,i])/N, SE[k] = Var(mu[k,1:N])/N
### Bd: Variance bretween, Bd = sum((mu[k] - mean(mu))^2)/(K-1)
### Total Variance = Wd + (1 + 1/K)*Bd

summ = function(x){
  return(c(mean(x), var(x)/length(x)))
}

tmp = apply(dtl,2,summ)

mu = mean(tmp[1,])
Wd = mean(tmp[2,]^2)
Bd = sum((tmp[1,] - mu)^2)/(numLoop - 1)
TV = Wd + (1+1/numLoop)*Bd

round(c(mu, Wd, Bd ,TV),4)
# 0.3055 0.0001 0.8101 0.8183
########################## summarize by Time0 (subgroup)

summ_bygroup = function(d3, numLoop, X){
  mul = sel = matrix(NA, length(unique(floor(d3[, get(X)]))), numLoop) 
  val = sort(unique(floor(d3[, get(X)]))) # get unique values of X
  
  for(h in 1:numLoop){
    tmp = data.frame(X = d3[idl[,h], get(X)], dlt = dtl[,h], A0 = d3[idl[,h], A0]) # add subgroup info here
    setDT(tmp)
    tmp[, X := floor(X)]
    muh = tmp[A0==0, mean(dlt), by = X]            # add subgroup info here
    seh = tmp[A0==0, var(dlt)/length(dlt), by = X] # add subgroup info here
    muh = muh[order(X)]; seh = seh[order(X)]
    ind = which(val %in% muh$X)
    mul[ind,h] = muh$V1; sel[ind,h] = seh$V1
  }
  
  Wd = apply(sel, 1, function(u) mean(u^2,na.rm=T))
  Bd = apply(mul, 1, function(u) var(u,na.rm=T))
  mu = apply(mul, 1, function(u) mean(u,na.rm=T))
  TV = Wd + (1+1/numLoop)*Bd
  
  return(list(val = val, mu = mu, TV = TV))
}

p1 = summ_bygroup(d3, numLoop, "Time0")
p2 = summ_bygroup(d3, numLoop, "YTime")
par(mfrow=c(2,2))
plot(p1$val, p1$mu, xlab = "Years since Enrollment", ylab = "Mean Delta",pch = 18, cex = 0.5, type="l", ylim = range(c(p1$mu,0))); abline(h=0,col="red")
plot(p1$val, p1$TV, xlab = "Years since Enrollment", ylab = "Total Variance", pch = 18, cex = 0.5, type="l")
maxt = 25
plot(p2$val[1:maxt], p2$mu[1:maxt], xlab = "Years since Onset", ylab = "Mean Delta",pch = 18, cex = 0.5, type="l", ylim = range(c(p2$mu[1:maxt],0))); abline(h=0,col="red")
plot(p2$val[1:maxt], p2$TV[1:maxt], xlab = "Years since Onset", ylab = "Total Variance", pch = 18, cex = 0.5, type="l")

##########################

tmp = simdatFit[, c("Time0", "dlt"), with=F]; tmp = tmp[order(Time0),]
plot(tmp$Time0, tmp$dlt, pch = 18, cex = 0.5, col = alpha("grey",0.3))
lo <- loess(dlt ~ Time0, data=tmp, span=0.30)
lines(tmp$Time0,predict(lo), col='red', lwd=1)
abline(h = c(0,-5), col = "blue")
#################################################################################################
##FIND THE OVERLAPPING REGIONS AT FIRST TIME POINT #########################
####looking of probability of getting treated for both control and treated groups

#simdatFit = copy(d3)
probTreat=predict(model2a, simdatFit, type="response")
summary(probTreat)

overlapTreat=c(max(min(probTreat[which(simdatFit$A1==1)]), min(probTreat[which(simdatFit$A1==0)])),
               min(max(probTreat[which(simdatFit$A1==1)]), max(probTreat[which(simdatFit$A1==0)])))
# 0.00138165 0.98462309

# the input of the following function MUST be a numeric list
plot(density(probTreat[simdatFit$A1==0], from=0, 
             to=1), lty=1, lwd=2, col="black", xlab="Propensity Score", main="")
lines(density(probTreat[simdatFit$A1==1],from=0, to=1), lty=2, lwd=2, col="red")
legend("topright", c("control","treated"), cex=1.2, lty=1:2, col=c("black", "red"))

##############################
####looking of probability of getting treated for both control and treated groups
probControl=1-predict(model2a, simdatFit, type="response")
overlapControl=c(max(min(probControl[which(simdatFit$A1==1)]), min(probControl[which(simdatFit$A1==0)])),
                 min(max(probControl[which(simdatFit$A1==1)]), max(probControl[which(simdatFit$A1==0)])))
#0.01537691 0.99861835

plot(density(probControl[simdatFit$A1==0], from=0, 
             to=1), lty=1, lwd=2, col="black", xlab="Propensity Score", main="")
lines(density(probControl[simdatFit$A1==1],from=0, to=1),
      lty=2, lwd=2, col="red")
legend("topleft", c("control","treated"), cex=1.2, lty=1:2, col=c("black", "red"))


simdatFit$includedT1=(probTreat >= overlapTreat[1] &  probTreat <= overlapTreat[2]) *
  (probControl >= overlapControl[1] &  probControl <= overlapControl[2])
sum(simdatFit$includedT1)  

