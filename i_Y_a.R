load("Data Processing V1.RData") #load d

### Plot of (pFVC, -mRSS, EF) above and (Pred, MTX, MMF, CTX)
mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
dat = d[, c("Patient.ID","YTime","FVC", "DLCO", "RVSP", "EF", "mRSS", mv), with=F]
colnames(dat)[7] = "nmRSS"; dat[,7] = -dat[,7]
dat = dat[order(Patient.ID, YTime),]

#length(unique(dat$Patient.ID)) # 1211 unique patients

dat$nonMMF = apply(dat[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)

if(trt_type == 1){
  nl = 1:4
  levlist = c("no trt","MMF","MMF+","others")
  dat[, type := 1] # no trt
  dat[MMF == 1 & nonMMF == 0, type := 2] # 2 MMF only
  dat[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
  dat[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
}

if(trt_type == 2){
  nl = 1:5
  levlist = c("no trt","MMF","MMF+","others", "Pred")
  dat$nonPred = apply(dat[,mv[-1],with=F],1,function(x) 1*(sum(x)!=0) )
  dat[, type := 1]
  dat[MMF == 1 & nonMMF == 0, type := 2] # 2 MMF only
  dat[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF and others
  dat[Pred == 1 & nonPred == 0, type := 5] # 5 Pred only
  dat[MMF == 0 & Pred == 0 & nonPred+nonMMF >= 1 , type := 4] # 4 only others (non Pred or MMF)
  
}

# ind is 1 if FVC or mRSS  or EF does not have any measurement
dat[, ind := 1*(sum(!is.na(FVC[1:.N]))*sum(!is.na(nmRSS[1:.N]))*sum(!is.na(EF[1:.N])) ==0), by = Patient.ID]
dat = dat[ind==0,]

layout(matrix(c(rep(1,3),2,rep(3,3),4,rep(5,3),6,rep(7,3),8), nrow = 8, ncol = 2))

for(id in idlist){
  pd = dat[Patient.ID == id,]
  ind1 = which(!is.na(pd$FVC)); timeFVC = pd$YTime[ind1]; FVC = pd$FVC[ind1]
  ind2 = which(!is.na(pd$nmRSS)); timeRSS = pd$YTime[ind2]; RSS = pd$nmRSS[ind2]
  ind3 = which(!is.na(pd$EF)); timeEF = pd$YTime[ind3]; EF = pd$EF[ind3]
  
  rgt = range(c(timeFVC,timeRSS))
  
  ## add extra space to right margin of plot within frame
  par(mar=c(0,4,3,4))
  
  ## Plot first set of data and draw its axis
  plot(timeFVC, FVC, pch=16, axes=FALSE, ylim = range(FVC), xlab="", ylab="", xlim = rgt,
       type="b",col="black", main=paste0("ID ", id))
  axis(2, ylim=range(FVC),col="black",las=1)  ## las=1 makes horizontal labels
  mtext("FVC",side=2,line=2.5)
  box()
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot the second plot and put axis scale on right
  plot(timeRSS, RSS, pch=15,  xlab="", ylab="", ylim=range(RSS), xlim = rgt,
       axes=FALSE, type="b", col="red")
  ## a little farther out (line=4) to make room for labels
  mtext("-mRSS",side=4,col="red",line=2.5) 
  axis(4, ylim=range(RSS), col="red",col.axis="red",las=1)
  
  ## Allow a third plot on the same graph
  par(new=TRUE)
  
  ## Plot the third plot and put axis scale on right
  plot(timeEF, EF, pch=15,  xlab="", ylab="", ylim=range(EF), xlim = rgt,
       axes=FALSE, type="b", col="blue")
  
  ## Add Legend
  #legend("topright",legend=c("pFVC","mRSS"),
  #       text.col=c("black","red"),pch=c(16,15),col=c("black","red"),cex=1.5)
  
  par(mar=c(4,4,0.5,4))
  ## Medication
  
  ind = which(!is.na(pd$type)); t = pd$YTime[ind]; M = pd$type[ind]
  
  plot(t, rep(0.5, length(t)), xlab="", ylab="",xlim = rgt, ylim = c(0,1),
       main="", col = M, pch = 19, axes=FALSE)
  
  axis(1,pretty(rgt,10))
  mtext("Time (Years since onset)",side=1,col="black",line=2.5)  
}