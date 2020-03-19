load("Data.RData") #load d

# Define length of follow-up
if(onset == 1){
  d$YTime = d$YTime1
} else {
  d$YTime = d$YTime0
}

d[, fup := YTime[.N], by = Patient.ID]
d[, sup := YTime[1], by = Patient.ID]

mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
pd = d[!is.na(med) & med == 1, c("Patient.ID", "YTime", "fup", "sup", mv),with=F]

pd$nonMMF = apply(pd[,mv[-3],with=F],1,function(x) 1*(sum(x)!=0) )# indicator of on trt (excluding MMF)

if(trt_type == 1){
  nl = 1:4
  levlist = c("no trt","MMF","MMF+","others")
  pd[, type := 1] # no trt
  pd[MMF == 1 & nonMMF == 0, type := 2] # 2 MMF only
  pd[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF+
  pd[MMF == 0 & nonMMF == 1 , type := 4] # 4 only others
}

if(trt_type == 2){
  nl = 1:5
  levlist = c("no trt","MMF","MMF+","others", "Pred")
  pd$nonPred = apply(pd[,mv[-1],with=F],1,function(x) 1*(sum(x)!=0) )
  pd[, type := 1]
  pd[MMF == 1 & nonMMF == 0, type := 2] # 2 MMF only
  pd[MMF == 1 & nonMMF == 1, type := 3] # 3 MMF and others
  pd[Pred == 1 & nonPred == 0, type := 5] # 5 Pred only
  pd[MMF == 0 & Pred == 0 & nonPred+nonMMF >= 1 , type := 4] # 4 only others (non Pred or MMF)
  
}

pd = pd[, c("Patient.ID", "YTime", "fup", "sup", "type"), with=F]

# nid number of patients into two columns of plots
nid = 100; pncol = 2
# Randomize nid number of patients and reassign IDs as 1 to nid
idlist = sample(unique(pd$Patient.ID), nid)
pdall = pd[Patient.ID %in% idlist, ]
len =  pdall[,.N, by = Patient.ID]$N; pdall[,id:= rep(1:nid, len)]

m = matrix(c(1,2,3,3),nrow = 2,ncol = pncol,byrow = TRUE)
layout(mat = m,heights = c(0.9,0.1))
for( i in 1:pncol ){
  if(i==1) par(mai=c(0,0,0,0.2))
  if(i==2) par(mai=c(0,0,0,0.2))
  nn =  floor(nid/pncol)-1
  ln = floor(nid/pncol)*(i-1)+1
  rg = ln:(ln+nn)
  pd = pdall[id %in% rg, ]; pdfup = unique(pd[, c("id", "fup", "sup", "Patient.ID"), with=F])
  
  plot(pd$YTime, pd$id, xlab="", ylab="",  yaxt='n', axes = F, xlim = c(0,max(pd$YTime,na.rm=T)+4),
       main="", col = pd$type, pch = 19, panel.first=for(k in 1:nrow(pdfup))  lines(c(pdfup[k,3], pdfup[k,2]), rep(pdfup[k,1],2), type = "l", lty = 1,col = rgb(red = 190/255, green =190/255, blue = 190/255)))
  #axis(2, at = unique(pd$id),labels =as.character(unique(pd$Patient.ID)), las = 2, cex.axis = 0.5)
  axis(1)
  abline(v = 0)
  text(pdfup$fup + 3, pdfup$id, pdfup$Patient.ID)
}
#par(mai=c(0,0,0,0))
par(mar=c(0,0,1.3,0))
plot.new()
legend(x = "top",inset = 0, levlist,
       text.col=nl,pch=19,col= nl,horiz=TRUE, bty="n", cex = 1.2)
