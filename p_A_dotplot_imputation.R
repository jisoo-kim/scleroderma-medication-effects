# Posterior mode prediction of medication type starting from the 2nd medicationr record
load("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/schlr/Code/res.RData")

nid = 50
set.seed(342020)
set.seed(100)
idlist = sample(unique(res$Patient.ID), nid)
res1 = res[Patient.ID %in% idlist, ]
len =  res1[,.N, by = Patient.ID]$N; res1[,id:= rep(1:nid, len)]
res1[, tmax := YTime[.N], by = Patient.ID]
#pdfup = unique(res1[, c("Patient.ID", "id" , "tmax"), with=F])

m = matrix(c(1,2,3,3),nrow = 2,ncol = pncol,byrow = TRUE)
layout(mat = m,heights = c(0.9,0.1))
par(mai=c(0,0,0,0))
plot(res1$YTime, res1$id, xlab="", ylab="",  yaxt='n',, xlim = c(0, max(res1$YTime) + 5),
     main="", col = res1$type, pch = 19)
#text(pdfup$tmax + 3, pdfup$id, pdfup$Patient.ID)
par(mai=c(0,0,0,0))
plot(res1$YTime, res1$id, xlab="", ylab="",  yaxt='n',
     main="", col = res1$pred, pch = 19)
#par(mai=c(0,0,0,0))
par(mar=c(0,0,1.1,0))
plot.new()
cl = sort(unique(res1$type))
legend(x = "top",inset = 0, c("no trt","MMF","MMF and others","others"),
       text.col=cl,pch=19,col= cl,horiz=TRUE, bty="n", cex = 1.2)
