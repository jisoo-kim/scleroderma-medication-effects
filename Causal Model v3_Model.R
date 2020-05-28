# Fit PENCOMP
library(nlme)
library(RcppArmadillo);library(Rcpp)
sourceCpp("cfunctions.cpp")

Vn = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1", "YTime")

load("tmp.RData")

### Restrict to only those taking either no trt or MMF
d3 = d3[type0 %in% c(1,2) & type1 %in% c(1,2) & type2 %in% c(1,2),]
#4365 windows from 552 ppl

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


#################################################################
#################################################################
#################################################################

d4 = d3[!is.na(mRSS1q),]

V = d4[, c(Vn,"A0"), with = F]; V = as.matrix(V)

Y = d4[, mRSS1q]

A = d4[, A1]  


# (Va,A) for logistic regression
# (Xv,Y) for linear regression, Y ~ X + l
FitWeight = function(v, Va, A, Xv, Y){
  # take one col out for logistic regression, select by AIC
  res = rep(NA,ncol(Va))
  for(j in 1:ncol(Va)){
    tmp = glm(A ~ Va[,-j], family = "binomial")
    res[j] = tmp$aic
  }
  j = which.min(res); colnames(Va)[j]
  
  model1a = glm(A ~ Va[,-j], family="binomial") 
  #P = model1a$fitted.values
  dlt = model1a$coefficients
  
  # logit P for fitting Y ~ X
  Vt = cbind(1,Xv[,-j]); l = Vt %*% dlt #log(P/(1-P))
  l = matrix(l, ncol = 1)
  
  X = cbind(Xv,l)
  x = cbind(v, cbind(1, v[,-j]) %*% dlt)
  return(x %*% solve(t(X)%*%X)%*%t(X))
}

W1 = FitWeight(V[A==0,], V, A, V[A==1,], Y[A==1])
W0 = FitWeight(V[A==1,], V, A, V[A==0,], Y[A==0])

#library('plot.matrix')
#W = cbind(W0, matrix(0, nrow(W0), ncol(W1)))
#W = rbind(W,  cbind(matrix(0, nrow(W1), ncol(W0)), W1)  )
#plot(W)

plot(W1[1:60,1:200],cex.axis = 0.4, col=topo.colors, 
     xlab = "D1",ylab = "Y(1) for D0", main = "Weights of Observed Y's for Calculating Potential Outcomes")

plot(matrix(W1[1,1:31^2],31,31),cex.axis = 0.35,  col=topo.colors,
     xlab = "D1",ylab = "D1", main = "Weights on Observed Y|D1 for Potential Outcome Y(1) of one person in D0")

plot(matrix(W0[1,1:57^2],57,57),cex.axis = 0.35, col=topo.colors, breaks=seq(-0.01,0.02,0.005),
     xlab = "D0",ylab = "D0", main = "Weights on Observed Y|D0 for Potential Outcome Y(0) of one person in D1")
##################################################################


# take one col out for logistic regression, select by AIC
res = rep(NA,ncol(V))
for(j in 1:ncol(V)){
  tmp = glm(A ~ V[,-j], family = "binomial")
  res[j] = tmp$aic
}
j = which.min(res); colnames(V)[j]

model1a = glm(A ~ V[,-j], family="binomial") 
#P = model1a$fitted.values
dlt = model1a$coefficients

foo = function(dlt, v, Xv, Y){
  Y = matrix(Y, ncol = 1)
  
  # logit P for fitting Y ~ X; log(P/(1-P))
  l  = matrix(cbind(1,Xv[,-j]) %*% dlt, ncol = 1) 
  lx = matrix(cbind(1, v[,-j]) %*% dlt, ncol = 1)
  
  m1 = solve(t(Xv)%*%Xv) %*% t(Xv)
  alpha = m1 %*% Y
  tau   = m1 %*% l 
    
  denom = as.numeric(t(l) %*% (l - Xv %*% tau))
  Ya = v %*% alpha + (lx - v %*% tau) %*% t(l) %*% (Y - Xv %*% alpha) /denom
  
  return(Ya)
}


#W1 = FitWeight(V[A==0,], V, A, V[A==1,], Y[A==1])
#Y1 = W1%*%Y[A==1] 
#W0 = FitWeight(V[A==1,], V, A, V[A==0,], Y[A==0])
#Y0 = W0%*%Y[A==1] 
Y0 = foo(dlt, V[A==1,], V[A==0,], Y[A==0])
Y1 = foo(dlt, V[A==0,], V[A==1,], Y[A==1])

Vl_ord = rbind( cbind(1,V[A==0,-j]), cbind(1,V[A==1,-j]) )
A_ord = matrix(c(rep(0, sum(A==0)), rep(1, sum(A==1))),ncol=1)
dlt_inf = infl_logistic(Vl_ord,matrix(dlt,ncol=1), A_ord)

Y0k = inflY0(V[A==0,], cbind(1,V[A==0,-j]),  matrix(Y[A==0],ncol=1), 
             matrix(dlt,ncol=1), dlt_inf,
             V[A==1,], cbind(1,V[A==1,-j]))

Y1k = inflY1(V[A==1,], cbind(1,V[A==1,-j]),  matrix(Y[A==1],ncol=1), 
             matrix(dlt,ncol=1), dlt_inf,
             V[A==0,], cbind(1,V[A==0,-j]))

dltY0 = matrix(rep(Y0,ncol(Y0k)), ncol = ncol(Y0k)) - Y0k
dltY1 = matrix(rep(Y1,ncol(Y1k)), ncol = ncol(Y1k)) - Y1k

plot(dltY1[1:60,3258+(1:200)],cex.axis = 0.4, col=topo.colors,
     xlab = "D1",ylab = "Y(1) - Y(1)[-j] for D0", main = "Influence of Observed (A,X,Y)'s on Potential Outcomes")

plot(matrix(dltY1[1,1:65^2],65,65),cex.axis = 0.35,  col=topo.colors,
     xlab = "(D0,D1)",ylab = "(D0,D1)", main = "Influence of Observed (A,X,Y)'s on Potential Outcome Y(1) of one person in D0")

plot(matrix(dltY0[1,1:65^2],65,65),cex.axis = 0.35,  col=topo.colors,
     xlab = "(D0,D1)",ylab = "(D0,D1)", main = "Influence of Observed (A,X,Y)'s on Potential Outcome Y(0) of one person in D1")

############################################
a = rbind(dltY1,dltY0); K = 10
#a = dltY1; K = 10

gap = (max(a) - min(a))/K
# anchor interval centers at 0, (-gap/2, gap/2]
ln = ceiling((-gap/2 - min(a))/gap) # number of intervals to the left of the anchor
rn = ceiling((max(a) - gap/2)/gap) # number of intervals to the right of the anchor
bl = -gap/2 - ln*gap; bu = gap/2 + rn*gap
grid = bl + gap*(0:(ln+rn+1))

e = ceiling((c(a) - bl)/gap)
e[e==0] = 1; sort(unique(e))

# match neutral color max(rn,ln)+1 to the anchor interval ln+1
rbPalu <- colorRampPalette(c('yellow','red'))
rbPall <- colorRampPalette(c('blue','green'))
m = max(rn,ln); Col = c(rbPall(ln), NA,rbPalu(rn))
#rbPal <- colorRampPalette(c('blue','red'))
#Col = rbPal(2*m+1)
#if(rn > ln){
#  Col = Col[-(1:(rn - ln))]
#} else if(rn < ln){
#  Col = Col[1:(ln+1+rn)]
#}

n=nrow(a)
dat = data.frame(x = rep(1:ncol(a), n), y = rep(1:n,each = ncol(a)), col = c(t(matrix(e, nrow = nrow(a))[1:n,])))
sum(dat$col == ln + 1)

# only plot those not in the anchor interval (which centers at zero) 
#99.5% of entries are in the anchor interval
pd = dat[dat$col != ln+1, ] 

# Color Demo
plot(1:length(Col), 1:length(Col), col = Col) 

png("InfluencePlot.png", width = 20, height = 20, units = 'in', res = 600)
plot(pd$x, pd$y, cex = 0.04,pch = 15,col = Col[pd$col], ylim = c(0,n+1), xlim = c(1,ncol(a)))
dev.off()

tiff("Plot3.tiff", width = 20, height = 20, units = 'in', res = 1000)
plot(pd$x, pd$y, cex = 0.1,pch = 15,col = Col[pd$col], ylim = c(0,n+1), xlim = c(1,ncol(a)))
dev.off()


if(0){
  
  plot(matrix(dltY1[1,1:57^2],57,57),cex.axis = 0.35,  col=topo.colors,
       xlab = "D0",ylab = "D0", main = "Influence of Observed (A,X,Y)|D0 on Potential Outcome Y(1) of one person in D0")
  
  plot(matrix(dltY1[1,3258+(1:31^2)],31,31),cex.axis = 0.35,  col=topo.colors,
       xlab = "D1",ylab = "D1", main = "Influence of Observed (A,X,Y)|D1 on Potential Outcome Y(1) of one person in D0")
  
  
  
  X = cbind(V,l)
  model1b = lm(Y ~ -1+X)
  model1b$coefficients
  
  m1 = solve(t(V)%*%V); H = V %*% m1 %*% t(V)
  m2 = H %*% Y
  
  b1denom = as.numeric(t(l) %*% (l - H%*%l))
  b1 = m1 %*% t(V) %*% Y - m1 %*% t(V) %*%( l%*% t(l) %*% (Y - m2)) / b1denom 
  
  b2 = t(l) %*% (Y - m2)
  b2 = b2 / b1denom
  
  cbind(model1b$coefficients, c(b1,b2), solve(t(X)%*%X)%*%t(X)%*%Y)
  
}