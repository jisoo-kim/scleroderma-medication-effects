library(RcppArmadillo);library(Rcpp)
sourceCpp("cfunctions.cpp")

### Function sweep and adjust pass matrix as reference!!!
X = matrix(c(1,1,1,1,1,1,
             1:3,1:3,
             1,1,1,-1,-1,-1), ncol=3)
Y = c(1,3,3,2,2,1)

XY = cbind(X,Y); colnames(XY) = NULL
M = t(XY) %*% XY

M1 = M
for(k in 1:ncol(X)){
  M = sweep(M, k)
  #sweep_web(M,k)
}
M
multsweep(M1, 1,ncol(X))

lm(Y~-1+X)
M[1:3,4]
sum((Y - X%*%M[1:3,4])^2) == M[4,4]

M[1:ncol(X),1:ncol(X)] #solve(t(X)%*%X)
M[1:ncol(X), ncol(X)+1] #solve(t(X)%*%X) %*% t(X) %*% Y

### Influence function of logistic regression
set.seed(1)
N = 10000
X = matrix(c( rnorm(3*N)), ncol = 3)
logit = function(p) return(log(p/(1-p)))
invlgt = function(xb) return(1/(1+exp(-xb)))
p = invlgt(X %*% matrix(c(2,0.5,5),ncol = 1))
A = rbinom(N, 1,p)
Y = rnorm(N, 0.2*X[,1] - A*X[,2] + X[,3] -X[,2]*X[,3] , 0.4)
dat = data.frame(X1 = X[,1], X2 = X[,2], X3 = X[,3], A = A, Y = Y)

dat = dat[order(dat$A),]

mod = glm(A ~ -1 + X1 + X2 , data = dat, family = "binomial")
dlt = as.numeric(mod$coefficients)
dlt_inf = infl_logistic(X[,-3],matrix(dlt,ncol=1), matrix(A,ncol=1))

phat = fitted(mod)
dat$pscore = log(phat/(1-phat))

modY = lm(Y ~ -1+ X1 + X2 + X3 + pscore, data = dat[A==0,])
yhat = fitted(modY, dat[A==1])
#apply(dlt_inf, 1, sum); dlt
a = infl(as.matrix(dat[dat$A==0,1:3]), as.matrix(dat[dat$A==0,1:2]), as.matrix(dat[dat$A==0,5]),
         matrix(dlt,ncol=1), dlt_inf,
         as.matrix(dat[dat$A==1,1:3]), as.matrix(dat[dat$A==1,1:2]))









tmp = proc.time()
a = infl_logistic(X,matrix(dlt,ncol=1), matrix(Y,ncol=1))
proc.time() - tmp
tmp = proc.time()
b = infl_dlt(X,dlt,Y)
proc.time() - tmp

if(0){
  mod = lm(Y ~ -1 + X1 + X2 + X3 + A, data = dat)
  #mod = glm(Y ~ -1 + X1 + X2 + X3, data = dat, family = "binomial")
  dlt = as.numeric(mod$coefficients)
  
  k = 9999 # take kth row out
  modk = lm(Y ~ -1 + X1 + X2 + X3 , data = dat[-k,])
  #modk = glm(Y ~ -1 + X1 + X2 + X3, data = dat[-k,], family = "binomial")
  dlt - as.numeric(modk$coefficients)
}


# multsweep()
tmp = matrix(c(1,1,1,3,2,-1,2,3,5),3); 
# the following three are equivalent:
solve(tmp)
sweep(sweep(sweep(tmp,1),2),3)
multsweep(tmp, 1,ncol(tmp)) # inverse

x = tmp
for(i in 1:3){
  x = sweep(x,i)
}

test(tmp,3)
