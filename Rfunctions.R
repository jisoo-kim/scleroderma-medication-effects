Blocks = function(A11,A12,A21,A22){
  rbind(cbind(A11,A12), cbind(A21,A22))
}

# rewrite in C++?
infl_dlt = function(X, dlt,Y){
  
  ph = invlgt(X%*%dlt); p1 = ph*(1-ph)
  sqV = diag(sqrt(c(p1)))
  Xt = sqV %*% X
  toswp = Blocks(t(Xt)%*%Xt, t(Xt), Xt, matrix(0,nrow(X),nrow(X))) # A22 = zero matrix
  H = multsweep(toswp, 1, ncol(Xt)) # - A21*inv(A11)*A12
  ind = (ncol(X)+1) : (ncol(H))
  H = -H[ind, ind] # A21*inv(A11)*A12
  
  #Wk = diag(N); Wk[k,k] = 0
  #resk = t(Xt) %*% Wk %*% Xt
  res = matrix(0,ncol(X),ncol(X))
  resl = vector("list", nrow(X))
  for(j in 1:nrow(X)){
    resl[[j]] = p1[j] * X[j,] %*% t(X[j,])
    res = res + resl[[j]]
  } 
  
  foo = function(k){
    resk = res - resl[[k]]
    pert = 1/(1-H[k,k]) * multsweep(resk, 1, ncol(res)) %*% X[k,] * (Y[k] - ph[k])
    return(pert)
  }
  
  pertl = lapply(1:N, foo)
  return(pertl)
}