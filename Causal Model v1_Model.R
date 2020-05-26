dat = copy(d4)

# Prepare data for modeling
## Create lag 1 outcomes
mv = c("bmit","FVCt","mRSSt","type")
for(v in mv){
    dat[, (paste0(v, "prev")) := c(NA, get(v)[-.N]), by = Patient.ID]
}

## Process baseline characteristics
### Sex
dat$Sex = as.numeric(dat$Sex)
### Race (1 white 2 black, use Race_1 and Race_2)
dat[is.na(Race), Race := 9]; d[Race == 9, Race := 0]
tmp = class.ind(dat[, Race]); colnames(tmp) = paste0("Race_",colnames(tmp))
dat = cbind(dat, tmp[, -ncol(tmp)]); dat[, Race := NULL]
## ethnic (use ethnic_0 and ethnic_1; indicator of latino)
dat[is.na(ethnic), ethnic := 9]
dat[, ethnic_0 := 1*(ethnic == 0)]; dat[, ethnic_1 := 1*(ethnic == 1)]; dat[, ethnic := NULL] 

# Model Y1 | Y0, X0, A0

V.R = diag(2); nu.R = 3 # uniform correlation
V.G = diag(c(1, 1, 0.005, 0.005)); nu.G = 4 # mean is V.G
b0 = c(matrix(c(0,0,0,0,0,2.17,2.17, 2.5,
                0,0,0,0,0,7.59,7.59, 8), byrow = T, nrow = 2))
V.B = diag(c(rep(1e+10,10), rep(c((1.65/1.96)^2,10.1^2),2 ) ,(2.6/1.96)^2, 8.3^2))

nitt = 10000; burnin = 1000
nitt = 10; burnin = 3
# pr=TRUE, store posterior distribution of random effects
colnames(dat)[1] = "pid"

m1 = MCMCglmm(
  fixed = cbind(FVCt, mRSSt) ~ trait:(-1 + ns(YTime, knots = c(10, 30), Boundary.knots= c(0, 40)) +
                                      age0 + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1 +
                                      medtype_n1 + medtype_0 + medtype_2 + medtype_3 + medtype_4),
  random = ~ us(trait + trait:YTime):pid, pr = TRUE,
  rcov = ~ us(trait):units,
  prior = list(R = list(V = V.R, nu = nu.R),
               G = list(G1 = list(V = V.G, nu = nu.G))),
  family = c("gaussian", "gaussian"), nitt = nitt, burnin = burnin,
  data = dat)

#modA = MCMCglmm (family=c("binomial"))
