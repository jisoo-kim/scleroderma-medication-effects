
##### Packages ##### 
Packages <- c("plyr","tidyr","ggplot2","data.table", "MCMCglmm","splines")
lapply(Packages, library, character.only = TRUE)

##### Data Read-in ##### 
load("tmp_JK1.Rdata")


##### Fit mcmcglmm to FVC, DLCO, MMF as outcome variables #####

Vn = c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1", "YTime")

V =  c("Time0", "MMFdos0", "age", "Sex", "ACA", "SCL70", "RNAPol", "cfcbFVC0q", "cfcbmRSS0q", "Race_1", "Race_2", "ethnic_0", "ethnic_1")

(FVC1q, mRSS1q, A1) ~ (FVC0, mRSS0, V)


V = "MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1"
Z = Time0 + YTime

dat <- dat %>% mutate(Time = YTime + Time1 - Time0, 
                      A1dup = A1,
                      FVC1qdup = FVC1q,
                      mRSS1qdup = mRSS1q)

Y1q ~  A1 + V (+ Zb)

A1  ~  Y0 + A0 + V (+ b)

Y0 = FVC0 + mRSS0
FVC0q + mRSS0q
cfcbFVC0 + cfcbmRSS0
cfcbFVC0q + cfcbmRSS0q

MCMCglmm(
  fixed = cbind(FVC1q, mRSS1q, A1) , data = dat)

colnames(dat)
str(dat)

burnin <- 500; nitt <- 5000

sum(is.na(dat$FVC0q))
dat$FVC0q

dat.comp %>% nrow
dat.comp <- dat[!is.na(FVC0q) & !is.na(mRSS0q), ]

# age -> age of onset
# time - ns(time, 3)

# conditional distributions of FVC1q, mRSS1q|A0, A1
# simulated values -> causal characteristics, plots counterfactual quantities
# categorize mmfdos0 - truncation at 5

# A1 expit(xbeta*) each iteration of the chain - check mixing
# evolving propensity score: p(A1 = 1|...knowledge of doctor, model)



fit1 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex  + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + (ACA + SCL70 + RNAPol)*ns(Time, 3)) +
                   at.level(trait, 2):(A1dup + (ACA + SCL70 + RNAPol)*ns(Time, 3)) +
                   at.level(trait, 3):(FVC0q + mRSS0q + A0 + ACA + SCL70 + RNAPol),
                 
                 random = ~ us(trait):Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = burnin, nitt = nitt, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp)



dat.comp$MMFdos0 %>% hist
test <- predict.MCMCglmm(fit1)
dat.comp$Race_1 %>% sum
dat.comp$Race_2 %>% sum



datn <- dat.comp %>% nrow

pred <- test[(datn*2 + 1):nrow(test)]


summary(fit1)
plot(fit1)


# using carry backward & foward imputation of FVC0 and mRSS0

dat.comp2 <- dat[!is.na(cfcbFVC0q) & !is.na(cfcbmRSS0q), ]

fit2 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup) +
                   at.level(trait, 2):(A1dup) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~  us(at.level(trait, 1):(1 + YTime) +
                                at.level(trait, 2):(1 + YTime) +
                                at.level(trait, 3):(1)):Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = burnin, nitt = nitt, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp2)

summary(fit2)
plot(fit2)

#dat.comp2$A1dup = dat.comp2$A1

##### Data Read-in ##### 
load("tmp_JK1.Rdata")
dat$A1dup = dat$A1
dat.comp3 <- dat[!is.na(cfcbFVC0q) & !is.na(cfcbmRSS0q), ]
fit3 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 2):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~ us(at.level(trait, 1):(1 + Time0) +
                                 at.level(trait, 2):(1 + Time0) +
                                 at.level(trait, 3):(1)):Patient.ID,
                   
                 rcov = ~ us(trait):units,
                 
                 burnin = 1000, nitt = 10000, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp3)

save(fit3, file = "Joint Model v1_res.RData")

load("Joint Model v1_res.RData")
summary(fit3)

par("mar")
par(mar=c(1,1,3,1))
plot(fit3)

###########################################
### try shared random effect

dat$A1dup = dat$A1
dat.comp4 <- dat[!is.na(cfcbFVC0q) & !is.na(cfcbmRSS0q), ]
fit4 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 2):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~ us(at.level(trait, 1):(1 + Time0) +
                                 at.level(trait, 2):(1 + Time0) +
                                 at.level(trait, 3):(1)):Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = 100, nitt = 200, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp4)

sf4 = summary(fit4)
sf4$Gcovariances[,1]


fit5 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 2):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~ us(trait:(1) + 
                                 at.level(trait, 1):(Time0) +
                                 at.level(trait, 2):(Time0)):Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = 100, nitt = 200, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp4)

sf5 = summary(fit5)
cbind(sf4$Gcovariances[,1],sf5$Gcovariances[,1])
rbind(names(sf4$Gcovariances[,1]),names(sf5$Gcovariances[,1])) #fit4 equiv fit5

#same rand intercept across outcome
fit6 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 2):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~ trait:Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = 100, nitt = 200, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp4)
sf6 = summary(fit6)
sf6$Gcovariances
names(sf6$Gcovariances[,1])


#Doesn't work
fit6 <- MCMCglmm(cbind(FVC1q, mRSS1q, A1) ~ trait:(MMFdos0 + age + Sex + ACA + SCL70 + RNAPol + Race_1 + Race_2 + ethnic_0 + ethnic_1) +
                   
                   at.level(trait, 1):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 2):(A1dup + ns(Time0, knots = c(10, 30), Boundary.knots= c(0, 40)) ) +
                   at.level(trait, 3):(cfcbFVC0q + cfcbmRSS0q + A0),
                 
                 random = ~ us(at.level(trait,1)+at.level(trait,2)):Patient.ID + 
                   at.level(trait,3):Patient.ID,
                 
                 rcov = ~ us(trait):units,
                 
                 burnin = 100, nitt = 200, pr = T,
                 family = c("gaussian", "gaussian", "categorical"), 
                 
                 data = dat.comp4)
sf6 = summary(fit6)
sf6$Gcovariances
names(sf6$Gcovariances[,1])

