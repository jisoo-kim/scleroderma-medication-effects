name = "YX"

if(name == "YX"){
  filename = "C:\\Users\\Yizhen Xu\\Documents\\data and code\\Jisoo Data for Visualization 9-11-2018.xlsx"
}
if(name == "JK"){
  
}

##### Packages ##### 

library(plyr)
library(tidyr)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(gridExtra)
library(lubridate)
library(data.table)
library(VIM)
library(poLCA)
library(lme4)
library(splines)
library(MCMCglmm)
library(nnet) # for function class.ind(), create indicators for multinomial variables
library(zoo) # for function na.locf(), carry forward/backward in case of missingness 

# Data processing
## Combining ()
## output: generate dataframe dat with cols:
#("Patient.ID", "Date", "YTime", "age", "Sex", "ethnic", "Race", "Height", "Weight", 
# "FVC", "DLCO", "RVSP", "EF", "mRSS",
# "Pred", "MTX", "MMF", "CTX", "IVIG", "AZA", "Rituximab", "Tocilizumab", "HCQ", "TNF", "LEF",
# "pft", "mrss", "echo", "med")

source("Data Processing V1.R")

# p_A_heatmap
## diag_cmb = TRUE if diagonal includes medication combinations
##            FALSE if diagonal indicates one medication only
##            "all" drawn side by side

diag_cmb = "all"; source("p_A_heatmap.R")

# p_A_dotplot
## horizontal dots representing treatment trajectories for 100 patients (randomized)
### trt_type : treatment types
###            1 if no trt, MMF only, MMF+, others
###            2 if no trt, MMF only, Pred only, MMF+, others 

set.seed(1)
trt_type = 1; source("p_A_dotplot.R")

# i_Y_a
## Outcome trajectories and medication plot
### idlist: ID numbers (of four people)

idlist = c(16,3229,916,1642) # idlist = sample(unique(dat$Patient.ID),4)

source("i_Y_a.R")

# Y ~ A
load("Data Processing for Outcome Treatment Models V1.RData") # data: pd

#prior <- list(R = list(V = diag(2)/3, n = 2), G = list(
#  G1 = list(V = diag(2)/3, n = 2), G2 = list(V = diag(2)/3, n = 2)))
pd = pd[YTime >= 0 & YTime <= 40,]


V.R = diag(2); nu.R = 3 # uniform correlation
V.G = diag(c(1, 1, 0.005, 0.005)); nu.G = 4 # mean is V.G
b0 = c(matrix(c(0,0,0,0,0,2.17,2.17, 2.5,
       0,0,0,0,0,7.59,7.59, 8), byrow = T, nrow = 2))
V.B = diag(c(rep(1e+10,10), rep(c((1.65/1.96)^2,10.1^2),2 ) ,(2.6/1.96)^2, 8.3^2))

#literatures mostly 12 months
#FVC: MMF from SLSII; others by subjective esitmate
#mRSS:MMF from Le et al. , others taken from IVIG in Poelman et al.

nitt = 10000; burnin = 1000
nitt = 10; burnin = 3
# pr=TRUE, store posterior distribution of random effects
colnames(pd)[1] = "pid"

m1 = MCMCglmm(
  fixed = cbind(FVC, mRSS) ~ trait:(-1 + ns(YTime, knots = c(10, 30), Boundary.knots= c(0, 40)) +
    medtype_n1 + medtype_0 + medtype_2 + medtype_3 + medtype_4),
  random = ~ us(trait + trait:YTime):pid, pr = TRUE,
  rcov = ~ us(trait):units,
  prior = list(R = list(V = V.R, nu = nu.R),
               G = list(G1 = list(V = V.G, nu = nu.G))),
  family = c("gaussian", "gaussian"), nitt = nitt, burnin = burnin,
  data = pd)

m2 = MCMCglmm(
  fixed = cbind(FVC, mRSS) ~ trait:(-1 + ns(YTime, knots = c(10, 30), Boundary.knots= c(0, 40)) +
                                      medtype_n1 + medtype_0 + medtype_2 + medtype_3 + medtype_4),
  random = ~ us(trait + trait:YTime):pid, pr = TRUE,
  rcov = ~ us(trait):units,
  prior = list(R = list(V = V.R, nu = nu.R),
               G = list(G1 = list(V = V.G, nu = nu.G)),
               B = list(mu = b0, V = V.B)),
  family = c("gaussian", "gaussian"), nitt = nitt, burnin = burnin,
  data = pd)

summary(m1)
plot(m1$VCV)
plot(m1$Sol)
