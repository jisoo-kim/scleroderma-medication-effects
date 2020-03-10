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
# define natural spline basis of Year time
tmp = ns(pd$YTime,df = 2)
pd$nt1 = tmp[,1]; pd$nt2 = tmp[,2]

#prior <- list(R = list(V = diag(2)/3, n = 2), G = list(
#  G1 = list(V = diag(2)/3, n = 2), G2 = list(V = diag(2)/3, n = 2)))
m1 = MCMCglmm(
  fixed = cbind(FVC, RVSP, mRSS) ~ - 1 + trait:nt1 + trait:nt2 +
    trait:medtype_0 + trait:medtype_1 + trait:medtype_2 + trait:medtype_3,
  random = ~ idh(trait):nt1 + idh(trait):nt2,
  rcov = ~ us(trait):units, 
  family = c("gaussian", "gaussian", "gaussian"), nitt = 60, burnin = 10,
  thin = 1, data = pd)

plot(m1$VCV)
plot(m1$Sol)