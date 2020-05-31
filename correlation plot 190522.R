
# Description:
# This code produces empirical correlation plot of observations at fixed time points after 
# taking away the effect of variable-wise trend in time by subtracting the value estimated 
# by fixed effects

# Reference code: "180822 Emp Corr Matrix without Fixed Effects.R"
# modified - 191120

##### Packages #####

Packages <- c("ggplot2", "splines", "magrittr", "lme4", "dplyr")
lapply(Packages, library, character.only = TRUE)


##### Read in Data #####

# Setwd to current directory
setwd(paste0(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), "/Data & Preprocessing"))
load("minndat4.Rdata")


##### function computing nearest value for all patients #####

# assign dataset
dat <- minndat
pidname <- "Patient.ID"
# length(unique(c(t(dat %>% select(pidname))))) 
variablenames <- c("pFVCt", "pDLCOt", "EFt", "RVSPt", "mRSSt")
datename <- "YearsSinceOnset"

# remove row with missing time & order by time
dat <- dat %>% filter(complete.cases(dat[, datename])) %>% 
  dplyr::arrange(!!as.symbol(pidname), !!as.symbol(datename))

# param
yearvec <- seq(0, 10, by = 1)
mintol <- 1

# set patient id & time variables
commonid <- c(t(unique(dat %>% select(pidname))))
repcommonid <- rep(commonid, each = length(yearvec))
repyearvec <- rep(yearvec, length(commonid))


# data on grid
griddat <- matrix(NA, nrow = length(repcommonid), ncol = 2 + length(variablenames)) %>%
  set_colnames(c(pidname, "time", variablenames)) %>% as.data.frame()
griddat[, 1:2] <- cbind(repcommonid, repyearvec)


# function to put data on grid for all markers
for(k in 1:length(variablenames)){
  
  yearvar <- NULL
  
  for(i in 1:length(commonid)){
    
    piddat <- dat %>% filter(!!as.symbol(pidname) == commonid[i] & !is.na(!!as.symbol(variablenames[k])))
    
    pvardate <- piddat %>% select(datename)
    pvariable <- piddat %>% select(variablenames[k])
    
    pyearvar <- NULL

    for(j in 1:length(yearvec)){
      
      absval <- abs(pvardate - yearvec[j])
      
      if(nrow(absval) == 0){
        
        pyearvar[j]  <- NA # assign NA if all NA
      
        }else{
          
        if(min(absval) > mintol | is.na(min(absval) > mintol)){ 
          
          pyearvar[j] <- NA # assign NA if no observation within time +/- mintol
          
        }else{
          
          minindex <- which.min(c(t(absval)))
          pyearvar[j] <- c(t(pvariable))[minindex]}}
      }
    
    yearvar <- c(yearvar, pyearvar)}
  
  griddat[, k + 2] <- yearvar
  
}

View(griddat)



##### Fit linear mixed model and subtract estimated fixed effect from observed values #####
# only useful when using cross-validation or having different fixed effects estimated
# since subtracting estimated fixed effects from the same model will not affect
# the empirical covariance matrix

fitdata <- cbind(repcommonid, repyearvec) %>% set_colnames(c(pidname, "time"))
origdat <- dat %>% select(c(pidname, datename, variablenames)) %>% 
  set_colnames(c(pidname, "time", variablenames))


for(i in 1:length(variablenames)){

  fitframe <- fitdata
  
  # dataset wiht complete outcome variable
  compdat <- griddat %>% filter(complete.cases(griddat[, i + 2])) 
  regdat <- origdat %>% filter(complete.cases(origdat[, i + 2])) 
  colnames(compdat)[1] = colnames(regdat)[1] <- "pid"
  
  # fit mixed effects model on original data
  fit <- lmer(data = regdat, 
              eval(as.name(variablenames[i])) ~ ns(time, df = 2) + (1 + time|pid), 
              control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

  # #fit using just the fixed effects
  # fitted.df <- data.frame(compdat[, 1:2], predict(fit, re.form = NA, newdata = compdat)) %>%
  #   set_colnames(c(pidname, "time", variablenames[i]))
  
  #fit using random effects
  fitted.df <- data.frame(compdat[, 1:2],
                          predict(fit, re.form = NULL, newdata = compdat)) %>%
    set_colnames(c(pidname, "time", variablenames[i]))
  
  fitdata <- merge(fitframe, fitted.df, by = c(pidname, "time"), all = T)
  
}


View(fitdata)

# subtract fitted values from observed
diffdat <- data.frame(griddat[, 1:2], (griddat[, 3:(2 + length(variablenames))] - fitdata[, 3:(2 + length(variablenames))]))


##### Empirical Variance Covariance #####

# choose data to plot
data <- griddat

# remove rows with all NAs
widedat <- reshape(data, idvar = pidname, timevar = "time", direction = "wide")[, -1]

# calculate 
cordat <- cor(widedat, method = "pearson", use = "pairwise.complete.obs")

# reorder by column name
cnames <- colnames(cordat); orderedcnames <- NULL

for(i in 1:length(variablenames)){
  orderedcnames <- c(orderedcnames,  cnames[startsWith(cnames, variablenames[i])])
}


# order by marker
cordat <- cordat[orderedcnames, orderedcnames]
melted_mat <- reshape2::melt(cordat)

# plot
ggplot(data = melted_mat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white", aes(y = reorder(Var2, desc(Var2)))) + labs(x = "", y = "") +
  scale_fill_gradient2(low = "blue", high = "red", "Corr", mid = "white", midpoint = 0) +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 10), legend.position = "none") + 
  scale_x_discrete(labels= rep(0:10, 5)) + 
  scale_y_discrete(labels= rep(0:10, 5)) +
  labs(x = "pFVC               pDLCO               EF                RVSP              mRSS", 
       y = "pFVC               pDLCO               EF                RVSP              mRSS")


# plot with corr labels
ggplot(data = melted_mat, aes(x=Var1, y=Var2, fill=value, label = round(value, 2))) + 
  geom_tile(color = "white", aes(y = reorder(Var2, desc(Var2)))) + labs(x = "", y = "") +
  scale_fill_gradient2(low = "blue", high = "red", "Corr", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text=element_text(size=8)) + 
  geom_text(colour = "black", size = 1.5, aes(y = reorder(Var2, desc(Var2))))













