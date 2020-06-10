
# read the chunk starting from line 355


##### Packages & Read in Data #####

Packages <- c("MCMCglmm", "splines", "dplyr", "plyr", "magrittr", "lme4", 
              "ggplot2", "tidyr", "matrixcalc", "abind", "gtable")

lapply(Packages, library, character.only = TRUE)

# Setwd to current directory
setwd("~/Desktop/Git/Data & Preprocessing")

# load data
load("minndat0.Rdata")

##### Filter Data #####

# assign dataset
dat <- minndat %>% filter(YearsSinceOnset >= 0 & YearsSinceOnset <= 40)
pidname <- "Patient.ID"
pidall <- unique(c(t(dat %>% select(pidname))))

variablenames <- c("pFVCt", "pDLCOt", "EFt") 

baselinevarnames <- c("AgeOnset", "race", "Type", "Sex", "ACA", "RNAPol", "SCL70")
datename <- "YearsSinceOnset"
timename <- "YearsSinceFirstObs"

# remove row with missing time & order by time
regdat <- dat %>% filter(complete.cases(!!as.symbol(datename))) %>%
  dplyr::arrange(!!as.symbol(pidname), !!as.symbol(datename)) %>%
  select(c(pidname, datename, timename, variablenames, baselinevarnames)) %>%
  set_colnames(c("pid", "time", "ref.time", variablenames, baselinevarnames))


# select patients with more than th observations
npid <- ddply(regdat, .(pid), summarize,
              nfvc = sum(!is.na(pFVCt)), ndlco = sum(!is.na(pDLCOt)),
              nef =  sum(!is.na(EFt)))
nrow(npid)

# for m = 3
th <- 2
filterdat <- npid %>% filter(nfvc >= th & ndlco >= th & nef >= th)
filterid <- filterdat$pid

# for m = 5

#th <- 3
#filterdat <- npid %>% filter(nfvc >= th & ndlco >= th & nef >= th & nrvsp >= th & nmrss >= th)
#filterid <- filterdat$pid

regdat.d <- regdat %>% filter(pid %in% filterid)

# create dataset with dummy variables
dummyvarnames <- c("race_AA", "Type_Diffuse", "Type_Sine", "Sex_Male", "ACApos", "RNAPolpos", "SCL70pos")
regdat.d <- fastDummies::dummy_cols(regdat.d, remove_most_frequent_dummy = T) %>%
  select(-c(baselinevarnames[-1])) %>% set_colnames(c("pid", "time", "ref.time", variablenames, "AgeOnset", dummyvarnames))

allid <- regdat.d$pid %>% unique


# create spline variables

# most recent obs for time 0
# observation going into the past k = 1, 2, 3
# knot at -10, 3

# last: most recent observation time since the time of onset
# btime: year since onset - last

regdat.d <- regdat.d %>% merge(ddply(regdat.d, .(pid), summarize, last = last(time)), by = "pid") %>%
  mutate(btime = time - last, 
         sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
         sp3yr = ifelse(btime + 3 > 0, btime + 3, 0))



##### model #####

# Set prior for combined model when m = 3

nmeas <- length(variablenames)

p <- 12

Mode <- diag(c(1, 1, 1, rep(0.005, 6), rep(0.00005, 3)))
nu <- 12

phi <- Mode * (nu + p + 1)
V <- phi/nu

# prior for R
p.R <- 3
Mode.R <- diag(3)

nu.R <- 3

phi.R <- Mode.R * (nu.R + p.R + 1)
V.R <- phi.R/nu.R


# set prior for separated model
nu.sep <- nu - (nmeas - 1) * (p/nmeas)
phi.sep <- diag(c(phi[1, 1], phi[4, 4], phi[7, 7], phi[10, 10]))
V.sep <- phi.sep/nu.sep

# prior for R
nu.R.sep <- nu.R - (nmeas - 1)

phi.R.sep <- phi.R[1, 1]
V.R.sep <- phi.R.sep/nu.R.sep


##### fit model #####

nitt <- 15000; burnin <- 2000

combmcmc <- MCMCglmm(as.formula(paste("cbind(", paste(variablenames, collapse = ","),
                                      ") ~ -1 + trait + trait:(ns(time, knots = c(10, 30), Boundary.knots= c(0, 40)) +
                                      AgeOnset + race_AA + Type_Diffuse + Type_Sine + Sex_Male + ACApos + RNAPolpos + SCL70pos)")),
                     #random = ~ us(trait + trait:time + trait:sp1yr):pid, pr = TRUE,
                     random = ~ us(trait + trait:btime + trait:sp10yr + trait:sp3yr):pid, pr = TRUE,
                     rcov = ~ us(trait):units,
                     prior = list(R = list(V = V.R, nu = nu.R),
                                  G = list(G1 = list(V = V, nu = nu))),
                     family = rep("gaussian", nmeas),
                     nitt = nitt, burnin = burnin, data = regdat.d)
summary(combmcmc)
save(combmcmc, file = "combmcmc_fvcdlcoef.Rdata")

load("combmcmc_3meas.Rdata")
load("combmcmc_fvcdlcoef.Rdata")



##### fitted values #####

# Combined model

# predicted values

# predval <- predict.MCMCglmm(combmcmc, marginal = NULL, re.form = NULL, interval = "confidence", level = 0.95)
# 
# re.form = NULL: using random effects structure as in the model
# martinal = NULL: use random effect estimates to predict
# 
# nfit <- dim(predval)[1]/3
# regdat.d <- regdat.d %>% mutate(fvcfit = predval[1:nfit], 
#                                 dlcofit = predval[(nfit + 1):(2*nfit)],
#                                 effit = predval[(2*nfit + 1):(3*nfit)])
# 
# fitdat <- cbind(predval[1:nfit, ], 
#                 dlcofit = predval[(nfit + 1):(2*nfit), ],
#                 effit = predval[(2*nfit + 1):(3*nfit), ]) %>%
#   set_colnames(paste0(rep(c("fvc", "dlco", "ef"), each = 3), 
#                       rep(c("fit", "lwr", "upr"), 3)))
# 
# regdat.d <- cbind(regdat.d, fitdat)



##### data for estimating regression parameters #####

regdat.d <- data.frame(cbind(ns(regdat.d$time, knots = c(10, 30), Boundary.knots= c(0, 40))),
                       regdat.d, int = 1)
colnames(regdat.d)[1:3] <- c("ns1", "ns2", "ns3")

FEcolnames <- c("int", "ns1", "ns2", "ns3", "AgeOnset", "race_AA", "Type_Diffuse", "Type_Sine", "Sex_Male",
                "ACApos", "RNAPolpos", "SCL70pos")

REcolnames <- c("int", "btime", "sp10yr", "sp3yr")

datfvc <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[1])))
datdlco <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[2])))
datef <- regdat.d %>% filter(complete.cases(!!as.symbol(variablenames[3])))



##### Varcov matrix from mcmcglmm ##### 

covdat <- apply(combmcmc$VCV, 2, mean)
ordervec <- c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12)

tmat <- covdat[endsWith(names(covdat), "pid")]
G <- matrix(tmat, sqrt(length(tmat)))
G2 <- G[ordervec, ordervec]

tmatr <- covdat[endsWith(names(covdat), "units")]
R2 <- matrix(tmatr, sqrt(length(tmatr)))

##### function to predict coefficients & variance parameters #####
est_coeff_Yvar <- function(G2, R2){
  
  # Var Y
  
  Vinv = VarYi = Zi = Xi = Yi = precB2 = xvinvy = Sigma2i = nXid1 = nXid2 = nXid3 = 
    nZid1 = nZid2 = nZid3 = nYi1 = nYi2 = nYi3 <- list()
  
  nobsmat <- matrix(NA, nrow = length(allid), ncol = nmeas)
  
  for(n in 1:length(allid)){
    
    # Z
    
    nZid1[[n]] <- nzid1 <- datfvc %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
    nZid2[[n]] <- nzid2 <- datdlco %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
    nZid3[[n]] <- nzid3 <- datef %>% filter(pid == allid[n]) %>% select(REcolnames) %>% as.matrix
    
    datlist <- list(nzid1, nzid2, nzid3)
    Zi[[n]] <- Z <- Reduce(direct.sum, datlist)
    
    nobsmat[n, ] <- c(nrow(nzid1), nrow(nzid2), nrow(nzid3))
    
    # X
    
    nXid1[[n]] <- x1 <- datfvc %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
    nXid2[[n]] <- x2 <- datdlco %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
    nXid3[[n]] <- x3 <- datef %>% filter(pid == allid[n]) %>% select(FEcolnames) %>% as.matrix
    
    
    Xi[[n]] <-  Reduce(direct.sum, list(x1, x2, x3))
    
    # Y 
    
    nYi1[[n]] <- y1 <- datfvc %>% filter(pid == allid[n]) %>% select(variablenames[1]) %>% as.matrix %>% set_colnames("y")
    nYi2[[n]] <- y2 <- datdlco %>% filter(pid == allid[n]) %>% select(variablenames[2]) %>% as.matrix %>% set_colnames("y")
    nYi3[[n]] <- y3 <- datef %>% filter(pid == allid[n]) %>% select(variablenames[3]) %>% as.matrix %>% set_colnames("y")
    
    
    Yi[[n]] <-  rbind(y1, y2, y3)
    
    
    # Sigma2
    
    Sigma2 <- NULL
    
    for (i in 1:nmeas){
      
      Cmat <- NULL
      
      for (j in 1:nmeas){
        
        time1 <- datlist[[i]][, 2]; time2 <- datlist[[j]][, 2]
        
        if (i == j){
          
          # diagonal matrices
          cmat <- R2[i, j] * diag(1, length(time1))
          
        }else{
          
          # off diagonal matrices
          covmat <- outer(1:length(time1), 1:length(time2), FUN = "paste", sep = ",")
          
          comb <- cbind(1:length(time1), match(time1, time2, nomatch = NA))
          
          if(all(is.na(comb[,2]))){
            
            cmat <- matrix(0, nrow = length(time1), ncol = length(time2))
            
          }else{
            
            covmat[covmat %in% apply(comb, 1, paste, collapse = "," )] <- 1
            covmat[covmat != 1] <- 0
            
            cmat <- (covmat %>% as.data.frame %>% sapply(as.numeric) - 1) * R2[i, j]
            
          }
        }
        
        Cmat <- cbind(Cmat, cmat)
      }
      
      Sigma2 <- rbind(Sigma2, Cmat)
      
    }
    
    VarYi[[n]] <- Z%*%G2%*%t(Z) + Sigma2
    
    Vinv[[n]] <- solve(VarYi[[n]])
    
    Sigma2i[[n]] <- Sigma2
    
    precB2[[n]] <- t(Xi[[n]]) %*% Vinv[[n]] %*% Xi[[n]]
    
    xvinvy[[n]] <- t(Xi[[n]]) %*% Vinv[[n]] %*% Yi[[n]]
    
    
  }
  
  # 1. beta
  
  varBeta <- solve(Reduce("+", precB2))
  
  XvinvY <- Reduce("+", xvinvy)
  
  estbeta <- varBeta %*% XvinvY
  
  
  # 2. bi
  
  estbi <- list()
  
  for (n in 1:length(allid)){
    
    estbi[[n]] <- G2 %*% t(Zi[[n]]) %*% Vinv[[n]] %*% (Yi[[n]] - Xi[[n]] %*% estbeta)
    
  }
  
  return(list(Beta = estbeta, bi = estbi, Vinv = Vinv, Zi = Zi, Xi = Xi, Yi = Yi, nobsmat = nobsmat))
  
}


##### MCMC - posterior distributions #####

vcvmat <- combmcmc$VCV

# select random rows from vcvmat
B <- 500
set.seed(123)
vcvmat <- vcvmat[sample(1:nrow(vcvmat), B), ]

# G, R, Beta, Vinv
G = R = Beta = bi = Vinv <- list()

for(b in 1:nrow(vcvmat)){
  
  tmat <- vcvmat[, endsWith(colnames(vcvmat), "pid")]
  G[[b]] = G2 <- matrix(tmat[b, ], sqrt(ncol(tmat)))[ordervec, ordervec]
  
  tmatr <- vcvmat[, endsWith(colnames(vcvmat), "units")]
  R[[b]] = R2 <- matrix(tmatr[b, ], sqrt(ncol(tmatr)))
  
  est <- est_coeff_Yvar(G2, R2)
  
  Beta[[b]] <- est$Beta
  Vinv[[b]] <- est$Vinv
  bi[[b]] <- est$bi
  
}


# save(Beta, file = "Beta.rdata")
# save(Vinv, file = "Vinv.rdata")
# save(bi, file = "bi.rdata")
# save(G, file = "G.rdata")
# save(R, file = "R.rdata")

Xi <- est$Xi; Zi <- est$Zi; Yi <- est$Yi; nobsmat <- est$nobsmat

load("Beta.rdata"); load("Vinv.rdata"); load("bi.rdata"); load("R.rdata"); load("G.rdata")


##### fitted values #####


vcvmat <- combmcmc$VCV
solmat <- combmcmc$Sol

# select random rows from vcvmat
B <- 500

nfe <- 36

FEsolmat <- solmat[sample(1:nrow(solmat), B), 1:nfe]
REsolmat <- solmat[sample(1:nrow(solmat), B), (nfe + 1):ncol(solmat)]
REpid <- REsolmat[, endsWith(colnames(REsolmat), paste0("pid.", allid[1]))]

# reformatting solutions by variables
colorder <- c(colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitpFVCt")], 
              colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitpDLCOt")],
              colnames(FEsolmat)[startsWith(colnames(FEsolmat), "traitEFt")])

ordFEsolmat <- FEsolmat[, order(match(colnames(FEsolmat), colorder))] 



preddat <- NULL

for(n in 1:length(allid)){
  
  Xn <- Xi[[n]]; Zn <- Zi[[n]]
  
  REpid <- REsolmat[, endsWith(colnames(REsolmat), paste0("pid.", allid[n]))]
  
  REcolorder <- c(colnames(REpid)[startsWith(colnames(REpid), "traitpFVCt")], 
                  colnames(REpid)[startsWith(colnames(REpid), "traitpDLCOt")],
                  colnames(REpid)[startsWith(colnames(REpid), "traitEFt")])
  
  ordREpid <- REpid[, order(match(colnames(REpid), REcolorder))] 
  
  trajpid <- Xn %*% t(ordFEsolmat) + Zn %*% t(ordREpid)
  
  btime <- apply(Zn[, c(2, 6, 10)], 1, min)
  
  newpreddat <- data.frame(allid[n], btime, apply(trajpid, 1, mean), 
                           apply(trajpid, 1, function(x){quantile(x, probs = c(0.025, 0.975))}) %>% t,
                           c(rep(variablenames[1], nobsmat[n, 1]), rep(variablenames[2], nobsmat[n, 2]),
                             rep(variablenames[3], nobsmat[n, 3]))) %>% 
    set_colnames(c("pid", "btime", "mean", "lwr", "upr", "measure"))
  
  print(n)
  preddat <- rbind(preddat, newpreddat)
  
}





##### distribution of predicted values #####

# function to generate design matrix for future value & calculate mean and var

predprob <- function(pid, predtime, G2, R2, vinv, beta){
  
  n <- which(allid == pid)
  
  # Zi+: random effects
  Zip <- data.frame(int = 1, btime = predtime) %>%
    mutate(sp10yr = ifelse(btime + 10 > 0, btime + 10, 0), 
           sp3yr = ifelse(btime + 3 > 0, btime + 3, 0)) %>% as.matrix
  
  Zp <- Reduce(direct.sum, list(Zip, Zip, Zip))
  
  # Xip: fixed effects
  mdat <- regdat.d %>% filter(pid == allid[n]) %>% head(1)
  
  Xip <- cbind(mdat$int, matrix(ns(mdat[, "last"] + predtime, knots = c(10, 30),
                                   Boundary.knots= c(0, 40))[1, ], nrow = 1),
               mdat %>% select(FEcolnames[-c(1:4)])) %>% as.matrix
  
  Xp <- Reduce(direct.sum, list(Xip, Xip, Xip))
  
  
  # Vp, Cp
  Vp <- Zp %*% G2 %*% t(Zp) + R2
  Cp <- Zi[[n]] %*% G2 %*% t(Zp)
  
  # predicted values & confidence interval
  mu <- (Xp %*% beta) + (t(Cp) %*% vinv[[n]] %*% (Yi[[n]] - Xi[[n]] %*% beta))
  sigma <- Vp - (t(Cp) %*% vinv[[n]] %*% Cp)
  
  return(list(mu = mu, sigma = sigma))
  
}








##### plot estimated health trajectories & predictions with uncertainty #####

# select patient
selectid <- 531

# set time of prediction = 1 year from the last fvc
predtime <- 1

# set threshold
fvc70 <- -0.41; fvc60 <- -0.91
ef50 <- -1.75; ef35 <- -2.43

#plot_fvc_pidest <- function(fvcres, selectid, predtime){


# distribution of predicted value
preddist <- predprob(pid = selectid, predtime = predtime, G2 = G2, R2 = R2, vinv = vinv, beta = ordered.postmean)

#fvcres = preddist$mu[1]; selectid = selectid; predtime = predtime

# estimation at predtime
predmean_fvc = fvcres<- preddist$mu[1]
predvar_fvc <- preddist$sigma[1, 1]

x <- seq(from = predmean_fvc - 2.5*sqrt(predvar_fvc), to = predmean_fvc + 2.5*sqrt(predvar_fvc), by = .01)
simdf <- data.frame(x = x, y = dnorm(x, mean = predmean_fvc, sd = sqrt(predvar_fvc)))

datpid <- regdat.d %>% filter(pid == selectid) 
preddat_meas <- preddat %>% filter(measure == variablenames[1]) %>% filter(pid == selectid)


# values for plotting
my_col <- "#00998a"
linecol <- "#0073C2FF"
axissize <- 9
ylimits <- c(-2.5, 2.5)


p1 <- ggplot() + 
  
  geom_ribbon(data = preddat_meas, aes(x = btime, ymin = lwr, ymax = upr), alpha = 0.2) +
  
  geom_point(data = preddat, mapping = aes(x = predtime, y = predmean_fvc), 
             size = 0.9, shape = 21, fill = "white") +
  
  geom_hline(yintercept = fvc70, color = my_col, size = 0.7, alpha = 0.3) + 
  geom_hline(yintercept = fvc60, color = my_col, size = 0.7, alpha = 0.5) + 
  geom_line(data = preddat_meas, aes(btime, mean), col = linecol, size = 0.7) +
  
  geom_point(data = datpid, aes(btime, pFVCt), size = 0.8) + 
  
  scale_y_continuous(limits = ylimits, expand = c(0, 0)) +
  scale_x_continuous(limits = c(min(datpid$btime[!is.na(datpid$pFVCt)]) - 1, predtime + 1.5), 
                     expand = c(0, 0)) +
  
  labs(x = "Time from current visit", y = "pFVC") +
  theme(axis.text.x = element_text(size = axissize),
        axis.text.y = element_text(size = axissize),  
        axis.title.x = element_text(size = axissize),
        axis.title.y = element_text(size = axissize))



if(subset(simdf, x >= -6 & x < fvc60) %>% nrow == 0){
  
  
  if(subset(simdf, x >= -6 & x < fvc70) %>% nrow == 0){
    
    a1 <- ggplot() + 
      
      geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
      geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                   size = 0.3, color = "grey") +
      
      coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
      
      scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$pFVCt)]) - 1, predtime + 1.5), 
                         expand = c(0, 0)) +    
      
      labs(x = "", y = "") +
      theme(axis.text.x = element_text(size = axissize),
            axis.text.y = element_text(size = axissize),  
            axis.title.x = element_text(size = axissize),
            axis.title.y = element_text(size = axissize))
  }else{
    
    a1 <- ggplot() + 
      
      geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
      geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                   size = 0.3, color = "grey") +
      
      geom_ribbon(data = subset(simdf, x >= -6 & x < fvc70), 
                  aes(x = x, ymin = predtime, ymax = y + predtime), 
                  fill = my_col, alpha = 0.3) +
      
      coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
      
      scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$pFVCt)]) - 1, predtime + 1.5), 
                         expand = c(0, 0)) +    
      
      labs(x = "", y = "") +
      theme(axis.text.x = element_text(size = axissize),
            axis.text.y = element_text(size = axissize),  
            axis.title.x = element_text(size = axissize),
            axis.title.y = element_text(size = axissize))
  }
  
}else{
  
  a1 <- ggplot() + 
    
    geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
    geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                 size = 0.3, color = "grey") +
    
    geom_ribbon(data = subset(simdf, x >= -6 & x < fvc70), 
                aes(x = x, ymin = predtime, ymax = y + predtime), 
                fill = my_col, alpha = 0.3) +
    geom_ribbon(data = subset(simdf, x >= -6 & x < fvc60), 
                aes(x = x, ymin = predtime, ymax = y + predtime), 
                fill = my_col, alpha = 0.5) + 
    
    coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
    
    scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$pFVCt)]) - 1, predtime + 1.5), 
                       expand = c(0, 0)) +    
    
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(size = axissize),
          axis.text.y = element_text(size = axissize),  
          axis.title.x = element_text(size = axissize),
          axis.title.y = element_text(size = axissize))
}



g1 <- ggplotGrob(p1)


#pdf(file = "test.pdf", width = 6, height = 4)
grid.newpage()
a1; grid.draw(g1)
#dev.off()


#}






#plot_ef_pidest <- function(efres, selectid, predtime){


#efres = efpidprob$postmean; selectid = selectid; predtime = predtime

# estimation at predtime
predmean_ef = efres <- preddist$mu[3]
predvar_ef <- preddist$sigma[3, 3]

x <- seq(from = predmean_ef - 3*sqrt(predvar_ef), to = predmean_ef + 3*sqrt(predvar_ef), by = .005)
simdf <- data.frame(x = x, y = dnorm(x, mean = predmean_ef, sd = sqrt(predvar_ef)))

datpid <- regdat.d %>% filter(pid == selectid) 
preddat_meas <- preddat %>% filter(measure == variablenames[3]) %>% filter(pid == selectid)


# values for plotting
my_col <- "red"
linecol <- "#0073C2FF"
axissize <- 9
ylimits <- c(-4, 3.5)


p3 <- ggplot() + 
  
  geom_ribbon(data = preddat_meas, aes(x = btime, ymin = lwr, ymax = upr), alpha = 0.2) +
  
  geom_point(data = preddat, mapping = aes(x = predtime, y = predmean_ef), 
             size = 0.9, shape = 21, fill = "white") +
  
  geom_hline(yintercept = ef50, color = my_col, size = 0.7, alpha = 0.3) + 
  geom_hline(yintercept = ef35, color = my_col, size = 0.7, alpha = 0.5) + 
  geom_line(data = preddat_meas, aes(btime, mean), col = linecol, size = 0.7) +
  
  geom_point(data = datpid, aes(btime, EFt), size = 0.8) + 
  
  scale_y_continuous(limits = ylimits, expand = c(0, 0)) +
  scale_x_continuous(limits = c(min(datpid$btime[!is.na(datpid$EFt)]) - 1, predtime + 1.5), 
                     expand = c(0, 0)) +
  
  labs(x = "Time from current visit", y = "EF") +
  theme(axis.text.x = element_text(size = axissize),
        axis.text.y = element_text(size = axissize),  
        axis.title.x = element_text(size = axissize),
        axis.title.y = element_text(size = axissize))


if(subset(simdf, x >= -6 & x < ef35) %>% nrow == 0){
  
  if(subset(simdf, x >= -6 & x < ef50) %>% nrow == 0){
    
    a3 <- ggplot() + 
      
      geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
      geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                   size = 0.3, color = "grey") +
      
      coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
      
      scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$EFt)]) - 1, predtime + 1.5), 
                         expand = c(0, 0)) +
      
      labs(x = "", y = "") +
      theme(axis.text.x = element_text(size = axissize),
            axis.text.y = element_text(size = axissize),  
            axis.title.x = element_text(size = axissize),
            axis.title.y = element_text(size = axissize))
    
  }else{
    
    a3 <- ggplot() + 
      
      geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
      geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                   size = 0.3, color = "grey") +
      
      geom_ribbon(data = subset(simdf, x >= -6 & x < ef50), 
                  aes(x = x, ymin = predtime, ymax = y + predtime), 
                  fill = my_col, alpha = 0.3) +
      
      coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
      
      scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$EFt)]) - 1, predtime + 1.5), 
                         expand = c(0, 0)) +
      
      labs(x = "", y = "") +
      theme(axis.text.x = element_text(size = axissize),
            axis.text.y = element_text(size = axissize),  
            axis.title.x = element_text(size = axissize),
            axis.title.y = element_text(size = axissize))
  }
  
}else{
  
  a3 <- ggplot() + 
    
    geom_line(data = simdf, aes(x = x, y = y + predtime), size = 0.3, color = "grey") + 
    geom_segment(aes(x = min(x), xend = max(x), y = predtime, yend = predtime), 
                 size = 0.3, color = "grey") +
    
    geom_ribbon(data = subset(simdf, x >= -6 & x < ef50), 
                aes(x = x, ymin = predtime, ymax = y + predtime), 
                fill = my_col, alpha = 0.3) +
    
    geom_ribbon(data = subset(simdf, x >= -6 & x < ef35), 
                aes(x = x, ymin = predtime, ymax = y + predtime), 
                fill = my_col, alpha = 0.5) + 
    
    coord_flip() + scale_x_continuous(limits = ylimits, expand = c(0, 0)) +
    
    scale_y_continuous(limits = c(min(datpid$btime[!is.na(datpid$EFt)]) - 1, predtime + 1.5), 
                       expand = c(0, 0)) +
    
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(size = axissize),
          axis.text.y = element_text(size = axissize),  
          axis.title.x = element_text(size = axissize),
          axis.title.y = element_text(size = axissize))
  
}



g3 <- ggplotGrob(p3)


#pdf(file = "test.pdf", width = 6, height = 4)
grid.newpage()
a3; grid.draw(g3)
#dev.off()


#}


















