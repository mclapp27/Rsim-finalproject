
## SETUP -------------------------------------------
# Source functions
source("overlapsim_func.R")
source("mutualsim_func.R")
source("rmvnorm.R")
source("cov_mat.R")
#load packages
library(beepr)
library(corrplot)
library(lavaan)
library(ggplot2)
library(plyr)
##RUN SIMULATION FUNCTIONS-------------------------
#mutualism parameters
p_x <- list(
  mu = as.matrix(rep(0, 9)),
  s  = as.matrix(rep(1, 9)),
  rho = 0.8
)
p_y <- list(
  b0 = 1,
  b = as.matrix(rep(1, 9)),
  r2 = 0.8
)
#run functions
r <- 100
M_Data <- list()
O_Data <- list()
for (i in 1:r) { 
  M_Data[[length(M_Data)+1]]<- mutual_sim(N = 1000, j = 100, p = 9, p_x, p_y)
  O_Data[[length(O_Data)+1]]<- overlap_sim(N = 1000, j = 100, p = 9, pleiotropy = c(0.7, 0.9))
  cat(i, " \r") #count function
  flush.console()
}

##CORRELATIONS-------------------
#format data
M_SP <- lapply(M_Data,"[[", 6)
O_SP <- lapply(O_Data,"[[", 6)
#Average corr matrix across iterations
Mcor <- lapply(M_SP, cor)#list of cor matrices per iterationns
m <- apply(simplify2array(Mcor), 1:2, mean) #simplify to one array of mean correlations across all iterations
csnps  <- m[101:109, 1:5] #matrix of snps x phenotypes correlation
cphens <- m[101:109, 101:109] #matrix of just phenotypes correlation

corrplot(as.matrix(cphens))
corrplot(as.matrix(csnps), is.corr = FALSE)

#Overlap: Average corr matrix across iterations
OVcor <- lapply(O_SP, cor) #create list of correlation matrices (1 per iteration)
ov <- apply(simplify2array(OVcor), 1:2, mean) #simplify to one array of mean correlation across all iterations
osnps  <- ov[101:109, 1:5] #matrix of snps x phenotypes correlation
ophens <- ov[101:109, 101:109] #matrix of just phenotypes correlation

corrplot(as.matrix(ophens))
corrplot(as.matrix(osnps), is.corr = FALSE) #remember: this won't look perfect because it's averaged across iterations, and every iteration samples differently 

###CFAS------------------------------------
#specify model
cfa1='F =~ pheno1+pheno2+pheno3+pheno4+pheno5+pheno6+pheno7+pheno8+pheno9'

#apply factor model to factor data
cfaM    <- lapply(M_SP, sem, model = cfa1) #fit CFA to data from each iteration: save all lavaan data in new list
Mutfit <- lapply(cfaM, fitMeasures, fit.measures = c("chisq", "pvalue", "cfi","rmsea","srmr")) #new list of only fit measures
Mutsol <- lapply(cfaM, standardizedsolution) #new list for standardized solultions only

#apply factor model to overlap data
cfaO    <- lapply(O_SP, sem, model = cfa1) #fit CFA to data from each iteration: save all lavaan data in new list
Overfit <- lapply(cfaO, fitMeasures, fit.measures = c("chisq", "pvalue", "cfi","rmsea","srmr")) #new list of only fit measures
Oversol <- lapply(cfaO, standardizedsolution) #new list for standardized solultions only

#Avg. Fit stats across iterations
MF <- apply(simplify2array(Mutfit), 1:2, mean) #take mean fit stats across iterations
OF <- apply(simplify2array(Overfit), 1:2, mean) #take mean fit stats across iterations
rowMeans(MF) #avg Factor generating model fits
rowMeans(OF) #avg Overlap generating model fits

#Avg. Solutions across iterations
#Factor
ms <- do.call(cbind, Mutsol) #combine all list dataframes into one array
avgMsol <- data.frame(path = paste(ms$lhs, ms$rhs)) #create dataframe of averages with one column of path labels
#subset into separate dataframes for each variable
jest<- ms[ , (names(ms) %in% c("est.std"))]
jse <- ms[ , (names(ms) %in% c("se"))]
jp  <- ms[ , (names(ms) %in% c("pvalue"))]
#take mean across rows to get average solution across all iterations
avgMsol$est  <- rowMeans(jest[ , ])
avgMsol$se   <- rowMeans(jse[ , ])
avgMsol$pval <- rowMeans(jp[ , ])

#Overlap
os <- do.call(cbind, Oversol) #combine all list dataframes into one array
avgOsol <- data.frame(path = paste(os$lhs, os$rhs)) #create dataframe of averages with one column of path labels
#subset into separate dataframes for each variable
jest<- os[ , (names(os) %in% c("est.std"))]
jse <- os[ , (names(os) %in% c("se"))]
jp  <- os[ , (names(os) %in% c("pvalue"))]
#take mean across rows to get average solution across all iterations
avgOsol$est  <- rowMeans(jest[ , ])
avgOsol$se   <- rowMeans(jse[ , ])
avgOsol$pval <- rowMeans(jp[ , ])

##QSNP ESTIMATION-------------------
#define function to generate model fits for each snp
q_calc <- function(SNPsPhenos){
  require(lavaan)
  #define models
  cp <- ' F=~pheno1+pheno2+pheno3+pheno4+pheno5+pheno6+pheno7+pheno8+pheno9
        F ~ SNP'
  ip <- ' F=~pheno1+pheno2+pheno3+pheno4+pheno5+pheno6+pheno7+pheno8+pheno9
        pheno1+pheno2+pheno3+pheno4+pheno5+pheno6+pheno7+pheno8+pheno9~SNP'
  SNPs <- SNPsPhenos[,1:100]
  Phens <- SNPsPhenos[,101:109]
  resultI <- list()
  resultC <- list()
  for(i in 1:ncol(SNPs)){
    resultC[[length(resultC)+1]] <- sem(model = cp, 
                                        data = cbind("SNP" = SNPs[,i], Phens))
    resultI[[length(resultI)+1]] <- sem(model = ip, 
                                        data = cbind("SNP" = SNPs[,i], Phens))
  }
  R <- list(
    CP = resultC, 
    IP = resultI
  )
  return(R)
}
#apply function to get all full models
MQ_Mods <- lapply(M_SP, q_calc)
OQ_Mods <- lapply(O_SP, q_calc)
#split into two lists for each: one CP list, one IP list
MCP_Full <- lapply(MQ_Mods,"[[", 1)
MIP_Full <- lapply(MQ_Mods,"[[", 2)
OCP_Full <- lapply(OQ_Mods,"[[", 1)
OIP_Full <- lapply(OQ_Mods,"[[", 2)
#functions for pulling fits and solutions from models
pull_fit <- function(FitList){
  result <- lapply(FitList, fitMeasures, fit.measures = c("chisq", "pvalue", "cfi","rmsea","srmr"))
  return(result)
}
pull_sols <- function(SolList){
  result <- lapply(SolList, standardizedsolution)
  return(result)
}
#pull model fits and solutions
MCP_Fit <- lapply(MCP_Full, pull_fit) 
MCP_Sol <- lapply(MCP_Full, pull_sols)
MIP_Fit <- lapply(MIP_Full, pull_fit)
MIP_Sol <- lapply(MIP_Full, pull_sols) 
OCP_Fit <- lapply(OCP_Full, pull_fit) 
OCP_Sol <- lapply(OCP_Full, pull_sols)
OIP_Fit <- lapply(OIP_Full, pull_fit)
OIP_Sol <- lapply(OIP_Full, pull_sols) 

#generate Qsnp by comparing models
##test line by line
q_comp <- function(ModList){
  result <- list()
  for (i in 1:length(ModList$CP)) { #loop that compares CP and IP to create Qsnp model comparison
    result[[length(result)+1]]<- lavTestLRT(ModList$CP[[i]], ModList$IP[[i]])
  }
  return(result)
}

MQ <- lapply(MQ_Mods, q_comp)
OQ <- lapply(OQ_Mods, q_comp)

###BUILD DATAFRAMES FOR GRAPHING------------------
#pull necessary data to format new lists
Mut_List <- list()
for(i in 1:r){
  Mut_List[[i]] <- list(
    CPSol = MCP_Sol[[i]], 
    Qdat = MQ[[i]]
  )
}

Over_List <- list()
for(i in 1:r){
  Over_List[[i]] <- list(
    CPSol = OCP_Sol[[i]], 
    Qdat = OQ[[i]]
  )
}
#function to build dataframe for each iteration of simulation
convert_data <- function(ModList){
  data <- data.frame(snp = rep(c(1:length(ModList$CPSol))))
  for (i in (1:length(ModList$CPSol))) {
    #Overlap
    data$Q[i]  <-    ModList$Qdat[[i]]$`Chisq diff`[2]
    data$QP[i] <-    ModList$Qdat[[i]]$`Pr(>Chisq)`[2]
    data$betaf[i] <- ModList$CPSol[[i]]$est.std[10]
    data$sef[i] <-   ModList$CPSol[[i]]$se[10]
    data$pf[i] <-    ModList$CPSol[[i]]$pvalue[10]
  }
  return(data)
}
#run function so each iteration has a dataframe including data for all snps
MutDF <- lapply(Mut_List, convert_data)
OverDF <- lapply(Over_List, convert_data)
#function to add remaining caluclated stats to dataframes
stat_convert <- function(data) {
  q2p<-pchisq(data$Q, 2, lower.tail = FALSE)
  data$Q2<-qchisq(q2p, 1, lower.tail = FALSE) 
  data$betachi <- (data$betaf/data$sef)^2
  data$ratio  <-  (data$betachi/data$Q2) 
  return(data)
}
#apply across dataframe lists
MutDF <- lapply(MutDF, stat_convert)
OverDF <- lapply(OverDF, stat_convert)

#function to get withiin iteration means
within_means <- function(data){
  QMean <- mean(data$Q2)
  AMean <- mean(data$betachi)
  return(c(QMean, AMean))
}

within_sds <- function(data){
  Qsd <- sd(data$Q2)
  Asd <- sd(data$betachi)
  return(c(Qsd, Asd))
}

#array of means by iteration for each generating model
Mmeans <- sapply(MutDF, within_means)
Omeans <- sapply(OverDF, within_means)

Msds <- sapply(MutDF, within_sds)
Osds <- sapply(OverDF, within_sds)

##distributions of Q and Association for visualization
chidist<-data.frame(group=c(rep("Mutualism", ncol(Mmeans)), 
                            rep("Overlap", ncol(Omeans))), 
                    chisquare=c(Mmeans[2, ], Omeans[2, ]))

#ggplot histogram of chisquare statistic 
chist<-ggplot(chidist, aes(x=chisquare, fill = group, color=group)) + 
  geom_histogram(binwidth = 0.1, position = "identity", alpha = 0.4) +
  geom_vline(aes(xintercept = mean(Mmeans[2, ])), linetype = "dashed", color = "#1B9E77", size = 0.5, alpha = 0.8) + #dist mean
  geom_vline(aes(xintercept = mean(Omeans[2, ])), linetype = "dashed", color = "#D95F02", size = 0.5, alpha = 0.8) + #dist mean
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+theme_minimal() +
  labs(title="Chisquare Distribution by Generating Model",x="Chisquare", y = "Density")
chist

#new dataframe to arrange and label data for graphing for Q
qdist<-data.frame(group = c(rep("Mutualism", ncol(Mmeans)), 
                          rep("Overlap", ncol(Omeans))), 
                  chisquare = c(Mmeans[1, ], Omeans[1, ]))
#ggplot histogram of Q statistic
qhist<-ggplot(qdist, aes(x=chisquare, fill = group, color=group)) + 
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.4) +
  geom_vline(aes(xintercept = mean(Mmeans[1, ])), linetype = "dashed", color="#1B9E77", size = 0.5, alpha = 0.8) +
  geom_vline(aes(xintercept = mean(Omeans[1, ])), linetype = "dashed", color="#D95F02", size = 0.5, alpha = 0.8) +
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+theme_minimal() +
  labs(title="Qsnp Distribution by Generating Model",x="Qsnp", y = "Density")
qhist

##RUN JACKKNIFE-----------------------------
#if Q >= 8 or association >= 1.3 then models detected as different
JK_Data <- vector(length = 100)
for(i in 1:r){
if(abs(mean(MutDF[[i]]$Q2) - mean(OverDF[[i]]$Q2)) >= 8
   || abs(mean(MutDF[[i]]$betachi) - mean(OverDF[[i]]$betachi)) >= 1.3){
  JK_Data[i] <- 1
}else{
  JK_Data[i] <- 0
}
}
#proportion success/fail
table(JK_Data)
#73% success
#27% #fail

#create jackknife function
jackknife <- function(data, func){
  # Obtain jackknife samples by using func and leaving one observation out
  jack <- sapply(seq_along(data), \(i) func(data[-i]))
  # Compute variance
  v <- ((length(jack)-1)/length(jack)) * sum((jack - mean(jack))^2)
  return(sqrt(v)) # Return SE
}


jackknife(JK_Data, mean)
#0.0446196

#results: percent that different is 73% - standard error of this probability is 5%
#what's the percent chance that they look different
#1 different, 0 the same - mean of that across all replications is a proportion
#jacknife: standard errors around that probability that represent the uncertainty of the probability estimate #empirical sampling variablilty of the probability that there's a detectable difference
