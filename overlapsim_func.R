#' Generates SNP and phenotype data based on overlap model
#'
#' @param N number of data elements to generate
#' @param j number of SNPs to generate
#' @param p number of phenotypes to generate
#' @param lamO proportion of variance explained by overlapping endophenotypes
#' @param lamU proportion of variance explained by error/the environment
#' @param pleiotropy a vector of 2 values, the minimum and maximum proportion of endophenotypes sampled for any given phenotype
#' @param maxphen maximum number of phenotypes any single endophenotype can affect
#' @return returns a list containing 5 data frame elements: SNP data, EP data, beta data, Design matrix, and phenotype data

overlap_sim <- function(N = 10000, j, p, lamO = 0.95, lamU = 0.3, pleiotropy = c(0.5, 0.8), maxphen = 7){
  #error catches
  conditions <- c(N, j, p, lamO, lamU, pleiotropy, maxphen)
  lams <- c(lamO, lamU)
  stopifnot(is.numeric(conditions), sum(lams^2) <= 1, pleiotropy <= 1)
  
  #create data structures
  SNPs <- matrix(NA, nrow = N , ncol = j)
  betas <- matrix(NA, nrow = j, ncol = 1)
  EPs <- matrix(NA, nrow = N, ncol = j) #endophenotypes (or component processes or whatever)
  Design <- matrix(0, nrow = p,ncol = j) #design matrix
  Phenos <- matrix(NA, nrow = N, ncol = p) #phenotypes
  #col and row labels
  colnames(SNPs) <- paste0("SNP", 1:j)
  rownames(betas) <- paste0("beta", 1:j)
  colnames(betas) <- "parameter"
  colnames(EPs) <- paste0("EP", 1:j)
  colnames(Design) <- paste0("EP", 1:j)
  rownames(Design) <- paste0("pheno", 1:p)
  colnames(Phenos) <- paste0("pheno", 1:p)
  #generate SNPs, EPs, betas
  for (i in 1:j) {
    SNPs[ ,i] <- rbinom(N, 2, runif(1, .05, .5)) #create SNPs with MAF=2pq
    betas[i, ] <- runif(1, .30, .50) * (2 * (rbinom(1, 1, .5) - .5))#SNP effect sizes in per standardized genotype metric (note that SNP is scaled below)
    EPs[, i] <- betas[i, ] * scale(SNPs[, i]) + (1 - (betas[i, ]^2))^.5 * rnorm(N, 0, 1) #endophenotype (or cognitive process or whatever) is equal to summation of 1 SNP effect and many other (unspecified genetic and environmental) effects that together form a normal distribution
  }
  #set parameters
  allSNP <- 1:j
  mn <- pleiotropy[1] #minimum % of snps to sample
  mx <- pleiotropy[2] #maximum % of snps to sample
  #generate phenotypes
  for (i in 1:p) { 
    design <- sample(allSNP[(colSums(Design) < maxphen) == TRUE], sample(round(mn*j):round(mx*j), 1), 1)
    Design[i, design] <- 1
    Phenos[, i] <- scale(rowMeans(EPs[, design])) * lamO + rnorm(N, 0, 1) * lamU
  }
  #store results
  extra_data <- list(
    betas = betas,
    EPs = EPs,
    Design = Design
  )
  #format for CFA eval
  SNPsPhenos <- cbind(SNPs, Phenos)
  #format SNP and phenotype data for Q comparison
  SNPlist <- lapply(seq_len(ncol(SNPs)), function(i) SNPs[,i]) #make every list element data from 1 snp
  Bind <- lapply(SNPlist, cbind, Phenos) #bind single SNP data to all pheno data
  Bind <- lapply(Bind, as.matrix) #convert data within list to matrices
  Bind <- lapply(Bind, function(df) { #change first col name to SNP for model read
    colnames(df)[1] <- "SNP"
    df})
  #output results as 3 lists
  results <- list(
    effects = extra_data, 
    CFAMatrix = SNPsPhenos, 
    QMatrix = Bind)
  return(results)
}
