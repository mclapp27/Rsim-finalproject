#' Generates SNP and phenotype data based on mutualism model 
#' 
#' goal: generate correlated residuals, intermediate phenotype = residual + snp effect
#'
#' @param N number of data elements to generate
#' @param j number of SNPs and phenotypes to generate
#' @param p_x list of x parameters including means, standard deviations, and correlations
#' @param p_y list of y parameters including the intercept (b0), vector of regression effects for each x (b), and r2
#' @return returns a list containing 5 data frame elements: SNP data, EP data, beta data, Design matrix, and phenotype data

mutual_sim <- function(N = 1000, j = 100, p, p_x, p_y, lamO = 0.95, lamU = 0.3){
  #create data structures
  SNPs <- matrix(NA, nrow = N , ncol = j)
  betas <- matrix(NA, nrow = j, ncol = 1)
  EPs <- matrix(NA, nrow = N, ncol = j)
  Y <- matrix(NA, nrow = N , ncol = p)
  Design <- matrix(0, nrow = p,ncol = j) #design matrix
  Phenos <- matrix(NA, nrow = N, ncol = p) #phenotypes
  colnames(Phenos) <- paste0("pheno", 1:p)
  
  #naming
  colnames(SNPs) <- paste0("SNP", 1:j)
  rownames(betas) <- paste0("beta", 1:j)
  colnames(EPs) <- paste0("EP", 1:j)
  colnames(Y) <- paste0("pheno", 1:p)
  colnames(Design) <- paste0("EP", 1:j)
  rownames(Design) <- paste0("pheno", 1:p)
  #generate SNPs and betas
  for (i in 1:j) {
    SNPs[ ,i] <- rbinom(N, 2, runif(1, .05, .5)) #create SNPs with MAF=2pq
    betas[i, ] <- runif(1, .30, .50) * (2 * (rbinom(1, 1, .5) - .5))#SNP effect sizes in per standardized genotype metric (note that SNP is scaled below)
    EPs[, i] <- betas[i, ] * scale(SNPs[, i]) + (1 - (betas[i, ]^2))^.5 * rnorm(N, 0, 1) #endophenotype is equal to summation of 1 SNP effect and many other (unspecified genetic and environmental) effects that together form a normal distribution
  }
  allSNP <- 1:j
  for (i in 1:p) { 
    design <- sample(allSNP, round(j/p))
    Design[i, design] <- 1
    allSNP <- allSNP[! allSNP %in% intersect(allSNP, design)]
    Phenos[, i] <- scale(rowMeans(EPs[, design])) * lamO + rnorm(N, 0, 1) * lamU
  }
  Sigma_x <- covmat_x(p_x) #create covariance matrix
  X <- Phenos + rmvnorm(N, p_x$mu, Sigma_x) #create X
  var_e <- (t(p_y$b) %*% cov(X) %*% p_y$b) * ((1/p_y$r2) - 1) #variance term
  SigY <- matrix(0, ncol = p, nrow = p)
  diag(SigY) <- var_e
  Y <- X + (lamU * rmvnorm(N, p_x$mu, SigY))
  #store results
  results <- list(
    SNPs = SNPs, 
    betas = betas, 
    EPs = EPs, 
    EP2 = X,
    Phenos = Y,
    SNPsPhenos = cbind(SNPs, Y))
  return(results)
}



