#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 4: Introduction to Monte Carlo Simulation
#
#' Generates multivariate normally distributed data
#'
#' @param n number of data elements to generate
#' @param mu a matrix of population means
#' @param Sigma a population covariance matrix 
#' @return returns a matrix of multivariate normally distributed data based on mu and Sigma
rmvnorm <- function(n, mu, Sigma){
  stopifnot(is.matrix(mu),
            is.matrix(Sigma),
            is.numeric(n),
            dim(t(mu))[2] == dim(Sigma)[1]
  )
  p <- dim(t(mu))[2]
  Z <- matrix(rnorm(p*n), nrow = n, ncol = p)
  M1 <- matrix(1, nrow = n, ncol = 1)
  X <- M1 %*% t(mu) + Z %*% chol(Sigma)
  return(X) 
}
