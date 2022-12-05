#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 4: Introduction to Monte Carlo Simulation
#
#' Generates a covariance matrix for data generation
#'
#' @param p_list list of parameters including standard deviation and correlation
#' @return returns a covariance matrix for data generation
#' 
covmat_x <- function(p_list) {
  with(p_list, {
    p <- NROW(s)
    # Create Sigma
    R <- matrix(rho, p, p)
    diag(R) <- 1
    D <- diag(c(s), p, p)
    D %*% R %*% D
  })
}

