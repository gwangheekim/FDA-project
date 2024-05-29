#' @export
spfda <- function(y, ...) {
  UseMethod("spfda") # must not equal to package name(or dynlib.. error occur!)
}

#' Spatial MCMC for SPFDA
#'
#' A function to perform spatial MCMC for spatial functional data analysis.
#' @param y A numeric vector of observations.
#' @param X A numeric matrix of predictors.
#' @param W A spatial weights matrix.
#' @param K Number of basis functions.
#' @param jump_lamb Jump size for lambda in the MCMC.
#' @param pr_b_sd Prior standard deviation for beta.
#' @param niter Number of MCMC iterations.
#' @param nburn Number of burn-in iterations.
#' @param nthin Thinning interval for MCMC.
#' @param type Type of model ("lag" or "error").
#' @param sc Logical indicating if the predictors should be scaled.
#' @return A list containing MCMC samples and estimates.
#' @examples
#' \donttest{
#'n = 100;  grid = seq(0, 1, length.out = 101)
#'beta1 = sin(grid * 2 * pi)
#'beta2 = -dnorm(grid, mean=.2, sd=.03) +
#'       3*dnorm(grid, mean=.5, sd=.04)+
#'       dnorm(grid, mean=.75, sd=.05)
#'
#'X <- matrix(0, nrow=n, ncol=length(grid))
#'for(i2 in 1:n){
#'  X[i2,]=X[i2,]+rnorm(length(grid), 0, 1)
#'  X[i2,]=X[i2,]+runif(1, 0, 5)
#'  X[i2,]=X[i2,]+rnorm(1, 1, 0.2)*grid
#'  for(j2 in 1:10){
#'    e=rnorm(2, 0, 1/j2^(2))
#'    X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
#'    X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
#'  } }
#'
#'W = matrix(0,n,n)
#'W[upper.tri(W)] = W[lower.tri(W)] = runif(n*(n-1)/2)
#'
#'lamb = 0.7
#'K = 30
#'
#'Y1_err = X %*% beta1*0.01 + solve(diag(n) - lamb*W)%*%rnorm(n, 0, 3)
#'Y2_err = X %*% beta2*0.01 + solve(diag(n) - lamb*W)%*%rnorm(n, 0, 3)
#'
#'Y1_lag = solve(diag(n) - lamb*W)%*%(X %*% beta1)*0.01 + rnorm(n, 0, 3)
#'Y2_lag = solve(diag(n) - lamb*W)%*%(X %*% beta2)*0.01 + rnorm(n, 0, 3)
#'
#'res_err1 = spfreg(y = Y1_err,X = X,W = W,K = K,
#'                  jump_lamb = 0.1,pr_b_sd = 3^2,
#'                  niter = 15000, nburn = 2500, nthin = 10,
#'                  type = "error", sc = F)
#'summary(res_err1)
#'diagnostic(res_err1, which = "beta", bidx = 3)
#'plot(res_err1, q = 0.1)
#'diagnostic(res_err1, which = "lambda")
#' }
#'
#' @export
spfreg <- function(y, type = "lag", sc = FALSE, ...){
  output <- spatial_mcmc(..., y)
  class(output) = "spfda"
  return(output)
}

