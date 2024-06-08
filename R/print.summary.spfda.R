#' Summary the result of spatial scalar on function regression.
#'
#' @description \link{summary.spfda} is used to summary the result of Spatial scalar on function regression.
#' @export
print.summary.spfda <- function(x, ...){
  cat("==========================","\n")
  cat("Summary of model","\n")
  cat("==========================","\n\n")
  cat(sprintf("MCMC sample of size %i, after burnin of %i iteration",
              x$mcmc_opt$nsamp, x$mcmc_opt$nburn),"\n\n")
  printCoefmat(x$coef)
  cat("\n---------------------------","\n\n")
  cat("Overall AIC (Smaller is better) :",x$aic,"\n")
  cat("Overall BIC (Smaller is better) :",x$bic,"\n")
}
