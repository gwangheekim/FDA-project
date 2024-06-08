#' Summary the result of spatial scalar on function regression.
#'
#' @description \link{summary.spfda} is used to summary the result of Spatial scalar on function regression.
#' @export
summary.spfda <- function(object, CI = 0.95){
  est = object$beta_est
  quant = t(apply(object$beta_mat, 2, function(x) quantile(x, probs = c((1 - CI)/2, 1-(1 - CI)/2))))
  ci.temp = c((1-CI)/2, (1+CI)/2)
  beta.summary = data.frame(cbind(est,quant))
  cname = paste("beta.coef", 1:ncol(object$beta_mat), sep=" ")
  colnames(beta.summary) = c("Posterior mean",
                             paste0(ci.temp[1]*100, '%'),
                             paste0(ci.temp[2]*100, '%'))
  rownames(beta.summary) = cname
  res = list(coef = beta.summary,
             mcmc_opt = object$mcmc_opt,
             aic = object$aic,
             bic = object$bic)
  class(res) <- "summary.spfda"
  res
}
