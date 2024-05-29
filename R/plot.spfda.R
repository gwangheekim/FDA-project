#' Plotting the result of Spatial scalar on functional regression
#' @description
#' \link{plot} is used to plot the estimated beta.
#'
#' @param object mcmc result
#' @param q credible interval parameter. default value is 0.05
#'
#' @export
plot.spfda <- function(object, q = 0.05){
  res = object
  Bf = fd(res$beta_est,res$my_basis)
  q1 = apply(res$beta_mat,2,function(x){quantile(x,c(q/2))})
  q3 = apply(res$beta_mat,2,function(x){quantile(x,c(1-q/2))})
  Bf_q1 = fd(q1,res$my_basis)
  Bf_q3 = fd(q3,res$my_basis)
  plot(NULL, ylim=c(min(Bf_q1$coefs),max(Bf_q3$coefs)), xlim = c(0,1),
       ylab = "beta", xlab = "time")
  points(res$grid,eval.fd(res$grid,Bf), type="l")
  lines(res$grid,eval.fd(res$grid,Bf_q1), lty=2, col="red")
  lines(res$grid,eval.fd(res$grid,Bf_q3), lty=2, col="red")
  abline(h=0,lty=2,lwd=2)
}
