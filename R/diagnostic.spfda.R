#' Diagnostic the result of Spatial scalar on functional regression
#'
#' @description \code{diagnostic} is used to diagnostic the result of Spatial scalar on functional regression.
#'
#' @param object object of class \code{spfda}.
#' @param which Select parameter for mcmc diagnostic. "beta", "sigma2", "lambda" are available.
#' @param bidx vector of integer, Index of beta for checking MCMC chain convergence.
#' @return \code{diagnostic} returns plots for checking MCMC convergence.
#'
#' @export diagnostic
diagnostic <- function(object, which = "beta", bidx = 1){
  UseMethod("diagnostic")
}

#' @export
diagnostic.spfda <- function(object, which = "beta", bidx = 1){
  if(which == "beta"){
    if(length(bidx) == 1){

      vals = object$beta_mat[,bidx]
      ACF = acf(vals, plot=F)
      acf.val = ACF
      acf.print = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
      par(mfrow = c(1,3))
      ts.plot(vals, ylab=paste0("Beta", bidx), main = "Trace Plot")
      plot(acf.val, main = "")
      title("ACF", line = 1.6)
      plot(density(vals), main = paste0("Density"))
      rug(jitter(vals))
      par(mfrow = c(1,1))

    }else{
      for(i in 1:length(bidx)){
        vals = object$beta_mat[,bidx[i]]
        ACF = acf(vals, plot=F)
        acf.val = ACF
        acf.print = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
        par(mfrow = c(1,3))
        ts.plot(vals, ylab=paste0("Beta", bidx[i]), main = "Trace Plot")
        plot(acf.val, main = "")
        title("ACF", line = 1.6)
        plot(density(vals), main = paste0("Density"))
        rug(jitter(vals))
        par(mfrow = c(1,1))
        if(i != length(bidx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

  }else if(which == "lambda"){

    ACF = acf(object$lamb_vec, plot=F)
    acf.val = ACF
    acf.print = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
    par(mfrow = c(1,3))
    ts.plot(object$lamb_vec, ylab="Lambda", main = "Trace Plot")
    plot(acf.val, main = "")
    title("ACF", line = 1.6)
    plot(density(object$lamb_vec), main = paste0("Density"))
    rug(jitter(object$lamb_vec))
    par(mfrow = c(1,1))

  }else if(which == "sigma2"){

    ACF = acf(object$s2_vec, plot=F)
    acf.val = ACF
    acf.print = setNames(drop(ACF$acf), format(drop(ACF$lag), digits = 3L))
    par(mfrow = c(1,3))
    ts.plot(object$s2_vec, ylab="Lambda", main = "Trace Plot")
    plot(acf.val, main = "")
    title("ACF", line = 1.6)
    plot(density(object$s2_vec), main = paste0("Density"))
    rug(jitter(object$s2_vec))
    par(mfrow = c(1,1))

  }
}
