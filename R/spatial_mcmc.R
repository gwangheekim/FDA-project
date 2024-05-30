spatial_mcmc <- function(y, X, W, K, jump_lamb, pr_b_sd, niter, nburn, nthin, type = "lag", sc = FALSE) {
  time_v <- seq(0, 1, length.out = ncol(X))
  my_basis <- create.bspline.basis(c(0, 1), nbasis = K)

  if (sc) {
    X <- apply(X, 2, scale)
  }
  N <- length(y)
  X1_F <- Data2fd(time_v, t(X), my_basis)
  A <- t(coef(X1_F))
  a = b = 1e-03

  beta_mat <- matrix(0, nrow = niter, ncol = K)
  s2_vec <- rep(0, niter)
  lamb_vec <- rep(0, niter)
  lik_vec <- rep(0, niter)
  accept_lamb <- 0

  beta_o <- rep(1, K)
  s2_o <- 1
  lamb_o <- 0.5
  mK <- rep(0,K)

  N = nrow(X)
  Z = A
  mK = rep(0,K)
  beta_s2 = pr_b_sd

  pb <- progress_bar$new(
    format = " Progress: [:bar] :percent, Estimated completion time: :eta",
    total = niter,
    clear = FALSE,
    width = 80
  )

  if(type == "lag"){
    print("spatial lag model")
    for(i in 1:niter){
      pb$tick()
      B = diag(N) - lamb_o*W
      SigK = pr_b_sd * diag(K)
      mean_b = solve(t(Z)%*%Z + s2_o*solve(SigK))%*%(t(Z)%*%B%*%y + s2_o*solve(SigK)%*%mK)
      Sigma_b = solve(t(Z)%*%Z + s2_o*solve(SigK))

      beta_o = MASS::mvrnorm(1,mean_b,Sigma_b)

      param1 = (N/2) + a
      param2 = (t(B%*%y - Z%*%beta_o)%*%(B%*%y - Z%*%beta_o) + 2*b)/2
      s2_o = rinvgamma(1,param1,param2)
      while(T){
        lamb_n = lamb_o + rnorm(1,0,jump_lamb)
        if((lamb_n <= 1)&(lamb_n >= -1)){
          break
        }
      }

      num = t((diag(N) - lamb_n*W)%*%y - Z%*%beta_o)%*%((diag(N) - lamb_n*W)%*%y - Z%*%beta_o)*(-1/(2*s2_o)) + log(abs(det(diag(N) - lamb_n*W)))
      den = t((diag(N) - lamb_o*W)%*%y - Z%*%beta_o)%*%((diag(N) - lamb_o*W)%*%y - Z%*%beta_o)*(-1/(2*s2_o)) + log(abs(det(diag(N) - lamb_o*W)))
      ratio = num - den
      u = log(runif(1))

      if(u < ratio){
        lamb_o = lamb_n
        accept_lamb = accept_lamb + 1/niter
      }else{
        lamb_o = lamb_o
      }

      beta_mat[i,] = beta_o
      s2_vec[i] = s2_o
      lamb_vec[i] = lamb_o
      lik_vec[i] = c(-(N/2)*log(2*pi) - (N/2)*log(s2_o) - (1/(2*s2_o))*t(B%*%y - Z%*%beta_o)%*%(B%*%y - Z%*%beta_o) + log(abs(det(B)))) + dmvnorm(beta_o, mK, SigK, log=T) + log(dinvgamma(s2_o,a,b)) + log(0.5)
    }

    inds = seq(nburn,niter,nthin)
    lamb_est = mean(lamb_vec[inds])
    s2_est = mean(s2_vec[inds])
    beta_est = colMeans(beta_mat[inds,])
    Best = diag(N) - lamb_est*W
    lik = c(-(N/2)*log(2*pi) - (N/2)*log(s2_est) - (1/(2*s2_est))*t(Best%*%y - Z%*%beta_est)%*%(Best%*%y - Z%*%beta_est) + log(abs(det(Best)))) + dmvnorm(beta_est, mK, SigK, log=T) + log(dinvgamma(s2_est,a,b)) + log(0.5)
    bic_res = -2*lik + (K+2)*log(N)
    aic_res = -2*lik + 2*(K+2)

  }else if(type == "error"){
    print("spatial error model")
    for(i in 1:niter){
      pb$tick()
      B = diag(N) - lamb_o*W
      SigK = beta_s2 * diag(K)
      mean_b = solve(t(Z)%*%t(B)%*%B%*%Z + s2_o*solve(SigK)) %*% (t(Z)%*%t(B)%*%B%*%y + s2_o*solve(SigK)%*%mK)
      Sigma_b = s2_o*solve(t(Z)%*%t(B)%*%B%*%Z + s2_o*diag(K))

      beta_o = MASS::mvrnorm(1,mean_b,Sigma_b)

      param1 = (N/2) + a
      param2 = (t(y - Z%*%beta_o)%*%t(B)%*%B%*%(y - Z%*%beta_o) + 2*b)/2
      s2_o = rinvgamma(1,param1,param2)
      while(T){
        lamb_n = lamb_o + rnorm(1,0,jump_lamb)
        if((lamb_n <= 1)&(lamb_n >= -1)){
          break
        }
      }

      num = t(y-Z%*%beta_o)%*%t(diag(N) - lamb_n*W)%*%(diag(N) - lamb_n*W)%*%(y-Z%*%beta_o)*(-1/(2*s2_o)) + log(abs(det(diag(N) - lamb_n*W)))
      den = t(y-Z%*%beta_o)%*%t(diag(N) - lamb_o*W)%*%(diag(N) - lamb_o*W)%*%(y-Z%*%beta_o)*(-1/(2*s2_o)) + log(abs(det(diag(N) - lamb_o*W)))
      ratio = num - den
      u = log(runif(1))

      if(u < ratio){
        lamb_o = lamb_n
        accept_lamb = accept_lamb + 1/niter
      }else{
        lamb_o = lamb_o
      }

      beta_mat[i,] = beta_o
      s2_vec[i] = s2_o
      lamb_vec[i] = lamb_o
      lik_vec[i] = c(-(N/2)*log(2*pi) - (N/2)*log(s2_o) - (1/(2*s2_o))*t(y - A%*%beta_o)%*%t(B)%*%B%*%(y-A%*%beta_o) + log(abs(det(B)))) + dmvnorm(beta_o, mK, SigK, log=T) + log(dinvgamma(s2_o,a,b)) + log(0.5)
    }

    inds = seq(nburn,niter,nthin)
    lamb_est = mean(lamb_vec[inds])
    s2_est = mean(s2_vec[inds])
    beta_est = colMeans(beta_mat[inds,])
    Best = diag(N) - lamb_est*W
    lik = c(-(N/2)*log(2*pi) - (N/2)*log(s2_est) - (1/(2*s2_est))*t(y - A%*%beta_est)%*%t(Best)%*%Best%*%(y-A%*%beta_est) + log(abs(det(Best)))) + dmvnorm(beta_est, mK, SigK, log=T) + log(dinvgamma(s2_est,a,b)) + log(0.5)
    bic_res = -2*lik + (K+2)*log(N)
    aic_res = -2*lik + 2*(K+2)
  }

  return(list(
    beta_mat = beta_mat[inds, ],
    s2_vec = s2_vec[inds],
    lamb_vec = lamb_vec[inds],
    lik_vec = lik_vec[inds],
    s2_est = s2_est,
    beta_est = beta_est,
    lamb_est = lamb_est,
    accept_lamb = accept_lamb,
    grid = time_v,
    A = Z,
    my_basis = my_basis,
    lik = lik,
    aic = aic_res,
    bic = bic_res,
    mcmc_opt = list(niter=niter,
                    nburn=nburn,
                    nthin=nthin)
  ))
}
