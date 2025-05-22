data_generation <- function(n = 1000, j = 101, p = 4, gamma0 = rep(0.5, p),
                            h0t = 1, cenrate = 0.5, tau = 1.55, mev = 0.5,
                            out.p = 0)
{
  
  # n: number of observations
  # j: number of grid points for the functional predictor
  # p: number of scalar predictors
  # gamma0: parameters for the scalar predictors
  # h0t: baseline hazard function
  # cenrate: censoring rate
  # upper bound of the uniform distribution which is used to simulate censoring time
  # mev: measurement error variance
  
  grid_points <- seq(0, 1, length.out = j)
  base_mat <- matrix(, nrow = 22, ncol = j)
  for(i in 1:10)
    base_mat[i,] <- sqrt(2) * sin(pi * (2*i-1) * grid_points)
  for(i in 11:20)
    base_mat[i,] <- sqrt(2) * cos(pi * (2*i-21) * grid_points)
  base_mat[21,] <- 1
  base_mat[22,] <- grid_points
  
  Sigma <- matrix(0, 22, 22)
  Sigma[1,1] <- 1
  Sigma[11,11] <- 1
  for(i in 2:10)
    Sigma[i,i] <- (i-1)^(-2)
  for(i in 12:20)
    Sigma[i,i] <- (i-11)^(-2)
  Sigma[21,21] <- 1
  Sigma[22,22] <- 1
  
  beta_fun <- function(t){
    sigma <- 0.3
    out <- 0.3 * (sin(pi*t) - cos(pi*t) + sin(3*pi*t) - cos(3*pi*t) + sin(5*pi*t)/9 - cos(5*pi*t)/9 +
                    sin(7*pi*t)/16 - cos(7*pi*t)/16 + sin(9*pi*t)/25 - cos(9*pi*t)/25 +
                    (1/sqrt(2*pi)/sigma) * exp(-(t-0.5)^2/2/sigma^2))
    return(out)
  }
  
  beta_t <- beta_fun(grid_points)
  
  corr_mat <- matrix(0, nrow(Sigma), p)
  corr_mat[1, 1:p] <- 0.1
  corr_mat[1:p, 1] <- 0.1
  
  SigmaZ <- 0.5^t(sapply(1:p, function(i, k) abs(i-k), 1:p))
  Bigsigma <-rbind(cbind(Sigma, corr_mat), cbind(t(corr_mat), SigmaZ))
  mu <- rep(0, nrow(Bigsigma))
  data <- mvrnorm(n, mu, Bigsigma)
  us <- data[, 1:nrow(Sigma)]
  Xt_smooth <- us %*% base_mat
  options(warn = -1)
  msm_err <- matrix(rnorm(n*j, 0, sqrt(mev)), n, j)
  Xt_me <- Xt_smooth + msm_err
  Z <- data[, (nrow(Sigma) + 1):(nrow(Sigma) + p)]
  
  params <- exp(Xt_me %*% beta_t * (grid_points[2]-grid_points[1]) + Z %*% gamma0)
  failure_time <- rexp(n, rate = params)
  cens_time <- runif(n, 0, tau)
  event <- rep(0, n)
  event <- as.numeric(failure_time < cens_time)
  time <- failure_time * event + cens_time * (rep(1, n) - event)
  
  ### --- OUTLIERS HANDLING ---
  if(out.p > 0){
    n_out <- ceiling(out.p * n)
    out_ids <- sample(1:n, n_out)
    
    mag_err <- r_ou(n=n, t = grid_points, mu=10, alpha = 1, sigma = 1,
                    x0=rnorm(n=n, mean = 10, sd = 1/sqrt(2*1)))$data
    
    ## 1. Functional Magnitude Outliers: moderate scale + shift
    for (id in out_ids) {
      Xt_me[id, ] <- Xt_me[id, ] + mag_err[id, ]
    }
    
    ## 3. Response Time Outliers
    half <- floor(n_out / 2)
    early_ids <- out_ids[1:half]
    late_ids <- out_ids[(half + 1):n_out]
    
    # Early failures
    time[early_ids] <- runif(length(early_ids), min = 0.001, max = 0.05)
    event[early_ids] <- 1
    
    # Late failures
    time[late_ids] <- runif(length(late_ids), min = tau + 0.3, max = tau + 1)
    event[late_ids] <- 0

  }
  
  return(list(time=time, event=event, Xt = Xt_me, Xt_true=Xt_smooth, Z=Z, 
              beta=beta_t, gamma=gamma0, gp=grid_points))
}
