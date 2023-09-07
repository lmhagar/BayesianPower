## code to approximate power curves for equivalence tests with the gamma model
## this procedure leverages Algorithm 1 and Algorithm 2 in two separate functions

## load required libraries
require(foreach)
require(doParallel)
require(doSNOW)
require(numDeriv)
require(qrng)
require(nleqslv)

## function to approximate power curve using approximation method in Algorithm 1
powerCurveGamma1 <- function(conviction, power, delta_L, delta_U, gamma_alpha.1, gamma_beta.1, gamma_alpha.2, 
                             gamma_beta.2, threshold = NULL, m = 1024, q = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for an equivalence or noninferiority test
  ## delta_U: upper interval endpoint for an equivalence or noninferiority test
  ## alpha_j has a Gamma(mu0.j, tau0.j) prior, where tau0.j is a rate for group j = 1, 2
  ## beta_j has a Gamma(kappa0.j, lambda0.j) prior, where lambda0.j is a rate for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  
  ## find good upper bound for the sample size to generate the asymptotic power curve
  if (is.finite(delta_L) & is.finite(delta_U)){
    int_len <- delta_U - delta_L
  }
  else if (is.finite(delta_L)){
    int_len <- 2*(1- delta_L)    
  }
  else{
    int_len <- 2*(delta_U - 1)
  }
  
  alpha1 <- gamma_alpha.1; alpha2 <- gamma_alpha.2; beta1 <- gamma_beta.1; beta2 <- gamma_beta.2
  ## define theta metrics in terms of design values
  theta1 <- 1 - pgamma(threshold, alpha1, beta1)
  theta2 <- 1 - pgamma(threshold, alpha2, beta2)
  
  ## compute partial derivatives of theta with respect to alpha and beta for each group
  int_1 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha1, beta = beta1)$value
  
  d_alpha1 <- digamma(alpha1)*pgamma(threshold,alpha1, beta1) - int_1
  
  d_beta1 <- (alpha1/beta1)*(pgamma(threshold, alpha1+1, beta1) - pgamma(threshold, alpha1, beta1))
  
  int_2 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha2, beta = beta2)$value
  
  d_alpha2 <- digamma(alpha2)*pgamma(threshold, alpha2, beta2) - int_2
  
  d_beta2 <- (alpha2/beta2)*(pgamma(threshold, alpha2+1, beta2) - pgamma(threshold, alpha2,beta2))
  
  ## compute Fisher information for each gamma model
  Fish1inv <- matrix((1/(trigamma(alpha1)*alpha1 - 1))*c(alpha1, beta1, beta1, trigamma(alpha1)*beta1^2), nrow = 2)
  Fish2inv <-  matrix((1/(trigamma(alpha2)*alpha2 - 1))*c(alpha2, beta2, beta2, trigamma(alpha2)*beta2^2), nrow = 2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_alpha1, d_beta1))%*%Fish1inv%*%c(d_alpha1, d_beta1)
  
  avar2 <- t(c(d_alpha2, d_beta2))%*%Fish2inv%*%c(d_alpha2, d_beta2)/q
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2
  
  ## return initial upper bound for the root finding algorithm
  mu_start <- (4*qnorm(1 - (1 - conviction)/2)*sqrt(Fish_ratio_mu)/int_len)^2
  
  ## generate a Sobol' sequence to find asymptotic power curve
  sob_pwr <- sobol(m, d = 4, randomize = "digital.shift")
  
  ## helper function used to compute an integral in the partial derivative
  ## of theta_j with respect to alpha_j
  fn_int <- function(x, alpha){
    x^(alpha-1)*exp(-x)*log(x)/gamma(alpha)
  }
  
  ## function to input into simplified uniroot function
  ## u is point from Sobol' sequence, targetP is the target power
  ## for the root-finding algorithm to obtain, params are the design
  ## values, deltas is the interval c(delta_L, delta_U), n_val is the 
  ## sample size presently explored by the root-finding algorithm, and 
  ## q is the constant for imbalanced sample size determination
  targetPower <- function(u, targetP, params, deltas, n_val, q){
    
    ## return negative power if sample size is not positive
    if (n_val <= 0){return(-1.5)}
    
    gamma_alpha.1 <- params[1]
    gamma_beta.1 <- params[2]
    gamma_alpha.2 <- params[3]
    gamma_beta.2 <- params[4]
    threshold <- params[5]
    
    ## generate approximately normal MLEs for group 1 using delta method
    rho1 <- 1/sqrt(gamma_alpha.1*trigamma(gamma_alpha.1))
    mat1 <- matrix(c(1/gamma_alpha.1, 1/gamma_alpha.1, 1/gamma_alpha.1, 
                     trigamma(gamma_alpha.1))/(trigamma(gamma_alpha.1)*gamma_alpha.1 - 1), nrow = 2)
    a1 <- qnorm(u[1], log(gamma_alpha.1), sqrt(mat1[1,1]/n_val))
    b1 <- qnorm(u[2], log(gamma_beta.1) + rho1*(a1 - log(gamma_alpha.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    ## generate approximately normal MLEs for group 2 using delta method
    rho2 <- 1/sqrt(gamma_alpha.2*trigamma(gamma_alpha.2))
    mat2 <- matrix(c(1/gamma_alpha.2, 1/gamma_alpha.2, 1/gamma_alpha.2, 
                     trigamma(gamma_alpha.2))/(trigamma(gamma_alpha.2)*gamma_alpha.2 - 1), nrow = 2)
    a2 <- qnorm(u[3], log(gamma_alpha.2), sqrt(mat2[1,1]/(q*n_val)))
    b2 <- qnorm(u[4], log(gamma_beta.2) + rho2*(a2 - log(gamma_alpha.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/(q*n_val)))
    
    ## exponentiate MLEs
    a1 <- exp(a1); b1 <- exp(b1); a2 <- exp(a2); b2 <- exp(b2)
    
    ## ensure no MLEs underflow to 0
    if(max(a1 <= 0, b1 <= 0, a2 <= 0, b2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pgamma(threshold, a1, b1)
    theta2 <- 1 - pgamma(threshold, a2, b2)
    
    ## compute partial derivatives of theta with respect to alpha and beta for each group
    int_1 <- tryCatch(integrate(fn_int,
                                lower = 0, upper = threshold*b1, alpha = a1)$value,
                      error = function(e) {-10001})
    
    ## say power criterion is not satisfied if the integral does not converge
    if (int_1 < -10000){return(-1.5)}
    
    d_alpha1 <- digamma(a1)*pgamma(threshold,a1, b1) - int_1
    
    d_beta1 <- (a1/b1)*(pgamma(threshold, a1+1, b1) - pgamma(threshold, a1, b1))
    
    int_2 <- tryCatch(integrate(fn_int,
                                lower = 0, upper = threshold*b2, alpha = a2)$value,
                      error = function(e) {-10001})
    
    if (int_2 < -10000){return(-1.5)}
    
    d_alpha2 <- digamma(a2)*pgamma(threshold, a2, b2) - int_2
    
    d_beta2 <- (a2/b2)*(pgamma(threshold, a2+1, b2) - pgamma(threshold, a2,b2))
    
    ## compute Fisher information for each gamma model
    Fish1inv <- matrix((1/(trigamma(a1)*a1 - 1))*c(a1, b1, b1, trigamma(a1)*b1^2), nrow = 2)
    Fish2inv <-  matrix((1/(trigamma(a2)*a2 - 1))*c(a2, b2, b2, trigamma(a2)*b2^2), nrow = 2)
    
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_alpha1, d_beta1))%*%Fish1inv%*%c(d_alpha1, d_beta1)
    
    avar2 <- t(c(d_alpha2, d_beta2))%*%Fish2inv%*%c(d_alpha2, d_beta2)/q
    
    ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
    Fish_ratio_mu <- tryCatch(avar1/theta1^2 + avar2/theta2^2, error = function(e) {-10001}) 
    
    ## return negative power if division causes NA/Inf values
    if (is.na(Fish_ratio_mu)){return(-1.5)}
    if (Fish_ratio_mu < -10000){return(-1.5)}
    
    ## return power based on normal approximation induced by Bernstein-von Mises
    realP <- pnorm(deltas[2], log(theta1/theta2), sqrt(Fish_ratio_mu)/sqrt(n_val)) - pnorm(deltas[1], log(theta1/theta2), sqrt(Fish_ratio_mu)/sqrt(n_val))
    
    ## return estimated power less target power (for root-finding algorithm)
    return(realP - targetP)
  }
  
  ## uniroot function without the error checking
  uu <- function(f, lower, upper, tol = 1e-4, maxiter =1000L, ...) {
    f.lower <- f(lower, ...)
    f.upper <- f(upper, ...)
    val <- .External2(stats:::C_zeroin2, function(arg) f(arg, ...),
                      lower, upper, f.lower, f.upper, tol, as.integer(maxiter))
    return(val[1])
  }
  
  ## implement the root-finding algorithm with suitable endpoints. We compute
  ## power at several sample sizes (i.e., power1, power2, ...) to ensure that 
  ## the solution for n for each Sobol' sequence point is between the endpoints
  ## chosen for the root-finding algorithm. We use a target power of 0 here (that
  ## is subtracted from the power at that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(0.5*mu_start),
                          params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(0.5*mu_start), targetP = conviction,
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_start),
                            params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(0.5*mu_start), upper = ceiling(mu_start), targetP = conviction, 
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        power3 <- targetPower(targetP = 0, n_val = ceiling(2*mu_start),
                              params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                              deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        if (power3 >= conviction){
          samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(2*mu_start), targetP = conviction, 
                             params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        }
        else{
          endpoints_vec[i] <- 1 
        }
      }
    }
  }
  
  ## find larger upper bound if the initial one is not sufficient
  ## (i.e., if some points still do not satisfy power criterion at 2*mu_start)
  last_group <- which(endpoints_vec == 1)
  if (length(last_group) == 0){
    upper_c <- 2
  } else{
    upper_c <- 2
    while(length(last_group) > 0){
      if (upper_c > 32){
        last_group <- NULL
      }
      upper_c <- 2*upper_c
      endpoints1_vec <- NULL
      for (i in 1:length(last_group)){
        endpoints1_vec[i] <- targetPower(targetP = 0, n_val = ceiling(upper_c*mu_start),
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      keep_vec <- ifelse(endpoints1_vec >= conviction, FALSE, TRUE)
      ## only keep points that still do not satisfy power criterion after increasing
      ## the upper bound for the sample size
      last_group <- last_group[keep_vec]
    }
  }
  
  ## implement the root-finding algorithm for each point in the Sobol' sequence
  ## that required a large upper bound (i.e., those in last_group)
  for (i in 1:nrow(sob_pwr)){
    if (endpoints_vec[i] == 1){
      samps_pwr[i] <- uu(targetPower, lower =ceiling(2*mu_start), upper = ceiling(upper_c*mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at recommended sample size and one
  ## sample size smaller. This can be commented out and the return statement
  ## could be modified to output only samps_pwr if desired.
  n_star <- quantile(samps_pwr, power)
  power_ceil <- NULL
  power_floor <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_ceil[i] <- targetPower(targetP = 0, n_val = ceiling(n_star),
                                 params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                 deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    power_floor[i] <- targetPower(targetP = 0, n_val = floor(n_star),
                                  params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                  deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
  }
  
  return(list(samps_pwr, c(n_star, mean(power_ceil >= conviction), mean(power_floor >= conviction))))
  
}

## function to approximate power curve using approximation method in Algorithm 2
powerCurveGamma2 <- function(conviction, power, delta_L, delta_U, gamma_alpha.1, gamma_beta.1, gamma_alpha.2, 
                             gamma_beta.2,  mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, mu0.2 = 2, 
                             tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25, threshold = NULL, m = 1024, q = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for an equivalence or noninferiority test
  ## delta_U: upper interval endpoint for an equivalence or noninferiority test
  ## gamma_alpha.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## gamma_beta.j: the design value for the rate parameter of the gamma distribution for group j = 1, 2
  ## alpha_j has a Gamma(mu0.j, tau0.j) prior, where tau0.j is a rate for group j = 1, 2
  ## beta_j has a Gamma(kappa0.j, lambda0.j) prior, where lambda0.j is a rate for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  
  
  ## concatenate all hyperparameters in matrix form (used in targetPower() function)
  hyper_mat <- rbind(c(mu0.1, tau0.1), c(kappa0.1, lambda0.1), c(mu0.2, tau0.2), c(kappa0.2, lambda0.2))
  
  ## find good upper bound for the sample size to generate the asymptotic power curve
  if (is.finite(delta_L) & is.finite(delta_U)){
    int_len <- delta_U - delta_L
  }
  else if (is.finite(delta_L)){
    int_len <- 2*(1- delta_L)    
  }
  else{
    int_len <- 2*(delta_U - 1)
  }
  
  alpha1 <- gamma_alpha.1; alpha2 <- gamma_alpha.2; beta1 <- gamma_beta.1; beta2 <- gamma_beta.2
  ## define theta metrics in terms of design values
  theta1 <- 1 - pgamma(threshold, alpha1, beta1)
  theta2 <- 1 - pgamma(threshold, alpha2, beta2)
  
  ## compute partial derivatives of theta with respect to alpha and beta for each group
  int_1 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha1, beta = beta1)$value
  
  d_alpha1 <- digamma(alpha1)*pgamma(threshold,alpha1, beta1) - int_1
  
  d_beta1 <- (alpha1/beta1)*(pgamma(threshold, alpha1+1, beta1) -  pgamma(threshold, alpha1, beta1))
  
  int_2 <- integrate(function(alpha, beta, x){dgamma(x, alpha, beta)*log(beta*x)},
                     lower = 0, upper = threshold, alpha = alpha2, beta = beta2)$value
  
  d_alpha2 <- digamma(alpha2)*pgamma(threshold, alpha2, beta2) - int_2
  
  d_beta2 <- (alpha2/beta2)*(pgamma(threshold, alpha2+1, beta2) - pgamma(threshold, alpha2,beta2))
  
  ## compute Fisher information for each gamma model
  Fish1inv <- matrix((1/(trigamma(alpha1)*alpha1 - 1))*c(alpha1, beta1, beta1, trigamma(alpha1)*beta1^2), nrow = 2)
  Fish2inv <-  matrix((1/(trigamma(alpha2)*alpha2 - 1))*c(alpha2, beta2, beta2, trigamma(alpha2)*beta2^2), nrow = 2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_alpha1, d_beta1))%*%Fish1inv%*%c(d_alpha1, d_beta1)
  
  avar2 <- t(c(d_alpha2, d_beta2))%*%Fish2inv%*%c(d_alpha2, d_beta2)/q
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2
  
  ## return initial upper bound for the root finding algorithm
  mu_start <- (4*qnorm(1 - (1 - conviction)/2)*sqrt(Fish_ratio_mu)/int_len)^2
  
  ## generate a Sobol' sequence to find asymptotic power curve
  sob_pwr <- sobol(m, d = 4, randomize = "digital.shift")
  
  ## function to solve for the posterior mode of each group
  ## x is the point (logalpha, logbeta) at which this function is evaluated at
  ## y = c(length(data), sum(log(data)), sum(data))
  fn_grad <- function(x, y, mu, tau, kappa, lambda) {
    
    res1 <- exp(x[1])*y[1]*x[2] - y[1]*digamma(exp(x[1]))*exp(x[1]) + exp(x[1])*(y[2] - tau) + mu
    
    res2 <- y[1]*exp(x[1]) - (y[3] + lambda)*exp(x[2]) + kappa
    
    return(c(res1, res2))
  }
  
  ## function to calculate the covariance matrix for the
  ## approximate posterior using Laplace's method. u is log(alpha_j), and
  ## v is log(beta_j). yy_star is c(length(data), sum(log(data)), sum(data)),
  ## and hyper = c(tau_j, lambda_j) are the relevant hyperparameters.
  calc_covar <- function(u, v, yy_star, hyper){
    a <- exp(u); b <- exp(v); n <- yy_star[1]
    tau <- hyper[1]; lambda <- hyper[2]
    d11 <- a*(n*digamma(a) + tau - n*v - yy_star[2]) + n*a^2*trigamma(a)
    d12 <- -1*a*n
    d22 <- b*(yy_star[3] + lambda)
    mat <- rbind(c(d11, d12), c(d12, d22))
    return(solve(mat))
  }
  
  ## helper function used to compute an integral in the partial derivative
  ## of theta_j with respect to alpha_j
  fn_int <- function(x, alpha){
    x^(alpha-1)*exp(-x)*log(x)/gamma(alpha)
  }
  
  ## function to input into simplified uniroot function
  ## u is point from Sobol' sequence, targetP is the target power
  ## for the root-finding algorithm to obtain, params are the design
  ## values, deltas is the interval c(delta_L, delta_U), n_val is the 
  ## sample size presently explored by the root-finding algorithm, hyper
  ## is a matrix of the hyperparameters, and q is the constant for 
  ## imbalanced sample size determination
  targetPower <- function(u, targetP, params, deltas, n_val, hyper, q){
    
    ## return negative power if sample size is not positive
    if (n_val <= 0){return(-1.5)}
    
    gamma_alpha.1 <- params[1]
    gamma_beta.1 <- params[2]
    gamma_alpha.2 <- params[3]
    gamma_beta.2 <- params[4]
    threshold <- params[5]
    
    ## generate approximately normal MLEs for group 1 using delta method
    rho1 <- 1/sqrt(gamma_alpha.1*trigamma(gamma_alpha.1))
    mat1 <- matrix(c(1/gamma_alpha.1, 1/gamma_alpha.1, 1/gamma_alpha.1, 
                     trigamma(gamma_alpha.1))/(trigamma(gamma_alpha.1)*gamma_alpha.1 - 1), nrow = 2)
    a1 <- qnorm(u[1], log(gamma_alpha.1), sqrt(mat1[1,1]/n_val))
    b1 <- qnorm(u[2], log(gamma_beta.1) + rho1*(a1 - log(gamma_alpha.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    ## generate approximately normal MLEs for group 2 using delta method
    rho2 <- 1/sqrt(gamma_alpha.2*trigamma(gamma_alpha.2))
    mat2 <- matrix(c(1/gamma_alpha.2, 1/gamma_alpha.2, 1/gamma_alpha.2, 
                     trigamma(gamma_alpha.2))/(trigamma(gamma_alpha.2)*gamma_alpha.2 - 1), nrow = 2)
    a2 <- qnorm(u[3], log(gamma_alpha.2), sqrt(mat2[1,1]/(q*n_val)))
    b2 <- qnorm(u[4], log(gamma_beta.2) + rho2*(a2 - log(gamma_alpha.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/(q*n_val)))
    
    ## exponentiate MLEs
    a1 <- exp(a1); b1 <- exp(b1); a2 <- exp(a2); b2 <- exp(b2)
    
    gamma_alpha.1 <- a1; gamma_beta.1 <- b1; gamma_alpha.2 <- a2; gamma_beta.2 <- b2
    
    ## summarize information from first group of data (faster computation)
    yy_star1 <- c(n_val, n_val*(digamma(gamma_alpha.1) - log(gamma_beta.1)), n_val*gamma_alpha.1/gamma_beta.1)
    ## find posterior modes for the first group (logalpha and logbeta)
    modes <- nleqslv(log(c(gamma_alpha.1, gamma_beta.1)), fn_grad, y = yy_star1, mu = hyper[1,1], tau = hyper[1,2],
                     kappa = hyper[2,1], lambda = hyper[2,2] )$x
    
    
    mat1_new <- calc_covar(modes[1], modes[2], yy_star1, c(hyper[1,2], hyper[2,2]))
    ## exponentiate modes to return to standard scale
    modes1 <- exp(modes)
    
    ## repeat all steps for the second group
    yy_star2 <- c(q*n_val, q*n_val*(digamma(gamma_alpha.2) - log(gamma_beta.2)), q*n_val*gamma_alpha.2/gamma_beta.2)
    modes <- nleqslv(log(c(gamma_alpha.2, gamma_beta.2)), fn_grad, y = yy_star2, mu = hyper[3,1], tau = hyper[3,2],
                     kappa = hyper[4,1], lambda = hyper[4,2] )$x
    
    mat2_new <- calc_covar(modes[1], modes[2], yy_star2, c(hyper[3,2], hyper[4,2]))
    modes2 <- exp(modes)
    a1 <- modes1[1]; b1 <- modes1[2]; a2 <- modes2[1]; b2 <- modes2[2]
    
    ## ensure no modes are 0 due to underflow errors
    if(max(a1 <= 0, b1 <= 0, a2 <= 0, b2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pgamma(threshold, a1, b1)
    theta2 <- 1 - pgamma(threshold, a2, b2)
    
    ## compute partial derivatives of theta with respect to logalpha and logbeta for each group
    ## this is different from the previous function that computes the derivatives with respect
    ## to alpha and beta
    int_1 <- tryCatch(integrate(fn_int,
                                lower = 0, upper = threshold*b1, alpha = a1)$value,
                      error = function(e) {-10001})
    
    if (int_1 < -10000){return(-1.5)}
    
    d_alpha1 <- (digamma(a1)*pgamma(threshold,a1, b1) - int_1)*a1
    
    d_beta1 <- ((a1/b1)*(pgamma(threshold, a1+1, b1) - pgamma(threshold, a1, b1)))*b1
    
    int_2 <- tryCatch(integrate(fn_int,
                                lower = 0, upper = threshold*b2, alpha = a2)$value,
                      error = function(e) {-10001})
    
    if (int_2 < -10000){return(-1.5)}
    
    d_alpha2 <- (digamma(a2)*pgamma(threshold, a2, b2) - int_2)*a2
    
    d_beta2 <- ((a2/b2)*(pgamma(threshold, a2+1, b2) - pgamma(threshold, a2,b2)))*b2
    
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_alpha1, d_beta1))%*%mat1_new%*%c(d_alpha1, d_beta1)
    
    avar2 <- t(c(d_alpha2, d_beta2))%*%mat2_new%*%c(d_alpha2, d_beta2)
    
    ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
    Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2 
    
    ## return negative power if division causes NA/Inf values
    if (is.na(Fish_ratio_mu)){return(-1.5)}
    if (Fish_ratio_mu < -10000){return(-1.5)}
    
    ## return power based on normal approximation induced by Bernstein-von Mises
    realP <- pnorm(deltas[2], log(theta1/theta2), sqrt(Fish_ratio_mu)) - pnorm(deltas[1], log(theta1/theta2), sqrt(Fish_ratio_mu))
    
    ## return estimated power less target power (for root-finding algorithm)
    return(realP - targetP)
  }
  
  ## uniroot function without the error checking
  uu <- function(f, lower, upper, tol = 1e-4, maxiter =1000L, ...) {
    f.lower <- f(lower, ...)
    f.upper <- f(upper, ...)
    val <- .External2(stats:::C_zeroin2, function(arg) f(arg, ...),
                      lower, upper, f.lower, f.upper, tol, as.integer(maxiter))
    return(val[1])
  }
  
  ## implement the root-finding algorithm with suitable endpoints. We compute
  ## power at several sample sizes (i.e., power1, power2, ...) to ensure that 
  ## the solution for n for each Sobol' sequence point is between the endpoints
  ## chosen for the root-finding algorithm. We use a target power of 0 here (that
  ## is subtracted from the power at that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(0.5*mu_start), hyper = hyper_mat,
                          params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(0.5*mu_start), targetP = conviction, hyper = hyper_mat,
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_start), hyper = hyper_mat,
                            params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(0.5*mu_start), upper = ceiling(mu_start), targetP = conviction, 
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        power3 <- targetPower(targetP = 0, n_val = ceiling(2*mu_start), hyper = hyper_mat,
                              params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                              deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        if (power3 >= conviction){
          samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(2*mu_start), targetP = conviction, 
                             params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        }
        else{
          endpoints_vec[i] <- 1 
        }
      }
    }
  }
  
  ## find larger upper bound if the initial one is not sufficient
  ## (i.e., if some points still do not satisfy power criterion at 2*mu_start)
  last_group <- which(endpoints_vec == 1)
  if (length(last_group) == 0){
    upper_c <- 2
  } else{
    upper_c <- 2
    while(length(last_group) > 0){
      if (upper_c > 32){
        last_group <- NULL
      }
      upper_c <- 2*upper_c
      endpoints1_vec <- NULL
      for (i in 1:length(last_group)){
        endpoints1_vec[i] <- targetPower(targetP = 0, n_val = ceiling(upper_c*mu_start), hyper = hyper_mat,
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      keep_vec <- ifelse(endpoints1_vec >= conviction, FALSE, TRUE)
      ## only keep points that still do not satisfy power criterion after increasing
      ## the upper bound for the sample size
      last_group <- last_group[keep_vec]
    }
  }
  
  ## implement the root-finding algorithm for each point in the Sobol' sequence
  ## that required a large upper bound (i.e., those in last_group)
  for (i in 1:nrow(sob_pwr)){
    if (endpoints_vec[i] == 1){
      samps_pwr[i] <- uu(targetPower, lower =ceiling(2*mu_start), upper = ceiling(upper_c*mu_start), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at recommended sample size and one
  ## sample size smaller. This can be commented out and the return statement
  ## could be modified to output only samps_pwr if desired.
  n_star <- quantile(samps_pwr, power)
  power_ceil <- NULL
  power_floor <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_ceil[i] <- targetPower(targetP = 0, n_val = ceiling(n_star), hyper = hyper_mat,
                                 params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                 deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    power_floor[i] <- targetPower(targetP = 0, n_val = floor(n_star), hyper = hyper_mat,
                                  params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                  deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
  }
  
  return(list(samps_pwr, c(n_star, mean(power_ceil >= conviction), mean(power_floor >= conviction))))
}

## read informative hyperparameters from .csv file for gamma distribution
informs <- read.csv("informs_gamma.csv")

## define the conviction threshold, target power, and delta_* for the 
## three scenarios as in the paper
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
equivalences <- c(0.25, 0.3, 0.15)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## Algorithm 1: Uninformative Priors for Figure 3
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                   .options.snow=opts) %dopar% {
                     unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                               delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                               gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024))
                   }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors for Figure 3
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors for Figure 3
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure 3
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Analysis with Bayes Factors in Section C.1 of the supplement (uses same functions)
## Only consider first setting; the conviction threshold must be redefined.
convictions <- c(0.5652)
powers <- c(0.6)
equivalences <- c(0.25)

## Algorithm 1: Uninformative Priors for Figure C.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bf_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bf_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors for Figure C.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bf_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bf_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## conviction threshold must be redefined for informative priors.
convictions <- c(0.5428)
## Algorithm 2: Uninformative Priors for Figure C.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bf_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bf_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure C.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bf_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bf_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Imbalanced Sample Size Determination for Section C.2 of the supplement
## conviction threshold set back to value for numerical study in the main paper
convictions <- c(0.5)

## Algorithm 1: Uninformative Priors for Figure C.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_q2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_q2_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors for Figure C.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_q2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_q2_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors for Figure C.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_q2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_q2_bpc_gamma_root_finding_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure C.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = log(1/(1 + equivalences[k])),
                                                delta_U = log(1 + equivalences[k]), gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_q2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_q2_bpc_gamma_root_finding_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

