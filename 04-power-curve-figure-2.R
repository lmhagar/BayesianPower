## code to approximate power curves for hypothesis tests with the gamma model
## this procedure leverages Algorithm 1 and Algorithm 2 in two separate functions

## these functions were used to produce Figure 2 of the main text along with 
## Figures D.1 and D.2 in the supplement

## load required libraries
require(foreach)
require(doParallel)
require(doSNOW)
require(numDeriv)
require(qrng)
require(nleqslv)
require(ggplot2)
require(cowplot)
require(ggpubr)
require(rjags)
require(coda)

## function to approximate power curve using approximation method in Algorithm 1
powerCurveGamma1 <- function(conviction, power, delta_L, delta_U, gamma_alpha.1, gamma_beta.1, gamma_alpha.2, 
                             gamma_beta.2, threshold = NULL, m = 1024, q = 1, seed = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for hypothesis test
  ## delta_U: upper interval endpoint for hypothesis test
  ## gamma_alpha.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## gamma_beta.j: the design value for the rate parameter of the gamma distribution for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  ## seed: seed used to generated randomized Sobol' sequence for reproducibility (default is 1)
  
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
  Fish_ratio_mu <- as.numeric(avar1/theta1^2 + avar2/theta2^2)
  
  ## uniroot without error checking
  uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.01, ...)
  {
    f <- function(x) fun(x, ...)
    x1 <- lower
    if (!is.null(f_lower)){
      f1 <- f_lower
    }
    else{
      f1 <- f(x1)
    }
    if (f1 > 0){return(x1)}
    x2 <- upper
    if (!is.null(f_upper)){
      f2 <- f_upper
    }
    else{
      f2 <- f(x2)
    }
    f2 <- f(x2)
    if (f2 < 0){return(x2)}
    x3 <- 0.5 * (lower + upper)
    niter <- 1
    while (niter <= maxiter) {
      f3 <- f(x3)
      if (abs(f3) < tol) {
        x0 <- x3
        return(x0)
      }
      if (f1 * f3 < 0) {
        upper <- x3}
      else {lower <- x3}
      if ((upper - lower) < tol2 * max(abs(upper), 1)) {
        x0 <- 0.5 * (lower + upper)
        return(x0)
      }
      denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
      numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
        (f2 - f3) + f1 * x2 * (f3 - f1)
      if (denom == 0) {
        dx <- upper - lower
      }
      else {
        dx <- f3 * numer/denom
      }
      x <- x3 + dx
      if ((upper - x) * (x - lower) < 0) {
        dx <- 0.5 * (upper - lower)
        x <- lower + dx
      }
      if (x1 < x3) {
        x2 <- x3
        f2 <- f3
      }
      else {
        x1 <- x3
        f1 <- f3
      }
      niter <- niter + 1
      if (abs(x - x3) < tol2) {
        x0 <- x
        return(x0)
      }
      x3 <- x
    }
    return(x0)
  }
  
  ## return initial upper bound for the root finding algorithm
  if (!is.finite(delta_U)){
    mu_start <- ((qnorm(power) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2
    mu_low <- min(((qnorm(1 - conviction + 0.05) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2, 0.5*mu_start)
    mu_high <- max(((qnorm((3 + power)/4) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2, 2*mu_start)
  }
  else if (!is.finite(delta_L)){
    mu_start <- ((qnorm(power) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2))))^2
    mu_low <- min((qnorm(1 - conviction + 0.05) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2)))^2, 0.5*mu_start)
    mu_high <- max((qnorm((3 + power)/4) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2)))^2, 2*mu_start)
  }
  else{
    ## find more conservative upper bound using criteria for credible intervals
    theta0 <- log(theta1) - log(theta2)
    a_cons <- (delta_U - theta0)/sqrt(Fish_ratio_mu)
    b_cons <- (delta_L - theta0)/sqrt(Fish_ratio_mu)
    c_cons <- qnorm((1-conviction)/2)
    ## lower bound for root-finding algorithm
    lower_cons <- -2*c_cons*sqrt(Fish_ratio_mu)/(delta_U - delta_L)
    upper_cons <- lower_cons
    
    fn_ci = function(n_sq, a, b, c, pwr){
      return(pnorm(a*n_sq + c) - pnorm(b*n_sq - c) - pwr)}
    
    upper_large <- FALSE
    while(upper_large == FALSE){
      upper_cons <- 10*upper_cons
      upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = power)
      if (upper_check > 0){
        upper_large <- TRUE
      }
    }
    
    upper_n <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = power, 
                   lower = lower_cons, upper = upper_cons))^2

    fn_start = function(n_trans, a, b, u, gam)
    {
      return(gam + pnorm(b, n_trans*qnorm(u), n_trans) - pnorm(a, n_trans*qnorm(u), n_trans))
    }
    uuu <- as.numeric(sobol(128, d = 1, randomize = "digital.shift", seed = seed))
    start_rep <- NULL
    for (i in 1:length(uuu)){
      check_upper_n <- fn_start(n_trans = sqrt(Fish_ratio_mu)/sqrt(upper_n), 
                                a = delta_U - theta0, u = uuu[i],
                                b = delta_L - theta0, gam = conviction)
      if (check_upper_n > 0){
        start_rep[i] <- upper_n
      }
      else{
        start_rep[i] <- (uu(fn_start, a = delta_U - theta0, u = uuu[i],
                        b = delta_L - theta0, gam = conviction, 
                        lower = (sqrt(Fish_ratio_mu)/sqrt(upper_n)), upper = 1000)/sqrt(Fish_ratio_mu))^(-2)
      }
    }
    mu_start <- quantile(start_rep, power)
    mu_low <- min(0.5*mu_start, quantile(start_rep, 1 - conviction + 0.05))
    mu_high <- max(2*mu_start, quantile(start_rep, (3 + power)/4))
  }
  
  ## generate a Sobol' sequence to find asymptotic power curve
  sob_pwr <- sobol(m, d = 4, randomize = "digital.shift", seed = seed)
  
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
    
    if(max(theta1 <= 0, theta1 >= 1, theta2 <= 0, theta2 >= 1)){return(-1.5)}
    
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
  
  ## implement the root-finding algorithm with suitable endpoints. We compute
  ## power at several sample sizes (i.e., power1, power2, ...) to ensure that 
  ## the solution for n for each Sobol' sequence point is between the endpoints
  ## chosen for the root-finding algorithm. We use a target power of 0 here (that
  ## is subtracted from the power at that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(mu_start),
                          params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      power1b <- targetPower(targetP = 0, n_val = ceiling(mu_low),
                             params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power1b >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(mu_low),  
                           f_upper = power1b - conviction, targetP = conviction,
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_low), upper = ceiling(mu_start), 
                           f_lower = power1b - conviction, f_upper = power1 - conviction, targetP = conviction,
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_high),
                            params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(mu_high), 
                           f_lower = power1 - conviction, f_upper = power2 - conviction, targetP = conviction, 
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        endpoints_vec[i] <- 1 
      }
    }
  }
  
  ## find larger upper bound if the initial one is not sufficient
  ## (i.e., if some points still do not satisfy power criterion at 2*mu_start)
  last_group <- which(endpoints_vec == 1)
  if (length(last_group) == 0){
    upper_c <- 1
  } else{
    upper_c <- 1
    while(length(last_group) > 0){
      if (upper_c > 32){
        last_group <- NULL
      }
      upper_c <- 2*upper_c
      endpoints1_vec <- NULL
      for (i in 1:length(last_group)){
        endpoints1_vec[i] <- targetPower(targetP = 0, n_val = ceiling(upper_c*mu_high),
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[last_group[i],], q = q)
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
      samps_pwr[i] <- uu(targetPower, lower = ceiling(mu_high), upper = ceiling(upper_c*mu_high), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at n_* (noninteger)
  n_star <- quantile(samps_pwr, power)
  power_star <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_star[i] <- targetPower(targetP = 0, n_val = n_star,
                                 params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                 deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
  }
  
  ## confirm that posterior probabilities are less than or at least gamma as expected based on root-finding
  ## consider all probabilities less than 0.001 from the target power as satisfying criterion
  consistency <- ifelse(samps_pwr < n_star, round(power_star,3) >= conviction, round(power_star,3) <= conviction)
  consistency <- ifelse(consistency == 1, 1, as.numeric(abs(power_star - conviction) < 0.001))
  power_star_round <- ifelse(consistency == 1, as.numeric(samps_pwr < n_star), as.numeric(samps_pwr > n_star))
  
  ## for any points where the root-finding algorithm has caused issues, 
  ## re-run the root-finding algorithm starting at n_*
  inconsistent <- which(consistency == 0)
  samps_pre <- NULL
  samps_post <- NULL
  if (length(inconsistent) > 0){
    for (i in 1:length(inconsistent)){
      samps_pre[i] <- samps_pwr[inconsistent[i]] 
      if (power_star[inconsistent[i]] < conviction){
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = n_star, upper = ceiling(upper_c*mu_start), 
                                         f_lower = power_star[inconsistent[i]] - conviction, targetP = conviction, 
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q)
      }
      else{
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = 2, 
                                         upper = n_star, f_upper = power_star[inconsistent[i]] - conviction, targetP = conviction, 
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q)
      }
      samps_post[i] <- samps_pwr[inconsistent[i]] 
    }
    
  }
  n_star <- quantile(samps_pwr, power)
  
  adjust <- c(samps_pre, samps_post, samps_pre - samps_post)
  adjust <- c(adjust, rep(0, 150 - length(adjust)))
  
  ## the first part are the sample sizes used to approximate the power curve; the last two parts
  ## are just to investigate how often the root-finding algorithm presents issues
  return(list(samps_pwr, c(n_star, mean(power_star_round >= conviction), mean(consistency)),
              adjust))
  
}

## function to approximate power curve using approximation method in Algorithm 2
powerCurveGamma2 <- function(conviction, power, delta_L, delta_U, gamma_alpha.1, gamma_beta.1, gamma_alpha.2, 
                             gamma_beta.2,  mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, mu0.2 = 2, 
                             tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25, threshold = NULL, m = 1024, q = 1,
                             seed = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for hypothesis test
  ## delta_U: upper interval endpoint for hypothesis test
  ## gamma_alpha.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## gamma_beta.j: the design value for the rate parameter of the gamma distribution for group j = 1, 2
  ## alpha_j has a Gamma(mu0.j, tau0.j) prior, where tau0.j is a rate for group j = 1, 2
  ## beta_j has a Gamma(kappa0.j, lambda0.j) prior, where lambda0.j is a rate for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  ## seed: seed used to generated randomized Sobol' sequence for reproducibility (default is 1)
  
  
  ## concatenate all hyperparameters in matrix form (used in targetPower() function)
  hyper_mat <- rbind(c(mu0.1, tau0.1), c(kappa0.1, lambda0.1), c(mu0.2, tau0.2), c(kappa0.2, lambda0.2))
  
  ## find good upper bound for the sample size to generate the asymptotic power curve
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
  Fish_ratio_mu <- as.numeric(avar1/theta1^2 + avar2/theta2^2)
  
  ## uniroot without error checking
  uu <- function (fun, lower, upper, f_lower = NULL, f_upper = NULL, maxiter = 1000, tol = 1e-4, tol2 = 0.01, ...)
  {
    f <- function(x) fun(x, ...)
    x1 <- lower
    if (!is.null(f_lower)){
      f1 <- f_lower
    }
    else{
      f1 <- f(x1)
    }
    if (f1 > 0){return(x1)}
    x2 <- upper
    if (!is.null(f_upper)){
      f2 <- f_upper
    }
    else{
      f2 <- f(x2)
    }
    f2 <- f(x2)
    if (f2 < 0){return(x2)}
    x3 <- 0.5 * (lower + upper)
    niter <- 1
    while (niter <= maxiter) {
      f3 <- f(x3)
      if (abs(f3) < tol) {
        x0 <- x3
        return(x0)
      }
      if (f1 * f3 < 0) {
        upper <- x3}
      else {lower <- x3}
      if ((upper - lower) < tol2 * max(abs(upper), 1)) {
        x0 <- 0.5 * (lower + upper)
        return(x0)
      }
      denom <- (f2 - f1) * (f3 - f1) * (f2 - f3)
      numer <- x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 *
        (f2 - f3) + f1 * x2 * (f3 - f1)
      if (denom == 0) {
        dx <- upper - lower
      }
      else {
        dx <- f3 * numer/denom
      }
      x <- x3 + dx
      if ((upper - x) * (x - lower) < 0) {
        dx <- 0.5 * (upper - lower)
        x <- lower + dx
      }
      if (x1 < x3) {
        x2 <- x3
        f2 <- f3
      }
      else {
        x1 <- x3
        f1 <- f3
      }
      niter <- niter + 1
      if (abs(x - x3) < tol2) {
        x0 <- x
        return(x0)
      }
      x3 <- x
    }
    return(x0)
  }
  
  ## return initial upper bound for the root finding algorithm
  if (!is.finite(delta_U)){
    mu_start <- ((qnorm(power) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2
    mu_low <- min(((qnorm(1 - conviction + 0.05) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2, 0.5*mu_start)
    mu_high <- max(((qnorm((3 + power)/4) + qnorm(conviction))*sqrt(Fish_ratio_mu)/((log(theta1) - log(theta2)) - delta_L))^2, 2*mu_start)
  }
  else if (!is.finite(delta_L)){
    mu_start <- ((qnorm(power) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2))))^2
    mu_low <- min((qnorm(1 - conviction + 0.05) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2)))^2, 0.5*mu_start)
    mu_high <- max((qnorm((3 + power)/4) + qnorm(conviction))*sqrt(Fish_ratio_mu)/(delta_U - (log(theta1) - log(theta2)))^2, 2*mu_start)
  }
  else{
    ## find more conservative upper bound using criteria for credible intervals
    theta0 <- log(theta1) - log(theta2)
    a_cons <- (delta_U - theta0)/sqrt(Fish_ratio_mu)
    b_cons <- (delta_L - theta0)/sqrt(Fish_ratio_mu)
    c_cons <- qnorm((1-conviction)/2)
    ## lower bound for root-finding algorithm
    lower_cons <- -2*c_cons*sqrt(Fish_ratio_mu)/(delta_U - delta_L)
    upper_cons <- lower_cons
    
    fn_ci = function(n_sq, a, b, c, pwr){
      return(pnorm(a*n_sq + c) - pnorm(b*n_sq - c) - pwr)}
    
    upper_large <- FALSE
    while(upper_large == FALSE){
      upper_cons <- 10*upper_cons
      upper_check <- fn_ci(n_sq = sqrt(upper_cons), a = a_cons, b = b_cons, c = c_cons, pwr = power)
      if (upper_check > 0){
        upper_large <- TRUE
      }
    }
    
    upper_n <- (uu(fn_ci, a = a_cons, b = b_cons, c = c_cons, pwr = power, 
                   lower = lower_cons, upper = upper_cons))^2
    
    fn_start = function(n_trans, a, b, u, gam)
    {
      return(gam + pnorm(b, n_trans*qnorm(u), n_trans) - pnorm(a, n_trans*qnorm(u), n_trans))
    }
    uuu <- as.numeric(sobol(128, d = 1, randomize = "digital.shift", seed = seed))
    start_rep <- NULL
    for (i in 1:length(uuu)){
      check_upper_n <- fn_start(n_trans = sqrt(Fish_ratio_mu)/sqrt(upper_n), 
                                a = delta_U - theta0, u = uuu[i],
                                b = delta_L - theta0, gam = conviction)
      if (check_upper_n > 0){
        start_rep[i] <- upper_n
      }
      else{
        start_rep[i] <- (uu(fn_start, a = delta_U - theta0, u = uuu[i],
                            b = delta_L - theta0, gam = conviction, 
                            lower = (sqrt(Fish_ratio_mu)/sqrt(upper_n)), upper = 1000)/sqrt(Fish_ratio_mu))^(-2)
      }
    }
    mu_start <- quantile(start_rep, power)
    mu_low <- min(0.5*mu_start, quantile(start_rep, 1 - conviction + 0.05))
    mu_high <- max(2*mu_start, quantile(start_rep, (3 + power)/4))
  }
  
  ## generate a Sobol' sequence to find asymptotic power curve
  sob_pwr <- sobol(m, d = 4, randomize = "digital.shift", seed = seed)
  
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
    
    if(max(theta1 <= 0, theta1 >= 1, theta2 <= 0, theta2 >= 1)){return(-1.5)}
    
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
  
  ## implement the root-finding algorithm with suitable endpoints. We compute
  ## power at several sample sizes (i.e., power1, power2, ...) to ensure that 
  ## the solution for n for each Sobol' sequence point is between the endpoints
  ## chosen for the root-finding algorithm. We use a target power of 0 here (that
  ## is subtracted from the power at that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(mu_start), hyper = hyper_mat,
                          params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      power1b <- targetPower(targetP = 0, n_val = ceiling(mu_low), hyper = hyper_mat,
                             params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power1b >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(mu_low), f_upper = power1b - conviction,
                           targetP = conviction, hyper = hyper_mat,
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_low), upper = ceiling(mu_start), 
                           f_lower = power1b - conviction, f_upper = power1 - conviction, targetP = conviction,
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_high), hyper = hyper_mat,
                            params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(mu_high), 
                           f_lower = power1 - conviction, f_upper = power2 - conviction, targetP = conviction, 
                           params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        endpoints_vec[i] <- 1 
      }
    }
  }
  
  ## find larger upper bound if the initial one is not sufficient
  ## (i.e., if some points still do not satisfy power criterion at 2*mu_start)
  last_group <- which(endpoints_vec == 1)
  if (length(last_group) == 0){
    upper_c <- 1
  } else{
    upper_c <- 1
    while(length(last_group) > 0){
      if (upper_c > 32){
        last_group <- NULL
      }
      upper_c <- 2*upper_c
      endpoints1_vec <- NULL
      for (i in 1:length(last_group)){
        endpoints1_vec[i] <- targetPower(targetP = 0, n_val = ceiling(upper_c*mu_high), hyper = hyper_mat,
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), sob_pwr[last_group[i],], q = q)
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
      samps_pwr[i] <- uu(targetPower, lower = ceiling(mu_high), upper = ceiling(upper_c*mu_high), targetP = conviction, 
                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold), hyper = hyper_mat,
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at n_* (noninteger)
  n_star <- quantile(samps_pwr, power)
  power_star <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_star[i] <- targetPower(targetP = 0, n_val = n_star,
                                 params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                 deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q, hyper = hyper_mat)
  }
  ## confirm that posterior probabilities are less than or at least gamma as expected based on root-finding
  ## consider all probabilities less than 0.001 from the target power as satisfying criterion
  consistency <- ifelse(samps_pwr < n_star, round(power_star,3) >= conviction, round(power_star,3) <= conviction)
  consistency <- ifelse(consistency == 1, 1, as.numeric(abs(power_star - conviction) < 0.001))
  power_star_round <- ifelse(consistency == 1, as.numeric(samps_pwr < n_star), as.numeric(samps_pwr > n_star))
  
  ## for any points where the root-finding algorithm has caused issues, 
  ## re-run the root-finding algorithm starting at n_*
  inconsistent <- which(consistency == 0)
  samps_pre <- NULL
  samps_post <- NULL
  if (length(inconsistent) > 0){
    for (i in 1:length(inconsistent)){
      samps_pre[i] <- samps_pwr[inconsistent[i]]
      if (power_star[inconsistent[i]] < conviction){
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = n_star, upper = ceiling(upper_c*mu_start), 
                                         f_lower = power_star[inconsistent[i]] - conviction, targetP = conviction, 
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q, hyper = hyper_mat)
      }
      else{
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = 2, upper = n_star, 
                                         f_upper = power_star[inconsistent[i]] - conviction, targetP = conviction, 
                                         params = c(gamma_alpha.1, gamma_beta.1, gamma_alpha.2, gamma_beta.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q, hyper = hyper_mat)
      }
      samps_post[i] <- samps_pwr[inconsistent[i]]
    }
    
  }
  n_star <- quantile(samps_pwr, power)
  
  adjust <- c(samps_pre, samps_post, samps_pre - samps_post)
  adjust <- c(adjust, rep(0, 150 - length(adjust)))

  ## the first part are the sample sizes used to approximate the power curve; the last two parts
  ## are just to investigate how often the root-finding algorithm presents issues
  return(list(samps_pwr, c(n_star, mean(power_star_round >= conviction), mean(consistency)),
              adjust))
  
}

## read informative hyperparameters from .csv file for gamma distribution
informs <- read.csv("informs_gamma.csv")

## define the conviction threshold, target power, and delta_* for the 
## three scenarios as in the paper
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
delta_Ls <- c(log(1/(1+0.25)), log(1/(1+0.3)), log(1/(1+0.15)))
delta_Us <- c(log(1 + 0.25), Inf, log(1 + 0.15))

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## Algorithm 1: Uninformative Priors for Figure 2
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                   .options.snow=opts) %dopar% {
                     unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                               delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                               gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k-1)*100 + j))
                   }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors for Figure 2
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+2)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors for Figure 2
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+5)*100 + j,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure 2
tic <- Sys.time()
for (k in 1:3){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+8)*100 + j,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Analysis with Bayes Factors in Section D.1 of the supplement (uses same functions)
## Only consider first setting; the conviction threshold must be redefined.
convictions <- c(0.5652)
powers <- c(0.6)
delta_Ls <- c(log(1/(1 + 0.25)))
delta_Us <- c(log(1 + 0.25))

## Algorithm 1: Uninformative Priors for Figure D.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+11)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bf_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bf_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bf_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors for Figure D.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25, seed = (k+12)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bf_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bf_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bf_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## conviction threshold must be redefined for informative priors.
convictions <- c(0.5428)

## Algorithm 1: Informative Priors for Figure D.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+13)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bf_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bf_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bf_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure D.1
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, seed = (k+14)*100 + j,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bf_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bf_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bf_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Imbalanced Sample Size Determination for Section D.2 of the supplement
## conviction threshold set back to value for numerical study in the main paper
convictions <- c(0.5)

## Algorithm 1: Uninformative Priors for Figure D.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2, seed = (k+15)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_q2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_q2_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_q2_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors for Figure D.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2, seed = (k+16)*100 + j))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_q2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_q2_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_q2_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors for Figure D.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2,
                                                mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, seed = (k+17)*100 + j,
                                                mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_q2_bpc_gamma_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_q2_bpc_gamma_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_q2_bpc_gamma_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors for Figure D.2
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveGamma2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                delta_U = delta_Us[k], gamma_alpha.1 = 2.11, gamma_beta.1 = 0.69,
                                                gamma_alpha.2 = 2.43, gamma_beta.2 = 0.79, threshold = 4.29, m = 1024, q = 2, seed = (k+18)*100 + j,
                                                mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_q2_bpc_gamma_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_q2_bpc_gamma_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_q2_bpc_gamma_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

closeAllConnections()

## The next section of code estimates the power curves for Figure 2 by simulating
## data (the red curves). This process is very computationally intensive. This code 
## also returns the red power curves for Figure 3.

## choose an array of sample sizes at which to approximate the power curve
rrr_lower <- c(50,100,100)
rrr_upper <- c(450,500,1750)
rrr_inc <- c(10, 20, 50)

## settings for three gamma scenarios
convictions <- c(0.5, 0.9, 0.8)
powers <- c(0.6, 0.7, 0.8)
delta_Ls <- c(log(1/(1+0.25)), log(1/(1+0.3)), log(1/(1+0.15)))
delta_Us <- c(log(1 + 0.25), Inf, log(1 + 0.15))

## design values
gamma_alpha.1 <- 2.11
gamma_beta.1 <- 0.69
gamma_alpha.2 <- 2.43
gamma_beta.2 <- 0.79
threshold <- 4.29

informs <- read.csv("informs_gamma.csv")

## function to do MCMC for the posterior of theta1/theta2
GammaPost <- function(y1, y2, mu1, tau1, kappa1, lambda1,
                      mu2, tau2, kappa2, lambda2, tau = 4.29, burnin = 1000, nchains = 1,
                      nthin = 2, ndraws = 200000){
  
  ## y1 and y2 are the food expenditure observations in each group
  ## alpha_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate
  ## beta_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate
  ## tau is the threshold for the tail probability
  ## burnin is the number of MCMC iterations to discard at the start of each chain
  ## nchains is the number of chains to generate
  ## nthin is the thinning parameter for the MCMC process
  ## ndraws is the number of draws to generate (excluding burnin but including thinned draws)
  
  ## precompute summary statistics to make JAGS faster
  sum_y1 <- sum(y1)
  sum_logy1 <- sum(log(y1))
  n1 <- length(y1)
  model1.fit <- jags.model(file="JAGS_gamma.txt",
                           data=list(n=n1, sum_y = sum_y1, sum_logy = sum_logy1, 
                                     tau0 = tau1, mu0 = mu1, zero = 0,
                                     kappa0 = kappa1, lambda0 = lambda1), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model1.fit, burnin, progress.bar = "none")
  model1.samples <- coda.samples(model1.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  ## get posterior draws for alpha_j and beta_j
  alpha.1 <- unlist(model1.samples[,1])
  beta.1 <- unlist(model1.samples[,2])
  
  sum_y2 <- sum(y2)
  sum_logy2 <- sum(log(y2))
  n2 <- length(y2)
  model2.fit <- jags.model(file="JAGS_gamma.txt",
                           data=list(n=n2, sum_y = sum_y2, sum_logy = sum_logy2, 
                                     tau0 = tau2, mu0 = mu2, zero = 0,
                                     kappa0 = kappa2, lambda0 = lambda2), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model2.fit, burnin, progress.bar = "none")
  model2.samples <- coda.samples(model2.fit, c("alpha", "beta"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  alpha.2 <- unlist(model2.samples[,1])
  beta.2 <- unlist(model2.samples[,2])
  
  ## obtain posterior draws for theta
  theta1 <- 1 - pgamma(tau, alpha.1, beta.1)
  theta2 <- 1 - pgamma(tau, alpha.2, beta.2)
  theta <- theta1/theta2
  theta <- ifelse(is.na(theta), Inf, theta)
  return(theta)
}

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## output power curve for informative priors
for (k in 1:3){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         ## generate data from design distributions
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         ## approximate posterior
                         theta.diff <- GammaPost(y_star1, y_star2, informs[1,1], informs[1,2],
                                                 informs[2,1], informs[2,2], informs[3,1], 
                                                 informs[3,2], informs[4,1], informs[4,2],
                                                 threshold)
                         
                         ## determine whether 100 x gamma % of posterior is in (delta_L, delta_U) using t3
                         ## also check this criterion for the equal-tailed credible is satisfied using t1
                         ## t2 (this is to repurpose the same posteriors for Figure 3)
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    ## results_rrr is the estimate for the power criterion in (1.4) of the main text
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    ## results_rrr_ci is the estimate for the power criterion in (1.7) of the main text
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("confirm_results_2", k, ".csv"), row.names = FALSE)
  }
}

## repeat this process to output power curve for the settings with uninformative priors
for (k in 1:3){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         theta.diff <- GammaPost(y_star1, y_star2, 2, 0.25,
                                                 2, 0.25, 2, 0.25, 2, 0.25,
                                                 threshold)
                         
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("confirm_results_1", k, ".csv"), row.names = FALSE)
  }
}

closeAllConnections()

## code to produce Figure 2
## define plotting parameters for each setting
powers <- c(0.6, 0.7, 0.8)
settings <- c("a", "b", "c")
jj <- 0
for (k in c(1:3)){
  jj <- jj + 1
  for (ii in 1:2){
    ## power curve approximations from Algorithm 4 with Algorithms 1 and 2
    res_1 <- read.csv(paste0("alg1_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    res_2 <- read.csv(paste0("alg2_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    ## power curve approximations from simulating data from design distributions
    rrr <- read.csv(paste0("confirm_results_",ii, k, ".csv"))$rrr
    results_rrr <- read.csv(paste0("confirm_results_", ii, k, ".csv"))$results_rrr
    power <- powers[jj]
    y_mat <- NULL
    z_mat <- NULL
    for (i in 1:100){
      ## use the empirical CDF to get coordinates to plot from each power curve
      funecdf <- ecdf(as.numeric(res_1[i,]))
      xxx <- seq(0, max(res_1[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 99, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
      funecdf <- ecdf(as.numeric(res_2[i,]))
      xxx <- seq(0, max(res_2[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 100, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
    }
    
    ## use assign() and get() functions to automate this process for all settings
    assign(paste0("power", ii, k), power)
    
    assign(paste0("data_full",ii, k), data.frame(n = y_mat$n, power = z_mat$power, curve = y_mat$curve))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$n > 0))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$power > 0))
    
    assign(paste0("data_sim",ii, k), data.frame(n = rrr, power = results_rrr))
    assign(paste0("title", ii, k), paste0("Setting ", ii, settings[jj]))
    
    ## Power curves obtained by simulating data are red, power curves with Algorithm 1 are yellow,
    ## and those with Algorithm 2 are blue
    assign(paste0("plot",ii, k), ggplot(get(paste0("data_full",ii, k)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = power, color=as.factor(curve), alpha = 0.25), size = 1) +
             labs(title=get(paste0("title", ii, k))) +
             labs(x= bquote(italic(n)), y= bquote('Power')) +
             theme(plot.title = element_text(size=20,face="bold",
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c(rep(c("#E6C032", "#56B4E9"), 100), "firebrick")) +
             theme(legend.text=element_text(size=18)) +
             ylim(0,1) + 
             xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
             geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
             geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
    )
  }
}

## create common legend for the bottom of the plot grid
assign(paste0("plot",ii, k, "legend"), ggplot(subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$curve < 103), aes(x=n)) + theme_bw() +
         geom_line(aes(y = power, color=as.factor(curve)), size = 1) +
         labs(title=get(paste0("title", ii, k))) +
         labs(x= bquote(italic(n)), y= bquote('Power')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ",
                            labels = c(rep(c("Algorithm 1  ", "Algorithm 2  "), 1), "Data"), 
                            values = c(rep(c("#E6C032", "#56B4E9"), 1), "firebrick")) +
         theme(legend.text=element_text(size=18)) +
         ylim(0,1) + 
         xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
         geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
         geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
)

## arrange in grid form
figp.row1 <- plot_grid(plot11 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot21 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row2 <- plot_grid(plot12 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot22 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp.row3 <- plot_grid(plot13 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot23 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))
figp <- plot_grid(figp.row1, figp.row2, figp.row3, nrow = 3)

fig_final <- plot_grid(figp, get_legend(get(paste0("plot",ii, k, "legend"))), ncol = 1, rel_heights = c(2, .1))

# output as .pdf file for the article
pdf(file = "Fig2BA.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 8.8) # The height of the plot in inches

fig_final

dev.off()

## The next section of code estimates the power curves for Figure D.1 by simulating
## data (the red curves). This process is very computationally intensive.

## choose an array of sample sizes at which to approximate the power curve
rrr_lower <- c(50,100,100)
rrr_upper <- c(250,500,1750)
rrr_inc <- c(10, 20, 50)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

convictions <- c(0.5428)
## output power curve for informative priors
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         ## generate data from design distributions
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         ## approximate posterior
                         theta.diff <- GammaPost(y_star1, y_star2, informs[1,1], informs[1,2],
                                                 informs[2,1], informs[2,2], informs[3,1], 
                                                 informs[3,2], informs[4,1], informs[4,2],
                                                 threshold)
                         
                         ## determine whether 100 x gamma % of posterior is in (delta_L, delta_U) using t3
                         ## also check this criterion for the equal-tailed credible is satisfied using t1
                         ## t2 (this is to repurpose the same posteriors for Figure 3)
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    ## results_rrr is the estimate for the power criterion in (1.4) of the main text
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    ## results_rrr_ci is the estimate for the power criterion in (1.7) of the main text
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("bf_confirm_results_2", k, ".csv"), row.names = FALSE)
  }
}

convictions <- c(0.5652)
## repeat this process to output power curve for the settings with uninformative priors
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         theta.diff <- GammaPost(y_star1, y_star2, 2, 0.25,
                                                 2, 0.25, 2, 0.25, 2, 0.25,
                                                 threshold)
                         
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("bf_confirm_results_1", k, ".csv"), row.names = FALSE)
  }
}

closeAllConnections()

## code to produce Figure D.1
## define plotting parameters for each setting
powers <- c(0.6, 0.7, 0.8)
settings <- c("a", "b", "c")
jj <- 0
for (k in c(1)){
  jj <- jj + 1
  for (ii in 1:2){
    ## power curve approximations from Algorithm 4 with Algorithms 1 and 2
    res_1 <- read.csv(paste0("alg1_bf_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    res_2 <- read.csv(paste0("alg2_bf_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    ## power curve approximations from simulating data from design distributions
    rrr <- read.csv(paste0("bf_confirm_results_",ii, k, ".csv"))$rrr
    results_rrr <- read.csv(paste0("bf_confirm_results_", ii, k, ".csv"))$results_rrr
    power <- powers[jj]
    y_mat <- NULL
    z_mat <- NULL
    for (i in 1:100){
      ## use the empirical CDF to get coordinates to plot from each power curve
      funecdf <- ecdf(as.numeric(res_1[i,]))
      xxx <- seq(0, max(res_1[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 99, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
      funecdf <- ecdf(as.numeric(res_2[i,]))
      xxx <- seq(0, max(res_2[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 100, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
    }
    
    ## use assign() and get() functions to automate this process for all settings
    assign(paste0("power", ii, k), power)
    
    assign(paste0("data_full",ii, k), data.frame(n = y_mat$n, power = z_mat$power, curve = y_mat$curve))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$n > 0))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$power > 0))
    
    assign(paste0("data_sim",ii, k), data.frame(n = rrr, power = results_rrr))
    assign(paste0("title", ii, k), paste0("Setting ", ii, settings[jj]))
    
    ## Power curves obtained by simulating data are red, power curves with Algorithm 1 are yellow,
    ## and those with Algorithm 2 are blue
    assign(paste0("plot",ii, k), ggplot(get(paste0("data_full",ii, k)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = power, color=as.factor(curve), alpha = 0.25), size = 1) +
             labs(title=get(paste0("title", ii, k))) +
             labs(x= bquote(italic(n)), y= bquote('Power')) +
             theme(plot.title = element_text(size=20,face="bold",
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c(rep(c("#E6C032", "#56B4E9"), 100), "firebrick")) +
             theme(legend.text=element_text(size=18)) +
             ylim(0,1) + 
             xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
             geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
             geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
    )
  }
}

## create common legend for the bottom of the plot grid
assign(paste0("plot",ii, k, "legend"), ggplot(subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$curve < 103), aes(x=n)) + theme_bw() +
         geom_line(aes(y = power, color=as.factor(curve)), size = 1) +
         labs(title=get(paste0("title", ii, k))) +
         labs(x= bquote(italic(n)), y= bquote('Power')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ",
                            labels = c(rep(c("Algorithm 1  ", "Algorithm 2  "), 1), "Data"), 
                            values = c(rep(c("#E6C032", "#56B4E9"), 1), "firebrick")) +
         theme(legend.text=element_text(size=18)) +
         ylim(0,1) + 
         xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
         geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
         geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
)

## arrange in grid form
figp.row1 <- plot_grid(plot11 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot21 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))

fig_final <- plot_grid(figp.row1, get_legend(get(paste0("plot",ii, k, "legend"))), ncol = 1, rel_heights = c(2/3, .1))

# output as .pdf file for the article
pdf(file = "FigD1BA.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 3.1) # The height of the plot in inches

fig_final

dev.off()

## The next section of code estimates the power curves for Figure D.2 by simulating
## data (the red curves). This process is very computationally intensive.

## choose an array of sample sizes at which to approximate the power curve
rrr_lower <- c(50,100,100)
rrr_upper <- c(200,500,1750)
rrr_inc <- c(10, 20, 50)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 10000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

convictions <- c(0.5)
## output power curve for informative priors
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,50,4), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         ## generate data from design distributions (sample size is twice
                         ## as large in second group)
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(2*rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         ## approximate posterior
                         theta.diff <- GammaPost(y_star1, y_star2, informs[1,1], informs[1,2],
                                                 informs[2,1], informs[2,2], informs[3,1], 
                                                 informs[3,2], informs[4,1], informs[4,2],
                                                 threshold)
                         
                         ## determine whether 100 x gamma % of posterior is in (delta_L, delta_U) using t3
                         ## also check this criterion for the equal-tailed credible is satisfied using t1
                         ## t2 (this is to repurpose the same posteriors for Figure 3)
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    ## results_rrr is the estimate for the power criterion in (1.4) of the main text
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    ## results_rrr_ci is the estimate for the power criterion in (1.7) of the main text
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("q2_confirm_results_2", k, ".csv"), row.names = FALSE)
  }
}

closeAllConnections()

## repeat this process to output power curve for the settings with uninformative priors
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,50,4), rrr)
  }
  if (k == 2){
    ## add more granularity for second sample size setting
    rrr <- c(5, seq(10,90,10), rrr)
  }
  results_rrr <- rep(0, length(rrr))
  results_rrr_ci <- rep(0, length(rrr))
  
  dL <- delta_Ls[k]
  dU <- delta_Us[k]
  alpha <- 1 - convictions[k]
  pwer <- convictions[k]
  
  for (i in 1:length(rrr)){
    print(rrr[i])
    
    pwr_rep <- foreach(j=1:10000, .combine='rbind', .packages = c("rjags", "coda"),
                       .options.snow=opts) %dopar% {
                         
                         ## generate data from design distributions (sample size is twice
                         ## as large in second group)
                         y_star1 <- rgamma(rrr[i], gamma_alpha.1, gamma_beta.1)
                         y_star2 <- rgamma(2*rrr[i], gamma_alpha.2, gamma_beta.2)
                         
                         theta.diff <- GammaPost(y_star1, y_star2, 2, 0.25,
                                                 2, 0.25, 2, 0.25, 2, 0.25,
                                                 threshold)
                         
                         t1 <- mean(ifelse(theta.diff > dU,1,0))
                         t2 <- mean(ifelse(theta.diff < dL,1,0))
                         t3 <- mean(ifelse(theta.diff > dL,
                                           ifelse(theta.diff < dU,1,0), 0))
                         c(t3 >= pwer, ifelse(t1 < alpha/2, t2 < alpha/2, 0))
                         
                       }
    
    results_rrr[i] <- mean(as.numeric(pwr_rep[,1]))
    results_rrr_ci[i] <- mean(as.numeric(pwr_rep[,2]))
    write.csv(data.frame(rrr = rrr , results_rrr = results_rrr, results_rrr_ci), 
              paste0("q2_confirm_results_1", k, ".csv"), row.names = FALSE)
  }
}

## code to produce Figure D.2
## define plotting parameters for each setting
powers <- c(0.6, 0.7, 0.8)
settings <- c("a", "b", "c")
jj <- 0
for (k in c(1)){
  jj <- jj + 1
  for (ii in 1:2){
    ## power curve approximations from Algorithm 4 with Algorithms 1 and 2
    res_1 <- read.csv(paste0("alg1_q2_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    res_2 <- read.csv(paste0("alg2_q2_bpc_gamma_samps_pwr_", ii, k, ".csv"))
    ## power curve approximations from simulating data from design distributions
    rrr <- read.csv(paste0("q2_confirm_results_",ii, k, ".csv"))$rrr
    results_rrr <- read.csv(paste0("q2_confirm_results_", ii, k, ".csv"))$results_rrr
    power <- powers[jj]
    y_mat <- NULL
    z_mat <- NULL
    for (i in 1:100){
      ## use the empirical CDF to get coordinates to plot from each power curve
      funecdf <- ecdf(as.numeric(res_1[i,]))
      xxx <- seq(0, max(res_1[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 99, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
      funecdf <- ecdf(as.numeric(res_2[i,]))
      xxx <- seq(0, max(res_2[i,]), length.out = 600)
      y_mat <- rbind(y_mat, data.frame(n = matrix(xxx, ncol = 1), curve = rep(2*i + 100, length(xxx))))
      z_mat <- rbind(z_mat, data.frame(power = matrix(funecdf(xxx), ncol = 1)))
      
    }
    
    ## use assign() and get() functions to automate this process for all settings
    assign(paste0("power", ii, k), power)
    
    assign(paste0("data_full",ii, k), data.frame(n = y_mat$n, power = z_mat$power, curve = y_mat$curve))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$n > 0))
    assign(paste0("data_full",ii, k), subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$power > 0))
    
    assign(paste0("data_sim",ii, k), data.frame(n = rrr, power = results_rrr))
    assign(paste0("title", ii, k), paste0("Setting ", ii, settings[jj]))
    
    ## Power curves obtained by simulating data are red, power curves with Algorithm 1 are yellow,
    ## and those with Algorithm 2 are blue
    assign(paste0("plot",ii, k), ggplot(get(paste0("data_full",ii, k)), aes(x=n)) + theme_bw() +
             geom_line(aes(y = power, color=as.factor(curve), alpha = 0.25), size = 1) +
             labs(title=get(paste0("title", ii, k))) +
             labs(x= bquote(italic(n)), y= bquote('Power')) +
             theme(plot.title = element_text(size=20,face="bold",
                                             margin = margin(t = 0, 0, 5, 0))) +
             theme(axis.text=element_text(size=16),
                   axis.title=element_text(size=18)) +
             theme(legend.position="none") +
             scale_color_manual(name = " ", 
                                values = c(rep(c("#E6C032", "#56B4E9"), 100), "firebrick")) +
             theme(legend.text=element_text(size=18)) +
             ylim(0,1) + 
             xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
             theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
             theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
             geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
             geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
    )
  }
}

## create common legend for the bottom of the plot grid
assign(paste0("plot",ii, k, "legend"), ggplot(subset(get(paste0("data_full",ii, k)), get(paste0("data_full",ii, k))$curve < 103), aes(x=n)) + theme_bw() +
         geom_line(aes(y = power, color=as.factor(curve)), size = 1) +
         labs(title=get(paste0("title", ii, k))) +
         labs(x= bquote(italic(n)), y= bquote('Power')) +
         theme(plot.title = element_text(size=20,face="bold",
                                         margin = margin(t = 0, 0, 5, 0))) +
         theme(axis.text=element_text(size=16),
               axis.title=element_text(size=18)) +
         theme(legend.position="bottom") +
         scale_color_manual(name = " ",
                            labels = c(rep(c("Algorithm 1  ", "Algorithm 2  "), 1), "Data"), 
                            values = c(rep(c("#E6C032", "#56B4E9"), 1), "firebrick")) +
         theme(legend.text=element_text(size=18)) +
         ylim(0,1) + 
         xlim(0, max(get(paste0("data_sim",ii, k))$n)) + 
         theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
         theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
         geom_line(data = get(paste0("data_sim",ii, k)), aes(x = n, y = power, color="firebrick"), size = 1) +
         geom_hline(yintercept=get(paste0("power", ii, k)), lty = 2)
)

## arrange in grid form
figp.row1 <- plot_grid(plot11 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")), 
                       plot21 + theme(plot.margin=unit(c(0.25,0.5,0.25,0.5),"cm")),
                       rel_widths = c(1, 1))

fig_final <- plot_grid(figp.row1, get_legend(get(paste0("plot",ii, k, "legend"))), ncol = 1, rel_heights = c(2/3, .1))

# output as .pdf file for the article
pdf(file = "FigD2BA.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 3.1) # The height of the plot in inches

fig_final

dev.off()