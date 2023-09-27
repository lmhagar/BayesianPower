## code to approximate power curves for hypothesis tests with the Weibull model
## this procedure leverages Algorithm 1 and Algorithm 2 in two separate functions

## these functions were used to produce Figure 4 of the main text

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

powerCurveWeibull1 <- function(conviction, power, delta_L, delta_U, wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2,  
                             threshold = NULL, m = 1024, q = 1, seed = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for hypothesis test
  ## delta_U: upper interval endpoint for hypothesis test
  ## wei_lambda.j: the design value for the scale parameter of the Weibull distribution for group j = 1, 2
  ## wei_k.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  ## seed: seed used to generated randomized Sobol' sequence for reproducibility (default is 1)
  
  lambda1 <- wei_lambda.1; k1 <- wei_k.2; lambda2 <- wei_lambda.2; k2 <- wei_k.2
  ## define theta metrics in terms of design values
  theta1 <- 1 - pweibull(threshold, k1, lambda1)
  theta2 <- 1 - pweibull(threshold, k2, lambda2)
  
  ## compute partial derivatives of theta with respect to lambda and k for each group
  d_lambda1 <- (k1/lambda1)*exp(-(threshold/lambda1)^(k1))*(threshold/lambda1)^(k1)
  
  d_k1 <- -1*exp(-(threshold/lambda1)^(k1))*(threshold/lambda1)^(k1)*log(threshold/lambda1)
  
  d_lambda2 <- (k2/lambda2)*exp(-(threshold/lambda2)^(k2))*(threshold/lambda2)^(k2)
  
  d_k2 <- -1*exp(-(threshold/lambda2)^(k2))*(threshold/lambda2)^(k2)*log(threshold/lambda2)
  
  emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
  ## compute Fisher information for each gamma model
  Fish1inv <- matrix((sqrt(6)/pi)^2*c((lambda1/k1)^2*(1 + emc), lambda1*(1 + digamma(1)), lambda1*(1 + digamma(1)), k1^2), nrow = 2)
  
  Fish2inv <- matrix((sqrt(6)/pi)^2*c((lambda2/k2)^2*(1 + emc), lambda2*(1 + digamma(1)), lambda2*(1 + digamma(1)), k2^2), nrow = 2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_lambda1, d_k1))%*%Fish1inv%*%c(d_lambda1, d_k1)
  
  avar2 <- t(c(d_lambda2, d_k2))%*%Fish2inv%*%c(d_lambda2, d_k2)/q
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- as.numeric(avar1/theta1^2 + avar2/theta2^2)
  
  ## simplified uniroot
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
  
  ## function to input into simplified uniroot function
  ## u is point from Sobol' sequence, targetP is the target power
  ## for the root-finding algorithm to obtain, and n_val is the 
  ## sample size presently explored by the root-finding algorithm
  targetPower <- function(u, targetP, params, deltas, n_val, q){
    
    ## return negative power if sample size is not positive
    if (n_val <= 0){return(-1.5)}
    
    wei_lambda.1 <- params[1]
    wei_k.1 <- params[2]
    wei_lambda.2 <- params[3]
    wei_k.2 <- params[4]
    threshold <- params[5]
    
    emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
    ## generate approximately normal MLEs for group 1 using delta method
    rho1 <- (1 + digamma(1))/sqrt(1 + emc)
    mat1 <- matrix((sqrt(6)/(pi*wei_k.1))^2*c(1 + emc, wei_k.1*(1 + digamma(1)), wei_k.1*(1 + digamma(1)), wei_k.1^2), nrow = 2)
    l1 <- qnorm(u[1], log(wei_lambda.1), sqrt(mat1[1,1]/n_val))
    k1 <- qnorm(u[2], log(wei_k.1) + rho1*(l1 - log(wei_lambda.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    ## generate approximately normal MLEs for group 2 using delta method
    rho2 <- (1 + digamma(1))/sqrt(1 + emc)
    mat2 <- matrix((sqrt(6)/(pi*wei_k.2))^2*c(1 + emc, wei_k.2*(1 + digamma(1)), wei_k.2*(1 + digamma(1)), wei_k.2^2), nrow = 2)
    l2 <- qnorm(u[3], log(wei_lambda.2), sqrt(mat2[1,1]/(q*n_val)))
    k2 <- qnorm(u[4], log(wei_k.2) + rho2*(l2 - log(wei_lambda.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/(q*n_val)))
    
    ## exponentiate MLEs
    l1 <- exp(l1); k1 <- exp(k1); l2 <- exp(l2); k2 <- exp(k2)
    
    ## ensure no MLEs underflow to 0
    if(max(l1 <= 0, k1 <= 0, l2 <= 0, k2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pweibull(threshold, k1, l1)
    theta2 <- 1 - pweibull(threshold, k2, l2)
    
    if(max(theta1 <= 0, theta1 >= 1, theta2 <= 0, theta2 >= 1)){return(-1.5)}
    
    ## compute partial derivatives of theta with respect to alpha and beta for each group
    d_lambda1 <- (k1/l1)*exp(-(threshold/l1)^(k1))*(threshold/l1)^(k1)
    
    d_k1 <- -1*exp(-(threshold/l1)^(k1))*(threshold/l1)^(k1)*log(threshold/l1)
    
    d_lambda2 <- (k2/l2)*exp(-(threshold/l2)^(k2))*(threshold/l2)^(k2)
    
    d_k2 <- -1*exp(-(threshold/l2)^(k2))*(threshold/l2)^(k2)*log(threshold/l2)
    
    ## compute Fisher information for each gamma model
    Fish1inv <- matrix((sqrt(6)/pi)^2*c((l1/k1)^2*(1 + emc), l1*(1 + digamma(1)), l1*(1 + digamma(1)), k1^2), nrow = 2)
    
    Fish2inv <- matrix((sqrt(6)/pi)^2*c((l2/k2)^2*(1 + emc), l2*(1 + digamma(1)), l2*(1 + digamma(1)), k2^2), nrow = 2)
    
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_lambda1, d_k1))%*%Fish1inv%*%c(d_lambda1, d_k1)
    
    avar2 <- t(c(d_lambda2, d_k2))%*%Fish2inv%*%c(d_lambda2, d_k2)/q
    
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
  
  ## find some better endpoints for the root-finding algorithm (closer to the 
  ## sample size that should be returned by the root-finding algorithm).
  ## we use a target power of 0 here (that is subtracted from the power at
  ## that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(mu_start),
                          params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      power1b <- targetPower(targetP = 0, n_val = ceiling(mu_low),
                             params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power1b >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(mu_low),
                           f_upper = power1b - conviction, targetP = conviction,
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_low), upper = ceiling(mu_start),
                           f_lower = power1b - conviction, f_upper = power1 - conviction, targetP = conviction,
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_high),
                            params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(mu_high), 
                           f_lower = power1 - conviction, f_upper = power2 - conviction, targetP = conviction, 
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
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
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                                         deltas = c(delta_L, delta_U), sob_pwr[last_group[i],], q = q)
      }
      keep_vec <- ifelse(endpoints1_vec >= conviction, FALSE, TRUE)
      ## only keep points that still do not satisfy power criterion after increasing
      ## the upper bound for the sample size
      last_group <- last_group[keep_vec]
    }
  }
  
  ## implement the root-finding algorithm for the last points in the Sobol' sequence
  ## with appropriate endpoints
  for (i in 1:nrow(sob_pwr)){
    if (endpoints_vec[i] == 1){
      samps_pwr[i] <- uu(targetPower, lower = ceiling(mu_high), upper = ceiling(upper_c*mu_high), targetP = conviction, 
                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at n_* (noninteger)
  n_star <- quantile(samps_pwr, power)
  power_star <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_star[i] <- targetPower(targetP = 0, n_val = n_star,
                                 params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
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
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q)
      }
      else{
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = 2, upper = n_star, 
                                         f_upper = power_star[inconsistent[i]] - conviction, targetP = conviction, 
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
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

## Algorithm 2
powerCurveWeibull2 <- function(conviction, power, delta_L, delta_U, wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, 
                             mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, mu0.2 = 2, 
                             tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25, threshold = NULL, m = 1024, q = 1, seed = 1){
  
  ## conviction: the conviction threshold in [0.5, 1), denoted gamma in the paper
  ## power: the target power, denoted capital Gamma in the paper
  ## delta_L: lower interval endpoint for hypothesis test
  ## delta_U: upper interval endpoint for hypothesis test
  ## wei_lambda.j: the design value for the scale parameter of the Weibull distribution for group j = 1, 2
  ## wei_k.j: the design value for the shape parameter of the gamma distribution for group j = 1, 2
  ## lambda_j has a Gamma(mu0.j, tau0.j) prior, where tau0.j is a rate for group j = 1, 2
  ## k_j has a Gamma(kappa0.j, lambda0.j) prior, where lambda0.j is a rate for group j = 1, 2
  ## threshold: the threshold for the tail probability in each group, denoted tau in the paper
  ## m: length of Sobol' sequence to approximate power curve
  ## q: constant for imbalanced sample size determination, where n_2 = q*n_1
  ## seed: seed used to generated randomized Sobol' sequence for reproducibility (default is 1)
  
  ## concatenate all hyperparameters in matrix form (used in getNormal() function)
  hyper_mat <- rbind(c(mu0.1, tau0.1), c(kappa0.1, lambda0.1), c(mu0.2, tau0.2), c(kappa0.2, lambda0.2))
  
  lambda1 <- wei_lambda.1; k1 <- wei_k.2; lambda2 <- wei_lambda.2; k2 <- wei_k.2
  ## define theta metrics in terms of design values
  theta1 <- 1 - pweibull(threshold, k1, lambda1)
  theta2 <- 1 - pweibull(threshold, k2, lambda2)
  
  ## compute partial derivatives of theta with respect to lambda and k for each group
  d_lambda1 <- (k1/lambda1)*exp(-(threshold/lambda1)^(k1))*(threshold/lambda1)^(k1)
  
  d_k1 <- -1*exp(-(threshold/lambda1)^(k1))*(threshold/lambda1)^(k1)*log(threshold/lambda1)
  
  d_lambda2 <- (k2/lambda2)*exp(-(threshold/lambda2)^(k2))*(threshold/lambda2)^(k2)
  
  d_k2 <- -1*exp(-(threshold/lambda2)^(k2))*(threshold/lambda2)^(k2)*log(threshold/lambda2)
  
  emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
  ## compute Fisher information for each gamma model
  Fish1inv <- matrix((sqrt(6)/pi)^2*c((lambda1/k1)^2*(1 + emc), lambda1*(1 + digamma(1)), lambda1*(1 + digamma(1)), k1^2), nrow = 2)
  
  Fish2inv <- matrix((sqrt(6)/pi)^2*c((lambda2/k2)^2*(1 + emc), lambda2*(1 + digamma(1)), lambda2*(1 + digamma(1)), k2^2), nrow = 2)
  
  ## apply the delta method to get the limiting variance for each theta_j metric
  avar1 <- t(c(d_lambda1, d_k1))%*%Fish1inv%*%c(d_lambda1, d_k1)
  
  avar2 <- t(c(d_lambda2, d_k2))%*%Fish2inv%*%c(d_lambda2, d_k2)/q
  
  ## apply the delta method to get the limiting variance for theta = logtheta1 - logtheta2
  Fish_ratio_mu <- avar1/theta1^2 + avar2/theta2^2
  
  ## simplified uniroot
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
                   lower = lower_cons, upper = upper_cons, tol2 = 0.1))^2
    
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
  ## x and y are parameterized as for the previous function
  fn_grad <- function(x, y, mu, tau, kappa, lambda) {
    
    res1 <- -y[3]*(x[1] - y[1])*exp(y[2])^2 + y[3]*(x[2] - y[2])*exp(y[2])*(1 + digamma(1)) - exp(x[1])*tau + mu
    
    res2 <- y[3]*(x[1] - y[1])*exp(y[2])*(1 + digamma(1)) - y[3]*(x[2] - y[2])*y[4] - exp(x[2])*lambda + kappa
    
    return(c(res1, res2))
  }
  
  calc_covar <- function(u, v, yy_star, hyper){
    l <- exp(u); k <- exp(v); n <- yy_star[3]
    tau <- hyper[1]; lambda <- hyper[2]
    d11 <- n*k^2 + l*tau
    d12 <- -n*k*(1 + digamma(1))
    d22 <- n*yy_star[4] + k*lambda
    mat <- rbind(c(d11, d12), c(d12, d22))
    return(solve(mat))
  }
  
  ## function to input into simplified uniroot function
  ## u is point from Sobol' sequence, targetP is the target power
  ## for the root-finding algorithm to obtain, and n_val is the 
  ## sample size presently explored by the root-finding algorithm
  targetPower <- function(u, targetP, params, deltas, n_val, hyper, q){
    
    ## return negative power if sample size is not positive
    if (n_val <= 0){return(-1.5)}
    
    wei_lambda.1 <- params[1]
    wei_k.1 <- params[2]
    wei_lambda.2 <- params[3]
    wei_k.2 <- params[4]
    threshold <- params[5]
    
    emc <- pi^2/6 + (digamma(1))^2 + 2*digamma(1)
    ## generate approximately normal MLEs for group 1 using delta method
    rho1 <- (1 + digamma(1))/sqrt(1 + emc)
    mat1 <- matrix((sqrt(6)/(pi*wei_k.1))^2*c(1 + emc, wei_k.1*(1 + digamma(1)), wei_k.1*(1 + digamma(1)), wei_k.1^2), nrow = 2)
    l1 <- qnorm(u[1], log(wei_lambda.1), sqrt(mat1[1,1]/n_val))
    k1 <- qnorm(u[2], log(wei_k.1) + rho1*(l1 - log(wei_lambda.1))*(sqrt(mat1[2,2])/sqrt(mat1[1,1])), sqrt(1- rho1^2)*sqrt(mat1[2,2]/n_val))
    
    ## generate approximately normal MLEs for group 2 using delta method
    rho2 <- (1 + digamma(1))/sqrt(1 + emc)
    mat2 <- matrix((sqrt(6)/(pi*wei_k.2))^2*c(1 + emc, wei_k.2*(1 + digamma(1)), wei_k.2*(1 + digamma(1)), wei_k.2^2), nrow = 2)
    l2 <- qnorm(u[3], log(wei_lambda.2), sqrt(mat2[1,1]/(q*n_val)))
    k2 <- qnorm(u[4], log(wei_k.2) + rho2*(l2 - log(wei_lambda.2))*(sqrt(mat2[2,2])/sqrt(mat2[1,1])), sqrt(1- rho2^2)*sqrt(mat2[2,2]/(q*n_val)))
    
    ## summarize information from first group of data (faster computation)
    yy_star1 <- c(l1, k1, n_val, 1 + emc)
    ## find posterior modes for the first group (logalpha and logbeta)
    modes <- nleqslv(c(l1, k1), fn_grad, y = yy_star1, mu = hyper[1,1], tau = hyper[1,2],
                     kappa = hyper[2,1], lambda = hyper[2,2] )$x
    
    mat1_new <- calc_covar(modes[1], modes[2], yy_star1, c(hyper[1,2], hyper[2,2]))
    ## exponentiate modes to return to standard scale
    modes1 <- exp(modes)
    
    ## repeat all steps for the second group
    yy_star2 <- c(l2, k2, q*n_val, 1 + emc)
    modes <- nleqslv(c(l2, k2), fn_grad, y = yy_star2, mu = hyper[3,1], tau = hyper[3,2],
                     kappa = hyper[4,1], lambda = hyper[4,2] )$x
    
    mat2_new <- calc_covar(modes[1], modes[2], yy_star2, c(hyper[3,2], hyper[4,2]))
    modes2 <- exp(modes)
    l1 <- modes1[1]; k1 <- modes1[2]; l2 <- modes2[1]; k2 <- modes2[2]
    
    ## ensure no modes are 0 due to underflow errors
    if(max(l1 <= 0, k1 <= 0, l2 <= 0, k2<= 0)){return(-1.5)}
    
    ## define theta metrics in terms of design values
    theta1 <- 1 - pweibull(threshold, k1, l1)
    theta2 <- 1 - pweibull(threshold, k2, l2)
    
    if(max(theta1 <= 0, theta1 >= 1, theta2 <= 0, theta2 >= 1)){return(-1.5)}
    
    ## compute partial derivatives of theta with respect to logalpha and logbeta for each group
    ## this is different from the previous function that computes the derivatives with respect
    ## to alpha and beta
    d_lambda1 <- k1*exp(-(threshold/l1)^(k1))*(threshold/l1)^(k1)
    
    d_k1 <- -k1*exp(-(threshold/l1)^(k1))*(threshold/l1)^(k1)*log(threshold/l1)
    
    d_lambda2 <- k2*exp(-(threshold/l2)^(k2))*(threshold/l2)^(k2)
    
    d_k2 <- -k2*exp(-(threshold/l2)^(k2))*(threshold/l2)^(k2)*log(threshold/l2)
    
    ## apply the delta method to get the limiting variance for each theta_j metric
    avar1 <- t(c(d_lambda1, d_k1))%*%mat1_new%*%c(d_lambda1, d_k1)
    
    avar2 <- t(c(d_lambda2, d_k2))%*%mat2_new%*%c(d_lambda2, d_k2)
    
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
  
  ## find some better endpoints for the root-finding algorithm (closer to the 
  ## sample size that should be returned by the root-finding algorithm).
  ## we use a target power of 0 here (that is subtracted from the power at
  ## that sample size) to obtain the actual power
  endpoints_vec <- rep(0, m)
  samps_pwr <- NULL
  for (i in 1:nrow(sob_pwr)){
    power1 <- targetPower(targetP = 0, n_val = ceiling(mu_start), hyper = hyper_mat,
                          params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                          deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    if (power1 >= conviction){
      tempi <- targetPower(targetP = 0, n_val = ceiling(0.5*mu_start), hyper = hyper_mat,
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (tempi < conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(0.5*mu_start), upper = ceiling(mu_start), 
                           f_lower = tempi - conviction, f_upper = power1 - conviction, 
                           targetP = conviction, hyper = hyper_mat,
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                           deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      }
      else{
        tempii <- targetPower(targetP = 0, n_val = ceiling(0.25*mu_start), hyper = hyper_mat,
                              params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                              deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        if (tempii < conviction){
          samps_pwr[i] <- uu(targetPower, lower =ceiling(0.25*mu_start), upper = ceiling(0.5*mu_start), 
                             f_lower = tempii - conviction, f_upper = tempi - conviction,
                             targetP = conviction, hyper = hyper_mat,
                             params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                             deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
        }
        else{
          tempiii <- targetPower(targetP = 0, n_val = ceiling(0.125*mu_start), hyper = hyper_mat,
                                 params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                                 deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
          if (tempiii < conviction){
            samps_pwr[i] <- uu(targetPower, lower =ceiling(0.125*mu_start), upper = ceiling(0.25*mu_start), 
                               f_lower = tempiii - conviction, f_upper = tempii - conviction,
                               targetP = conviction, hyper = hyper_mat,
                               params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                               deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
          }
          else{
            samps_pwr[i] <- uu(targetPower, lower =2, upper = ceiling(0.125*mu_start), 
                               f_lower = tempiii - conviction, targetP = conviction, hyper = hyper_mat,
                               params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                               deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
          }
        }
      }
    }
    else{
      power2 <- targetPower(targetP = 0, n_val = ceiling(mu_high), hyper = hyper_mat,
                            params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
                            deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
      if (power2 >= conviction){
        samps_pwr[i] <- uu(targetPower, lower =ceiling(mu_start), upper = ceiling(mu_high), 
                           f_lower = power1 - conviction, f_upper = power2 - conviction, targetP = conviction, 
                           params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold), hyper = hyper_mat,
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
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
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
                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold), hyper = hyper_mat,
                         deltas = c(delta_L, delta_U), u = sob_pwr[i,], q = q)
    }
  }
  
  ## confirmation power estimates at n_* (noninteger)
  n_star <- quantile(samps_pwr, power)
  power_star <- NULL
  for (i in 1:nrow(sob_pwr)){
    power_star[i] <- targetPower(targetP = 0, n_val = n_star, hyper = hyper_mat,
                                 params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold),
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
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold), hyper = hyper_mat,
                                         deltas = c(delta_L, delta_U), u = sob_pwr[inconsistent[i],], q = q)
      }
      else{
        samps_pwr[inconsistent[i]] <- uu(targetPower, lower = 2, upper = n_star, 
                                         f_upper = power_star[inconsistent[i]] - conviction, targetP = conviction, hyper = hyper_mat,
                                         params = c(wei_lambda.1, wei_k.1, wei_lambda.2, wei_k.2, threshold), 
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

## read informative hyperparameters from .csv file
informs <- read.csv("informs_weibull.csv")
## reorder matrix rows since they were output with the order shape -> scale
## and the function takes the reverse order
informs <- informs[c(2,1,4,3),]

## define the conviction threshold, target power, and delta_* for the 
## three scenarios as in the paper
convictions <- c(0.5)
powers <- c(0.6)
delta_Ls <- c(log(1/(1+0.25)))
delta_Us <- c(log(1 + 0.25))

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 100, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## Algorithm 1: Uninformative Priors
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                   .options.snow=opts) %dopar% {
                     unlist(powerCurveWeibull1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                               delta_U = delta_Us[k], wei_lambda.1 = 3.39, wei_k.1 = 1.41, seed = (k+23)*100 + j,
                                               wei_lambda.2 = 3.42, wei_k.2 = 1.49, threshold = 4.29, m = 1024))
                   }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_weibull_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_weibull_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bpc_weibull_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 1: Informative Priors
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveWeibull1(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                  delta_U = delta_Us[k], wei_lambda.1 = 3.39, wei_k.1 = 1.41, seed = (k+24)*100 + j,
                                                  wei_lambda.2 = 3.42, wei_k.2 = 1.49, threshold = 4.29, m = 1024))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg1_bpc_weibull_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg1_bpc_weibull_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg1_bpc_weibull_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Uninformative Priors
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveWeibull2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                  delta_U = delta_Us[k], wei_lambda.1 = 3.39, wei_k.1 = 1.41, 
                                                  wei_lambda.2 = 3.42, wei_k.2 = 1.49, threshold = 4.29, m = 1024, seed = (k+25)*100 + j,
                                                  mu0.1 = 2, tau0.1 = 0.25, kappa0.1 = 2, lambda0.1 = 0.25, 
                                                  mu0.2 = 2, tau0.2 = 0.25, kappa0.2 = 2, lambda0.2 = 0.25))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_weibull_samps_pwr_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_weibull_consistency_1", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bpc_weibull_adjust_1", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

## Algorithm 2: Informative Priors
tic <- Sys.time()
for (k in 1){
  res_temp <- foreach(j=1:100, .combine='rbind', .packages = c("numDeriv", "qrng", "nleqslv"),
                      .options.snow=opts) %dopar% {
                        unlist(powerCurveWeibull2(conviction = convictions[k], power = powers[k], delta_L = delta_Ls[k],
                                                  delta_U = delta_Us[k], wei_lambda.1 = 3.39, wei_k.1 = 1.41, 
                                                  wei_lambda.2 = 3.42, wei_k.2 = 1.49, threshold = 4.29, m = 1024, seed = (k+26)*100 + j,
                                                  mu0.1 = informs[1,1], tau0.1 = informs[1,2], kappa0.1 = informs[2,1], lambda0.1 = informs[2,2], 
                                                  mu0.2 = informs[3,1], tau0.2 = informs[3,2], kappa0.2 = informs[4,1], lambda0.2 = informs[4,2]))
                      }
  
  write.csv(res_temp[,1:1024], paste0("alg2_bpc_weibull_samps_pwr_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1025:1027], paste0("alg2_bpc_weibull_consistency_2", k, ".csv"), row.names = FALSE)
  write.csv(res_temp[,1028:1177], paste0("alg2_bpc_weibull_adjust_2", k, ".csv"), row.names = FALSE)
}
toc <- Sys.time()
toc - tic

closeAllConnections()

## power curve for Weibull data
## settings for simulation (values to compute power)
rrr_lower <- c(50,100,100)
rrr_upper <- c(250,600,1750)
rrr_inc <- c(10, 20, 50)

wei_lambda.1 = 3.39; wei_k.1 = 1.41 
wei_lambda.2 = 3.42; wei_k.2 = 1.49
threshold <- 4.29

informs <- read.csv("informs_weibull.csv")

## function to do MCMC for the posterior of theta1/theta2
WeibullPost <- function(y1, y2, mu1, tau1, kappa1, lambda1,
                        mu2, tau2, kappa2, lambda2, tau = 4.29, burnin = 1000, nchains = 1,
                        nthin = 1, ndraws = 50000){
  
  if (length(y1) <= 100){
    ndraws <- 75000
  }
  else if (length(y1) >= 200){
    ndraws <- 25000
  }
  
  ## y1 and y2 are the food expenditure observations in each group
  ## lambda_j has a Gamma(mu_j, tau_j) prior, where tau_j is a rate
  ## k_j has a Gamma(kappa_j, lambda_j) prior, where lambda_j is a rate (nu_j in paper)
  ## tau is the threshold for the tail probability
  ## burnin is the number of MCMC iterations to discard at the start of each chain
  ## nchains is the number of chains to generate
  ## nthin is the thinning parameter for the MCMC process
  ## ndraws is the number of draws to generate (excluding burnin but including thinned draws)
  
  n1 <- length(y1)
  model1.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n1, y = y1, 
                                     tau0 = tau1, mu0 = mu1,
                                     kappa0 = kappa1, lambda0 = lambda1), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model1.fit, burnin, progress.bar = "none")
  model1.samples <- coda.samples(model1.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.1 <- unlist(model1.samples[,2])
  lambda.1 <- unlist(model1.samples[,1])
  
  n2 <- length(y2)
  model2.fit <- jags.model(file="JAGS_weibull.txt",
                           data=list(n=n2, y = y2, 
                                     tau0 = tau2, mu0 = mu2,
                                     kappa0 = kappa2, lambda0 = lambda2), 
                           n.chains = nchains, quiet = TRUE)
  
  update(model2.fit, burnin, progress.bar = "none")
  model2.samples <- coda.samples(model2.fit, c("nu", "l"), n.iter=ndraws, thin=nthin, progress.bar = "none")
  
  nu.2 <- unlist(model2.samples[,2])
  lambda.2 <- unlist(model2.samples[,1])
  
  theta1 <- 1 - pweibull(tau, nu.1, lambda.1)
  theta2 <- 1 - pweibull(tau, nu.2, lambda.2)
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
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
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
                         
                         y_star1 <- rweibull(rrr[i], wei_k.1, wei_lambda.1)
                         y_star2 <- rweibull(rrr[i], wei_k.2, wei_lambda.2)
                         
                         theta.diff <- WeibullPost(y_star1, y_star2, informs[1,1], informs[1,2],
                                                   informs[2,1], informs[2,2], informs[3,1], informs[3,2],
                                                   informs[4,1], informs[4,2],
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
              paste0("wei_confirm_results_2", k, ".csv"), row.names = FALSE)
  }
}

## output power curve for uninformative priors
for (k in 1){
  rrr <- seq(rrr_lower[k], rrr_upper[k], rrr_inc[k])
  if (k == 1){
    ## add more granularity for smallest sample size setting
    rrr <- c(seq(10,45,5), rrr)
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
                         
                         y_star1 <- rweibull(rrr[i], wei_k.1, wei_lambda.1)
                         y_star2 <- rweibull(rrr[i], wei_k.2, wei_lambda.2)
                         
                         theta.diff <- WeibullPost(y_star1, y_star2, 2, 1,
                                                   2, 1, 2, 1, 2, 1,
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
              paste0("wei_confirm_results_1", k, ".csv"), row.names = FALSE)
  }
}

## code to produce Figure 4
## define plotting parameters for each setting
powers <- c(0.6, 0.7, 0.8)
settings <- c("a", "b", "c")
jj <- 0
for (k in c(1)){
  jj <- jj + 1
  for (ii in 1:2){
    ## power curve approximations from Algorithm 4 with Algorithms 1 and 3
    res_1 <- read.csv(paste0("alg1_bpc_weibull_samps_pwr_", ii, k, ".csv"))
    res_2 <- read.csv(paste0("alg2_bpc_weibull_samps_pwr_", ii, k, ".csv"))
    ## power curve approximations from simulating data from design distributions
    rrr <- read.csv(paste0("wei_confirm_results_",ii, k, ".csv"))$rrr
    results_rrr <- read.csv(paste0("wei_confirm_results_", ii, k, ".csv"))$results_rrr
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
    ## and those with Algorithm 3 are blue
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
                            labels = c(rep(c("Algorithm 1  ", "Algorithm 3  "), 1), "Data"), 
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
pdf(file = "Fig4BA.pdf",   # The directory you want to save the file in
    width = 10.5, # The width of the plot in inches
    height = 3.1) # The height of the plot in inches

fig_final

dev.off()