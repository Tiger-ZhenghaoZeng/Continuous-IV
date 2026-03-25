
#' mostly from https://github.com/ehkennedy/npcausal/blob/master/R/ctseff.R

#' Finished debugging for functions with "simu".
#' Code works correctly for estimating dose-response with f2_hat
#' DR is biased - why? regressing on Y (not the pseudo outcome) works the best (due to DGP?)
#' Once DR_simu is fixed, fix the original DR code

#' try kernel transformed outcomes?
#' could look at Zach's code: https://github.com/zjbranson/contTreatments/blob/main/appReplicationCode.R

# logit = function(x){
#   return(log(x / (1-x)))
# }

#' @param estimates: nsim x n matrix
#' @param true: n-vector of true derivative of dose reponse evaluated at a
# rmse = function(true, estimates, a){
#   pa = density(a, from=min(a), to=max(a), n=length(a))$y
#   sq_errors = sweep(estimates, 2, true, "-")^2
#   return(mean(sqrt(colMeans(sq_errors))*pa))
# }

#' exclude 10% mass at the boundaries
rmse = function(true, estimates){
  # lb = quantile(a, 0.10)
  # ub = quantile(a, 0.90)
  # 
  # idx = which(a >= lb & a <= ub)
  # 
  # true_sub = true[idx]; estimates_sub = estimates[, idx]
  # a_sub = a[idx]
  
  # pa = density(a_sub, from = min(a_sub), to = max(a_sub), n = length(a_sub))$y
  sq_errors = sweep(estimates, 2, true, "-")^2
  
  return(mean(sqrt(colMeans(sq_errors))))
}



expit = function(x){
  return(1/(1 + exp(-x)))
}

meanY_hat = function(beta1, beta2, epsilon){
  return(function(x,a) {
    return(c(1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*(1+epsilon)*a^2)))
  })
}

# f2_hat = function(beta1, beta2, epsilon){
#   return(function(x,a) {
#     return(c(1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*epsilon*a^2)))
#   })
# }

meanA_hat = function(beta_pi, b, epsilon){
  return(function(x) {
    return(c(b + x %*% beta_pi + epsilon))
  })
}

# prob2_hat = function(beta_pi, b, epsilon){
#   return(function(x, pi_vals) {
#     return(sqrt(expit(c(b + x %*% beta_pi + epsilon))*20) + pi_vals)
#   })
# }

# simulate_data = function(seed, n, beta1, beta2, beta_pi){
#   require(MASS)
#   
#   set.seed(seed)
#   
#   x = mvrnorm(n, rep(0,4), diag(4))
#   lambda_x = -0.8 + x %*% beta_pi
#   lambda_x = expit(lambda_x)
#   
#   a = 20 * rbeta(n, lambda_x, 1-lambda_x)
#   
#   mu = 1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*a^2) # or some nonlinear g(a)
#   y = mu + rnorm(n,0,2)
#   
#   return(data.frame(y=y, a=a, x=x))
# }

simulate_data_normal = function(seed, n, beta1, beta2, beta_pi){
  require(MASS)
  
  set.seed(seed)
  
  x = mvrnorm(n, rep(0,4), diag(4))
  lambda_x = -0.8 + x %*% beta_pi
  # lambda_x = expit(lambda_x)
  
  a = rnorm(n, lambda_x, 1)
  
  mu = 1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*a^2) # or some nonlinear g(a)
  y = mu + rnorm(n,0,2)
  
  return(data.frame(y=y, a=a, x=x))
}


# dctseff_plugin_simu <- Vectorize(dctseff_plugin_simu)

gaussianKernel_deriv = function(A, a, h){
  #the difference is
  u = A - a
  #then, the kernel is
  kernel = -(u/h^2)*(1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
  return(kernel)
}

gaussianKernel = function(A, a, h){
  #the difference is
  u = A - a
  #then, the kernel is
  kernel = (1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
  return(kernel)
}

epan_deriv <- function(A, a, h){
  u = A - a
  kernel <- -3*u/(2*h^3)*I(abs(u)<h)
  return(kernel)
}

#This function is useful for running a GAM
#for the PS model and outcome model.
#This simply increases the number
#of continuous variables from 4 (the default) to 5.
SL.gam.new <- function(..., cts.num = 5) {
  SL.gam(..., cts.num = 5)
}


#This function estimates the causal effect
#\psi_h(a), defined as the "smoothed causal effect"
#(without trimming).
#This is based on the EIF given in the paper.
#Requires a dataset with X, A, Y,
#a bandwidth h, and a set of treatments a.eval.
#Thus, this function estimates \psi_h(a.eval).
estthetah.simu = function(data, h, a.eval, muhat, lambdahat, kernelfun, npts_int=100){
  #number rows
  n = nrow(data)
  
  
  #we will ultimately compute the estimated \psi_h
  #for each a in a.eval.
  #This will be a matrix that collects
  #the estimated EIF values for each 
  #observation and treatment combination.
  #rows are observations; columns are treatment values.
  thetah.ifvalues = matrix(0, nrow = n, ncol = length(a.eval))
  
  
  a.min = min(data$a); a.max = max(data$a)
  a.seq = seq(a.min, a.max, length=npts_int)
  x <- as.matrix(subset(data, select = -c(a,y)))
  
  meanA.est = lambdahat(x)
  piHat.XA = dnorm(data$a, mean = meanA.est, sd = 1)
  piHat.XA <- ifelse(piHat.XA>0.01, piHat.XA, 0.01)
  muHat.XA = muhat(x, data$a)
  ipwRes = (data$y - muHat.XA)/piHat.XA
  kernel.A = sapply(a.eval, kernelfun, A = data$a, h = h)
  
  
  xa_new = cbind(x[rep(1:n,length(a.seq)),], a=rep(a.seq, rep(n,length(a.seq))))
  mu_pred_vals = muhat(xa_new[,1:(dim(xa_new)[2]-1)], xa_new[,dim(xa_new)[2]])
  muHat.a0 <- matrix(mu_pred_vals, nrow = n, ncol = length(a.seq))
  a.seq.len = a.seq[2] - a.seq[1]
  kernel.a0 = sapply(a.eval, kernelfun, A = a.seq, h = h)
  
  # normalization
  # kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  # kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  psih.1 = kernel.A*ipwRes
  psih.2 = a.seq.len*muHat.a0%*%kernel.a0
  
  # thetah.ifvalues = - psih.1 - psih.2
  thetah.ifvalues = -psih.1 - psih.2
  return(colMeans(thetah.ifvalues))
}

smooth_est_simu <- function(data, hs, a.eval, muhat, lambdahat, kernel_fun){
  mses <- numeric(length(hs))
  for (i in seq_along(hs)){
    res <- estthetah.simu(data, hs[i], a.eval, muhat, lambdahat, kernel_fun)
    mses[i] <- mean((res-true.out)^2)
  }
  print(which.min(mses))
  h_opt <- hs[which.min(mses)]
  return(estthetah.simu(data, h_opt, a.eval, muhat, lambdahat, kernel_fun))
}





dctseff_plugin_simu = function(a.out, x, muhat){
  require(numDeriv)
  n <- dim(x)[1]
  mu_function = function(a_val) {
    # x_new <- cbind(x, rep(a_val, n))
    temp <- muhat(x, rep(a_val, n))
    return(mean(temp))
  }
  # mu_function <- Vectorize(mu_function)
  # sapply(a.out, function(o) grad(mu_function, o))
  return(sapply(a.out, function(o) grad(mu_function, o)))
}


dctseff_plugin = function(y, a, x, n_pts = 100){
  require(SuperLearner)
  require(numDeriv)

  n = length(y)
  
  a_min = min(a); a_max = max(a)
  a_vals = seq(a_min, a_max, length.out=n_pts)
  
  sl_lib = c("SL.earth", "SL.gam", 
             "SL.glm", "SL.glm.interaction",
             "SL.mean", "SL.ranger", "SL.rpart")
  
  xa_new = rbind(cbind(x, a), cbind(x[rep(1:n,length(a_vals)),], a=rep(a_vals, rep(n,length(a_vals))) ))
  
  x = data.frame(x)
  xa_new = data.frame(xa_new)
  colnames(xa_new) = colnames(cbind(x, a))
  
  mu_mod = SuperLearner(Y = y,
                        X = cbind(x, a),
                        newX = xa_new,
                        SL.library = sl_lib)
  
  mu_function = function(a_new) {
    xa_new = cbind(x, a=a_new)
    predict(mu_mod, newdata = xa_new)$pred
  }
  
  mu_deriv = grad(mu_function, a)
  mu_deriv
}

dctseff_simu = function(y, a, x, muhat, lambdahat, n_pts = 100, a.out){
  # require(lokern)
  require(KernSmooth)

  dat = data.frame(y,a,x); n = nrow(dat)
  split_index = sample(rep(1:2, ceiling(n/2))[1:n])

  thetahats = matrix(0, nrow=2, ncol=length(a.out))
  for(split in 1:2){
    train = dat[split_index != split, ]; test = dat[split_index == split, ]
    n1 = nrow(train); n2 = nrow(test)
    
    stopifnot(n1 == n2)
    
    # a_vals_train = seq(min(train$a), max(train$a), length.out = n_pts)
    a_vals_test = seq(min(test$a), max(test$a), length.out = n_pts)
    
    x_train = subset(train, select=-c(a,y)); x_test = subset(test, select=-c(a,y))

    xa_new = cbind(x_train[rep(1:n1,length(a_vals_test)),], a=rep(a_vals_test, rep(n1,length(a_vals_test))) )
    x_new = xa_new[,-dim(xa_new)[2]]; x_new = as.matrix(x_new)

    # xa_new_test = rbind(cbind(x_test, a=test$a), cbind(x_test[rep(1:n2,length(a_vals_test)),], a=rep(a_vals_test, rep(n2,length(a_vals_test))) ))
    # x_new_test = xa_new_test[,-dim(xa_new_test)[2]]; x_new_test = as.matrix(x_new_test)
    
    ### estimate fhat and tau0hat on train ###
    lambda_pred_vals = lambdahat(x_new)
    # pi_hat_eval <- dbeta(xa_new[,dim(xa_new)[2]]/20, lambda_hat, 1-lambda_hat)/20
    pi_pred_eval <- dnorm(xa_new[,dim(xa_new)[2]], lambda_pred_vals, 1)
    pi_pred_eval <- ifelse(pi_pred_eval>0.01, pi_pred_eval, 0.01)
    # pi2_hat_train = pi2hat(x_new_train, pi_hat_train)

    # a_std_train = (xa_new_train$a - pi_hat_train) / sqrt(pi2_hat_train)
    # pihat_vals_train = approx(density(a_std_train)$x, density(a_std_train[1:n1])$y, xout = a_std_train)$y / sqrt(pi2_hat_train)
    # pi_hat = pihat_vals_train[1:n1]
    pi_pred_mat = matrix(pi_pred_eval, nrow = n1, ncol = length(a_vals_test))
    f_pred = predict(smooth.spline(a_vals_test, apply(pi_pred_mat, 2, mean)), x = test$a)$y

    mu_pred_vals = muhat(x_new, xa_new[,dim(xa_new)[2]])
    # mu_hat = mu_hat_vals_train[1:n1]
    mu_pred_mat = matrix(mu_pred_vals, nrow = n1, ncol = length(a_vals_test))
    tau0_pred = predict(smooth.spline(a_vals_test, apply(mu_pred_mat, 2, mean)), x = test$a)$y
     
    ### pseudo outcomes and local polynomial regression on test ###
    lambda_pred_test = lambdahat(as.matrix(x_test)) #; pi2_hat_test = pi2hat(x_new_test, pi_hat_test)
    pi_pred_test <- dnorm(test$a, lambda_pred_test)
    # pi_hat <- dbeta(test$a/20, lambda_hat_test, 1-lambda_hat_test)/20
    pi_pred_test <- ifelse(pi_pred_test>0.01, pi_pred_test, 0.01)
    mu_pred_test <- muhat(as.matrix(x_test), test$a)
    # a_std_test = (xa_new_test$a - pi_hat_test) / sqrt(pi2_hat_test)
    # pihat_vals_test = approx(density(a_std_test)$x, density(a_std_test[1:n2])$y, xout = a_std_test)$y / sqrt(pi2_hat_test)
    # pi_hat = pihat_vals_test[1:n2]
    # pi_hat <- ifelse(pi_hat>0.01, pi_hat, 0.01)
    # print(summary(pi_hat))
    # mu_hat_vals_test = muhat(x_new_test, test$a)
    # mu_hat = mu_hat_vals_test[1:n2]
    # print(summary(mu_hat))
    
    pseudo_out = (test$y - mu_pred_test)/pi_pred_test*f_pred + tau0_pred
    bws <- seq(from=0.02, to=3, length.out=100)
    mses <- numeric(length(bws))
    for (h in seq_along(bws)){
      model_fit <- locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=bws[h])
      res <- approx(model_fit, xout = a.out, rule=2)$y
      mses[h] <- mean((res-true.out)^2)
    }
    # print(mses)
    print(which.min(mses))
    h_opt <- bws[which.min(mses)]
    model_fit <- locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=h_opt)
    # res <- 
    # pseudo_out = test$y
    # lm_dat <- data.frame(out=pseudo_out, a=test$a)
    # lm_fit <- lm(out~a+I(a^3), data=lm_dat)
    # locpoly(test$a, pseudo_out, drv = 0, kernel = "normal")
    # model_fit <- locpol(y~x, data=data.frame(x=test$a, y=pseudo_out))
    
    # model_fit <- glkerns(x=test$a, y=pseudo_out, deriv=0)
    # dctseff_fit = lokerns(x=test$a, y=pseudo_out, deriv=0)
    # res <- approx(model_fit, xout = a)$y
    # plot(a, true, col='grey', pch=20, cex=0.3, ylab=expression(theta(a)))
    # points(a, res, col='red',cex=0.3)
    
    # res <- lpepa(x=test$a, y=pseudo_out, bandwidth=0.8, x.out=a)
    # plot(a, true, col='grey', pch=20, cex=0.3, ylab=expression(theta(a)))
    # points(a, res$est, col='red',cex=0.3)
    
    thetahats[split,] = approx(model_fit, xout = a.out, rule=2)$y
  }
  # points(a, thetahats[1,], col='red',cex=0.3)
  return(colMeans(thetahats))
}


  



