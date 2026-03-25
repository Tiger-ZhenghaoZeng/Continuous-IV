
library(MASS)
library(KernSmooth)
library(splines)
library(utils)

rmse = function(true, estimates){
  sq_errors = sweep(estimates, 2, true, "-")^2
  sq_errors[sq_errors>10] <- NA
  return(mean(sqrt(colMeans(sq_errors, na.rm=T)), na.rm=T))
}

gaussianKernel_deriv = function(A, a, h){
  u = A - a
  kernel = -(u/h^2)*(1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
  return(kernel)
}

epan_deriv <- function(A, a, h){
  u = A - a
  kernel <- -3*u/(2*h^3)*I(abs(u)<h)
  return(kernel)
}


simulate_data = function(seed, beta1, beta2, beta_pi1, beta_pi2, n) {
  set.seed(seed)

  x = mvrnorm(n, rep(0,4), diag(4))
  lambda_x = 2 + x %*% beta_pi1

  z = rnorm(n, lambda_x, 1)

  # shouldn't a be binary?
  a = 1 + x %*% beta_pi2 + 0.1*z + rnorm(n)

  y = 1 + x %*% beta1 + z*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*z^2) + rnorm(n,0,1)

  return(data.frame(y=y,a=a,x=x,z=z))
}

meanY_hat = function(beta1, beta2, epsilon){
  return(function(x,z) {
    return(c(1 + x %*% beta1 + z*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*(1+epsilon)*z^2)))
  })
}

meanZ_hat = function(beta_pi1, epsilon){
  return(function(x) {
    return(c(2 + x %*% beta_pi1 + epsilon))
  })
}

meanA_hat = function(beta_pi2, epsilon){
  return(function(x,z) {
    return(c(1 + x %*% beta_pi2 + 0.1*z + epsilon))
  })
}

estthetah_numer_simu = function(y, z, x, h, z.eval, muhat, lambdahat, kernelfun, npts_int=100){

  n = length(y)
  thetah.ifvalues = matrix(0, nrow = n, ncol = length(z.eval))
  
  z.min = min(z); z.max = max(z)
  z.seq = seq(z.min, z.max, length=npts_int)
  x = as.matrix(x)
  
  meanZ.est = lambdahat(x)
  piHat.XZ = dnorm(z, mean = meanZ.est, sd = 1)
  piHat.XZ = ifelse(piHat.XZ>0.01, piHat.XZ, 0.01)
  muHat.XZ = muhat(x,z)
  ipwRes = (y - muHat.XZ)/piHat.XZ
  kernel.Z = sapply(z.eval, kernelfun, A = z, h = h)
  
  xz_new = cbind(x[rep(1:n,length(z.seq)),], z=rep(z.seq, rep(n,length(z.seq))))
  mu_pred_vals = muhat(xz_new[,1:(dim(xz_new)[2]-1)], xz_new[,dim(xz_new)[2]])
  muHat.z0 = matrix(mu_pred_vals, nrow = n, ncol = length(z.seq))
  z.seq.len = z.seq[2] - z.seq[1]
  kernel.z0 = sapply(z.eval, kernelfun, A = z.seq, h = h)
  
  # normalization
  # kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  # kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  psih.1 = kernel.Z*ipwRes
  psih.2 = z.seq.len*muHat.z0%*%kernel.z0
  
  # thetah.ifvalues = - psih.1 - psih.2
  thetah.ifvalues = -psih.1 - psih.2
  return(colMeans(thetah.ifvalues))
}

estthetah_denom_simu = function(a, z, x, h, z.eval, etahat, lambdahat, kernelfun, npts_int=100){
  
  n = length(y)
  #we will ultimately compute the estimated \psi_h
  #for each a in a.eval.
  #This will be a matrix that collects
  #the estimated EIF values for each
  #observation and treatment combination.
  #rows are observations; columns are treatment values.
  thetah.ifvalues = matrix(0, nrow = n, ncol = length(z.eval))
  
  z.min = min(z); z.max = max(z)
  z.seq = seq(z.min, z.max, length=npts_int)
  x = as.matrix(x)
  
  meanZ.est = lambdahat(x)
  piHat.XZ = dnorm(z, mean = meanZ.est, sd = 1)
  piHat.XZ = ifelse(piHat.XZ>0.01, piHat.XZ, 0.01)
  etaHat.XZ = etahat(x,z)
  ipwRes = (a - etaHat.XZ)/piHat.XZ
  kernel.Z = sapply(z.eval, kernelfun, A = z, h = h)
  
  xz_new = cbind(x[rep(1:n,length(z.seq)),], z=rep(z.seq, rep(n,length(z.seq))))
  eta_pred_vals = etahat(xz_new[,1:(dim(xz_new)[2]-1)], xz_new[,dim(xz_new)[2]])
  etaHat.z0 = matrix(eta_pred_vals, nrow = n, ncol = length(z.seq))
  z.seq.len = z.seq[2] - z.seq[1]
  kernel.z0 = sapply(z.eval, kernelfun, A = z.seq, h = h)
  
  # normalization
  # kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  # kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  psih.1 = kernel.Z*ipwRes
  psih.2 = z.seq.len*etaHat.z0%*%kernel.z0
  
  thetah.ifvalues = -psih.1 - psih.2
  return(colMeans(thetah.ifvalues))
}

dctseff_numer_simu = function(y, z, x, muhat, lambdahat, n_pts = 100, z.out, h){
  
  dat = data.frame(y,z,x); n = nrow(dat)
  split_index = sample(rep(1:2, ceiling(n/2))[1:n])
  
  thetahats = matrix(0, nrow=2, ncol=length(z.out))
  for(split in 1:2){
    train = dat[split_index != split, ]; test = dat[split_index == split, ]
    n1 = nrow(train); n2 = nrow(test)
    
    stopifnot(n1 == n2)
    
    z_vals_test = seq(min(test$z), max(test$z), length.out = n_pts)
    
    x_train = subset(train, select=-c(z,y)); x_test = subset(test, select=-c(z,y))
    
    xz_new = cbind(x_train[rep(1:n1,length(z_vals_test)),], z=rep(z_vals_test, rep(n1,length(z_vals_test))) )
    x_new = xz_new[,-dim(xz_new)[2]]; x_new = as.matrix(x_new)
    
    ### estimate fhat and tau0hat on train ###
    lambda_pred_vals = lambdahat(x_new)
    pi_pred_eval = dnorm(xz_new[,dim(xz_new)[2]], lambda_pred_vals, 1)
    pi_pred_eval = ifelse(pi_pred_eval>0.01, pi_pred_eval, 0.01)
    pi_pred_mat = matrix(pi_pred_eval, nrow = n1, ncol = length(z_vals_test))
    f_pred = predict(smooth.spline(z_vals_test, apply(pi_pred_mat, 2, mean)), x = test$z)$y
    
    mu_pred_vals = muhat(x_new, xz_new[,dim(xz_new)[2]])
    mu_pred_mat = matrix(mu_pred_vals, nrow = n1, ncol = length(z_vals_test))
    tau0_pred = predict(smooth.spline(z_vals_test, apply(mu_pred_mat, 2, mean)), x = test$z)$y
    
    lambda_pred_test = lambdahat(as.matrix(x_test)) 
    pi_pred_test = dnorm(test$z, lambda_pred_test)
    pi_pred_test = ifelse(pi_pred_test>0.01, pi_pred_test, 0.01)
    mu_pred_test = muhat(as.matrix(x_test), test$z)
    
    pseudo_out = (test$y - mu_pred_test)/pi_pred_test*f_pred + tau0_pred
    model_fit = locpoly(x=test$z, y=pseudo_out, drv=1, bandwidth=h)
    thetahats[split,] = approx(model_fit, xout = z.out)$y
  }
  return(colMeans(thetahats))
}

dctseff_denom_simu = function(a, z, x, etahat, lambdahat, n_pts=100, z.out, h){
  
  dat = data.frame(a,z,x); n = nrow(dat)
  split_index = sample(rep(1:2, ceiling(n/2))[1:n])
  
  thetahats = matrix(0, nrow=2, ncol=length(z.out))
  for(split in 1:2){
    train = dat[split_index != split, ]; test = dat[split_index == split, ]
    n1 = nrow(train); n2 = nrow(test)
    
    stopifnot(n1 == n2)
    
    z_vals_test = seq(min(test$z), max(test$z), length.out = n_pts)
    
    x_train = subset(train, select=-c(z,a)); x_test = subset(test, select=-c(z,a))
    
    xz_new = cbind(x_train[rep(1:n1,length(z_vals_test)),], z=rep(z_vals_test, rep(n1,length(z_vals_test))) )
    x_new = xz_new[,-dim(xz_new)[2]]; x_new = as.matrix(x_new)
    
    ### estimate fhat and tau0hat on train ###
    lambda_pred_vals = lambdahat(x_new)
    pi_pred_eval = dnorm(xz_new[,dim(xz_new)[2]], lambda_pred_vals, 1)
    pi_pred_eval = ifelse(pi_pred_eval>0.01, pi_pred_eval, 0.01)
    pi_pred_mat = matrix(pi_pred_eval, nrow = n1, ncol = length(z_vals_test))
    f_pred = predict(smooth.spline(z_vals_test, apply(pi_pred_mat, 2, mean)), x = test$z)$y
    
    eta_pred_vals = etahat(x_new, xz_new[,dim(xz_new)[2]])
    eta_pred_mat = matrix(eta_pred_vals, nrow = n1, ncol = length(z_vals_test))
    tau0_pred = predict(smooth.spline(z_vals_test, apply(eta_pred_mat, 2, mean)), x = test$z)$y

    lambda_pred_test = lambdahat(as.matrix(x_test)) 
    pi_pred_test = dnorm(test$z, lambda_pred_test)
    pi_pred_test = ifelse(pi_pred_test>0.01, pi_pred_test, 0.01)
    eta_pred_test = etahat(as.matrix(x_test), test$z)
    
    pseudo_out = (test$a - eta_pred_test)/pi_pred_test*f_pred + tau0_pred
  
    model_fit = locpoly(x=test$z, y=pseudo_out, drv=1, bandwidth=h)
    thetahats[split,] = approx(model_fit, xout = z.out)$y
  }
  return(colMeans(thetahats))
}

npliv_dr_simu = function(y, z, a, x, n_pts=100, z.out, muhat, lambdahat, etahat, true.out){

  # bws = expand.grid(
  #   bws.y = seq(from=1, to=3, length.out=10),
  #   bws.a = seq(from=1, to=3, length.out=10)
  # )
  # 
  # mses = numeric(nrow(bws))
  # 
  # for (i in 1:nrow(bws)){
  #   print(paste0('bandwidths for y, a: ', bws[i,1], ', ', bws[i,2]))
  #   numer = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i,1])
  #   denom = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i,2])
  #   mses[i] = mean((numer/denom-true.out)^2, na.rm=T)
  # }
  # 
  # h_opt = bws[which.min(mses),]
  # print('optimal bandwidths:')
  # print(paste(h_opt))
  # 
  # numer_opt = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=h_opt[1]$bws.y)
  # denom_opt = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=h_opt[2]$bws.a)
  # 
  # return(numer_opt/denom_opt)
  bws <- seq(from=0.02, to=1.5, by=0.01)
  mses_num = numeric(length(bws))
  for (i in 1:length(bws)){
    # print(paste0('bandwidths for y ', bws[i,1], ', ', bws[i,2]))
    numer = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i])
    # denom = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i,2])
    mses_num[i] = mean((numer-true_num.out)^2)
  }
  hy_opt <- bws[which.min(mses_num)]
  print(paste('optimal bandwidth for y:', hy_opt))
  numer_opt = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=hy_opt)
  mses_denom = numeric(length(bws))
  for (i in 1:length(bws)){
    # print(paste0('bandwidths for y ', bws[i,1], ', ', bws[i,2]))
    # numer = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i])
    denom = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i])
    mses_denom[i] = mean((denom-true_den.out)^2)
  }
  ha_opt <- bws[which.min(mses_denom)]
  print(paste('optimal bandwidth for a:', ha_opt))
  denom_opt = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=ha_opt)
  return(numer_opt/denom_opt)
}

npliv_smooth_simu = function(y, z, a, x, n_pts=100, z.out, muhat, lambdahat, etahat, kernelfun, true.out){
  
  # bws = expand.grid(
  #   bws.y = seq(from=1, to=3, length.out=10),
  #   bws.a = seq(from=1, to=3, length.out=10)
  # )
  # 
  # mses = numeric(nrow(bws))
  
  # for (i in 1:nrow(bws)){
  #     print(paste0('bandwidths for y, a: ', bws[i,1], ', ', bws[i,2]))
  #     numer = estthetah_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=bws[i,1])
  #     denom = estthetah_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=bws[i,2])
  #     mses[i] = mean((numer/denom-true.out)^2, na.rm=T)
  # }
  # 
  # h_opt = bws[which.min(mses),]
  # print('optimal bandwidths:')
  # print(paste(h_opt))
  # 
  # numer_opt = estthetah_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=h_opt[1]$bws.y)
  # denom_opt = estthetah_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=h_opt[2]$bws.a)
  # 
  # return(numer_opt/denom_opt)
  # Choose optimal bandwidth for numerator and denominator separately to increase speed
  bws <- seq(from=0.02, to=3, by=0.02)
  mses_num = numeric(length(bws))
  for (i in 1:length(bws)){
    # print(paste0('bandwidths for y ', bws[i,1], ', ', bws[i,2]))
    numer = estthetah_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=bws[i])
    # denom = dctseff_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i,2])
    mses_num[i] = mean((numer-true_num.out)^2)
  }
  hy_opt <- bws[which.min(mses_num)]
  print(paste('optimal bandwidth for y (smooth approximation):', hy_opt))
  numer_opt = estthetah_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=hy_opt)
  mses_denom = numeric(length(bws))
  for (i in 1:length(bws)){
    # print(paste0('bandwidths for y ', bws[i,1], ', ', bws[i,2]))
    # numer = dctseff_numer_simu(y=y, z=z, x=x, muhat=muhat, lambdahat=lambdahat, n_pts=n_pts, z.out=z.out, h=bws[i])
    denom = estthetah_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=bws[i])
    mses_denom[i] = mean((denom-true_den.out)^2)
  }
  ha_opt <- bws[which.min(mses_denom)]
  print(paste('optimal bandwidth for a (smooth approximation)::', ha_opt))
  denom_opt = estthetah_denom_simu(a=a, z=z, x=x, etahat=etahat, lambdahat=lambdahat, kernelfun=kernelfun, npts_int=n_pts, z.eval=z.out, h=ha_opt)
  return(numer_opt/denom_opt)
}

##### need to work on this #####
spliv_dr_simu = function(y, a, z, x, lambdahat,
                         muhat,
                         etahat, n_pts=100) {
  
  # dat = data.frame(y=y, a=a, z=z, x)
  
  n <- length(y)
  w.fn = function(tval) {
    dnorm(tval, 2, sqrt(1.07))
  }
  w.dz.fn = function(tval) {
    - ((tval-2) / 1.07) * dnorm(tval, 2, sqrt(1.07))
  }
  
  w  = w.fn(z)
  wdz = w.dz.fn(z)
  
  mu_pred = muhat(x,z)
  eta_pred <- etahat(x,z)
  lambda_pred <- lambdahat(x)
  pi_pred <- dnorm(z, lambda_pred, 1)
  pi_pred <- ifelse(pi_pred>0.01, pi_pred, 0.01)
  
  z.grid = seq(min(z), max(z), length.out=n_pts)
  dz <- z.grid[2]-z.grid[1]
  
  w.grid  = w.fn(z.grid)
  wdz.grid = w.dz.fn(z.grid)
  
  # m.zvals = numeric(length(z.grid))
  # l.zvals = numeric(length(z.grid))

  # for(i in seq_along(z.grid)){
  #   zi = z.grid[i]
  #   newdata_i = as.data.frame(cbind(z=zi, x))
  #   names(newdata_i) = c("z", colnames(x))
  #   pred_m_i = muhat(zi, x)
  #   # pred_l_i = etahat(zi, x)
  #   # pred_m_i = predict(mfit, newdata=newdata_i)$pred
  #   # pred_l_i = predict(lamod, newdata=newdata_i)$pred
  #   m.zvals[i] = mean(pred_m_i)
  #   # l.zvals[i] = mean(pred_l_i)
  # }
  xz_new = cbind(x[rep(1:n,length(z.grid)),], a=rep(z.grid, rep(n,length(z.grid))) )
  x_new = xz_new[,-dim(xz_new)[2]]; x_new = as.matrix(x_new)
  
  mu_vals = muhat(x_new, xz_new[,dim(xz_new)[2]])
  mu_vals_mat = matrix(mu_vals, nrow = n, ncol = length(z.grid))
  m_vals <- apply(mu_vals_mat, 2, mean)
  meana_vals = etahat(x_new, xz_new[,dim(xz_new)[2]])
  meana_vals_mat <- matrix(meana_vals, nrow = n, ncol = length(z.grid))
  l_vals <- apply(meana_vals_mat, 2, mean)
  
  I1 <- sum((w.fn(z.grid)+z.grid*w.dz.fn(z.grid))*m_vals*dz)
  I2 <- sum((2*z.grid*w.fn(z.grid)+z.grid^2*w.dz.fn(z.grid))*l_vals*dz)
  # sm_m  = smooth.spline(z.grid, m.zvals, df=spline_df)
  # sm_l  = smooth.spline(z.grid, l.zvals, df=spline_df)
  # 
  # mhat.fn = function(tt) predict(sm_m, x=tt)$y
  # lhat.fn = function(tt) predict(sm_l, x=tt)$y
  # 
  # g2.fn = function(tt) (w.fn(tt) + tt*w.dz.fn(tt)) * mhat.fn(tt)
  # g1.fn = function(tt) tt*(2*w.fn(tt) + tt*w.dz.fn(tt)) * lhat.fn(tt)
  # 
  # lims = c(min(z), max(z))
  # reg.num = integrate(function(u) sapply(u, g2.fn), lower=lims[1], upper=lims[2])$value
  # reg.den = integrate(function(u) sapply(u, g1.fn), lower=lims[1], upper=lims[2])$value
  
  y_resid = y - mu_pred
  a_resid = a - eta_pred
  
  ipw_num = mean((w + z*wdz) * y_resid / pi_pred)
  ipw_den = mean(z*(2*w + z*wdz) * a_resid / pi_pred)
  
  dr_est = (ipw_num + I1) / (ipw_den + I2)
  
  return(dr_est)
}
