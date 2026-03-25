
library(SuperLearner)
library(foreach)
library(doParallel)

epsilon <- 0.05


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


SL.gam.new <- function(..., cts.num = 5) {
  SL.gam(..., cts.num = 5)
}


#This function estimates the derivative of the dose-response curve using the smooth approach
smooth_est = function(y, a, x, a.out, h, nsplits = 2,
                         sl.lib.ps = c("SL.gam", 
                                       "SL.glm", "SL.glm.interaction",
                                       "SL.mean", "SL.ranger"),
                         sl.lib.outcome = c("SL.gam", 
                                            "SL.glm", "SL.glm.interaction",
                                            "SL.mean", "SL.ranger"), npts_int=100,
                         family="gaussian", ifs=T){
  # family = 'gaussian' if the outcome y is continuous, family = 'binomial' if the outcome is binary
  # ifs=T if influence functions are required in the outputs, ifs are used in variance calculations
  n = length(y)
  
  #first, divide the data into
  #a train-test split. Define the split indexes:
  splitIndex = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  #we will ultimately compute the estimated \psi_h
  #for each a in a.out.
  #This will be a matrix that collects
  #the estimated EIF values for each 
  #observation and treatment combination.
  #rows are observations; columns are treatment values.
  thetah.ifvalues = matrix(0, nrow = n, ncol = length(a.out))
  
  
  #now we will repeat this process for
  #each sample split...
  for(split in 1:nsplits){
    #divide into training and testing
    #which splitIndex is equal to this split?
    train <- splitIndex != split
    test <- splitIndex == split
    if (nsplits == 1) {
      train <- test
    }
    #then, the datasets are
    y.train = y[train]
    y.test = y[test]
    a.train <- a[train]
    a.test <- a[test]
    
    a.min = min(a.test); a.max = max(a.test)
    a.seq = seq(a.min, a.max, length=npts_int)
    #it'll also be useful to have a covariate matrix
    #for the training data and for the test data.
    X.train = as.data.frame(x[train,])
    X.test = as.data.frame(x[test,])
    
    #first, estimate pi
    psMod.sl = SuperLearner(Y = a.train,
                            X = X.train, newX = X.test,
                            SL.library = sl.lib.ps)
    #thus, the mean estimate is (on the test data):
    meanA.est = psMod.sl$SL.predict
    #To flexibly estimate the variance, it'll be helpful
    #to first compute squared residuals from this model:
    meanA.est.train = as.numeric(predict(psMod.sl, newdata = X.train)$pred)
    piRes2 = (a.train - meanA.est.train)^2
    #Then, a flexible model for the variance is:
    pi2mod = SuperLearner(Y = piRes2,
                          X = X.train, newX = X.test,
                          SL.library = sl.lib.ps)
    varA.est = pi2mod$SL.predict
    varA.est.train = as.numeric(predict(pi2mod, newdata = X.train)$pred)
    #Then, \hat{\pi}(X, A) is:
    kde_eps <- density((a.train-meanA.est.train)/sqrt(varA.est.train), bw="SJ")
    a_std = (a.test - meanA.est) / sqrt(varA.est)
    piHat.XA_test <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(varA.est)
    piHat.XA_test <- ifelse(piHat.XA_test<=epsilon|is.na(piHat.XA_test), epsilon, piHat.XA_test)
    # piHat.XA = dnorm(data.test$a, mean = meanA.est, sd = sqrt(varA.est))
    
    #meanwhile, the estimated mu is...
    XA.train <- data.frame(x=x[train,],a=a[train])
    XA.test <- data.frame(x=x[test,],a=a[test])   
    if (family=='gaussian'){
      muModel = SuperLearner(Y = y[train],
                             X = XA.train,
                             newX = XA.test,
                             SL.library = sl.lib.outcome)
    }else{
      muModel = SuperLearner(Y = y[train],
                             X = XA.train,
                             newX = XA.test,
                             SL.library = sl.lib.outcome, family=binomial())
    }
    
    
    #First, compute \hat{\mu}(X, A): 
    muHat.XA_test = as.numeric(predict(muModel)$pred)
    
    #computing the IPW residuals
    ipwRes = (y[test] - muHat.XA_test)/piHat.XA_test
    
    #It'll be useful to put this into an n * length(a.eval) matrix,
    #where each column is just ipwRes
    # ipwRes.mat = matrix(rep(ipwRes, length(a.eval)), 
    #                     ncol = length(a.eval))
    
    #now produce a matrix of estimates from the mu model,
    #where the rows correspond to subjects and the columns
    #correspond to the a.seq.
    #Thus, we must create a new data.frame, where
    #each row corresponds to a subject/a.seq combination.
    #first, create the covariate data
    mu.preddata = X.test[rep(1:length(y.test),length(a.seq)),]
    #now add the treatment data, corresponding to each a.seq.
    mu.preddata$a = rep(a.seq, rep(length(y.test),length(a.seq)))
    colnames(mu.preddata) <- colnames(XA.train)
    #now predict across these a.seq
    #this will create a matrix, where
    #each row corresponds to a subject
    #and each column corresponds to a a.seq
    muHat.a0 = matrix(as.numeric(predict(muModel, newdata = mu.preddata)$pred),
                      ncol = length(a.seq))
    #Note that this is \hat{mu}(X, a_0) in the paper,
    #which is used within the integral of the estimator.
    
    #when computing the integral manually,
    #we need to know the "gap" between the a0s
    #(in order to compute a Riemann sum)
    a.seq.len = a.seq[2] - a.seq[1]
    
    #Now we just need to compute the kernel terms,
    #which we'll store in matrices.
    #kernel.A will be n * length(a.eval), where each
    #column contains K_h(A) for a particular a.
    kernel.A = sapply(a.out, gaussianKernel_deriv, A = a[test], h = h)    
    #kernel.a0 will be length(a.seq) * length(a.eval), where each
    #column contains K_h(a0) for a particular a.
    kernel.a0 = sapply(a.out, gaussianKernel_deriv, A = a.seq, h = h)
    #normalize the kernels
    # kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
    # kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
    
    #the estimator consists of two terms.
    #The first term is:
    psih.1 = kernel.A*c(ipwRes)
    
    #The second term then is:
    psih.2 = a.seq.len*muHat.a0%*%kernel.a0
    
    #then, the IF values are
    thetah.ifvalues[test,] = - psih.1 - psih.2
  }
  if(ifs==F){
    return(colMeans(thetah.ifvalues))
  }
  else{
    return(list(est=colMeans(thetah.ifvalues), ifs=thetah.ifvalues))
  }
}

#' localized kernel regression estimator with pseudo outcomes
#' with cross-fitting
dctseff = function(y, a, x, a.out, h, n_pts = 100, sl_lib = c("SL.gam", 
                                                              "SL.glm", 
                                                              "SL.mean", "SL.ranger"), family="gaussian", ifs=T, parallel=T){
  # family = 'gaussian' if the outcome y is continuous, family = 'binomial' if the outcome is binary
  # ifs=T if influence functions are required in the outputs, ifs are used in variance calculations
  
  require(SuperLearner)
  # require(lokern)
  require(KernSmooth)
  if(parallel){
    num_cores <- detectCores() - 2  # Use all but one core
    cl <- makeCluster(num_cores)  # Create cluster
    # on.exit(stopCluster(cl), add = TRUE)
    registerDoParallel(cl)  # Register the parallel backend
    on.exit(stopCluster(cl), add = TRUE)  # Ensure cleanup
  }
  
  
  n = length(y); thetahats = matrix(NA, nrow=3, ncol=length(a.out))
  # colnamex <- colnames(x)
  dat = data.frame(y=y, a=a, x)
  
  if (ifs){
    theta_ifs <- matrix(NA, nrow=n, ncol=length(a.out))
  }
  ### 1. split samples into train and test ###
  split_index = sample(1:3, size=n, replace = T)
  for(split in 1:3){
    
    train1 = dat[split_index %% 3== split%%3, ]
    train2 <- dat[split_index %% 3== (split+1)%%3,]
    test <- dat[split_index %% 3== (split+2)%%3,]
    # n_train = nrow(train)
    # shuffled_train = train[sample(n_train), ]
    # 
    # half_point = floor(n_train / 2)
    # d1 = shuffled_train[1:half_point, ]
    # d2 = shuffled_train[(half_point + 1):n_train, ]
    
    x_train1 = subset(train1, select=-c(a,y))
    x_train2 = subset(train2, select=-c(a,y)); x_test = subset(test, select=-c(a,y))
    
    n1 = nrow(train1); n2 = nrow(x_train2); n3 = nrow(test)
    
    ### 2. nuisance function training ###
    
    # set up evaluation points and matrices for predictions
    a_min = min(test$a); a_max = max(test$a)
    a_vals = seq(a_min, a_max, length.out=n_pts)
    
    # xa_new = rbind(cbind(x_test, a=test$a), 
    #                cbind(x_test[rep(1:n2,length(a_vals)),], 
    #                      a=rep(a_vals, rep(n2,length(a_vals))) ))
    
    # x_new = xa_new[,-dim(xa_new)[2]]
    # x_new = data.frame(x_new)
    
    xa_new = cbind(x_train2[rep(1:n2,length(a_vals)),], a=rep(a_vals, rep(n2,length(a_vals))) )
    
    x_new = xa_new[,-dim(xa_new)[2]]; x_new = as.matrix(x_new)
    
    
    x_train1 = data.frame(x_train1)
    xa_train1 = subset(train1, select=-c(y))
    xa_test = subset(test, select=-c(y))
    
    # colnames(x_train1) = colnames(x_new) # colnames need to match
    xa_new = data.frame(xa_new)
    # colnames(xa_train1) = colnames(xa_new)
    
    # estimate nuisance functions via super learner
    meanA_mod = SuperLearner(Y = train1$a,
                             X = x_train1, newX = rbind(x_train1, x_train2, x_test),
                             SL.library = sl_lib)
    print(paste("Number of NA in mean of treatment", sum(is.na(meanA_mod$SL.predict))))
    meanA_pred1 <- meanA_mod$SL.predict[1:n1]
    meanA_pred2 <- meanA_mod$SL.predict[(n1+1):(n1+n2)]
    meanA_pred_test <- meanA_mod$SL.predict[(n1+n2+1):(n1+n2+n3)]
    # pimod_vals = pi_mod$SL.predict
    
    
    varA_mod = SuperLearner(Y = (train1$a - meanA_pred1)^2, 
                            X = x_train1, newX = rbind(x_train1, x_train2, x_test), 
                            SL.library = sl_lib)
    print(paste("Number of NA in variance of treatment", sum(is.na(varA_mod$SL.predict))))
    varA_pred1 <- varA_mod$SL.predict[1:n1]
    varA_pred2 <- varA_mod$SL.predict[(n1+1):(n1+n2)]
    varA_pred_test <- varA_mod$SL.predict[(n1+n2+1):(n1+n2+n3)]
    
    if (family=='gaussian'){
      meanY_mod = SuperLearner(Y = train1$y,
                               X = xa_train1, newX = rbind(xa_test, xa_new),
                               SL.library = sl_lib, family=gaussian())
    }else{
      meanY_mod = SuperLearner(Y = train1$y,
                               X = xa_train1, newX = rbind(xa_test, xa_new),
                               SL.library = sl_lib, family=binomial())
    }
    
    meanY_pred = meanY_mod$SL.predict[1:nrow(xa_new)]
    meanY_pred_test <- meanY_mod$SL.predict[(1+nrow(xa_new)):(nrow(test)+nrow(xa_new))]
    print(paste("Number of NA in mean of outcome", sum(is.na(meanY_mod$SL.predict))))
    
    
    # construct estimated pi/varpi/mu values
    kde_eps <- density((train1$a-meanA_pred1)/sqrt(varA_pred1), bw="SJ")
    
    a_std = (xa_new$a - rep(meanA_pred2, n_pts)) / sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- ifelse(pi_pred_vals<=epsilon|is.na(pi_pred_vals), epsilon, pi_pred_vals)
    # pihat_vals = approx(density(a_std)$x, density(a_std[1:n2])$y, xout = a_std)$y / sqrt(pi2mod_vals)
    pi_pred_test <- approx(kde_eps$x, kde_eps$y, xout=(test$a-meanA_pred_test)/sqrt(varA_pred_test), rule = 2)$y/ sqrt(varA_pred_test)
    pi_pred_test <- ifelse(pi_pred_test<=epsilon|is.na(pi_pred_test), epsilon, pi_pred_test)
    # pihat = pihat_vals[1:n2]
    
    
    pi_pred_mat = matrix(pi_pred_vals, nrow = n2, ncol = n_pts)
    # print(dim(pi_pred_mat))
    f_pred = predict(smooth.spline(a_vals, apply(pi_pred_mat, 2, mean, na.rm=T)), x = test$a)$y
    # muhat = muhat_vals[1:n2]
    
    mu_pred_mat = matrix(meanY_pred, nrow = n2, ncol = n_pts)
    
    tau0_pred = predict(smooth.spline(a_vals, apply(mu_pred_mat, 2, mean, na.rm=T)), x = test$a)$y
    
    # mhat.mat = matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
    
    # fhat = mean(pihat/varpihat); tau0hat = mean(muhat)
    # print(paste('fhat:', fhat))
    # print(paste('tau0hat:', tau0hat))
    
    ### 3. form adjusted/pseudo outcome xi ###
    pseudo_out = (test$y - meanY_pred_test) / (pi_pred_test / f_pred) + tau0_pred
    print(paste("Number of NA in pseudo-outcomes", sum(is.na(pseudo_out))))
    
    ### 4. local polynomial regression ###
    # bws <- seq(from=0.1,to=4,by=0.06)
    # mses <- numeric(length(bws))
    # for (h in seq_along(bws)){
    #   model_fit <- locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=bws[h])
    #   res <- approx(model_fit, xout = a.out, rule = 2)$y
    #   mses[h] <- mean((res-true.out)^2)
    # }
    # h_opt <- bws[which.min(mses)]
    model_fit = locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=h) 
    thetahats[split,] = approx(model_fit, xout = a.out, rule = 2)$y
    if(ifs){
      a_points <- train2[sample(1:n2,min(100, n2), replace=F),]$a
      xa_new_var <- cbind(x_test[rep(1:n3, rep(length(a_points), n3)),], a=rep(a_points, n3))
      mu_evals <- predict(meanY_mod, newdata = xa_new_var)$pred
      mu_evals_mat <- matrix(mu_evals, nrow=length(a_points))
      if (parallel){
        ifs_test <- foreach(j=seq_along(a.out), .combine='cbind', .export='gaussianKernel') %dopar%{
          g <- matrix(0, nrow=n3, ncol=3)
          g[,1] <- 1
          g[,2] <- (test$a-a.out[j])/h
          g[,3] <- g[,2]^2
          K <- gaussianKernel(test$a, a.out[j], h)
          D_hat <- t(g)%*%diag(K)%*%g/n3
          cross_moment <- t(g) %*% (K*pseudo_out)/n3
          beta_hat <- solve(D_hat)%*%cross_moment
          term1 <- c((solve(D_hat)%*%t(g)%*%diag(K*(pseudo_out-c(g%*%beta_hat))))[2,])
          K2 <- gaussianKernel(a_points, a.out[j], h)
          g_eval <- matrix(0, nrow=length(a_points), ncol=3)
          g_eval[,1] <- 1
          g_eval[,2] <- (a_points-a.out[j])/h
          g_eval[,3] <- g_eval[,2]^2
          term2 <- c((solve(D_hat) %*% (t(g_eval) %*% (mu_evals_mat*K2)/length(a_points)))[2,])
          term3 <- c(beta_hat)[2]
          (term1+term2-term3)/h
        }
        theta_ifs[split_index %% 3== (split+2)%%3,] <- ifs_test
      }
      else{
        for (j in seq_along(a.out)){
          g <- matrix(0, nrow=n3, ncol=3)
          g[,1] <- 1
          g[,2] <- (test$a-a.out[j])/h
          g[,3] <- g[,2]^2
          K <- gaussianKernel(test$a, a.out[j], h)
          D_hat <- t(g)%*%diag(K)%*%g/n3
          cross_moment <- t(g) %*% (K*pseudo_out)/n3
          beta_hat <- solve(D_hat)%*%cross_moment
          term1 <- c((solve(D_hat)%*%t(g)%*%diag(K*(pseudo_out-c(g%*%beta_hat))))[2,])
          # a_points <- train2[sample(1:n2,min(100, n2), replace=F),]$a
          # xa_new_var <- cbind(x_test[rep(1:n3, rep(length(a_points), n3)),], a=rep(a_points, n3))
          # mu_evals <- predict(meanY_mod, newdata = xa_new_var)$pred
          # mu_evals_mat <- matrix(mu_evals, nrow=length(a_points))
          K2 <- gaussianKernel(a_points, a.out[j], h)
          g_eval <- matrix(0, nrow=length(a_points), ncol=3)
          g_eval[,1] <- 1
          g_eval[,2] <- (a_points-a.out[j])/h
          g_eval[,3] <- g_eval[,2]^2
          term2 <- c((solve(D_hat) %*% (t(g_eval) %*% (mu_evals_mat*K2)/length(a_points)))[2,])
          term3 <- c(beta_hat)[2]
          theta_ifs[split_index %% 3== (split+2)%%3,j] <- (term1+term2-term3)/h
        }
      }
      rm(meanA_mod, varA_mod, meanY_mod, mu_evals, mu_evals_mat, xa_new_var)
      gc()
    }
    
  }
  # if (parallel) stopCluster(cl)
  if (ifs){
    return(list(est=colMeans(thetahats), ifs=theta_ifs))
  }else{
    return(colMeans(thetahats))
  }
}

npliv_lp = function(y, a, x, z, z.out, h1, h2, family.y, family.a='binomial'){
  # family.y, family.a specify the type of y,a. By default a is binary in LIV curves.
  numer = dctseff(y, z, x, a.out=z.out, h1, family=family.y, ifs=T); denom = dctseff(a, z, x, a.out=z.out, h2, family=family.a, ifs=T)
  ifs1 <- t(numer$ifs)
  ifs2 <- t(denom$ifs)
  est1 <- numer$est
  est2 <- denom$est
  ifs_ratio <- ifs1/est2-ifs2*est1/(est2^2)
  sd_est <- apply(ifs_ratio, 1, sd)/sqrt(length(y))
  return(list(est=est1/est2, sd_error=sd_est))
}

npliv_smooth <- function(y, a, x, z, z.out, h1, h2, family.y, family.a='binomial'){
  # family.y, family.a specify the type of y,a. By default a is binary in LIV curves.
  numer = smooth_est(y, z, x, a.out=z.out, h1, family=family.y, ifs=T); denom = smooth_est(a, z, x, a.out=z.out, h2, family=family.a,ifs=T)
  ifs1 <- t(numer$ifs)
  ifs2 <- t(denom$ifs)
  est1 <- numer$est
  est2 <- denom$est
  ifs_ratio <- ifs1/est2-ifs2*est1/(est2^2)
  sd_est <- apply(ifs_ratio, 1, sd)/sqrt(length(y))
  return(list(est=est1/est2, sd_error=sd_est))
}


cv = function(y, a, x, hs, method, sl_lib= c("SL.gam", 
                                             "SL.glm", "SL.glm.interaction",
                                             "SL.mean", "SL.ranger"), n_pts=100, family.y='gaussian')
{
  # method can be "lp" or "smooth"
  # family.y specifies the type of the outcome y. 
  # When y is the treatment and a is the IV, set family.y='binomial' in LIV settings.
  # When y is the outcome and a is the IV, set family.y='gaussian' in LIV settings.
  dat = data.frame(y,a,x); n = nrow(dat)
  
  #sl_lib = c("SL.earth", "SL.gam", 
  #           "SL.glm", "SL.glm.interaction",
  #           "SL.mean", "SL.ranger", "SL.rpart")
  
  risks = matrix(0, nrow=length(hs), ncol=2)
  n_split = 2
  split_index = sample(rep(1:n_split, ceiling(n/n_split))[1:n])
  for (split in 1:n_split){
    train = dat[split_index != split, ]; test = dat[split_index == split, ]
    n1 = nrow(train); n2 = nrow(test)
    ## set up data for training ##
    x_train = subset(train, select=-c(a,y)); x_test = subset(test, select=-c(a,y))
    a_seq = seq(min(test$a), max(test$a), length.out = n_pts)
    ax_train = subset(train, select=-c(y)); ax_test = subset(test, select=-c(y))
    ax_new = cbind(a=rep(a_seq, rep(n2,length(a_seq))), x_test[rep(1:n2,length(a_seq)),])
    
    ## mu function used for residual term and integral 
    if(family.y=='gaussian'){
      meanY_mod = SuperLearner(Y = train$y,
                               X = ax_train, newX = rbind(ax_new, ax_test),
                               SL.library = sl_lib)
    }else{
      meanY_mod = SuperLearner(Y = train$y,
                               X = ax_train, newX = rbind(ax_new, ax_test),
                               SL.library = sl_lib, family=binomial())
    }
    
    meanY_pred = meanY_mod$SL.predict[1:nrow(ax_new)]
    meanY_pred_mat <- matrix(meanY_pred, nrow=n2) #n_2 * 100 matrix for evaluating the integral
    meanY_pred_test <- meanY_mod$SL.predict[(1+nrow(ax_new)):(nrow(test)+nrow(ax_new))]
    
    # pi function used for the residual term
    meanA_mod = SuperLearner(Y = train$a,
                             X = x_train, newX = rbind(x_train, x_test),
                             SL.library = sl_lib)
    meanA_pred <- meanA_mod$SL.predict[1:n1]
    meanA_pred_test <- meanA_mod$SL.predict[(n1+1):(n1+n2)]
    
    varA_mod = SuperLearner(Y = (train$a - meanA_pred)^2, 
                            X = x_train, newX = rbind(x_train, x_test), 
                            SL.library = sl_lib)
    varA_pred <- varA_mod$SL.predict[1:n1]
    varA_pred_test <- varA_mod$SL.predict[(n1+1):(n1+n2)]
    
    kde_eps <- density((train$a-meanA_pred)/sqrt(varA_pred), bw="SJ")
    a_std = (test$a - meanA_pred_test) / sqrt(varA_pred_test)
    pi_pred_test <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(varA_pred_test)
    pi_pred_test <- ifelse(pi_pred_test<=epsilon|is.na(pi_pred_test), epsilon, pi_pred_test)
    
    # marginal density f
    kde_a <- density(train$a, bw="SJ")
    density_a_seq <- approx(kde_a$x, kde_a$y, xout=a_seq, rule = 2)$y
    
    
    for (i in 1:length(hs)){
      ## Fit theta
      if (method=='lp'){
        theta_hat <- dctseff(y=train$y, train$a, x_train, a.out=c(test$a, a_seq), hs[i], n_pts = 100, family=family.y, ifs=F, parallel=F)
        theta_hat_test <- theta_hat[1:n2]
        theta_hat_seq <- theta_hat[(n2+1):(n2+n_pts)]
      }
      else if(method=='smooth'){
        theta_hat <- smooth_est(y=train$y, train$a, x_train, a.out=c(test$a, a_seq), hs[i], family=family.y, ifs=F)
        theta_hat_test <- theta_hat[1:n2]
        theta_hat_seq <- theta_hat[(n2+1):(n2+n_pts)]
      }
      
      
      ## Fit the derivative for product of theta and f
      spline_fit <- smooth.spline(a_seq, density_a_seq*theta_hat_seq)
      spline_deri <- predict(spline_fit, c(test$a, a_seq), deriv = 1)$y
      deriv_test <- spline_deri[1:n2]
      deriv_seq <- spline_deri[(n2+1):(n2+n_pts)]
      
      # Calculate the risk
      integral <- c(meanY_pred_mat%*%deriv_seq*(a_seq[2]-a_seq[1]))
      risks[i,split] = mean(theta_hat_test^2 + 2*(integral +  deriv_test* (test$y-meanY_pred_test)/pi_pred_test))
    }
  }
  return(list(hs=hs, risks=apply(risks, 1, mean)))
}

#' dose-response estimator with pseudo outcomes
#' with cross-fitting
drc = function(y, a, x, a.out, n_pts = 100, sl_lib = c("SL.gam", 
                                                              "SL.glm", 
                                                              "SL.mean", "SL.ranger"), family="gaussian", ifs=T){
  # family = 'gaussian' if the outcome y is continuous, family = 'binomial' if the outcome is binary
  # ifs=T if influence functions are required in the outputs, ifs are used in variance calculations
  
  require(SuperLearner)
  # require(lokern)
  require(KernSmooth)
  
  
  n = length(y); thetahats = matrix(NA, nrow=3, ncol=length(a.out))
  # colnamex <- colnames(x)
  dat = data.frame(y=y, a=a, x)
  
  if (ifs){
    theta_ifs <- matrix(NA, nrow=n, ncol=length(a.out))
  }
  ### 1. split samples into train and test ###
  split_index = sample(1:3, size=n, replace = T)
  for(split in 1:3){
    
    train1 = dat[split_index %% 3== split%%3, ]
    train2 <- dat[split_index %% 3== (split+1)%%3,]
    test <- dat[split_index %% 3== (split+2)%%3,]
    # n_train = nrow(train)
    # shuffled_train = train[sample(n_train), ]
    # 
    # half_point = floor(n_train / 2)
    # d1 = shuffled_train[1:half_point, ]
    # d2 = shuffled_train[(half_point + 1):n_train, ]
    
    x_train1 = subset(train1, select=-c(a,y))
    x_train2 = subset(train2, select=-c(a,y)); x_test = subset(test, select=-c(a,y))
    
    n1 = nrow(train1); n2 = nrow(x_train2); n3 = nrow(test)
    
    ### 2. nuisance function training ###
    
    # set up evaluation points and matrices for predictions
    a_min = min(test$a); a_max = max(test$a)
    a_vals = seq(a_min, a_max, length.out=n_pts)
    
    # xa_new = rbind(cbind(x_test, a=test$a), 
    #                cbind(x_test[rep(1:n2,length(a_vals)),], 
    #                      a=rep(a_vals, rep(n2,length(a_vals))) ))
    
    # x_new = xa_new[,-dim(xa_new)[2]]
    # x_new = data.frame(x_new)
    
    xa_new = cbind(x_train2[rep(1:n2,length(a_vals)),], a=rep(a_vals, rep(n2,length(a_vals))) )
    
    x_new = xa_new[,-dim(xa_new)[2]]; x_new = as.matrix(x_new)
    
    
    x_train1 = data.frame(x_train1)
    xa_train1 = subset(train1, select=-c(y))
    xa_test = subset(test, select=-c(y))
    
    # colnames(x_train1) = colnames(x_new) # colnames need to match
    xa_new = data.frame(xa_new)
    # colnames(xa_train1) = colnames(xa_new)
    
    # estimate nuisance functions via super learner
    meanA_mod = SuperLearner(Y = train1$a,
                             X = x_train1, newX = rbind(x_train1, x_train2, x_test),
                             SL.library = sl_lib)
    print(paste("Number of NA in mean of treatment", sum(is.na(meanA_mod$SL.predict))))
    meanA_pred1 <- meanA_mod$SL.predict[1:n1]
    meanA_pred2 <- meanA_mod$SL.predict[(n1+1):(n1+n2)]
    meanA_pred_test <- meanA_mod$SL.predict[(n1+n2+1):(n1+n2+n3)]
    # pimod_vals = pi_mod$SL.predict
    
    
    varA_mod = SuperLearner(Y = (train1$a - meanA_pred1)^2, 
                            X = x_train1, newX = rbind(x_train1, x_train2, x_test), 
                            SL.library = sl_lib)
    print(paste("Number of NA in variance of treatment", sum(is.na(varA_mod$SL.predict))))
    varA_pred1 <- varA_mod$SL.predict[1:n1]
    varA_pred2 <- varA_mod$SL.predict[(n1+1):(n1+n2)]
    varA_pred_test <- varA_mod$SL.predict[(n1+n2+1):(n1+n2+n3)]
    
    if (family=='gaussian'){
      meanY_mod = SuperLearner(Y = train1$y,
                               X = xa_train1, newX = rbind(xa_test, xa_new),
                               SL.library = sl_lib, family=gaussian())
    }else{
      meanY_mod = SuperLearner(Y = train1$y,
                               X = xa_train1, newX = rbind(xa_test, xa_new),
                               SL.library = sl_lib, family=binomial())
    }
    
    meanY_pred = meanY_mod$SL.predict[1:nrow(xa_new)]
    meanY_pred_test <- meanY_mod$SL.predict[(1+nrow(xa_new)):(nrow(test)+nrow(xa_new))]
    print(paste("Number of NA in mean of outcome", sum(is.na(meanY_mod$SL.predict))))
    
    
    # construct estimated pi/varpi/mu values
    kde_eps <- density((train1$a-meanA_pred1)/sqrt(varA_pred1), bw="SJ")
    
    a_std = (xa_new$a - rep(meanA_pred2, n_pts)) / sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- ifelse(pi_pred_vals<=epsilon|is.na(pi_pred_vals), epsilon, pi_pred_vals)
    # pihat_vals = approx(density(a_std)$x, density(a_std[1:n2])$y, xout = a_std)$y / sqrt(pi2mod_vals)
    pi_pred_test <- approx(kde_eps$x, kde_eps$y, xout=(test$a-meanA_pred_test)/sqrt(varA_pred_test), rule = 2)$y/ sqrt(varA_pred_test)
    pi_pred_test <- ifelse(pi_pred_test<=epsilon|is.na(pi_pred_test), epsilon, pi_pred_test)
    # pihat = pihat_vals[1:n2]
    
    
    pi_pred_mat = matrix(pi_pred_vals, nrow = n2, ncol = n_pts)
    # print(dim(pi_pred_mat))
    f_pred = predict(smooth.spline(a_vals, apply(pi_pred_mat, 2, mean, na.rm=T)), x = test$a)$y
    # muhat = muhat_vals[1:n2]
    
    mu_pred_mat = matrix(meanY_pred, nrow = n2, ncol = n_pts)
    
    tau0_pred = predict(smooth.spline(a_vals, apply(mu_pred_mat, 2, mean, na.rm=T)), x = test$a)$y
    
    # mhat.mat = matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
    
    # fhat = mean(pihat/varpihat); tau0hat = mean(muhat)
    # print(paste('fhat:', fhat))
    # print(paste('tau0hat:', tau0hat))
    
    ### 3. form adjusted/pseudo outcome xi ###
    pseudo_out = (test$y - meanY_pred_test) / (pi_pred_test / f_pred) + tau0_pred
    print(paste("Number of NA in pseudo-outcomes", sum(is.na(pseudo_out))))
    
    ### 4. local polynomial regression ###
    # bws <- seq(from=0.1,to=4,by=0.06)
    # mses <- numeric(length(bws))
    # for (h in seq_along(bws)){
    #   model_fit <- locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=bws[h])
    #   res <- approx(model_fit, xout = a.out, rule = 2)$y
    #   mses[h] <- mean((res-true.out)^2)
    # }
    # h_opt <- bws[which.min(mses)]
    h <- dpill(test$a, pseudo_out)
    model_fit = locpoly(x=test$a, y=pseudo_out, bandwidth=h) 
    thetahats[split,] = approx(model_fit, xout = a.out, rule = 2)$y
    if(ifs){
        for (j in seq_along(a.out)){
          g <- matrix(0, nrow=n3, ncol=2)
          g[,1] <- 1
          g[,2] <- (test$a-a.out[j])/h
          K <- gaussianKernel(test$a, a.out[j], h)
          D_hat <- t(g)%*%diag(K)%*%g/n3
          cross_moment <- t(g) %*% (K*pseudo_out)/n3
          beta_hat <- solve(D_hat)%*%cross_moment
          term1 <- c((solve(D_hat)%*%t(g)%*%diag(K*(pseudo_out-c(g%*%beta_hat))))[1,])
          a_points <- train2[sample(1:n2,min(500, n2), replace=F),]$a
          xa_new_var <- cbind(x_test[rep(1:n3, rep(length(a_points), n3)),], a=rep(a_points, n3))
          mu_evals <- predict(meanY_mod, newdata = xa_new_var)$pred
          mu_evals_mat <- matrix(mu_evals, nrow=length(a_points))
          K2 <- gaussianKernel(a_points, a.out[j], h)
          g_eval <- matrix(0, nrow=length(a_points), ncol=2)
          g_eval[,1] <- 1
          g_eval[,2] <- (a_points-a.out[j])/h
          term2 <- c((solve(D_hat) %*% (t(g_eval) %*% (mu_evals_mat*K2)/length(a_points)))[1,])
          term3 <- c(beta_hat)[1]
          theta_ifs[split_index %% 3== (split+2)%%3,j] <- term1+term2-term3
        }
    }
  }
  if (ifs){
    return(list(est=colMeans(thetahats), ifs=theta_ifs))
  }else{
    return(colMeans(thetahats))
  }
}


comp_class = function(y, a, x, z, z.out, family.y, family.a='binomial'){
  # Compute a 2 by 2 matrix for treatment effects among the maximal complier class
  # (1,1) element is the estimated proportion of the complier
  # (2,1) element is the estimated variance for the proportion of the complier
  # (1,2) element is the estimated effects among the complier class
  # (2,2) element is the estimated variance for the estimated effects
  # family.y, family.a specify the type of y,a. By default a is binary.
  
  numer = drc(y,z,x,z.out, ifs=T, family=family.y); denom = drc(a,z,x,z.out, ifs=T, family=family.a)
  res <- matrix(NA, nrow=2, ncol=2)
  colnames(res) <- c('proportion', 'effect')
  rownames(res) <- c('est','sd')
  ifs1 <- c(numer$ifs[,2] - numer$ifs[,1]) # IF of the numerator
  ifs2 <- c(denom$ifs[,2] - denom$ifs[,1]) # IF of the denominator
  est1 <- numer$est[2]-numer$est[1]
  est2 <- denom$est[2]-denom$est[1]
  res[1,1] <- est2
  res[2,1] <- sd(ifs2)/sqrt(length(y))
  res[1,2] <- est1/est2
  ifs_ratio <- ifs1/est2-ifs2*est1/(est2^2)
  res[2,2] <- sd(ifs_ratio)/sqrt(length(y))
  return(res)
}


# Example for confidence intervals for the derivative function

simulate_data_normal = function(seed, n, beta1, beta2, beta_pi){
  require(MASS)
  
  set.seed(seed)
  
  x = mvrnorm(n, rep(0,4), diag(4))
  lambda_x = -0.8 + x %*% beta_pi
  # lambda_x = expit(lambda_x)
  
  a = rnorm(n, lambda_x, 1)
  # mu = 1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3]) # or some nonlinear g(a)
  mu = 1 + x %*% beta1 + a*(beta2[1] + beta2[2]*x[,1] + beta2[3]*x[,3] - beta2[4]^2*a^2) # or some nonlinear g(a)
  y = mu + rnorm(n,0,2)
  
  return(data.frame(y=y, a=a, x=x))
}


beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0.1, -0.1, 0.1, 0.13); beta_pi = c(0.1, 0.1, -0.1, 0.2)

n <- 10000
thre <- 0.1
# Set up some evaluation data 
a.out <- simulate_data_normal(521, 125, beta1, beta2, beta_pi)$a
lb = quantile(a.out, thre)
ub = quantile(a.out, 1-thre)
idx = which(a.out >= lb & a.out <= ub)
a.out <- a.out[idx]
a.out <- a.out[order(a.out)]
# true.out = 1 +  a.out*(beta2[1]  - beta2[4]^2*a.out^2)
true.out = beta2[1] - 3*beta2[4]^2*a.out^2

# DGP and estimation
dat = simulate_data_normal(2024, n, beta1, beta2, beta_pi)
y = dat$y; a = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4)
order_a = order(a)

a = a[order_a]
y = y[order_a]
x = x[order_a, ]


bws <- c(0.1,0.5,1,2,3,4,5,8,10)
mses_lp <- numeric(length(bws))
for (ii in seq_along(bws)){
  temp_lp <- dctseff(y,a,x,a.out, bws[ii],ifs=F)
  mses_lp[ii] <- mean((temp_lp-true.out)^2)
}
# print(mses)
print(which.min(mses_lp))
h_opt <- bws[which.min(mses_lp)]
res_lp <- dctseff(y,a,x,a.out, h_opt,ifs=T)

library(ggplot2)

# Create a data frame with three curves
data_plot <- data.frame(
  a = rep(a.out, 3),  # Repeat x values for each curve
  y = c(res_lp$est, res_lp$est-1.96*apply(res_lp$ifs, 2, sd)/sqrt(n), res_lp$est+1.96*apply(res_lp$ifs, 2, sd)/sqrt(n)), 
  Curve = rep(c("Estimates", "Lower CI", "Upper CI"), each = length(a.out))  # Label curves
)


custom_colors <- c("Estimates" = "black", "Lower CI" = "blue", "Upper CI" = "blue")


ggplot(data_plot, aes(x = a, y = y, color = Curve)) +
  geom_line(linewidth = 1.2) +  
  scale_color_manual(values = custom_colors) +  
  labs(
    title = "Plot of Estimates with 95% CI",
    x = "Instrument",
    y = "Outcome",
    color = "Lines"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


bws <- c(0.1,0.5,1,2,3,4,5,8,10)
mses_sm <- numeric(length(bws))
for (ii in seq_along(bws)){
  temp_sm <- smooth_est(y,a,x,a.out, bws[ii],ifs=F)
  mses_sm[ii] <- mean((temp_sm-true.out)^2)
}
# print(mses)
print(which.min(mses_sm))
h_opt_sm <- bws[which.min(mses_sm)]

res_sm <- smooth_est(y,a,x,a.out, h_opt_sm ,ifs=T)

data_plot <- data.frame(
  a = rep(a.out, 3),  # Repeat x values for each curve
  y = c(res_sm$est, res_sm$est-1.96*apply(res_sm$ifs, 2, sd)/sqrt(n), res_sm$est+1.96*apply(res_sm$ifs, 2, sd)/sqrt(n)), 
  Curve = rep(c("Estimates", "Lower CI", "Upper CI"), each = length(a.out))  # Label curves
)


custom_colors <- c("Estimates" = "black", "Lower CI" = "blue", "Upper CI" = "blue")


ggplot(data_plot, aes(x = a, y = y, color = Curve)) +
  geom_line(linewidth = 1.2) +  
  scale_color_manual(values = custom_colors) +  
  labs(
    title = "Plot of Estimates with 95% CI",
    x = "Instrument",
    y = "Outcome",
    color = "Lines"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


# Estimate the LIV curve
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


beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0, -0.1, 0.1, 0.13)
beta_pi1 = c(0.1, 0.1, -0.1, 0.2); beta_pi2 = c(0.1, -0.2, 0.3, 0.1)

# Set up evaluation data
thre <- 0.1
z.out <- simulate_data(521, beta1, beta2, beta_pi1, beta_pi2, 125)$z
lb = quantile(z.out, thre)
ub = quantile(z.out, 1-thre)
idx = which(z.out >= lb & z.out <= ub)
z.out <- z.out[idx]
z.out <- z.out[order(z.out)]
true_num.out <- -3*0.13^2*z.out^2
true_den.out <- rep(0.1,length(z.out))
true.out = -0.507*z.out^2

# DGP and estimation
dat = simulate_data(666, beta1, beta2, beta_pi1, beta_pi2, 10000)
y = dat$y; a = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4); x = as.matrix(x); z = dat$z
order_z = order(z)
a = a[order_z]; y = y[order_z]; x = x[order_z, ]; z = z[order_z]


liv_sm_est <- npliv_smooth(y, a, x, z, z.out, 1, 0.3, family.y='gaussian', family.a='gaussian')

data_plot <- data.frame(
  z = rep(z.out, 3),  # Repeat x values for each curve
  y = c(liv_sm_est$est, liv_sm_est$est-1.96*liv_sm_est$sd_error/sqrt(n), liv_sm_est$est+1.96*liv_sm_est$sd_error/sqrt(n)), 
  Curve = rep(c("Estimates", "Lower CI", "Upper CI"), each = length(z.out))  # Label curves
)


custom_colors <- c("Estimates" = "black", "Lower CI" = "blue", "Upper CI" = "blue")


ggplot(data_plot, aes(x = z, y = y, color = Curve)) +
  geom_line(linewidth = 1.2) +  
  scale_color_manual(values = custom_colors) +  
  labs(
    title = "Plot of Estimates with 95% CI",
    x = "Instrument",
    y = "Outcome",
    color = "Lines"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Centered title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


# Estimate the maximal complier class and the treatment effects
beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0.1, -0.1, 0.1, 0.13); beta_pi = c(0.1, 0.1, -0.1, 0.2)

n <- 10000

dat = simulate_data_normal(2024, n, beta1, beta2, beta_pi)
y = dat$y; z = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4)
order_z = order(z)

z = z[order_z]
y = y[order_z]
x = x[order_z, ]

z.out <- c(-2,0.5) # Set the evaluation points to be the minimum and maximum of the range
res_complier <- comp_class(y=y,a=y,x,z=z,z.out=z.out, family.y='gaussian', family.a='gaussian') # Also set a as y so the true estimand should be 1

print(res_complier)
c(res_complier[1,1]-1.96*res_complier[2,1], res_complier[1,1]+1.96*res_complier[2,1]) # CI for the proportion of compliers
c(res_complier[1,2]-1.96*res_complier[2,2], res_complier[1,2]+1.96*res_complier[2,2]) # CI for the treatment effects








