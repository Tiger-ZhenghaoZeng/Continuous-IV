
library(MASS)
library(KernSmooth)
library(splines)
library(SuperLearner)
library(parallel)
library(foreach)
library(doParallel)

pos <- function(v, eps = 1e-6) {
  v[!is.finite(v)] <- eps
  pmax(v, eps)
}

gaussianKernel <- function(A, a, h){
  u <- A - a
  (1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
}

gaussianKernel_deriv <- function(A, a, h){
  u <- A - a
  -(u/h^2)*(1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
}

var_from_ifs <- function(ifs_mat, n_eff = nrow(ifs_mat)) {
  apply(ifs_mat, 2, var, na.rm = TRUE) / n_eff
}

simulate_data_liv <- function(seed, n, alpha = c(0.2, 0.2, 0.3, -0.1)) {
  set.seed(seed)
  X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
  
  Z <- runif(n, 0, 1)
  U <- runif(n, 0, 1)
  A <- as.integer(U <= Z)
  
  Y <- as.numeric(X %*% alpha) + A * (1 + U^2 + X[,1] - 0.5 * X[,2]^2 * U) +  rnorm(n)
  
  data.frame(y = Y, a = A, z = Z, x1 = X[,1], x2 = X[,2], x3 = X[,3], x4 = X[,4])
}

true_gamma <- function(z) {
  1 + z^2 - 0.5*z
}

true_denom <- function(z) {
  rep(1, length(z))
}


dctseff <- function(y, a, x, a.out, h,
                    n_pts = 100,
                    sl_lib = c("SL.glm", "SL.ranger"),
                    family = "gaussian",
                    ifs = TRUE,
                    parallel = TRUE,
                    epsilon = 0.05) {
  
  require(SuperLearner)
  require(KernSmooth)
  require(foreach)
  require(doParallel)
  require(parallel)
  
  safe_solve <- function(A, b = NULL, ridge = 1e-10) {
    tryCatch({
      if (is.null(b)) solve(A) else solve(A, b)
    }, error = function(e) {
      A2 <- A + diag(ridge, nrow(A))
      if (is.null(b)) solve(A2) else solve(A2, b)
    })
  }
  
  created_cluster <- FALSE
  if (isTRUE(parallel) && isTRUE(ifs)) {
    # if nothing registered (or only doSEQ), create cluster
    if (!foreach::getDoParRegistered() || foreach::getDoParName() == "doSEQ") {
      num_cores <- max(1, parallel::detectCores() - 2)
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      created_cluster <- TRUE
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
  }
  
  n <- length(y)
  thetahats <- matrix(NA_real_, nrow = 3, ncol = length(a.out))
  dat <- data.frame(y = y, a = a, x)
  
  if (ifs) {
    theta_ifs <- matrix(NA_real_, nrow = n, ncol = length(a.out))
  }
  
  split_index <- sample(1:3, size = n, replace = TRUE)
  
  for (split in 1:3) {
    train1 <- dat[split_index %% 3 == split %% 3, ]
    train2 <- dat[split_index %% 3 == (split + 1) %% 3, ]
    test   <- dat[split_index %% 3 == (split + 2) %% 3, ]
    
    x_train1 <- subset(train1, select = -c(a, y))
    x_train2 <- subset(train2, select = -c(a, y))
    x_test   <- subset(test,   select = -c(a, y))
    
    n1 <- nrow(train1)
    n2 <- nrow(x_train2)
    n3 <- nrow(test)
    
    # evaluation grid for nuisance integration
    a_min <- min(test$a); a_max <- max(test$a)
    a_vals <- seq(a_min, a_max, length.out = n_pts)
    
    xa_new <- cbind(x_train2[rep(1:n2, length(a_vals)), , drop = FALSE],
                    a = rep(a_vals, rep(n2, length(a_vals))))
    x_train1 <- data.frame(x_train1)
    xa_train1 <- subset(train1, select = -c(y))
    xa_test   <- subset(test,   select = -c(y))
    xa_new_df <- data.frame(xa_new)
    
    # nuisance: E[a|x], Var(a|x)
    meanA_mod <- SuperLearner(
      Y = train1$a,
      X = x_train1,
      newX = rbind(x_train1, x_train2, x_test),
      SL.library = sl_lib
    )
    meanA_pred1 <- meanA_mod$SL.predict[1:n1]
    meanA_pred2 <- meanA_mod$SL.predict[(n1 + 1):(n1 + n2)]
    meanA_pred_test <- meanA_mod$SL.predict[(n1 + n2 + 1):(n1 + n2 + n3)]
    
    varA_mod <- SuperLearner(
      Y = (train1$a - meanA_pred1)^2,
      X = x_train1,
      newX = rbind(x_train1, x_train2, x_test),
      SL.library = sl_lib
    )
    varA_pred1 <- pos(varA_mod$SL.predict[1:n1])
    varA_pred2 <- pos(varA_mod$SL.predict[(n1 + 1):(n1 + n2)])
    varA_pred_test <- pos(varA_mod$SL.predict[(n1 + n2 + 1):(n1 + n2 + n3)])
    
    std1 <- (train1$a - meanA_pred1) / sqrt(varA_pred1)
    std1 <- std1[is.finite(std1)]
    if (length(std1) < 20) std1 <- rnorm(50)
    
    kde_eps <- density(std1, bw = "SJ")
    
    # nuisance: E[y|x,a]
    if (family == "gaussian") {
      meanY_mod <- SuperLearner(
        Y = train1$y,
        X = xa_train1,
        newX = rbind(xa_test, xa_new_df),
        SL.library = sl_lib,
        family = gaussian()
      )
    } else {
      meanY_mod <- SuperLearner(
        Y = train1$y,
        X = xa_train1,
        newX = rbind(xa_test, xa_new_df),
        SL.library = sl_lib,
        family = binomial()
      )
    }
    
    pred <- meanY_mod$SL.predict
    meanY_pred_test <- pred[1:nrow(xa_test)]
    meanY_pred      <- pred[(nrow(xa_test) + 1):(nrow(xa_test) + nrow(xa_new_df))]
    
    # conditional density pi(a|x) evaluated for xa_new and test
    a_std <- (xa_new_df$a - rep(meanA_pred2, n_pts)) / sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- approx(kde_eps$x, kde_eps$y, xout = a_std, rule = 2)$y / sqrt(rep(varA_pred2, n_pts))
    pi_pred_vals <- pos(pi_pred_vals, eps = epsilon)
    
    pi_pred_test <- approx(kde_eps$x, kde_eps$y,
                           xout = (test$a - meanA_pred_test) / sqrt(varA_pred_test),
                           rule = 2)$y / sqrt(varA_pred_test)
    pi_pred_test <- pos(pi_pred_test, eps = epsilon)
    
    pi_pred_mat <- matrix(pi_pred_vals, nrow = n2, ncol = n_pts)
    f_pred <- predict(smooth.spline(a_vals, apply(pi_pred_mat, 2, mean, na.rm = TRUE)),
                      x = test$a)$y
    f_pred <- pos(f_pred, eps = epsilon)
    
    mu_pred_mat <- matrix(meanY_pred, nrow = n2, ncol = n_pts)
    tau0_pred <- predict(smooth.spline(a_vals, apply(mu_pred_mat, 2, mean, na.rm = TRUE)),
                         x = test$a)$y
    
    # pseudo outcome
    pseudo_out <- (test$y - meanY_pred_test) / (pi_pred_test / f_pred) + tau0_pred
    
    # point estimate via local polynomial derivative
    model_fit <- locpoly(x = test$a, y = pseudo_out, drv = 1, bandwidth = h, degree = 2)
    thetahats[split, ] <- approx(model_fit, xout = a.out, rule = 2)$y
    
    if (ifs) {
      idx_test <- which(split_index %% 3 == (split + 2) %% 3)
      
      a_points <- train2[sample(1:n2, min(200, n2), replace = FALSE), ]$a
      
      xa_new_var <- cbind(x_test[rep(1:n3, rep(length(a_points), n3)), , drop = FALSE],
                          a = rep(a_points, n3))
      
      mu_evals <- predict(meanY_mod, newdata = data.frame(xa_new_var))$pred
      mu_evals_mat <- matrix(mu_evals, nrow = length(a_points))  # (len(a_points) x n3)
      
      # helper to compute IFs at one a0 = a.out[j]
      compute_if_at_a0 <- function(a0) {
        g <- cbind(1, (test$a - a0) / h, ((test$a - a0) / h)^2)
        K <- gaussianKernel(test$a, a0, h)
        
        D_hat <- t(g) %*% (g * K) / n3
        
        # guard: if local design is degenerate, return NAs (so na.rm can drop)
        rc <- suppressWarnings(tryCatch(rcond(D_hat), error = function(e) NA_real_))
        if (!is.finite(rc) || rc < 1e-12) return(rep(NA_real_, n3))
        
        cross_moment <- t(g) %*% (K * pseudo_out) / n3
        beta_hat <- safe_solve(D_hat, cross_moment)
        resid    <- pseudo_out - c(g %*% beta_hat)
        
        Ginv  <- safe_solve(D_hat, t(g))         # 3 x n3
        term1 <- Ginv[2, ] * (K * resid)         # length n3
        
        K2 <- gaussianKernel(a_points, a0, h)
        g_eval <- cbind(1, (a_points - a0) / h, ((a_points - a0) / h)^2)
        
        Dinv  <- safe_solve(D_hat)
        term2 <- (Dinv %*% (t(g_eval) %*% (mu_evals_mat * K2) / length(a_points)))[2, ]
        
        term3 <- c(beta_hat)[2]
        
        (term1 + term2 - term3) / h
      }
      
      if (isTRUE(parallel)) {
        ifs_test <- foreach(j = seq_along(a.out), .combine = "cbind",
                            .export = c("gaussianKernel")) %dopar% {
                              # safe_solve + compute_if_at_a0 are in the parent environment of the foreach call,
                              # and will be captured; gaussianKernel is explicitly exported.
                              compute_if_at_a0(a.out[j])
                            }
        theta_ifs[idx_test, ] <- ifs_test
      } else {
        for (j in seq_along(a.out)) {
          theta_ifs[idx_test, j] <- compute_if_at_a0(a.out[j])
        }
      }
      
      rm(meanA_mod, varA_mod, meanY_mod, mu_evals, mu_evals_mat, xa_new_var)
      gc()
    }
  }
  
  if (ifs) {
    list(est = colMeans(thetahats), ifs = theta_ifs)
  } else {
    colMeans(thetahats)
  }
}

smooth_est <- function(y, a, x, a.out, h,
                       nsplits = 2,
                       sl.lib.ps = c("SL.glm", "SL.ranger"),
                       sl.lib.outcome = c("SL.glm", "SL.ranger"),
                       family = "gaussian",
                       ifs = TRUE,
                       epsilon = 0.05) {
  
  n <- length(y)
  splitIndex <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
  
  thetah.ifvalues <- matrix(0, nrow = n, ncol = length(a.out))
  
  for (split in 1:nsplits) {
    train <- splitIndex != split
    test <- splitIndex == split
    if (nsplits == 1) train <- test
    
    y.train <- y[train]
    y.test  <- y[test]
    a.train <- a[train]
    a.test  <- a[test]
    
    a.seq <- seq(min(a.out) - 3*h, max(a.out) + 3*h, length=100)
    
    X.train <- as.data.frame(x[train, , drop = FALSE])
    X.test  <- as.data.frame(x[test,  , drop = FALSE])
    
    # pi(a|x): mean + variance + KDE residual density
    psMod.sl <- SuperLearner(Y = a.train, X = X.train, newX = X.test, SL.library = sl.lib.ps)
    meanA.est <- psMod.sl$SL.predict
    meanA.est.train <- as.numeric(predict(psMod.sl, newdata = X.train)$pred)
    piRes2 <- (a.train - meanA.est.train)^2
    
    pi2mod <- SuperLearner(Y = piRes2, X = X.train, newX = X.test, SL.library = sl.lib.ps)
    varA.est <- pos(pi2mod$SL.predict)
    varA.est.train <- pos(as.numeric(predict(pi2mod, newdata = X.train)$pred))
    
    std_train <- (a.train - meanA.est.train) / sqrt(varA.est.train)
    std_train <- std_train[is.finite(std_train)]
    if (length(std_train) < 20) std_train <- rnorm(50)
    kde_eps <- density(std_train, bw = "SJ")
    
    a_std <- (a.test - meanA.est) / sqrt(varA.est)
    piHat.XA_test <- approx(kde_eps$x, kde_eps$y, xout = a_std, rule = 2)$y / sqrt(varA.est)
    piHat.XA_test <- pos(piHat.XA_test, eps = epsilon)
    
    # mu(x,a): outcome regression
    XA.train <- data.frame(x = x[train, , drop = FALSE], a = a[train])
    XA.test  <- data.frame(x = x[test,  , drop = FALSE], a = a[test])
    
    if (family == "gaussian") {
      muModel <- SuperLearner(Y = y.train, X = XA.train, newX = XA.test,
                              SL.library = sl.lib.outcome, family = gaussian())
    } else {
      muModel <- SuperLearner(Y = y.train, X = XA.train, newX = XA.test,
                              SL.library = sl.lib.outcome, family = binomial())
    }
    
    muHat.XA_test <- as.numeric(predict(muModel)$pred)
    ipwRes <- (y.test - muHat.XA_test) / piHat.XA_test
    
    # predictions mu(x,a0) over a.seq
    mu.preddata <- X.test[rep(1:length(y.test), length(a.seq)), , drop = FALSE]
    mu.preddata$a <- rep(a.seq, rep(length(y.test), length(a.seq)))
    colnames(mu.preddata) <- colnames(XA.train)
    
    muHat.a0 <- matrix(as.numeric(predict(muModel, newdata = mu.preddata)$pred),
                       ncol = length(a.seq))
    
    a.seq.len <- a.seq[2] - a.seq[1]
    
    kernel.A  <- sapply(a.out, gaussianKernel_deriv, A = a.test, h = h)
    kernel.a0 <- sapply(a.out, gaussianKernel_deriv, A = a.seq,  h = h)
    
    psih.1 <- kernel.A * c(ipwRes)
    psih.2 <- a.seq.len * muHat.a0 %*% kernel.a0
    
    thetah.ifvalues[test, ] <- -psih.1 - psih.2
  }
  
  if (!ifs) {
    colMeans(thetah.ifvalues)
  } else {
    list(est = colMeans(thetah.ifvalues), ifs = thetah.ifvalues)
  }
}

npliv_lp <- function(y, a, x, z, z.out, h1, h2, family.y, family.a = 'binomial', parallel = FALSE) {
  # y: outcome; a: treatment (binary); z: instrument (continuous)
  # h1: bandwidth for numerator derivative; h2: bandwidth for denominator derivative
  numer <- dctseff(y, z, x, a.out = z.out, h = h1, family = family.y, ifs = TRUE, parallel = parallel)
  denom <- dctseff(a, z, x, a.out = z.out, h = h2, family = family.a, ifs = TRUE, parallel = parallel)
  
  ifs1 <- t(numer$ifs) # D x n
  ifs2 <- t(denom$ifs) # D x n
  est1 <- numer$est
  est2 <- denom$est
  
  ifs_ratio <- ifs1 / est2 - ifs2 * est1 / (est2^2)
  sd_est <- apply(ifs_ratio, 1, sd, na.rm = TRUE) / sqrt(length(y))
  
  list(est = est1 / est2, sd_error = sd_est)
}

npliv_smooth <- function(y, a, x, z, z.out, h1, h2, family.y, family.a = 'binomial') {
  numer <- smooth_est(y, z, x, a.out = z.out, h = h1, family = family.y, ifs = TRUE)
  denom <- smooth_est(a, z, x, a.out = z.out, h = h2, family = family.a, ifs = TRUE)
  
  ifs1 <- t(numer$ifs)
  ifs2 <- t(denom$ifs)
  est1 <- numer$est
  est2 <- denom$est
  
  ifs_ratio <- ifs1 / est2 - ifs2 * est1 / (est2^2)
  sd_est <- apply(ifs_ratio, 1, sd, na.rm = TRUE) / sqrt(length(y))
  
  list(est = est1 / est2, sd_error = sd_est)
}
