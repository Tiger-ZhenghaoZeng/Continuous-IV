#'  Goal: evaluate performance of DR and plug-in estimators 
#'  when the nuisance functions are estimated with different convergence rate O(n^{\alpha})
#'  Key: guarantee that RMSE of muhat and pihat are O_P(n^\alpha)
#'  to remain agnostic the types of estimators
#'  as alpha gets closer to the parametric rate, plug-in should be doing better

rm(list=ls())

source('simu1_utils.R')
library(foreach)
library(doParallel)

nsim = 100; n = 20000
alphas = seq(from = 0.02, to = 0.3, by=0.01)  # Probably don't need this. We just use the nonparametric estimates.

beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0.1, -0.1, 0.1, 0.13); beta_pi = c(0.1, 0.1, -0.1, 0.2)


thre <- 0.1
a.out <- simulate_data_normal(521, 500, beta1, beta2, beta_pi)$a
lb = quantile(a.out, thre)
ub = quantile(a.out, 1-thre)
idx = which(a.out >= lb & a.out <= ub)
a.out <- a.out[idx]
a.out <- a.out[order(a.out)] # There are 400 points right now, maybe too many, can first go with it
# true = 1 +  a*(beta2[1]  - beta2[4]^2*a^2)
true.out = beta2[1] - 3*beta2[4]^2*a.out^2


# plugin_estimates = dr_estimates = smooth_estimates = matrix(0, nrow = nsim, ncol = length(a.out))
bws <- seq(from=0.02, to=3, length.out=100) # potential sets of bandwidth, should fix one before running repeated experiments

# Un-parallel version
# rmses = matrix(0, nrow = 3, ncol = length(alphas))
# for(i in seq_along(alphas)){
#   plugin_estimates = dr_estimates = smooth_estimates = matrix(0, nrow = nsim, ncol = length(a.out))
#   for (sim in 1:nsim){
#     dat = simulate_data_normal(sim*2024, n, beta1, beta2, beta_pi)
#     y = dat$y; a = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4)
#     order_a = order(a)
#     
#     a = a[order_a]
#     y = y[order_a]
#     x = x[order_a, ]
#     
#     # res = ctseff(y, a, x, bw.seq = seq(0.2, 2, 100))
#     
#     
#     # muhat = f_hat(beta1, beta2, 0)
#     # mu2hat = f2_hat(beta1, beta2, 1)
#     # pihat = prob_hat(beta_pi, -0.8, 0)
#     # pi2hat = prob2_hat(beta_pi, -0.8, 0)
#     
#     muhat = meanY_hat(beta1, beta2, rnorm(1, n^(-alphas[alpha]), n^(-alphas[alpha])))
#     # mu2hat = f2_hat(beta1, beta2, rnorm(1, 0, n^(-alphas[alpha])))
#     lambdahat = meanA_hat(beta_pi, -0.8, rnorm(1, n^(-alphas[alpha]), n^(-alphas[alpha])))
#     # pi2hat = prob2_hat(beta_pi, -0.8, rnorm(1, 0, n^(-alphas[alpha])))
#     
#     plugin_estimates[sim,] = dctseff_plugin_simu(a.out, x, muhat)
#     dr_estimates[sim,] = dctseff_simu(y, a, x, muhat, lambdahat, a.out=a.out)
#     smooth_estimates[sim,] <- smooth_est_simu(dat, hs, a.eval, muhat, lambdahat, gaussianKernel_deriv)
#     
#     # if (sim %% 20 == 0){
#     #   plot(a.out, true.out, col='grey', pch=20, cex=0.3, ylab=expression(theta(a)))
#     #   points(a.out, dr_estimates[sim,], col='black', pch=20, cex=0.3)
#     #   points(a.out, plugin_estimates[sim,], col='red', pch=20, cex=0.3)
#     #   legend("bottomleft", c("True", "Doubly Robust", "Plug-in"), lty=c(1,1,1), col=c("grey","black", "red"), cex=0.6)
#     # }
#   }
#   # c(rmse(true.out, plugin_estimates), rmse(true.out, dr_estimates), rmse(true.out, smooth_estimates))
#   rmses[1,alpha] = rmse(true.out, plugin_estimates)
#   rmses[2,alpha] = rmse(true.out, dr_estimates)
#   rmses[3,alpha] = rmse(true.out, smooth_estimates)
# }

# Parallel version
set.seed(521)
myCluster <- makeCluster(detectCores()-1, type="FORK")
registerDoParallel(myCluster)
rmses <- foreach(alpha = alphas, .combine = 'cbind') %dopar%{
  plugin_estimates = dr_estimates = smooth_estimates = matrix(0, nrow = nsim, ncol = length(a.out))
  for (sim in 1:nsim){
    dat = simulate_data_normal(sim*2024, n, beta1, beta2, beta_pi)
    y = dat$y; a = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4)
    order_a = order(a)
    
    a = a[order_a]
    y = y[order_a]
    x = x[order_a, ]
    
    # res = ctseff(y, a, x, bw.seq = seq(0.2, 2, 100))
    

    # muhat = f_hat(beta1, beta2, 0)
    # mu2hat = f2_hat(beta1, beta2, 1)
    # pihat = prob_hat(beta_pi, -0.8, 0)
    # pi2hat = prob2_hat(beta_pi, -0.8, 0)

    muhat = meanY_hat(beta1, beta2, 5*rnorm(1, n^(-alpha), n^(-alpha)))  # Replace with nonparametric estimates
    # mu2hat = f2_hat(beta1, beta2, rnorm(1, 0, n^(-alphas[alpha])))
    lambdahat = meanA_hat(beta_pi, -0.8, 5*rnorm(1, n^(-alpha), n^(-alpha))) # Replace with nonparametric estimates
    # pi2hat = prob2_hat(beta_pi, -0.8, rnorm(1, 0, n^(-alphas[alpha])))

    plugin_estimates[sim,] = dctseff_plugin_simu(a.out, x, muhat) # Don't need to compute it
    dr_estimates[sim,] = dctseff_simu(y, a, x, muhat, lambdahat, a.out=a.out) # Replace with functions in realdata_utils.R
    smooth_estimates[sim,] <- smooth_est_simu(dat, bws, a.out, muhat, lambdahat, gaussianKernel_deriv)

  }
  c(rmse(true.out, plugin_estimates), rmse(true.out, dr_estimates), rmse(true.out, smooth_estimates)) # Should return derivative estimates, variance estimates, coverage
  # rmses[1,alpha] = rmse(true.out, plugin_estimates)
  # rmses[2,alpha] = rmse(true.out, dr_estimates)
  # rmses[3,alpha] = rmse(true.out, smooth_estimates)
}
stopCluster(myCluster)
write.csv(as.data.frame(rmses), file = "rmse20000.csv", row.names = FALSE)


data <- data.frame(
  x = rep(alphas, 3),  
  y = c(t(as.matrix(rmse20000))),  
  line = rep(c("Plug-in", "LP", "Smooth"), each = length(alphas))  
)

# Step 2: Plot the three lines using ggplot2
ggplot(data, aes(x = x, y = y, color = line, group = line)) +
  geom_line(linewidth = 1) +  # Add lines
  geom_point(size = 1.5) +   # Add points to mark each line
  labs(x = expression(alpha),
       y = "RMSE",
       color = "Lines") +  # Label for the legend
  theme_minimal() +  # Apply a clean theme
  theme(
    legend.text = element_text(size = 16),    # Increase legend item text size
    legend.title = element_text(size = 16),  # Increase legend title text size
    legend.key.size = unit(1.5, "cm"),        # Increase legend key size
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size title
    axis.text = element_text(size = 16),  # Increase axis text size
    axis.title = element_text(size = 16)  # Increase axis label size
  )

# ns <- c(10000, 50000, 100000, 200000)
# hs <- seq(from=0.02, to=1.2, by=0.02)
# plot(a.out, true.out, col=1, ylim=c(-0.2,0.15), type='l')
# for (i in seq_along(ns)){
#   dat = simulate_data_normal(2024, ns[i], beta1, beta2, beta_pi)
#   res <- smooth_est_simu(dat, hs, a.out, muhat, lambdahat, gaussianKernel_deriv)
#   lines(a.out, res, col=i+1)
# }
# legend("topright", paste0('n=', ns), lty=c(1,1,1,1), pch=c(1,1,1,1), col=2:5, cex=0.6)

# plot(x = alphas, rmses[2,], main = paste0('n=',n), ylim = range(c(rmses[1,],rmses[2,], rmses[3,])), xlab=expression(alpha),
#      ylab = "RMSE", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
# lines(x = alphas, y = rmses[2,])
# points(x = alphas, y = rmses[1,], col = "red")
# lines(x = alphas, y = rmses[1,], col = "red", lty=2)
# points(x = alphas, y = rmses[3,], col = "blue")
# lines(x = alphas, y = rmses[3,], col="blue", lty=3)
# legend("topright", c("Doubly Robust", "Plug-in", "Smooth"), lty=1:3, pch=c(1,1,1), col=c("black", "red", "blue"), cex=0.6)

