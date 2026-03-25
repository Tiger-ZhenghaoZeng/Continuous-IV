
#' TODO()
#' everything blows up sometimes...
#' increase proportion of compliers 
#' try coding up the linear projection
#' try fitting two separate bandwidths? 

rm(list=ls())

library(foreach)
library(doParallel)
library(ggplot2)

source('simu2_utils.R')

alphas = seq(from = 0.01, to = 0.3, by=0.01)

nsim = 100; n = 2000
beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0, -0.1, 0.1, 0.13)
beta_pi1 = c(0.1, 0.1, -0.1, 0.2); beta_pi2 = c(0.1, -0.2, 0.3, 0.1)

thre <- 0.1
z.out <- simulate_data(521, beta1, beta2, beta_pi1, beta_pi2, 500)$z
lb = quantile(z.out, thre)
ub = quantile(z.out, 1-thre)
idx = which(z.out >= lb & z.out <= ub)
z.out <- z.out[idx]
z.out <- z.out[order(z.out)]
true_num.out <- -3*0.13^2*z.out^2
true_den.out <- rep(0.1,length(z.out))
true.out = -0.507*z.out^2

myCluster = makeCluster(detectCores()-1, type="PSOCK")
registerDoParallel(myCluster)

clusterEvalQ(myCluster, {
  library(foreach)
  library(doParallel)
  library(MASS) 
  library(KernSmooth) 
  #source('sims2_utils_updates.R')
})

t0 = proc.time()[3]
rmses = foreach(alpha = alphas, .combine = 'cbind') %dopar%{
  spliv_estimates = npliv_lp_estimates = npliv_smooth_estimates = matrix(0, nrow=nsim, ncol=n*(1-2*thre))
  
  for (i in 1:nsim){
    dat = simulate_data(i*666, beta1, beta2, beta_pi1, beta_pi2, n)
    y = dat$y; a = dat$a; x = cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4); x = as.matrix(x); z = dat$z
    order_z = order(z)
    a = a[order_z]; y = y[order_z]; x = x[order_z, ]; z = z[order_z]
    
    # z.out = z
    # lb = quantile(z.out, thre); ub = quantile(z.out, 1-thre)
    # idx = which(z.out >= lb & z.out <= ub)
    # z.out = z.out[idx]
    # y = y[idx]; a = a[idx]; x = x[idx,]; z = z[idx]
    # z.out = z.out[order(z.out)]
    
    # true.out = -0.507*z.out^2
    
    # muhat = meanY_hat(beta1, beta2, 0)
    # lambdahat = meanZ_hat(beta_pi1, 0)
    # etahat = meanA_hat(beta_pi2, 0)
    
    muhat = meanY_hat(beta1, beta2, rnorm(1, n^(-alpha), n^(-alpha)))
    lambdahat = meanZ_hat(beta_pi1, rnorm(1, n^(-alpha), n^(-alpha)))
    etahat = meanA_hat(beta_pi2, rnorm(1, n^(-alpha), n^(-alpha)))
    
    # # Check numerator
    # num_lp_est <- dctseff_simu(y, z, x, muhat, lambdahat, a.out=z.out)
    # # num_smooth_est <- smooth_est_simu(dat, seq(from=0.02, to=3, length.out=100), z.out, muhat, lambdahat, gaussianKernel_deriv)
    # plot(z.out, true_num.out, type='l')
    # lines(z.out, num_lp_est,col='red')
    # # lines(z.out, num_smooth_est,col='blue')
    # 
    # # Check denominator
    # denom_lp_est <- dctseff_simu(a, z, x, muhat, lambdahat, a.out=z.out)
    # plot(z.out, rep(0.1, length(z.out)), type='l')
    # lines(z.out, denom_lp_est,col='red')
    
    psi_hat = spliv_dr_simu(y=y, a=a, z=z, x=x, muhat=muhat, etahat=etahat, lambdahat=lambdahat)
    npliv_lp_estimates[i,] = npliv_dr_simu(y=y, a=a, z=z, x=x, z.out=z.out, muhat=muhat, lambdahat=lambdahat, etahat=etahat, true.out=true.out)
    npliv_smooth_estimates[i,] = npliv_smooth_simu(y=y, a=a, z=z, x=x, z.out=z.out, muhat=muhat, lambdahat=lambdahat, etahat=etahat, kernelfun=gaussianKernel_deriv, true.out=true.out)
    spliv_estimates[i,] <- psi_hat*z.out
    # plot(z.out, true.out, type='l')
    # lines(z.out, psi_hat*z.out, col=2)
    # lines(z.out, lp_est, col=3)
    # lines(z.out, smooth_est, col=4)
    # legend("topright", c("true", "linear", "lp", "smooth"), col=1:4, lty=rep(1,4))
  }
    c(rmse(true.out, spliv_estimates), rmse(true.out, npliv_lp_estimates), rmse(true.out, npliv_smooth_estimates))
    
}
stopCluster(myCluster)
t1 = proc.time()[3]
print(paste('time elapsed:', round(t1-t0,2), 'seconds'))

# save(rmses, file = paste0("rmse", n, ".RData"))
write.csv(as.data.frame(rmses), file = paste0("simu2rmse", n, ".csv"), row.names = FALSE)
xx <- as.data.frame(lapply(simu2rmse2000, function(x) if(is.numeric(x)) round(x, 6) else x))
xx <- as.data.frame(lapply(simu2rmse20000, function(x) if(is.numeric(x)) round(x, 6) else x))
# load
# load('rmse2000.RData')
save

data_plot = data.frame(
  x = rep(alphas, 3),  
  y = c(t(as.matrix(simu2rmse2000))),  
  line = rep(c("Linear", "LP", "Smooth"), each = length(alphas))  
)

# Step 2: Plot the three lines using ggplot2
ggplot(data_plot, aes(x = x, y = y, color = line, group = line)) +
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


