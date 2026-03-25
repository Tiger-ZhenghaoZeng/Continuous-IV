rm(list=ls())

library(numDeriv)
library(SuperLearner)
library(KernSmooth)
library(MASS)

source('utils_updates.R')

# cv = function(y, a, x, hs,  sl_lib= c("SL.earth", "SL.gam", 
#                                                 "SL.glm", "SL.glm.interaction",
#                                                 "SL.mean", "SL.ranger", "SL.rpart"), n_pts=100)
#   {
#   
#   dat = data.frame(y,a,x); n = nrow(dat)
#   
#   #sl_lib = c("SL.earth", "SL.gam", 
#   #           "SL.glm", "SL.glm.interaction",
#   #           "SL.mean", "SL.ranger", "SL.rpart")
#   
#   risks = numeric(length(hs))
#   
#   for (i in 1:length(hs)){
#     print(paste('h:', hs[i]))
#     n_split = 2
#     split_index = sample(rep(1:n_split, ceiling(n/n_split))[1:n])
#     
#     risk = numeric(n_split)
#     for (split in 1:n_split){
#       train = dat[split_index != split, ]; test = dat[split_index == split, ]
#       n1 = nrow(train); n2 = nrow(test)
#       
#       ## set up data for training ##
#       x_train = subset(train, select=-c(a,y)); x_test = subset(test, select=-c(a,y))
#       a_seq = seq(min(test$a), max(test$a), length.out = n_pts)
#       ax_train = subset(train, select=-c(y)); ax_test = subset(test, select=-c(y))
#       ax_new = cbind(a=rep(a_seq, rep(n2,length(a_seq))), x_test[rep(1:n2,length(a_seq)),])
#       
#       
#       ## mu function used for residual term and integral 
#       meanY_mod = SuperLearner(Y = train$y,
#                                X = ax_train, newX = rbind(ax_new, ax_test),
#                                SL.library = sl_lib)
#       meanY_pred = meanY_mod$SL.predict[1:nrow(ax_new)]
#       meanY_pred_mat <- matrix(meanY_pred, nrow=n2) #n_2 * 100 matrix for evaluating the integral
#       meanY_pred_test <- meanY_mod$SL.predict[(1+nrow(ax_new)):(nrow(test)+nrow(ax_new))]
#       
#       # pi function used for the residual term
#       meanA_mod = SuperLearner(Y = train$a,
#                                X = x_train, newX = rbind(x_train, x_test),
#                                SL.library = sl_lib)
#       meanA_pred <- meanA_mod$SL.predict[1:n1]
#       meanA_pred_test <- meanA_mod$SL.predict[(n1+1):(n1+n2)]
#       
#       varA_mod = SuperLearner(Y = (train$a - meanA_pred)^2, 
#                               X = x_train, newX = rbind(x_train, x_test), 
#                               SL.library = sl_lib)
#       varA_pred <- varA_mod$SL.predict[1:n1]
#       varA_pred_test <- varA_mod$SL.predict[(n1+1):(n1+n2)]
#       
#       kde_eps <- density((train$a-meanA_pred)/sqrt(varA_pred), bw="SJ")
#       a_std = (test$a - meanA_pred_test) / sqrt(varA_pred_test)
#       pi_pred_test <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(varA_pred_test)
#       pi_pred_test <- ifelse(pi_pred_test<=0.01|is.na(pi_pred_test), 0.01, pi_pred_test)
#       
#       ## product of theta and f
#       theta_hat <- dctseff(y=train$y, train$a, x_train, a.out=c(test$a, a_seq), hs[i], n_pts = 100)
#       theta_hat_test <- theta_hat[1:n2]
#       theta_hat_seq <- theta_hat[(n2+1):(n2+n_pts)]
#       
#       kde_a <- density(train$a, bw="SJ")
#       density_a_seq <- approx(kde_a$x, kde_a$y, xout=a_seq, rule = 2)$y
#       
#       ## Fit the derivative
#       spline_fit <- smooth.spline(a_seq, density_a_seq*theta_hat_seq)
#       spline_deri <- predict(spline_fit, c(test$a, a_seq), deriv = 1)$y
#       deriv_test <- spline_deri[1:n2]
#       deriv_seq <- spline_deri[(n2+1):(n2+n_pts)]
#  
#       #temp <- meanY_pred_mat*deriv_seq
#       integral <- meanY_pred_mat%*%deriv_seq*(a_seq[2]-a_seq[1])
#       risk[split] = mean(theta_hat_test^2 + 2*(integral +  deriv_test* (test$y-meanY_pred_test)/pi_pred_test))
#       
#       # x_train = data.frame(x_train)
#       
#       # xa_new = data.frame(xa_new)
#       
#       # colnames(x_train) = colnames(x_new)
#       # colnames(xa_train) = colnames(xa_new)
#         
#       ## estimate nuisance functions via super learner ##
#       # make sure these are correct
#       # pi_mod = SuperLearner(Y = train$a, X = x_train, newX = x_new, SL.library = sl_lib)
#       # pimod_vals = pi_mod$SL.predict
#       # pimod_vals <- ifelse(pimod_vals<=0.01|is.na(pimod_vals), 0.01, pimod_vals)
#       # pi_pred_mat = matrix(pimod_vals, nrow = n2, ncol = length(a_vals))
#       # f_pred = predict(smooth.spline(a_vals, apply(pi_pred_mat, 2, mean)), x = test$a)$y
#       
#       # mu_mod = SuperLearner(Y = train$y, X = xa_train, newX = xa_new, SL.library = sl_lib)
#       # muhat_vals = mu_mod$SL.predict
#       # muhat = muhat_vals[1:n2]
#       
#       # pi2_mod = SuperLearner(Y = (train$a - pimod_vals[1:n2])^2,X = x_train, newX = x_new, SL.library = sl_lib)
#       # pi2mod_vals = pi2_mod$SL.predict
#       # a_std = (xa_new$a - pimod_vals) / sqrt(pi2mod_vals)
#       # pihat_vals = approx(density(a_std)$x, density(a_std[1:n2])$y, xout = a_std)$y / sqrt(pi2mod_vals)
#       # pihat = pihat_vals[1:n2]
#       # 
#       # mu_pred_mat = matrix(muhat_vals, nrow = n2, ncol = length(a_vals))
#       # tau0_pred = predict(smooth.spline(a_vals, apply(mu_pred_mat, 2, mean)), x = test$a)$y
#       # 
#       # marginal_fit = density(train$a, bw=hs[i])
#       # 
#       # pseudo_out = (test$y - muhat)/pihat*f_pred + tau0_pred
#       # theta_fit = locpoly(x=test$a, y=pseudo_out, drv=1, bandwidth=hs[i]) 
#       # thetahat = approx(theta_fit, xout = test$a, rule = 2)$y ## extrapolation
#       
#       ## computes the d{fhat(z)thetahat(z)}/dz ##
#       # ftheta_deriv = function(a.new, h){
#       #   fhat = approx(marginal_fit$x, marginal_fit$y, xout = a.new, rule = 2)$y ## extrapolation
#       #   fhat * approx(theta_fit, xout = a.new, rule = 2)$y ## extrapolation
#       # }
#       
#       ## integrand ##
#       # int_func = function(a.new, h){
#       #   xa_new = cbind(x_test[1,], a=a.new); xa_new = data.frame(xa_new)
#       #   muhat = predict(mu_mod, xa_new)$pred
#       #   
#       #   grad(function(a.new) ftheta_deriv(a.new, h), a.new) * muhat
#       # }
# 
#       ## compute the integral and pseudo risk ##
#       # a.seq = seq(min(test$a), max(test$a), length.out = 100)
#       # a.seq.len = diff(a.seq)[1]
#       # int_val = sum(sapply(a.seq, function(a) int_func(a, h=hs[i])) * a.seq.len)
#       # 
#       # risk[split] = mean(thetahat^2 + 2*(int_val + grad(function(a) ftheta_deriv(a, h=hs[i]), a) * (test$y-muhat)/pihat))
#     }
#     
#     risks[i] = mean(risk)
#   }
#   
#   # hhat = hs[which.min(risks)]
#   
#   return(list(hs=hs, risks=risks))
# }


cv = function(y, a, x, hs, method, sl_lib= c("SL.earth", "SL.gam", 
                                      "SL.glm", "SL.glm.interaction",
                                      "SL.mean", "SL.ranger", "SL.rpart"), n_pts=100)
{
  # method can be "lp" or "smooth"
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
    meanY_mod = SuperLearner(Y = train$y,
                             X = ax_train, newX = rbind(ax_new, ax_test),
                             SL.library = sl_lib)
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
    pi_pred_test <- ifelse(pi_pred_test<=0.01|is.na(pi_pred_test), 0.01, pi_pred_test)
    
    # marginal density f
    kde_a <- density(train$a, bw="SJ")
    density_a_seq <- approx(kde_a$x, kde_a$y, xout=a_seq, rule = 2)$y
    
    
    for (i in 1:length(hs)){
      ## Fit theta
      if (method=='lp'){
        theta_hat <- dctseff(y=train$y, train$a, x_train, a.out=c(test$a, a_seq), hs[i], n_pts = 100)
        theta_hat_test <- theta_hat[1:n2]
        theta_hat_seq <- theta_hat[(n2+1):(n2+n_pts)]
      }
      else if(method=='smooth'){
        theta_hat <- estPsih.split(y=train$y, train$a, x_train, hs[i], a.out=c(test$a, a_seq))
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
#   for (i in 1:length(hs)){
#     print(paste('h:', hs[i]))
#     n_split = 2
#     split_index = sample(rep(1:n_split, ceiling(n/n_split))[1:n])
#     
#     risk = numeric(n_split)
#     for (split in 1:n_split){
#       train = dat[split_index != split, ]; test = dat[split_index == split, ]
#       n1 = nrow(train); n2 = nrow(test)
#       
#       ## set up data for training ##
#       x_train = subset(train, select=-c(a,y)); x_test = subset(test, select=-c(a,y))
#       a_seq = seq(min(test$a), max(test$a), length.out = n_pts)
#       ax_train = subset(train, select=-c(y)); ax_test = subset(test, select=-c(y))
#       ax_new = cbind(a=rep(a_seq, rep(n2,length(a_seq))), x_test[rep(1:n2,length(a_seq)),])
#       
#       
#       ## mu function used for residual term and integral 
#       meanY_mod = SuperLearner(Y = train$y,
#                                X = ax_train, newX = rbind(ax_new, ax_test),
#                                SL.library = sl_lib)
#       meanY_pred = meanY_mod$SL.predict[1:nrow(ax_new)]
#       meanY_pred_mat <- matrix(meanY_pred, nrow=n2) #n_2 * 100 matrix for evaluating the integral
#       meanY_pred_test <- meanY_mod$SL.predict[(1+nrow(ax_new)):(nrow(test)+nrow(ax_new))]
#       
#       # pi function used for the residual term
#       meanA_mod = SuperLearner(Y = train$a,
#                                X = x_train, newX = rbind(x_train, x_test),
#                                SL.library = sl_lib)
#       meanA_pred <- meanA_mod$SL.predict[1:n1]
#       meanA_pred_test <- meanA_mod$SL.predict[(n1+1):(n1+n2)]
#       
#       varA_mod = SuperLearner(Y = (train$a - meanA_pred)^2, 
#                               X = x_train, newX = rbind(x_train, x_test), 
#                               SL.library = sl_lib)
#       varA_pred <- varA_mod$SL.predict[1:n1]
#       varA_pred_test <- varA_mod$SL.predict[(n1+1):(n1+n2)]
#       
#       kde_eps <- density((train$a-meanA_pred)/sqrt(varA_pred), bw="SJ")
#       a_std = (test$a - meanA_pred_test) / sqrt(varA_pred_test)
#       pi_pred_test <- approx(kde_eps$x, kde_eps$y, xout=a_std, rule = 2)$y/ sqrt(varA_pred_test)
#       pi_pred_test <- ifelse(pi_pred_test<=0.01|is.na(pi_pred_test), 0.01, pi_pred_test)
#       
#       ## product of theta and f
#       theta_hat <- dctseff(y=train$y, train$a, x_train, a.out=c(test$a, a_seq), hs[i], n_pts = 100)
#       theta_hat_test <- theta_hat[1:n2]
#       theta_hat_seq <- theta_hat[(n2+1):(n2+n_pts)]
#       
#       kde_a <- density(train$a, bw="SJ")
#       density_a_seq <- approx(kde_a$x, kde_a$y, xout=a_seq, rule = 2)$y
#       
#       ## Fit the derivative
#       spline_fit <- smooth.spline(a_seq, density_a_seq*theta_hat_seq)
#       spline_deri <- predict(spline_fit, c(test$a, a_seq), deriv = 1)$y
#       deriv_test <- spline_deri[1:n2]
#       deriv_seq <- spline_deri[(n2+1):(n2+n_pts)]
#       
#       #temp <- meanY_pred_mat*deriv_seq
#       integral <- meanY_pred_mat%*%deriv_seq*(a_seq[2]-a_seq[1])
#       risk[split] = mean(theta_hat_test^2 + 2*(integral +  deriv_test* (test$y-meanY_pred_test)/pi_pred_test))
#       
#     }
#     
#     risks[i] = mean(risk)
#   }
#   
#   # hhat = hs[which.min(risks)]
#   
#   return(list(hs=hs, risks=risks))
# }


n = 500
beta1 = c(0.2, 0.2, 0.3, -0.1); beta2 = c(0.1, -0.1, 0.1, 0.13); beta_pi = c(0.1, 0.1, -0.1, 0.2)

data = simulate_data_normal(2024, n, beta1, beta2, beta_pi)
y = data$y; a = data$a; x = cbind(data$x.1, data$x.2, data$x.3, data$x.4)
a.out = data$a
thre <- 0.1
lb = quantile(a.out, thre)
ub = quantile(a.out, 1-thre)
idx = which(a.out >= lb & a.out <= ub)
a.out <- a.out[idx]
a.out <- a.out[order(a.out)]
true.out = beta2[1] - 3*beta2[4]^2*a.out^2

n <- 10000
data = simulate_data_normal(2024, n, beta1, beta2, beta_pi)
y = data$y; a = data$a; x = cbind(data$x.1, data$x.2, data$x.3, data$x.4)
order_a = order(a)
a = a[order_a]
y = y[order_a]
x = x[order_a, ]
hs = seq(0.5, 2.5, by=0.05)

cv_res = cv(y, a, x, hs, sl_lib=c("SL.gam", 
                           "SL.glm", "SL.glm.interaction",
                           "SL.mean", "SL.ranger"))
hhat = cv_res$hs[which.min(cv_res$risks)]
dr_res = dctseff(y=y, a=a, x=x, a.out=a.out, h=hhat)

results_df = data.frame(
  a.out = c(a.out, a.out),
  value = c(true.out, dr_res),
  category = factor(rep(c("True", "DR"), each = length(a.out)))
)

library(ggplot2)
library(ggthemes)

ggplot(results_df, aes(x = a.out, y = value, color = category)) +
  geom_point(aes(shape = category), size = 1) +
  scale_color_manual(values = c("blue", "grey")) +
  labs(
       x = "a", y = expression(theta(a)), title = "n=20000, h=4") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) #+  ylim(min(true.out), max(true.out))

