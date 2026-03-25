rm(list = ls())

library(future)
library(future.apply)
library(ggplot2)
library(grid)

source("simu3_utils.R")

seed <- 521

# true coefficients
beta_pi <- c(0.1, 0.1, -0.1, 0.2)
beta1   <- c(0.2, 0.2, 0.3, -0.1)
beta2   <- c(0.1, -0.1, 0.1, 0.13)

# grids / tuning
ns   <- seq(500, 5000, by = 500)
hs   <- seq(0.02, 2, by = 0.02)
thre <- 0.1
D    <- 80
M    <- 1000

# choose a fixed support using a big pilot sample
n0 <- 200000
dat0 <- simulate_data_normal(seed, n0, beta1, beta2, beta_pi)
lb0 <- quantile(dat0$a, thre)
ub0 <- quantile(dat0$a, 1 - thre)

# ---- fixed grid ----
a_grid0    <- seq(lb0, ub0, length.out = D)
true_grid0 <- beta2[1] - 3 * beta2[4]^2 * a_grid0^2

zcrit <- qnorm(0.975)

n_workers <- min(availableCores(), 8)
plan(multisession, workers = n_workers)


select_bandwidths <- function(n, seed_n, beta1, beta2, beta_pi,
                              hs, lb0, ub0, a_grid0, true_grid0) {
  
  dat <- simulate_data_normal(seed_n, n, beta1, beta2, beta_pi)
  y <- dat$y
  a <- dat$a
  x <- as.matrix(cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4))
  
  # use the same fixed support across all n
  lb <- lb0
  ub <- ub0
  keep <- which(a >= lb & a <= ub)
  
  y <- y[keep]
  a <- a[keep]
  x <- x[keep, , drop = FALSE]
  
  a.out    <- a_grid0
  true.out <- true_grid0
  
  mses_lp     <- numeric(length(hs))
  mses_smooth <- numeric(length(hs))
  
  for (i in seq_along(hs)) {
    cat("--- n =", n, "| h =", hs[i], "---\n")
    
    # avoid nested parallelism: keep dctseff parallel = FALSE here
    lp_fit <- dctseff(
      y = y, a = a, x = x, a.out = a.out, h = hs[i],
      ifs = FALSE, parallel = FALSE
    )
    
    smooth_fit <- smooth_est(
      y = y, a = a, x = x, a.out = a.out, h = hs[i],
      ifs = FALSE
    )
    
    mses_lp[i]     <- mean((lp_fit - true.out)^2, na.rm = TRUE)
    mses_smooth[i] <- mean((smooth_fit - true.out)^2, na.rm = TRUE)
  }
  
  h_lp <- hs[which.min(mses_lp)]
  h_sm <- hs[which.min(mses_smooth)]
  
  # undersmoothing
  h_lp <- h_lp / log(n)
  h_sm <- h_sm / log(n)
  
  list(h_lp = h_lp, h_sm = h_sm, lb = lb, ub = ub)
}

run_one_rep <- function(m, n, seed_n, lb, ub, a_grid,
                        beta1, beta2, beta_pi, h_lp, h_sm) {
  
  dat <- simulate_data_normal(
    seed = seed_n + m, n = n,
    beta1 = beta1, beta2 = beta2, beta_pi = beta_pi
  )
  
  y <- dat$y
  a <- dat$a
  x <- as.matrix(cbind(dat$x.1, dat$x.2, dat$x.3, dat$x.4))
  
  # truncate to fixed support
  keep <- which(a >= lb & a <= ub)
  y <- y[keep]
  a <- a[keep]
  x <- x[keep, , drop = FALSE]
  n_eff <- length(y)
  
  # LP
  fit_lp <- dctseff(
    y = y, a = a, x = x, a.out = a_grid, h = h_lp,
    ifs = TRUE, parallel = FALSE
  )
  est_lp <- fit_lp$est
  var_lp <- apply(fit_lp$ifs, 2, var, na.rm = TRUE) / n_eff
  
  # Smooth
  fit_sm <- smooth_est(
    y = y, a = a, x = x, a.out = a_grid, h = h_sm,
    ifs = TRUE, nsplits = 2
  )
  est_sm <- fit_sm$est
  var_sm <- apply(fit_sm$ifs, 2, var, na.rm = TRUE) / n_eff
  
  list(
    est_lp = est_lp, var_lp = var_lp,
    est_sm = est_sm, var_sm = var_sm
  )
}

run_for_n <- function(n, seed, beta1, beta2, beta_pi,
                      hs, M, zcrit, lb0, ub0, a_grid0, true_grid0) {
  
  seed_n <- seed + 100000L * n  # separate seeds across n
  
  # bandwidth selection on the SAME grid used below
  bw <- select_bandwidths(
    n = n,
    seed_n = seed_n,
    beta1 = beta1,
    beta2 = beta2,
    beta_pi = beta_pi,
    hs = hs,
    lb0 = lb0,
    ub0 = ub0,
    a_grid0 = a_grid0,
    true_grid0 = true_grid0
  )
  
  h_lp <- bw$h_lp
  h_sm <- bw$h_sm
  lb   <- bw$lb
  ub   <- bw$ub
  
  # same grid used in bw selection
  a_grid    <- a_grid0
  true_grid <- true_grid0
  
  # replicae
  res_list <- future_lapply(
    X = 1:M,
    FUN = run_one_rep,
    n = n, seed_n = seed_n,
    lb = lb, ub = ub, a_grid = a_grid,
    beta1 = beta1, beta2 = beta2, beta_pi = beta_pi,
    h_lp = h_lp, h_sm = h_sm,
    future.seed = TRUE
  )
  
  # combine
  est_lp <- do.call(rbind, lapply(res_list, `[[`, "est_lp"))
  var_lp <- do.call(rbind, lapply(res_list, `[[`, "var_lp"))
  est_sm <- do.call(rbind, lapply(res_list, `[[`, "est_sm"))
  var_sm <- do.call(rbind, lapply(res_list, `[[`, "var_sm"))
  
  # coverage
  cover_lp <- (true_grid >= est_lp - zcrit * sqrt(var_lp)) &
    (true_grid <= est_lp + zcrit * sqrt(var_lp))
  
  cover_sm <- (true_grid >= est_sm - zcrit * sqrt(var_sm)) &
    (true_grid <= est_sm + zcrit * sqrt(var_sm))
  
  cover_prob_lp <- colMeans(cover_lp, na.rm = TRUE)
  cover_prob_sm <- colMeans(cover_sm, na.rm = TRUE)
  avg_cover_lp  <- mean(cover_prob_lp, na.rm = TRUE)
  avg_cover_sm  <- mean(cover_prob_sm, na.rm = TRUE)
  
  # average interval length
  len_lp <- 2 * zcrit * sqrt(var_lp)
  len_sm <- 2 * zcrit * sqrt(var_sm)
  
  avg_len_lp <- mean(len_lp, na.rm = TRUE)
  avg_len_sm <- mean(len_sm, na.rm = TRUE)
  
  # compare mean(var_hat) to empirical var(est)
  var_emp_lp <- apply(est_lp, 2, var, na.rm = TRUE)
  var_emp_sm <- apply(est_sm, 2, var, na.rm = TRUE)
  var_hat_lp <- colMeans(var_lp, na.rm = TRUE)
  var_hat_sm <- colMeans(var_sm, na.rm = TRUE)
  
  den_lp <- pmax(var_emp_lp, 1e-12)
  den_sm <- pmax(var_emp_sm, 1e-12)
  
  relerr_var_lp <- median(abs(var_hat_lp - var_emp_lp) / den_lp, na.rm = TRUE)
  relerr_var_sm <- median(abs(var_hat_sm - var_emp_sm) / den_sm, na.rm = TRUE)
  
  saveRDS(
    list(
      n = n,
      a_grid = a_grid,
      truth = true_grid,
      bandwidths = list(lp = h_lp, smooth = h_sm),
      estimates = list(lp = est_lp, smooth = est_sm),
      var_hat_rep = list(lp = var_lp, smooth = var_sm),
      interval_length_rep = list(lp = len_lp, smooth = len_sm),
      cov_ind = list(lp = cover_lp, smooth = cover_sm),
      var_emp = list(lp = var_emp_lp, smooth = var_emp_sm),
      var_hat_mean = list(lp = var_hat_lp, smooth = var_hat_sm),
      avg_cover = list(lp = avg_cover_lp, smooth = avg_cover_sm),
      avg_interval_length = list(lp = avg_len_lp, smooth = avg_len_sm),
      relerr_var = list(lp = relerr_var_lp, smooth = relerr_var_sm)
    ),
    file = paste0("simu3_results_", n, ".rds")
  )
  
  data.frame(
    n = n,
    avg_cover_lp = avg_cover_lp,
    avg_cover_sm = avg_cover_sm,
    avg_len_lp = avg_len_lp,
    avg_len_sm = avg_len_sm,
    relerr_var_lp = relerr_var_lp,
    relerr_var_sm = relerr_var_sm
  )
}

summ <- do.call(
  rbind,
  lapply(
    ns,
    run_for_n,
    seed = seed,
    beta1 = beta1,
    beta2 = beta2,
    beta_pi = beta_pi,
    hs = hs,
    M = M,
    zcrit = zcrit,
    lb0 = lb0,
    ub0 = ub0,
    a_grid0 = a_grid0,
    true_grid0 = true_grid0
  )
)

plan(sequential)

data_cov <- rbind(
  data.frame(n = summ$n, y = summ$avg_cover_lp, line = "LP"),
  data.frame(n = summ$n, y = summ$avg_cover_sm, line = "Smooth")
)

p_cov <- ggplot(data_cov, aes(x = n, y = y, color = line, group = line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  labs(x = "n",
       y = "Average coverage",
       color = "Method") +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )

ggsave("avg_cover.pdf", p_cov, width = 7, height = 7)

data_len <- rbind(
  data.frame(n = summ$n, y = summ$avg_len_lp, line = "LP"),
  data.frame(n = summ$n, y = summ$avg_len_sm, line = "Smooth")
)

p_len <- ggplot(data_len, aes(x = n, y = y, color = line, group = line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(x = "n",
       y = "Average interval length",
       color = "Method") +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.size = unit(1.5, "cm"),
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16)
  )

ggsave("avg_len.pdf", p_len, width = 7, height = 7)

# p_var <- ggplot(data_var, aes(x = n, y = y, color = line, group = line)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 1.5) +
#   labs(x = "n",
#        y = "Relative variance error",
#        color = "Method") +
#   theme_minimal() +
#   theme(
#     legend.text = element_text(size = 16),
#     legend.title = element_text(size = 16),
#     legend.key.size = unit(1.5, "cm"),
#     plot.title = element_text(hjust = 0.5, size = 16),
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 16)
#   )
# 
# ggsave("simu3_avg_var.pdf", p_var, width = 7, height = 7)
