rm(list = ls())

library(future)
library(future.apply)
library(ggplot2)
library(grid)

source("simu4_avg_utils.R")

seed <- 521

## grids / tuning
ns   <- seq(500, 5000, by = 500)
hs   <- seq(0.3, 1.5, by = 0.02)
thre <- 0.10
D    <- 80
M    <- 1000

## fixed support + truths
lb0 <- thre
ub0 <- 1 - thre

z_grid0  <- seq(lb0, ub0, length.out = D)
truth0   <- true_gamma(z_grid0)
truth_d0 <- true_denom(z_grid0)

zcrit <- qnorm(0.975)

n_workers <- min(availableCores(), 8)
plan(multisession, workers = n_workers)

select_bandwidths <- function(n, seed_n, hs, z_grid, truth_num, truth_den) {
  
  dat <- simulate_data_liv(seed = seed_n, n = n)
  y <- dat$y
  a <- dat$a
  z <- dat$z
  x <- as.matrix(dat[, c("x1", "x2", "x3", "x4")])
  
  mses_num_lp <- numeric(length(hs))
  mses_den_lp <- numeric(length(hs))
  mses_num_sm <- numeric(length(hs))
  mses_den_sm <- numeric(length(hs))
  
  for (i in seq_along(hs)) {
    h <- hs[i]
    cat("--- n =", n, "| h =", h, "---\n")
    
    ## avoid nested parallelism + IF computation during bandwidth search
    fit_num_lp <- dctseff(y, z, x, a.out = z_grid, h = h,
                          family = "gaussian", ifs = FALSE, parallel = FALSE)
    fit_den_lp <- dctseff(a, z, x, a.out = z_grid, h = h,
                          family = "binomial", ifs = FALSE, parallel = FALSE)
    
    fit_num_sm <- smooth_est(y, z, x, a.out = z_grid, h = h,
                             family = "gaussian", ifs = FALSE)
    fit_den_sm <- smooth_est(a, z, x, a.out = z_grid, h = h,
                             family = "binomial", ifs = FALSE)
    
    mses_num_lp[i] <- mean((fit_num_lp - truth_num)^2, na.rm = TRUE)
    mses_den_lp[i] <- mean((fit_den_lp - truth_den)^2, na.rm = TRUE)
    mses_num_sm[i] <- mean((fit_num_sm - truth_num)^2, na.rm = TRUE)
    mses_den_sm[i] <- mean((fit_den_sm - truth_den)^2, na.rm = TRUE)
  }
  
  hnum_lp <- hs[which.min(mses_num_lp)]
  hden_lp <- hs[which.min(mses_den_lp)]
  hnum_sm <- hs[which.min(mses_num_sm)]
  hden_sm <- hs[which.min(mses_den_sm)]
  
  hnum_lp <- hnum_lp / log(n)
  hden_lp <- hden_lp / log(n)
  hnum_sm <- hnum_sm / log(n)
  hden_sm <- hden_sm / log(n)
  
  list(
    hnum_lp = hnum_lp,
    hden_lp = hden_lp,
    hnum_sm = hnum_sm,
    hden_sm = hden_sm,
    mses = list(
      hs = hs,
      lp = list(num = mses_num_lp, den = mses_den_lp),
      smooth = list(num = mses_num_sm, den = mses_den_sm)
    )
  )
}

run_one_rep <- function(m, n, seed_n, z_grid,
                        hnum_lp, hden_lp, hnum_sm, hden_sm) {
  
  dat <- simulate_data_liv(seed = seed_n + m, n = n)
  y <- dat$y
  a <- dat$a
  z <- dat$z
  x <- as.matrix(dat[, c("x1", "x2", "x3", "x4")])
  
  fit_lp <- npliv_lp(
    y = y, a = a, x = x, z = z, z.out = z_grid,
    h1 = hnum_lp, h2 = hden_lp,
    family.y = "gaussian", family.a = "binomial", parallel = FALSE
  )
  
  fit_sm <- npliv_smooth(
    y = y, a = a, x = x, z = z, z.out = z_grid,
    h1 = hnum_sm, h2 = hden_sm,
    family.y = "gaussian", family.a = "binomial"
  )
  
  list(
    est_lp = fit_lp$est,
    se_lp  = fit_lp$sd_error,
    est_sm = fit_sm$est,
    se_sm  = fit_sm$sd_error
  )
}

run_for_n <- function(n, seed, hs, M, zcrit, z_grid0, truth0, truth_d0) {
  
  seed_n <- seed + 100000L * n  
  
  # bandwidth selection
  bw <- select_bandwidths(
    n = n,
    seed_n = seed_n,
    hs = hs,
    z_grid = z_grid0,
    truth_num = truth0,
    truth_den = truth_d0
  )
  
  hnum_lp <- bw$hnum_lp
  hden_lp <- bw$hden_lp
  hnum_sm <- bw$hnum_sm
  hden_sm <- bw$hden_sm
  
  z_grid <- z_grid0
  truth  <- truth0
  
  # replicate
  res_list <- future_lapply(
    X = 1:M,
    FUN = run_one_rep,
    n = n, seed_n = seed_n,
    z_grid = z_grid,
    hnum_lp = hnum_lp, hden_lp = hden_lp,
    hnum_sm = hnum_sm, hden_sm = hden_sm,
    future.seed = TRUE
  )
  
  # combine
  est_lp <- do.call(rbind, lapply(res_list, `[[`, "est_lp"))
  se_lp  <- do.call(rbind, lapply(res_list, `[[`, "se_lp"))
  est_sm <- do.call(rbind, lapply(res_list, `[[`, "est_sm"))
  se_sm  <- do.call(rbind, lapply(res_list, `[[`, "se_sm"))
  
  # coverage 
  cov_lp <- (truth >= est_lp - zcrit * se_lp) & (truth <= est_lp + zcrit * se_lp)
  cov_sm <- (truth >= est_sm - zcrit * se_sm) & (truth <= est_sm + zcrit * se_sm)
  
  cover_prob_lp <- colMeans(cov_lp, na.rm = TRUE)
  cover_prob_sm <- colMeans(cov_sm, na.rm = TRUE)
  avg_cover_lp  <- mean(cover_prob_lp, na.rm = TRUE)
  avg_cover_sm  <- mean(cover_prob_sm, na.rm = TRUE)
  
  # interval lengths
  len_lp <- 2 * zcrit * se_lp
  len_sm <- 2 * zcrit * se_sm
  
  avg_len_lp <- mean(len_lp, na.rm = TRUE)
  avg_len_sm <- mean(len_sm, na.rm = TRUE)
  
  # variance estimation error
  var_emp_lp <- apply(est_lp, 2, var, na.rm = TRUE)
  var_emp_sm <- apply(est_sm, 2, var, na.rm = TRUE)
  
  var_hat_lp <- colMeans(se_lp^2, na.rm = TRUE)
  var_hat_sm <- colMeans(se_sm^2, na.rm = TRUE)
  
  den_lp <- pmax(var_emp_lp, 1e-12)
  den_sm <- pmax(var_emp_sm, 1e-12)
  
  relerr_var_lp <- median(abs(var_hat_lp - var_emp_lp) / den_lp, na.rm = TRUE)
  relerr_var_sm <- median(abs(var_hat_sm - var_emp_sm) / den_sm, na.rm = TRUE)
  
  saveRDS(list(
    n = n,
    z_grid = z_grid,
    truth = truth,
    truth_denom = truth_d0,
    bandwidths = list(
      lp = list(num = hnum_lp, den = hden_lp),
      smooth = list(num = hnum_sm, den = hden_sm)
    ),
    bandwidth_search = bw$mses,
    estimates = list(lp = est_lp, smooth = est_sm),
    se_rep = list(lp = se_lp, smooth = se_sm),
    interval_length_rep = list(lp = len_lp, smooth = len_sm),
    cov_ind = list(lp = cov_lp, smooth = cov_sm),
    cover_prob = list(lp = cover_prob_lp, smooth = cover_prob_sm),
    var_emp = list(lp = var_emp_lp, smooth = var_emp_sm),
    var_hat_mean = list(lp = var_hat_lp, smooth = var_hat_sm),
    avg_cover = list(lp = avg_cover_lp, smooth = avg_cover_sm),
    avg_interval_length = list(lp = avg_len_lp, smooth = avg_len_sm),
    relerr_var = list(lp = relerr_var_lp, smooth = relerr_var_sm)
  ), file = paste0("simu4_results_", n, ".rds"))
  
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
    ns, run_for_n,
    seed = seed,
    hs = hs,
    M = M,
    zcrit = zcrit,
    z_grid0 = z_grid0,
    truth0 = truth0,
    truth_d0 = truth_d0
  )
)

plan(sequential) 

print(summ)
saveRDS(summ, file = "avg_metrics_vs_n_simu4.rds")

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

ggsave("simu4_avg_cover.pdf", p_cov, width = 7, height = 7)

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

ggsave("simu4_avg_len.pdf", p_len, width = 7, height = 7)

data_var <- rbind(
  data.frame(n = summ$n, y = summ$relerr_var_lp, line = "LP"),
  data.frame(n = summ$n, y = summ$relerr_var_sm, line = "Smooth")
)

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
# ggsave("simu4_avg_var.pdf", p_var, width = 7, height = 7)