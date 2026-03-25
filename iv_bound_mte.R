library(gurobi)
# install.packages("ivmte")
library(ivmte)

# A template
res <- ivmte(
  data    = dat,
  target  = "ate",              # we want bounds on the ATE
  
  ## Regression-based criterion:
  treat   = D,                  # treatment indicator in the dataset dat
  outcome = Y,                  # outcome in the dataset dat
  
  ## MTR specifications (MTR(0,u,x), MTR(1,u,x)):
  ## Here: cubic B-splines in u + linear in covariates
  m0 = ~ uSplines(degree = 3,
                  knots  = c(0.25, 0.50, 0.75)) +
    X1 + X2, # Replace with our covariates, I think linear is enough...
  m1 = ~ uSplines(degree = 3,
                  knots  = c(0.25, 0.50, 0.75)) +
    X1 + X2, # Replace with our covariates, I think linear is enough...
  
  uname = u,                    # name of the latent U in the MTRs (symbolic)
  
  ## Propensity score model P(D=1 | Z, X):
  propensity = D ~ Z + X1 + X2, # Replace with our covariates, I think linear is enough...
  
  ## Optimization solver for the bounds:
  solver = "gurobi",        
  
  ## (optional) no bootstrap inference for now:
  bootstraps = 0,
  point = F
)

# ATE bounds:
res$bounds


# A toy example


set.seed(123)

## Simulate a simple binary-treatment IV dataset
n <- 2000
dat <- data.frame(
  X1 = rnorm(n),
  X2 = rbinom(n, 1, 0.5),
  Z  = rbinom(n, 1, 0.5)        # instrument
)

# Treatment assignment: D depends on Z and X
dat$D <- rbinom(
  n, 1,
  plogis(-0.3 + 0.8 * dat$Z + 0.4 * dat$X1 - 0.5 * dat$X2)
)

# Outcome with treatment effect heterogeneity (just for illustration)
tau <- 1 + 0.5 * dat$X1        # heterogeneous effect
mu0 <- 0.5 + 0.3 * dat$X1 - 0.2 * dat$X2
dat$Y <- mu0 + tau * dat$D + rnorm(n)

## Regression-based MTE bounds on the ATE
res <- ivmte(
  data    = dat,
  target  = "ate",
  
  treat   = D,
  outcome = Y,
  
  # MTR(0,u,x) and MTR(1,u,x): splines in u + covariates
  m0 = ~ uSplines(degree = 3,
                  knots  = c(0.25, 0.5, 0.75)) +
    X1 + X2,
  m1 = ~ uSplines(degree = 3,
                  knots  = c(0.25, 0.5, 0.75)) +
    X1 + X2,
  
  uname = u,
  
  # Propensity P(D=1 | Z, X)
  propensity = D ~ Z + X1 + X2,
  
  solver     = "gurobi",     # or "lpSolveAPI" etc.
  bootstraps = 0,
  point = F
)

res$bounds





