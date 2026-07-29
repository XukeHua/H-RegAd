rm(list = ls())
options(stringsAsFactors = FALSE)

load("HIV_data/sCD14wC.rda")

source("proj_linf_ball.R")
source("soft_thresh.R")
source("sGS2w.R")

library(robregcc)
library(MASS)
library(dplyr)
library(readr)

robust_sigma_hat <- function(x) {
  s <- mad(x, center = median(x), constant = 1.4826)
  if (!is.finite(s) || s <= 1e-12) s <- sd(x)
  if (!is.finite(s) || s <= 1e-12) s <- 1
  s
}

huber_loss_mean_c <- function(r, c) {
  a <- abs(r)
  mean(ifelse(a <= c, 0.5 * r^2, c * (a - 0.5 * c)))
}

count_selected <- function(beta, tol = 1e-8) {
  sum(abs(beta[-1]) > tol)
}

count_outliers <- function(gamma, tol = 1e-8) {
  sum(abs(gamma) > tol)
}

C <- matrix(1, nrow = 1, ncol = 60)

W <- sCD14[, 1:60]
X <- t(apply(W, 1, function(x) {
  indx <- x != 0
  x[!indx] <- 0.5
  x / sum(x)
}))

y <- sCD14[, 61]

exp_seed <- 123
set.seed(exp_seed)

Xt <- cbind(1, log(X))                               
Ct <- cbind(matrix(0, nrow = nrow(C), ncol = 1), C) 
C_sp <- t(Ct)                                        

n <- nrow(Xt)
p <- ncol(Xt)
bw <- c(0, rep(1, p - 1))   


control <- robregcc_option(maxiter = 5000, tol = 1e-12, lminfac = 1e-12)
fit.nr <- classo(Xt, y, Ct, we = bw, type = 1, control = control)

set.seed(exp_seed)
b1 <- 0.25
cc1 <- 2.937
control <- robregcc_option(maxiter = 5000, tol = 1e-4, lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y, alp = 0.4, cfac = 2, b1 = b1, cc1 = cc1, Ct, bw, 1, control)

set.seed(exp_seed)
control <- robregcc_option()
beta.wt <- fit.init$betaR
beta.wt[1] <- 0

control$gamma <- 1
control$spb <- 50 / p
control$outMiter <- 1000
control$inMiter <- 3000
control$nlam <- 100
control$lmaxfac <- 1
control$lminfac <- 1e-8
control$tol <- 1e-20
control$out.tol <- 1e-16
control$kfold <- 10
control$sigmafac <- 2

# [A]
fit.ada <- robregcc_sp(
  Xt, y, Ct,
  beta.init = beta.wt, cindex = 1,
  gamma.init = fit.init$residuals,
  control = control,
  penalty.index = 1, alpha = 0.95
)

# [E]
set.seed(exp_seed)
control$lmaxfac <- 1
control$lminfac <- 1e-8
fit.soft <- robregcc_sp(
  Xt, y, Ct, cindex = 1,
  control = control,
  penalty.index = 2, alpha = 0.95
)

# [H]
set.seed(exp_seed)
control$lmaxfac <- 1e1
control$lminfac <- 1e-2
control$sigmafac <- 2
fit.hard <- robregcc_sp(
  Xt, y, Ct,
  beta.init = fit.init$betaf,
  gamma.init = fit.init$residuals,
  cindex = 1,
  control = control,
  penalty.index = 3, alpha = 0.95
)

set.seed(exp_seed)
K <- 10
fold_id <- sample(rep(1:K, length.out = n))

nlam_sp <- 150
lam_max <- max(abs(as.numeric(crossprod(Xt, y)))) / n
lam_min <- lam_max * 1e-12
lambda_seq_sp <- exp(seq(log(lam_max), log(lam_min), length.out = nlam_sp))

cv_sp <- numeric(nlam_sp)

for (ii in seq_along(lambda_seq_sp)) {
  lam <- lambda_seq_sp[ii]
  score_k <- numeric(K)
  
  for (kk in 1:K) {
    idx_te <- which(fold_id == kk)
    idx_tr <- setdiff(seq_len(n), idx_te)
    
    X_tr <- Xt[idx_tr, , drop = FALSE]
    y_tr <- y[idx_tr]
    X_te <- Xt[idx_te, , drop = FALSE]
    y_te <- y[idx_te]
    
    fit.sp <- spadmm_lasso_C_with_stop(
      X = X_tr, y = y_tr, C = C_sp,
      lambda1 = lam,
      tau = 1.345, sigma = 1.0, rho = 1.0,
      bw = bw,
      maxiter = 3000, tol = 1e-6,
      verbose = FALSE
    )
    
    yhat_te <- as.numeric(X_te %*% fit.sp$beta)
    res_te  <- y_te - yhat_te
    
    sigma_hat <- robust_sigma_hat(y_tr)
    c_thr <- 1.345 * sigma_hat
    
    score_k[kk] <- huber_loss_mean_c(res_te, c = c_thr) / (c_thr^2)
  }
  
  cv_sp[ii] <- mean(score_k)
  if (ii %% 10 == 0) cat("[H-ADMM CV]", ii, "/", nlam_sp, "\n")
}

i_min_sp <- which.min(cv_sp)
lambda_use_sp <- lambda_seq_sp[i_min_sp]

cat("H-ADMM chosen lambda =", lambda_use_sp, "\n")
cat("lambda_use_sp / lam_max =", lambda_use_sp / lam_max, "\n")

set.seed(exp_seed)
fit.sp.full <- spadmm_lasso_C_with_stop(
  X = Xt, y = y, C = C_sp,
  lambda1 = lambda_use_sp,
  tau = 1.345, sigma = 1.0, rho = 1.0,
  bw = bw,
  maxiter = 10000, tol = 1e-10,
  verbose = TRUE
)

beta_sp <- fit.sp.full$beta
res_sp0 <- drop(y) - as.numeric(Xt %*% beta_sp)
c_thr_full <- 1.345 * robust_sigma_hat(y)

kappa_seq <- as.numeric(
  quantile(abs(res_sp0),
           probs = seq(0.80, 0.995, length.out = 120),
           names = FALSE)
)

obj_seq <- numeric(length(kappa_seq))

for (i in seq_along(kappa_seq)) {
  kappa <- kappa_seq[i]
  gamma_hat <- soft_thresh(res_sp0, kappa)
  res_corr <- res_sp0 - gamma_hat
  obj_seq[i] <- huber_loss_mean_c(res_corr, c = c_thr_full) / (c_thr_full^2)
}

kappa_best <- kappa_seq[which.min(obj_seq)]
gamma_sp <- soft_thresh(res_sp0, kappa_best)

cat("H-ADMM selected kappa =", kappa_best, "\n")

s_index <- 1   

beta_lam_cv <- cbind(
  NR       = fit.nr$beta,
  A        = fit.ada$betaE[, s_index],
  E        = fit.soft$betaE[, s_index],
  H        = fit.hard$betaE[, s_index],
  `H-ADMM` = beta_sp
)

shift_est_lam_cv <- cbind(
  NR       = rep(0, length(y)),
  A        = fit.ada$gamma0[, s_index],
  E        = fit.soft$gamma0[, s_index],
  H        = fit.hard$gamma0[, s_index],
  `H-ADMM` = gamma_sp
)

fitted_mat <- Xt %*% beta_lam_cv + shift_est_lam_cv
residual_mat <- matrix(y, nrow = length(y), ncol = ncol(fitted_mat)) - fitted_mat
colnames(residual_mat) <- colnames(beta_lam_cv)

R_square <- 1 - apply(residual_mat, 2, function(e) sum(e^2)) / sum((y - mean(y))^2)

summary_table <- lapply(colnames(beta_lam_cv), function(m) {
  res_m <- residual_mat[, m]
  beta_m <- beta_lam_cv[, m]
  gamma_m <- shift_est_lam_cv[, m]
  
  data.frame(
    Method = m,
    Selected_Taxa = count_selected(beta_m),
    Flagged_Outliers = count_outliers(gamma_m),
    R2 = as.numeric(R_square[m]),
    MAE = mean(abs(res_m)),
    MedAE = median(abs(res_m)),
    Scaled_Huber_Error = huber_loss_mean_c(res_m, c = c_thr_full) / (c_thr_full^2)
  )
}) %>% bind_rows()

summary_table <- summary_table %>%
  mutate(
    R2 = round(R2, 4),
    MAE = round(MAE, 4),
    MedAE = round(MedAE, 4),
    Scaled_Huber_Error = round(Scaled_Huber_Error, 4)
  )

print(summary_table)

predictor_name <- c("intercept", colnames(X))
ord_predictor <- order(predictor_name)

predictor_est_c1 <- data.frame(
  beta_lam_cv[ord_predictor, ],
  v.name = predictor_name[ord_predictor]
)

write_csv(summary_table, "hiv_model_summary.csv")

save(
  summary_table,
  predictor_est_c1,
  file = "hiv_model_summary.rda"
)

save(
  y, Xt, X,
  beta_lam_cv,
  shift_est_lam_cv,
  fitted_mat,
  residual_mat,
  R_square,
  predictor_est_c1,
  summary_table,
  file = "hiv_fit_objects.rda"
)

cat("\nSaved files:\n")
cat("  - hiv_model_summary.csv\n")
cat("  - hiv_model_summary.rda\n")
cat("  - hiv_fit_objects.rda\n")
