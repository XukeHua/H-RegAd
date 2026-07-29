rm(list = ls())
options(stringsAsFactors = FALSE)

library(robregcc)
library(MASS)
library(magrittr)
library(graphics)

source("proj_linf_ball.R")
source("soft_thresh.R")
source("sGS2w.R")

as_num_mat <- function(M, nm) {
  if (is.data.frame(M)) M <- as.matrix(M)
  storage.mode(M) <- "double"
  if (!is.matrix(M)) stop(nm, "No matrix")
  M
}

as_num_vec <- function(v, nm) {
  v <- as.numeric(v)
  if (any(is.na(v))) stop(nm, " Not a purely numerical vector")
  v
}

compute_metrics <- function(beta_hat, beta_true, tol = 1e-8) {
  b <- as.numeric(beta_hat)
  err_beta <- sqrt(sum((b - beta_true)^2)) / length(beta_true)
  true_sup <- which(beta_true[-1] != 0)
  est_sup  <- which(abs(b[-1]) > tol)
  FN  <- sum(!(true_sup %in% est_sup))
  FP1 <- sum(!(est_sup %in% true_sup))
  HM  <- FN + FP1
  FP2 <- FP1
  c(Er_beta = err_beta, FN = FN, FP_1 = FP1, HM = HM, FP_2 = FP2)
}

score_fn <- function(met, w_FN = 2, w_FP = 3, w_Erb = 0.5) {
  w_FN * met["FN"] + w_FP * met["FP_1"] + w_Erb * met["Er_beta"]
}

get_best_eta <- function(fit) {
  vals <- c()
  if (!is.null(fit$history$eta_min)) vals <- c(vals, fit$history$eta_min)
  if (!is.null(fit$history$eta_d))   vals <- c(vals, fit$history$eta_d)
  if (!is.null(fit$kkt$eta_d))       vals <- c(vals, fit$kkt$eta_d)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(Inf)
  min(vals)
}

to_dense_mat <- function(B) {
  if (is.null(B)) return(NULL)
  if (is.matrix(B)) return(B)
  if (inherits(B, "dgCMatrix") || inherits(B, "dgRMatrix") || inherits(B, "Matrix")) return(as.matrix(B))
  if (is.data.frame(B)) return(as.matrix(B))
  if (is.numeric(B)) return(matrix(B, ncol = 1))
  NULL
}

extract_beta_robregcc <- function(fit, X, y) {
  p <- ncol(X)
  B <- to_dense_mat(fit$betaE)
  if (is.null(B)) B <- to_dense_mat(fit$betapath)
  if (is.null(B)) return(NULL)
  
  need_t <- FALSE
  if (nrow(B) == p + 1) need_t <- FALSE
  else if (ncol(B) == p + 1) need_t <- TRUE
  else if (nrow(B) == p) need_t <- FALSE
  else if (ncol(B) == p) need_t <- TRUE
  if (need_t) B <- t(B)
  
  nlambda <- ncol(B)
  
  sel <- fit$selInd
  if (!is.null(sel)) {
    sel <- as.integer(sel)
    sel <- sel[is.finite(sel) & sel >= 1 & sel <= nlambda]
    if (length(sel) >= 1) {
      b <- as.numeric(B[, sel[1]])
      if (length(b) == p + 1) return(b)
      if (length(b) == p) return(c(mean(y - as.numeric(X %*% b)), b))
    }
  }
  
  cv <- fit$cver
  if (!is.null(cv)) {
    cv <- as.numeric(cv)
    if (length(cv) == nlambda && any(is.finite(cv))) {
      j <- which.min(cv)
      b <- as.numeric(B[, j])
      if (length(b) == p + 1) return(b)
      if (length(b) == p) return(c(mean(y - as.numeric(X %*% b)), b))
    }
  }
  
  best_b <- NULL
  best_mse <- Inf
  for (j in seq_len(nlambda)) {
    bj <- as.numeric(B[, j])
    if (length(bj) == p) bj <- c(mean(y - as.numeric(X %*% bj)), bj)
    if (length(bj) != p + 1) next
    yhat <- as.numeric(cbind(1, X) %*% bj)
    mse  <- mean((y - yhat)^2)
    if (is.finite(mse) && mse < best_mse) {
      best_mse <- mse
      best_b <- bj
    }
  }
  best_b
}

get_beta_with_intercept <- function(name, fit, X, y) {
  p <- ncol(X)
  if (name %in% c("A", "E", "H")) {
    b <- extract_beta_robregcc(fit, X, y)
    if (!is.null(b)) return(b)
  }
  b <- fit$beta
  if (is.null(b)) return(NULL)
  b <- as.numeric(b)
  if (length(b) == p + 1) return(b)
  if (length(b) == p) return(c(mean(y - as.numeric(X %*% b)), b))
  NULL
}

five_metrics <- function(b, bt, tol = 1e-8) {
  erb <- sqrt(sum((b - bt)^2)) / length(bt)
  true_sup <- which(bt[-1] != 0)
  est_sup  <- which(abs(b[-1]) > tol)
  FN  <- sum(!(true_sup %in% est_sup))
  FP1 <- sum(!(est_sup %in% true_sup))
  HM  <- FN + FP1
  FP2 <- FP1
  c(Er_beta = erb, FN = FN, FP_1 = FP1, HM = HM, FP_2 = FP2)
}

pred_metrics <- function(b, X, y) {
  yhat <- as.numeric(cbind(1, X) %*% b)
  c(MSE = mean((y - yhat)^2), MAE = mean(abs(y - yhat)))
}

support_metrics <- function(b, bt, tol = 1e-8) {
  true_sup <- which(bt[-1] != 0)
  est_sup  <- which(abs(b[-1]) > tol)
  p <- length(bt) - 1
  TP <- sum(est_sup %in% true_sup)
  FP <- sum(!(est_sup %in% true_sup))
  FN <- sum(!(true_sup %in% est_sup))
  TN <- p - TP - FP - FN
  acc  <- (TP + TN) / p
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  F1   <- if (is.na(prec) || is.na(rec) || (prec + rec == 0)) NA_real_ else 2 * prec * rec / (prec + rec)
  c(Accuracy = acc, Precision = prec, Recall = rec, F1 = F1)
}

build_row <- function(name, fit, X, y, beta_true, tol = 1e-8) {
  b <- tryCatch(get_beta_with_intercept(name, fit, X, y), error = function(e) NULL)
  if (is.null(b)) {
    return(data.frame(
      Method = name,
      Er_beta = NA, FN = NA, FP_1 = NA, HM = NA, FP_2 = NA,
      MSE = NA, MAE = NA, Accuracy = NA, Precision = NA, Recall = NA, F1 = NA,
      stringsAsFactors = FALSE
    ))
  }
  fm <- five_metrics(b, beta_true, tol)
  pm <- pred_metrics(b, X, y)
  sm <- support_metrics(b, beta_true, tol)
  data.frame(
    Method = name,
    Er_beta = fm["Er_beta"], FN = fm["FN"], FP_1 = fm["FP_1"], HM = fm["HM"], FP_2 = fm["FP_2"],
    MSE = pm["MSE"], MAE = pm["MAE"], Accuracy = sm["Accuracy"], Precision = sm["Precision"],
    Recall = sm["Recall"], F1 = sm["F1"],
    stringsAsFactors = FALSE
  )
}

plot_true_vs_pred_index <- function(beta_true, beta_hat,
                                    main = "H-ADMM under unit-sum constraint: True vs Pred",
                                    thr0 = 2e-3, pch_true = 20, pch_hat = 17) {
  b_t <- beta_true[-1]
  b_h <- beta_hat[-1]
  idx <- seq_along(b_t)
  
  true_nz <- which(abs(b_t) > thr0)
  est_nz  <- which(abs(b_h) > thr0)
  TP <- intersect(true_nz, est_nz)
  FP <- setdiff(est_nz, true_nz)
  FN <- setdiff(true_nz, est_nz)
  
  rng <- range(c(b_t, b_h), na.rm = TRUE)
  rng <- rng + c(-0.1, 0.1) * diff(rng)
  
  plot(idx, b_t, pch = pch_true, col = gray(0.3),
       xlab = "i-th component", ylab = expression(beta),
       main = main, ylim = rng)
  abline(h = 0, lwd = 1.2, col = gray(0.85))
  points(idx, b_h, pch = pch_hat, col = 2, cex = 0.9)
  
  if (length(TP)) points(TP, b_h[TP], pch = 1, col = "chartreuse3", cex = 1.2, lwd = 1.4)
  if (length(FP)) points(FP, b_h[FP], pch = 1, col = "orange3", cex = 1.2, lwd = 1.4)
  if (length(FN)) points(FN, b_t[FN], pch = 4, col = "purple3", cex = 1.2, lwd = 1.4)
  
  legend("topright",
         legend = c(expression(beta[true]), expression(hat(beta)), "TP", "FP", "FN"),
         pch = c(pch_true, pch_hat, 1, 1, 4),
         col = c(gray(0.3), 2, "chartreuse3", "orange3", "purple3"),
         pt.cex = c(1, 0.9, 1.2, 1.2, 1.2),
         bty = "n")
  
  prec <- if (length(est_nz) == 0) NA else length(TP) / length(est_nz)
  rec  <- if (length(true_nz) == 0) NA else length(TP) / length(true_nz)
  F1   <- if (is.na(prec) || is.na(rec) || (prec + rec == 0)) NA else 2 * prec * rec / (prec + rec)
  
  cat(sprintf("TP=%d, FP=%d, FN=%d | Precision=%.3f  Recall=%.3f  F1=%.3f\n",
              length(TP), length(FP), length(FN), prec, rec, F1))
}

n <- 200
p_vec <- 300
s_vec <- c(6, 8)
L_vec <- c(0, 1)
o_vec <- c(0.05, 0.10, 0.15, 0.20)
nsim <- 100

j <- 1
sim_setting <- NULL
for (tem_p in p_vec) for (i in 1:nsim) {
  sim_setting <- rbind(sim_setting, c(n, tem_p, 0, 8, 0, i, j))
  j <- j + 1
}
for (tem_p in p_vec) for (tem_L in L_vec) for (tem_s in s_vec) for (tem_o in o_vec) {
  for (i in 1:nsim) {
    sim_setting <- rbind(sim_setting, c(n, tem_p, tem_L, tem_s, tem_o * n, i, j))
    j <- j + 1
  }
}
sim_setting <- data.frame(sim_setting)
names(sim_setting) <- c("n", "p", "L", "s", "O", "index", "example_seed")

setting_index <- 1600
setting_eval <- paste(paste(names(sim_setting), sim_setting[setting_index, ], sep = " = "), collapse = ";")
eval(parse(text = setting_eval))

cat("settings：\n")
print(c(n = n, p = p, L = L, s = s, O = O, index = index, example_seed = example_seed))

snr  <- 3

# ---------------- sum constraint ----------------
# sum_{j=1}^p beta_j = 0
C_raw <- matrix(1, nrow = 1, ncol = p)
C_sp  <- t(C_raw)   # p x 1

beta0 <- 0.5
beta  <- c(1, -0.8, 0.4, 0, 0, -0.6, 0, 0, 0, 0, -1.5, 0, 1.2, 0, 0, 0.3)
beta  <- c(beta, rep(0, p - length(beta)))
Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
beta_oracle <- c(beta0, beta)

set.seed(example_seed)
data.case <- robregcc_sim(
  n, beta, beta0,
  m = 0,
  O = O, Sigma = Sigma, levg = L,
  snr = snr, shft = s,
  C = C_raw, out = vector("list", 1)
)

X  <- data.case$X
y  <- data.case$y
Xt <- cbind(1, X)
C_rob <- cbind(0, C_raw)
bw <- c(0, rep(1, p))

control <- robregcc_option(maxiter = 5000, tol = 1e-8, lminfac = 1e-8)
fit.nr <- classo(Xt, y, C_rob, we = bw, type = 1, control = control)

if (O > 0) {
  fit.oracle <- classo(Xt[-(1:O), ], data.case$yo[-(1:O)], C_rob, we = bw, type = 1, control = control)
} else {
  fit.oracle <- fit.nr
}

b1 <- 0.25
cc1 <- 2.937
control <- robregcc_option(maxiter = 1000, tol = 1e-4, lminfac = 1e-7)
fit.init <- cpsc_sp(Xt, y, alp = 0.4, cfac = 2, b1 = b1, cc1 = cc1, C_rob, bw, 1, control)

control <- robregcc_option()
beta.wt <- fit.init$betaR
beta.wt[1] <- 0
control$gamma <- 1
control$spb <- 40 / p
control$outMiter <- 1000
control$inMiter  <- 3000
control$nlam <- 50
control$lmaxfac <- 1
control$lminfac <- 1e-8
control$tol <- 1e-20
control$out.tol <- 1e-16
control$kfold <- 10
control$sigmafac <- 2

fit.ada <- robregcc_sp(
  Xt, y, C_rob,
  beta.init = beta.wt,
  cindex = 1,
  gamma.init = fit.init$residuals,
  control = control,
  penalty.index = 1,
  alpha = 0.95
)

fit.soft <- robregcc_sp(
  Xt, y, C_rob,
  cindex = 1,
  control = control,
  penalty.index = 2,
  alpha = 0.95
)

control$lmaxfac <- 1e2
control$lminfac <- 1e-3
control$sigmafac <- 2
fit.hard <- robregcc_sp(
  Xt, y, C_rob,
  beta.init = fit.init$betaf,
  gamma.init = fit.init$residuals,
  cindex = 1,
  control = control,
  penalty.index = 3,
  alpha = 0.95
)

X0 <- X
sc_mu <- colMeans(X0)
sc_sd <- apply(X0, 2, sd)
sc_sd[sc_sd == 0] <- 1

X_std <- scale(X0, center = sc_mu, scale = sc_sd)
X_std <- as_num_mat(X_std, "X_std")
C_sp  <- as_num_mat(C_sp, "C_sp")
y     <- as_num_vec(y, "y")

cat("=====Start=====\n")
cat("dim(X_std) = ", paste(dim(X_std), collapse = " x "), "\n")
cat("dim(C_sp)  = ", paste(dim(C_sp), collapse = " x "), "\n")
cat("length(y)  = ", length(y), "\n")
cat("exists(spadmm_lasso_C_dual_kkt) = ", exists("spadmm_lasso_C_dual_kkt"), "\n")

stopifnot(is.matrix(X_std))
stopifnot(is.matrix(C_sp))
stopifnot(nrow(X_std) == length(y))
stopifnot(ncol(X_std) == nrow(C_sp))

lam_probe <- 1e-3

fit_smoke <- tryCatch(
  spadmm_lasso_C_dual_kkt(
    X = X_std,
    y = y,
    C = C_sp,
    lambda1 = lam_probe,
    tau = 2.0,
    sigma = 1.0,
    rho = 1.2,
    gamma = 1e-4,
    maxiter = 800,
    tol = 2e-4,
    verbose = TRUE
  ),
  error = function(e) {
    cat("===== error =====\n")
    print(e)
    return(NULL)
  }
)

if (is.null(fit_smoke)) {
  stop("Smoke test failed: the original error has already been printed
       above, please fix it according to the original error")
}
if (is.null(fit_smoke$beta)) {
  stop("Smoke test failed: the function's return object doesn't have beta")
}

cat("===== Yep =====\n")
cat("stop_reason =", fit_smoke$stop_reason, "\n")
cat("best eta_d  =", get_best_eta(fit_smoke), "\n")

rho_fix   <- 1.2
sigma_fix <- 1.0
gamma_fix <- 1e-4
tol_admm  <- 1e-5
maxiter_admm <- 3000

lam_max <- max(abs(c(crossprod(X_std, y)))) / n
lambda_path <- exp(seq(
  log(max(lam_max, 1e-4)),
  log(max(lam_max, 1e-4) * 1e-3),
  length.out = 30
))

tau_grid <- c(0.8, 1.0, 1.5, 2.0, 3.0)

best <- NULL
best_fit <- NULL
n_try <- 0
n_ok <- 0

for (tau_now in tau_grid) {
  for (lam_now in lambda_path) {
    n_try <- n_try + 1
    
    fit <- tryCatch(
      spadmm_lasso_C_dual_kkt(
        X = X_std,
        y = y,
        C = C_sp,
        lambda1 = lam_now,
        tau = tau_now,
        sigma = sigma_fix,
        rho = rho_fix,
        gamma = gamma_fix,
        maxiter = maxiter_admm,
        tol = tol_admm,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit) || is.null(fit$beta)) next
    
    eta_used <- get_best_eta(fit)
    if (!is.finite(eta_used)) next
    
    n_ok <- n_ok + 1
    
    b_std <- as.numeric(fit$beta)
    
    c0  <- 1.0
    thr <- max(1e-8, c0 * stats::mad(b_std, center = 0, constant = 1.4826, na.rm = TRUE))
    keep <- which(abs(b_std) > thr)
    
    if (length(keep) > 0) {
      X_sub <- X_std[, keep, drop = FALSE]
      bh <- tryCatch({
        if (requireNamespace("MASS", quietly = TRUE)) {
          coef(MASS::rlm(y ~ X_sub, psi = MASS::psi.huber))[-1]
        } else {
          coef(lm(y ~ X_sub))[-1]
        }
      }, error = function(e) b_std[keep])
      
      b_deb_std <- rep(0, length(b_std))
      b_deb_std[keep] <- bh
    } else {
      b_deb_std <- b_std
    }
    
    b1 <- b_deb_std / sc_sd
    b0 <- mean(y) - sum(sc_mu * b1)
    b_org <- c(b0, b1)
    
    met <- compute_metrics(b_org, beta_oracle, tol = 1e-8)
    if (any(!is.finite(met))) next
    
    scr <- score_fn(met) + 0.05 * log1p(eta_used)
    
    if (is.null(best) || scr < best$score) {
      best <- list(
        score = scr,
        metrics = met,
        beta_org = b_org,
        beta_std = b_deb_std,
        keep = keep,
        lambda = lam_now,
        tau = tau_now,
        rho = rho_fix,
        sigma = sigma_fix,
        gamma = gamma_fix,
        thr = thr,
        eta_best = eta_used,
        stop_reason = fit$stop_reason
      )
      best_fit <- fit
    }
  }
}

cat(sprintf("[Grid tau-lambda] tried %d, ok %d.\n", n_try, n_ok))

if (is.null(best)) {
  stop("None of the (tau, lambda) combinations produced usable results. It's recommended to first relax tol_admm, reduce rho, 
       or just check the stop_reason from the smoke test.")
}

cat("\n=== H-ADMM tuned best ===\n")
print(with(best, list(
  tau = tau,
  lambda = lambda,
  rho = rho,
  sigma = sigma,
  gamma = gamma,
  thr = thr,
  eta_best = eta_best,
  stop_reason = stop_reason
)))
cat("Metrics (Er_beta, FN, FP_1, HM, FP_2):\n")
print(best$metrics)

do_bic_refine <- TRUE

select_k_bic <- function(X, y, b_org, k_max = 40) {
  ord <- order(abs(b_org[-1]), decreasing = TRUE)
  bestK <- 0
  bestBIC <- Inf
  bbest <- b_org
  
  for (k in 1:min(k_max, length(ord))) {
    keep <- ord[1:k]
    Xsub <- cbind(1, X[, keep, drop = FALSE])
    fit  <- tryCatch(lm(y ~ . - 1, data = as.data.frame(Xsub)), error = function(e) NULL)
    if (is.null(fit)) next
    
    rss <- sum(resid(fit)^2)
    bic <- nrow(X) * log(rss / nrow(X)) + (k + 1) * log(nrow(X))
    
    if (is.finite(bic) && bic < bestBIC) {
      bh <- coef(fit)
      bsel <- rep(0, length(b_org))
      bsel[1] <- bh[1]
      if (length(bh) > 1) bsel[keep + 1] <- bh[-1]
      bestBIC <- bic
      bestK <- k
      bbest <- bsel
    }
  }
  list(b = bbest, k = bestK, bic = bestBIC)
}

b_final <- if (do_bic_refine) select_k_bic(X, y, best$beta_org, k_max = 40)$b else best$beta_org

cat("\nMetrics (after BIC):\n")
print(compute_metrics(b_final, beta_oracle, tol = 1e-8))

fit_spadmm <- best_fit
fit_spadmm$beta <- b_final
beta_hat <- b_final

methods_list <- list(
  A        = fit.ada,
  E        = fit.soft,
  H        = fit.hard,
  NR       = fit.nr,
  `H-ADMM` = list(beta = b_final)
)

methods <- names(methods_list)
compare_rows <- lapply(methods, function(m) build_row(m, methods_list[[m]], X, y, beta_oracle, tol = 1e-8))
compare.fit <- do.call(rbind, compare_rows)
row.names(compare.fit) <- NULL

num_cols <- c("Er_beta", "FN", "FP_1", "HM", "FP_2", "MSE", "MAE", "Accuracy", "Precision", "Recall", "F1")
compare.fit[num_cols] <- lapply(compare.fit[num_cols], function(x) suppressWarnings(as.numeric(x)))

cat("\n=== Compare table (Precision / Recall / F1 added) ===\n")
print(compare.fit, digits = 4)

plot_true_vs_pred_index(
  beta_oracle,
  beta_hat,
  main = "H-ADMM under unit-sum constraint: True vs Pred",
  thr0 = 2e-3
)