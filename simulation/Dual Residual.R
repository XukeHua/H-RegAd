rm(list = ls())
options(stringsAsFactors = FALSE)

quick_test <- FALSE       
zoom_max_iter <- 300      
tol_kkt <- 1e-4          
maxiter_main <- 2000      

source("proj_linf_ball.R")
source("soft_thresh.R")
source("sGS2w.R")

library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

n <- 200
p <- 300

beta0 <- 0.5

beta <- c(
  1, -0.8, 0.4, 0, 0,
  -0.6, 0, 0, 0, 0,
  -1.5, 0, 1.2, 0, 0, 0.3
)

beta <- c(beta, rep(0, p - length(beta)))

Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))

snr  <- 3
sfix <- 6

if (quick_test) {
  O_vec <- c(0.20)
  L_vec <- c(0)
} else {
  O_vec <- c(0, 0.10, 0.20, 0.30, 0.40)
  L_vec <- c(0, 1)
}

# ---- zero-sum constraint ----
C_raw_zero <- matrix(1, nrow = 1, ncol = p)
C_sp_zero  <- t(C_raw_zero)

# ---- subcomposition constraint ----
ngrp <- 4
C1 <- matrix(0, nrow = ngrp, ncol = 23)
tind <- c(0, 10, 16, 20, 23)

for (ii in 1:ngrp) {
  C1[ii, (tind[ii] + 1):tind[ii + 1]] <- 1
}

C_raw_sub <- matrix(0, nrow = ngrp, ncol = p)
C_raw_sub[, 1:ncol(C1)] <- C1
C_sp_sub <- t(C_raw_sub)

generate_data_case <- function(n, p, beta, beta0, Sigma,
                               O_cnt, Lfix, snr, sfix, seed = 123) {
  set.seed(seed)
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  signal <- as.numeric(beta0 + X %*% beta)
  
  sig_var <- var(signal)
  if (!is.finite(sig_var) || sig_var <= 0) {
    sig_var <- 1
  }
  
  sigma_eps <- sqrt(sig_var / snr)
  eps <- rnorm(n, mean = 0, sd = sigma_eps)
  
  y <- signal + eps
  
  out_idx <- integer(0)
  
  # response contamination
  if (O_cnt > 0) {
    out_idx <- sample(seq_len(n), O_cnt, replace = FALSE)
    y[out_idx] <- y[out_idx] + sfix * sigma_eps
  }
  
  # leverage contamination
  if (Lfix == 1 && O_cnt > 0) {
    lev_cols <- 1:min(10, p)
    shift_vec <- apply(X[, lev_cols, drop = FALSE], 2, sd)
    shift_vec[shift_vec == 0] <- 1
    
    X[out_idx, lev_cols] <- sweep(
      X[out_idx, lev_cols, drop = FALSE],
      2,
      3 * shift_vec,
      FUN = "+"
    )
  }
  
  list(
    X = X,
    y = y,
    out_idx = out_idx,
    sigma_eps = sigma_eps
  )
}

trim_kkt_history <- function(hist) {
  if (is.null(hist$eta_d)) {
    stop("eta_d is missing in history. Please make sure your source is sGS2w.R.")
  }
  
  valid_idx <- which(is.finite(hist$eta_d))
  
  if (length(valid_idx) == 0) {
    stop("There are no valid records in history$eta_d.")
  }
  
  last_k <- max(valid_idx)
  
  out <- data.frame(
    iter  = seq_len(last_k),
    eta1  = if (!is.null(hist$eta1)) hist$eta1[1:last_k] else NA_real_,
    eta2  = if (!is.null(hist$eta2)) hist$eta2[1:last_k] else NA_real_,
    eta3  = if (!is.null(hist$eta3)) hist$eta3[1:last_k] else NA_real_,
    eta4  = if (!is.null(hist$eta4)) hist$eta4[1:last_k] else NA_real_,
    eta_d = hist$eta_d[1:last_k]
  )
  
  out$eta_d <- pmax(out$eta_d, 1e-12)
  
  out
}

tune_HRegAd_once <- function(X, y, C_sp,
                             tau_grid = c(1.0, 1.345, 1.5, 2.0),
                             nlam = 6,
                             rho = 1.0,
                             sigma = 1.0,
                             maxiter = 600,
                             tol = 1e-4) {
  
  sc_mu <- colMeans(X)
  sc_sd <- apply(X, 2, sd)
  sc_sd[sc_sd == 0] <- 1
  
  X_std <- scale(X, center = sc_mu, scale = sc_sd)
  
  lam_max <- max(abs(crossprod(X_std, y))) / nrow(X_std)
  lam_max <- max(lam_max, 1e-4)
  
  lambda_path <- exp(
    seq(log(lam_max), log(lam_max * 1e-3), length.out = nlam)
  )
  
  best <- list(
    score = Inf,
    lambda1 = NA_real_,
    tau = NA_real_
  )
  
  for (tau in tau_grid) {
    for (lam in lambda_path) {
      
      fit <- tryCatch(
        spadmm_lasso_C_with_stop(
          X = X_std,
          y = as.numeric(y),
          C = as.matrix(C_sp),
          lambda1 = lam,
          tau = tau,
          sigma = sigma,
          rho = rho,
          maxiter = maxiter,
          tol = tol,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit) || is.null(fit$history)) {
        next
      }
      
      hist_df <- tryCatch(
        trim_kkt_history(fit$history),
        error = function(e) NULL
      )
      
      if (is.null(hist_df)) {
        next
      }
      
      final_score <- tail(hist_df$eta_d, 1)
      
      if (is.finite(final_score) && final_score < best$score) {
        best$score   <- final_score
        best$lambda1 <- lam
        best$tau     <- tau
      }
    }
  }
  
  if (!is.finite(best$lambda1) || !is.finite(best$tau)) {
    stop("Parameter tuning failed: no valid (lambda1, tau) found")
  }
  
  best
}

run_HRegAd_with_eta_trace <- function(X, y, C_sp,
                                      lambda1, tau,
                                      rho = 1.0,
                                      sigma = 1.0,
                                      maxiter = 2000,
                                      tol = 1e-4,
                                      verbose = TRUE) {
  
  sc_mu <- colMeans(X)
  sc_sd <- apply(X, 2, sd)
  sc_sd[sc_sd == 0] <- 1
  
  X_std <- scale(X, center = sc_mu, scale = sc_sd)
  
  fit <- spadmm_lasso_C_with_stop(
    X = X_std,
    y = as.numeric(y),
    C = as.matrix(C_sp),
    lambda1 = lambda1,
    tau = tau,
    sigma = sigma,
    rho = rho,
    maxiter = maxiter,
    tol = tol,
    verbose = verbose
  )
  
  hist_df <- trim_kkt_history(fit$history)
  
  hist_df$iter_end <- fit$iter
  hist_df$stop_reason <- fit$stop_reason
  
  hist_df
}

trace_all <- list()
idx <- 1

for (Lfix in L_vec) {
  for (o_prop in O_vec) {
    
    O_cnt <- as.integer(round(o_prop * n))
    O_lab <- paste0(round(o_prop * 100), "%")
    
    dat_case <- generate_data_case(
      n = n,
      p = p,
      beta = beta,
      beta0 = beta0,
      Sigma = Sigma,
      O_cnt = O_cnt,
      Lfix = Lfix,
      snr = snr,
      sfix = sfix,
      seed = 1000 + 100 * Lfix + O_cnt
    )
    
    X <- dat_case$X
    y <- dat_case$y
    
    for (constraint_name in c("Zero-sum constraint", "Subcomposition constraint")) {
      
      if (constraint_name == "Zero-sum constraint") {
        C_sp_use <- C_sp_zero
      } else {
        C_sp_use <- C_sp_sub
      }
      
      tuned <- tune_HRegAd_once(
        X = X,
        y = y,
        C_sp = C_sp_use,
        tau_grid = c(1.0, 1.345, 1.5, 2.0),
        nlam = 6,
        rho = 1.0,
        sigma = 1.0,
        maxiter = 600,
        tol = tol_kkt
      )
      
      cat(sprintf(
        "\nScenario: L=%d | O=%s | Constraint=%s --> tau=%.3f, lambda1=%.3e\n",
        Lfix, O_lab, constraint_name, tuned$tau, tuned$lambda1
      ))
      
      hist_df <- run_HRegAd_with_eta_trace(
        X = X,
        y = y,
        C_sp = C_sp_use,
        lambda1 = tuned$lambda1,
        tau = tuned$tau,
        rho = 1.0,
        sigma = 1.0,
        maxiter = maxiter_main,
        tol = tol_kkt,
        verbose = TRUE
      )
      
      hist_df$O <- O_lab
      hist_df$L <- paste0("L = ", Lfix)
      hist_df$Constraint <- constraint_name
      
      trace_all[[idx]] <- hist_df
      idx <- idx + 1
    }
  }
}

trace_df <- bind_rows(trace_all)

trace_df$O <- factor(
  trace_df$O,
  levels = c("0%", "10%", "20%", "30%", "40%")
)

trace_df$L <- factor(
  trace_df$L,
  levels = c("L = 0", "L = 1")
)

trace_df$Constraint <- factor(
  trace_df$Constraint,
  levels = c("Zero-sum constraint", "Subcomposition constraint")
)

trace_df_zoom <- trace_df %>%
  filter(iter <= zoom_max_iter)

check_df <- trace_df %>%
  group_by(L, Constraint, O) %>%
  summarise(
    first_eta_d = first(eta_d),
    last_eta_d  = last(eta_d),
    min_eta_d   = min(eta_d, na.rm = TRUE),
    n_iter      = max(iter),
    .groups = "drop"
  )

print(check_df)

write.csv(
  check_df,
  file = "HRegAd_dualKKT_check_summary.csv",
  row.names = FALSE
)

my_cols <- c(
  "0%"  = "#1B9E77",
  "10%" = "#D95F02",
  "20%" = "#7570B3",
  "30%" = "#E7298A",
  "40%" = "#66A61E"
)

theme_kkt <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      strip.background = element_blank(),
    
      strip.text = element_text(face = "plain", size = 10),
      strip.text.x = element_text(face = "plain", size = 10, margin = margin(b = 4, t = 2)),
      strip.text.y = element_text(face = "plain", size = 10),
      strip.text.y.right = element_text(face = "plain", size = 10, angle = -90),
      
      legend.position = "bottom",
      legend.title = element_text(face = "plain", size = 10),
      legend.text = element_text(face = "plain", size = 9.5),
      
      axis.title = element_text(face = "plain", size = 11),
      axis.title.x = element_text(face = "plain", size = 11),
      axis.title.y = element_text(face = "plain", size = 11),
      axis.text = element_text(face = "plain", size = 9.5),
      
      plot.title = element_blank(),
      plot.margin = margin(8, 10, 8, 8)
    )
}

p_eta_d <- ggplot(trace_df, aes(x = iter, y = eta_d, color = O)) +
  geom_line(linewidth = 0.95) +
  geom_hline(
    yintercept = tol_kkt,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.6
  ) +
  facet_grid(L ~ Constraint) +
  scale_color_manual(values = my_cols) +
  scale_y_log10(labels = label_scientific(digits = 2)) +
  labs(
    x = "Iteration",
    y = expression(paste("Dual KKT residual ", eta[d])),
    color = "Outlier proportion"
  ) +
  theme_kkt()

print(p_eta_d)

p_eta_d_zoom <- ggplot(trace_df_zoom, aes(x = iter, y = eta_d, color = O)) +
  geom_line(linewidth = 0.95) +
  geom_hline(
    yintercept = tol_kkt,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.6
  ) +
  facet_grid(L ~ Constraint) +
  scale_color_manual(values = my_cols) +
  scale_y_log10(labels = label_scientific(digits = 2)) +
  labs(
    x = "Iteration",
    y = expression(paste("Dual KKT residual ", eta[d])),
    color = "Outlier proportion"
  ) +
  theme_kkt()

print(p_eta_d_zoom)

ggsave(
  filename = "HRegAd_dualKKT_eta_d.pdf",
  plot = p_eta_d,
  width = 9.2,
  height = 5.6,
  dpi = 300
)

ggsave(
  filename = "HRegAd_dualKKT_eta_d.png",
  plot = p_eta_d,
  width = 9.2,
  height = 5.6,
  dpi = 300
)

ggsave(
  filename = "HRegAd_dualKKT_eta_d_zoom.pdf",
  plot = p_eta_d_zoom,
  width = 9.2,
  height = 5.6,
  dpi = 300
)

ggsave(
  filename = "HRegAd_dualKKT_eta_d_zoom.png",
  plot = p_eta_d_zoom,
  width = 9.2,
  height = 5.6,
  dpi = 300
)

cat("\nSaved files:\n")
cat("  - HRegAd_dualKKT_eta_d.pdf / .png\n")
cat("  - HRegAd_dualKKT_eta_d_zoom.pdf / .png\n")
cat("  - HRegAd_dualKKT_check_summary.csv\n")