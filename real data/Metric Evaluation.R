###fitted model
###Strain selection
###Stability selection

rm(list = ls())
options(stringsAsFactors = FALSE)

load("hiv_fit_objects.rda")

library(ggplot2)
library(dplyr)
library(scales)

tol <- 1e-8

plot_df <- data.frame(
  y = as.numeric(y),
  yhat_NR = as.numeric(fitted_mat[, "NR"]),
  yhat_HADMM = as.numeric(fitted_mat[, "H-ADMM"]),
  gamma_HADMM = as.numeric(shift_est_lam_cv[, "H-ADMM"])
)

plot_df <- plot_df %>%
  mutate(
    outlier_HADMM = abs(gamma_HADMM) > tol,
    point_group = factor(
      ifelse(outlier_HADMM, "Flagged outlier", "Regular sample"),
      levels = c("Regular sample", "Flagged outlier")
    )
  )

fit_hadmm_line <- lm(yhat_HADMM ~ y, data = plot_df %>% filter(!outlier_HADMM))
fit_nr_line <- lm(yhat_NR ~ y, data = plot_df)

xgrid <- seq(min(plot_df$y), max(plot_df$y), length.out = 300)

line_df <- data.frame(
  y = xgrid,
  yhat_HADMM = predict(fit_hadmm_line, newdata = data.frame(y = xgrid)),
  yhat_NR = predict(fit_nr_line, newdata = data.frame(y = xgrid))
)

if (exists("summary_table")) {
  r2_hadmm <- summary_table$R2[summary_table$Method == "H-ADMM"]
  r2_nr    <- summary_table$R2[summary_table$Method == "NR"]
} else {
  r2_hadmm <- 1 - sum((plot_df$y - plot_df$yhat_HADMM)^2) / sum((plot_df$y - mean(plot_df$y))^2)
  r2_nr    <- 1 - sum((plot_df$y - plot_df$yhat_NR)^2) / sum((plot_df$y - mean(plot_df$y))^2)
}

lab_hadmm <- paste0("R²(H-RegAd) = ", sprintf("%.2f", r2_hadmm))
lab_nr    <- paste0("R²(RobRegCC-NR) = ", sprintf("%.2f", r2_nr))

xy_min <- min(c(plot_df$y, plot_df$yhat_HADMM, plot_df$yhat_NR), na.rm = TRUE)
xy_max <- max(c(plot_df$y, plot_df$yhat_HADMM, plot_df$yhat_NR), na.rm = TRUE)

pad <- 0.04 * (xy_max - xy_min)
xy_limits <- c(xy_min - pad, xy_max + pad)

x_text <- xy_limits[1] + 0.06 * diff(xy_limits)
y_text1 <- xy_limits[2] - 0.08 * diff(xy_limits)
y_text2 <- xy_limits[2] - 0.16 * diff(xy_limits)

p_fit <- ggplot(plot_df, aes(x = y, y = yhat_HADMM)) +

  geom_abline(
    intercept = 0, slope = 1,
    color = "grey70", linewidth = 0.8, linetype = "dotted"
  ) +
  geom_point(
    aes(shape = point_group, color = point_group),
    size = 2.5, alpha = 0.92
  ) +
  geom_line(
    data = line_df,
    aes(x = y, y = yhat_HADMM),
    color = "#C0392B", linewidth = 1.2
  ) +
  geom_line(
    data = line_df,
    aes(x = y, y = yhat_NR),
    color = "#2C7FB8", linewidth = 1.1, linetype = "dashed"
  ) +
  annotate(
    "text",
    x = x_text, y = y_text1,
    label = lab_hadmm,
    hjust = 0, size = 4.7
  ) +
  annotate(
    "text",
    x = x_text, y = y_text2,
    label = lab_nr,
    hjust = 0, size = 4.7
  ) +
  scale_color_manual(
    values = c(
      "Regular sample" = "black",
      "Flagged outlier" = "#C0392B"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Regular sample" = 16,
      "Flagged outlier" = 17
    )
  ) +
  coord_cartesian(xlim = xy_limits, ylim = xy_limits) +
  scale_x_continuous(labels = label_number(big.mark = ",")) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(
    x = "Observed sCD14",
    y = "Fitted sCD14",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

print(p_fit)

ggsave(
  filename = "hiv_fit_HADMM_vs_NR_refined2.pdf",
  plot = p_fit,
  width = 6.3,
  height = 5.4,
  dpi = 300
)

ggsave(
  filename = "hiv_fit_HADMM_vs_NR_refined2.png",
  plot = p_fit,
  width = 6.3,
  height = 5.4,
  dpi = 300
)

cat("\nSaved files:\n")
cat("  - hiv_fit_HADMM_vs_NR_refined2.pdf\n")
cat("  - hiv_fit_HADMM_vs_NR_refined2.png\n")









rm(list = ls())
options(stringsAsFactors = FALSE)
load("hiv_fit_objects.rda")

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

if (!exists("beta_lam_cv")) {
  stop("beta_lam_cv not found, please run hiv_realdata_metrics.R first")
}

if (exists("X") && !is.null(colnames(X))) {
  genus_names <- colnames(X)
} else if (exists("predictor_est_c1")) {
  genus_names <- predictor_est_c1$v.name[predictor_est_c1$v.name != "intercept"]
} else {
  genus_names <- paste0("Genus_", seq_len(nrow(beta_lam_cv) - 1))
}

beta_no_intercept <- beta_lam_cv[-1, , drop = FALSE]

if (length(genus_names) != nrow(beta_no_intercept)) {
  stop("The length of genus_names doesn’t match the number of r
       ows in beta_no_intercept, please check the object")
}

rownames(beta_no_intercept) <- genus_names

tol <- 1e-8

select_mat <- abs(beta_no_intercept) > tol

coef_hadmm <- beta_no_intercept[, "H-ADMM"]

selected_hadmm_tbl <- data.frame(
  Genus = genus_names,
  Coef_HADMM = as.numeric(coef_hadmm),
  AbsCoef = abs(as.numeric(coef_hadmm)),
  Direction = ifelse(coef_hadmm > 0, "Positive", "Negative"),
  Selected_NR = select_mat[, "NR"],
  Selected_A = select_mat[, "A"],
  Selected_E = select_mat[, "E"],
  Selected_H = select_mat[, "H"],
  Selected_HADMM = select_mat[, "H-ADMM"]
) %>%
  filter(Selected_HADMM) %>%
  mutate(
    Shared_Count = Selected_NR + Selected_A + Selected_E + Selected_H,
    Shared_Pattern = paste0(
      ifelse(Selected_NR, "NR;", ""),
      ifelse(Selected_A, "A;", ""),
      ifelse(Selected_E, "E;", ""),
      ifelse(Selected_H, "H;", "")
    ),
    Shared_Pattern = ifelse(Shared_Pattern == "", "H-ADMM only", Shared_Pattern)
  ) %>%
  arrange(desc(AbsCoef))

print(selected_hadmm_tbl)

comparison_tbl <- data.frame(
  Genus = genus_names,
  Coef_NR = beta_no_intercept[, "NR"],
  Coef_A  = beta_no_intercept[, "A"],
  Coef_E  = beta_no_intercept[, "E"],
  Coef_H  = beta_no_intercept[, "H"],
  Coef_HADMM = beta_no_intercept[, "H-ADMM"],
  Selected_NR = select_mat[, "NR"],
  Selected_A  = select_mat[, "A"],
  Selected_E  = select_mat[, "E"],
  Selected_H  = select_mat[, "H"],
  Selected_HADMM = select_mat[, "H-ADMM"]
) %>%
  mutate(
    N_Selected_Methods = Selected_NR + Selected_A + Selected_E + Selected_H + Selected_HADMM
  ) %>%
  arrange(desc(abs(Coef_HADMM)))

overlap_summary <- data.frame(
  Quantity = c(
    "Selected by NR",
    "Selected by A",
    "Selected by E",
    "Selected by H",
    "Selected by H-ADMM",
    "Selected by both H-ADMM and NR",
    "Selected by both H-ADMM and A",
    "Selected by both H-ADMM and E",
    "Selected by both H-ADMM and H",
    "Selected only by H-ADMM"
  ),
  Value = c(
    sum(select_mat[, "NR"]),
    sum(select_mat[, "A"]),
    sum(select_mat[, "E"]),
    sum(select_mat[, "H"]),
    sum(select_mat[, "H-ADMM"]),
    sum(select_mat[, "H-ADMM"] & select_mat[, "NR"]),
    sum(select_mat[, "H-ADMM"] & select_mat[, "A"]),
    sum(select_mat[, "H-ADMM"] & select_mat[, "E"]),
    sum(select_mat[, "H-ADMM"] & select_mat[, "H"]),
    sum(select_mat[, "H-ADMM"] &
          !select_mat[, "NR"] &
          !select_mat[, "A"] &
          !select_mat[, "E"] &
          !select_mat[, "H"])
  )
)

print(overlap_summary)

plot_tbl <- selected_hadmm_tbl %>%
  mutate(
    Genus = factor(Genus, levels = rev(Genus)),
    Direction = factor(Direction, levels = c("Positive", "Negative"))
  )

p_genus <- ggplot(plot_tbl, aes(x = Genus, y = Coef_HADMM, fill = Direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(
    values = c("Positive" = "#C0392B", "Negative" = "#2C7FB8")
  ) +
  labs(
    x = NULL,
    y = "Coefficient",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_genus)

write_csv(selected_hadmm_tbl, "hiv_selected_genera_HADMM.csv")
write_csv(comparison_tbl, "hiv_selected_genera_comparison.csv")
write_csv(overlap_summary, "hiv_selected_genera_overlap_summary.csv")

ggsave(
  filename = "hiv_selected_genera_HADMM.pdf",
  plot = p_genus,
  width = 7.2,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "hiv_selected_genera_HADMM.png",
  plot = p_genus,
  width = 7.2,
  height = 5.8,
  dpi = 300
)

cat("\nSaved files:\n")
cat("  - hiv_selected_genera_HADMM.csv\n")
cat("  - hiv_selected_genera_comparison.csv\n")
cat("  - hiv_selected_genera_overlap_summary.csv\n")
cat("  - hiv_selected_genera_HADMM.pdf\n")
cat("  - hiv_selected_genera_HADMM.png\n")



rm(list = ls())
options(stringsAsFactors = FALSE)

load("HIV_data/sCD14wC.rda")

source("proj_linf_ball.R")
source("soft_thresh.R")
source("sGS2uzhenshishuju.R")

library(dplyr)
library(ggplot2)
library(readr)
library(scales)
exp_seed <- 123

B <- 50
subsample_ratio <- 0.80

K_inner <- 5
nlam_sp <- 40
tau_sp <- 1.345
sigma_sp <- 1.0
rho_sp <- 1.0

tol <- 1e-8

lambda_use_manual <- NA_real_

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

C <- matrix(1, nrow = 1, ncol = 60)

W <- sCD14[, 1:60]
X <- t(apply(W, 1, function(x) {
  indx <- x != 0
  x[!indx] <- 0.5
  x / sum(x)
}))

y <- sCD14[, 61]

Xt <- cbind(1, log(X))
Ct <- cbind(matrix(0, nrow = nrow(C), ncol = 1), C)   
C_sp <- t(Ct)                                         

n <- nrow(Xt)
p <- ncol(Xt)
bw <- c(0, rep(1, p - 1))  

if (is.null(colnames(X))) {
  genus_names <- paste0("Genus_", 1:ncol(X))
} else {
  genus_names <- colnames(X)
}

tune_hadmm_lambda <- function(X_train, y_train, C_sp, bw,
                              K_inner = 5, nlam_sp = 40,
                              tau_sp = 1.345, sigma_sp = 1.0, rho_sp = 1.0,
                              seed_inner = 123) {
  n_tr <- nrow(X_train)
  
  lam_max <- max(abs(as.numeric(crossprod(X_train, y_train)))) / n_tr
  lam_min <- lam_max * 1e-8
  lambda_seq <- exp(seq(log(lam_max), log(lam_min), length.out = nlam_sp))
  
  set.seed(seed_inner)
  fold_id <- sample(rep(1:K_inner, length.out = n_tr))
  
  cv_score <- rep(NA_real_, nlam_sp)
  
  for (ii in seq_along(lambda_seq)) {
    lam <- lambda_seq[ii]
    score_k <- rep(NA_real_, K_inner)
    
    for (kk in 1:K_inner) {
      idx_val <- which(fold_id == kk)
      idx_fit <- setdiff(seq_len(n_tr), idx_val)
      
      X_fit <- X_train[idx_fit, , drop = FALSE]
      y_fit <- y_train[idx_fit]
      X_val <- X_train[idx_val, , drop = FALSE]
      y_val <- y_train[idx_val]
      
      fit_tmp <- tryCatch(
        spadmm_lasso_C_with_stop(
          X = X_fit, y = y_fit, C = C_sp,
          lambda1 = lam,
          tau = tau_sp, sigma = sigma_sp, rho = rho_sp,
          bw = bw,
          maxiter = 3000, tol = 1e-6,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_tmp)) next
      
      yhat_val <- as.numeric(X_val %*% fit_tmp$beta)
      res_val <- y_val - yhat_val
      
      sigma_hat <- robust_sigma_hat(y_fit)
      c_thr <- 1.345 * sigma_hat
      
      score_k[kk] <- huber_loss_mean_c(res_val, c = c_thr) / (c_thr^2)
    }
    
    cv_score[ii] <- mean(score_k, na.rm = TRUE)
    
    if (ii %% 10 == 0) {
      cat("[full-data tuning]", ii, "/", nlam_sp, "\n")
    }
  }
  
  lambda_seq[which.min(cv_score)]
}

if (is.na(lambda_use_manual)) {
  set.seed(exp_seed)
  lambda_use <- tune_hadmm_lambda(
    X_train = Xt, y_train = y, C_sp = C_sp, bw = bw,
    K_inner = K_inner, nlam_sp = nlam_sp,
    tau_sp = tau_sp, sigma_sp = sigma_sp, rho_sp = rho_sp,
    seed_inner = exp_seed
  )
} else {
  lambda_use <- lambda_use_manual
}

cat("Selected lambda for stability analysis =", lambda_use, "\n")

coef_mat <- matrix(0, nrow = p - 1, ncol = B)
rownames(coef_mat) <- genus_names
colnames(coef_mat) <- paste0("Rep_", 1:B)

for (b in 1:B) {
  set.seed(exp_seed + b)
  
  idx_sub <- sort(sample(seq_len(n), size = floor(subsample_ratio * n), replace = FALSE))
  X_sub <- Xt[idx_sub, , drop = FALSE]
  y_sub <- y[idx_sub]
  
  fit_sub <- tryCatch(
    spadmm_lasso_C_with_stop(
      X = X_sub, y = y_sub, C = C_sp,
      lambda1 = lambda_use,
      tau = tau_sp, sigma = sigma_sp, rho = rho_sp,
      bw = bw,
      maxiter = 5000, tol = 1e-8,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit_sub)) {
    cat("Subsample", b, "failed.\n")
    next
  }
  
  beta_sub <- as.numeric(fit_sub$beta)
  
  coef_mat[, b] <- beta_sub[-1]
  
  if (b %% 10 == 0) {
    cat("[stability subsampling]", b, "/", B, "\n")
  }
}

selected_mat <- abs(coef_mat) > tol

selection_freq <- rowMeans(selected_mat, na.rm = TRUE)

mean_coef_selected <- sapply(1:nrow(coef_mat), function(i) {
  vals <- coef_mat[i, selected_mat[i, ], drop = TRUE]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
})

median_coef_selected <- sapply(1:nrow(coef_mat), function(i) {
  vals <- coef_mat[i, selected_mat[i, ], drop = TRUE]
  if (length(vals) == 0) return(NA_real_)
  median(vals)
})

sign_consistency <- sapply(1:nrow(coef_mat), function(i) {
  vals <- coef_mat[i, selected_mat[i, ], drop = TRUE]
  if (length(vals) == 0) return(NA_real_)
  abs(mean(sign(vals)))
})

dominant_direction <- ifelse(mean_coef_selected > 0, "Positive", "Negative")

stability_tbl <- data.frame(
  Genus = genus_names,
  Selection_Frequency = selection_freq,
  Mean_Coef_When_Selected = mean_coef_selected,
  Median_Coef_When_Selected = median_coef_selected,
  Sign_Consistency = sign_consistency,
  Dominant_Direction = dominant_direction
) %>%
  arrange(desc(Selection_Frequency), desc(abs(Mean_Coef_When_Selected)))

top_tbl <- stability_tbl %>%
  filter(Selection_Frequency > 0) %>%
  slice(1:12)

print(top_tbl)

plot_tbl <- top_tbl %>%
  mutate(
    Genus = factor(Genus, levels = rev(Genus)),
    Dominant_Direction = factor(Dominant_Direction, levels = c("Positive", "Negative"))
  )

p_stab <- ggplot(plot_tbl, aes(x = Genus, y = Selection_Frequency, fill = Dominant_Direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = label_number(accuracy = 0.1)
  ) +
  scale_fill_manual(
    values = c("Positive" = "#C0392B", "Negative" = "#2C7FB8")
  ) +
  labs(
    x = NULL,
    y = "Selection frequency",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_stab)

write_csv(stability_tbl, "hiv_selected_genera_stability.csv")
write_csv(top_tbl, "hiv_selected_genera_stability_top.csv")

ggsave(
  filename = "hiv_selected_genera_stability.pdf",
  plot = p_stab,
  width = 7.0,
  height = 5.8,
  dpi = 300
)

ggsave(
  filename = "hiv_selected_genera_stability.png",
  plot = p_stab,
  width = 7.0,
  height = 5.8,
  dpi = 300
)

cat("\nSaved files:\n")
cat("  - hiv_selected_genera_stability.csv\n")
cat("  - hiv_selected_genera_stability_top.csv\n")
cat("  - hiv_selected_genera_stability.pdf\n")
cat("  - hiv_selected_genera_stability.png\n")