vec_norm <- function(x) {
  sqrt(sum(x^2))
}

safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(Inf)
  num / den
}

spadmm_lasso_C_dual_kkt <- function(X, y, C,
                                    lambda1 = 1e-3,
                                    tau = 1.345,
                                    sigma = 1.0,
                                    rho = 1.0,
                                    eta = NULL,
                                    gamma = 1e-4,
                                    maxiter = 1000,
                                    tol = 1e-4,
                                    verbose = TRUE,
                                    init = list()) {
  
  ## Input
  if (is.data.frame(X)) X <- as.matrix(X)
  if (is.data.frame(C)) C <- as.matrix(C)
  
  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  storage.mode(C) <- "double"
  
  if (!is.matrix(X)) stop("X must be a numeric matrix.")
  if (!is.numeric(y)) stop("y must be a numeric vector.")
  if (!is.matrix(C)) stop("C must be a numeric matrix.")
  
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (length(y) != n) stop("length(y) must be equal to nrow(X).")
  if (nrow(C) != p) stop("nrow(C) must be equal to ncol(X).")
  
  r <- ncol(C)
  
  if (r == 0) {
    C <- matrix(0, nrow = p, ncol = 0)
  }
  
  if (!is.finite(lambda1) || lambda1 < 0) stop("lambda1 must be nonnegative and finite.")
  if (!is.finite(tau) || tau <= 0) stop("tau must be positive and finite.")
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive and finite.")
  if (!is.finite(rho) || rho <= 0) stop("rho must be positive and finite.")
  if (!is.finite(gamma) || gamma <= 0) stop("gamma must be positive and finite.")
  if (!is.finite(maxiter) || maxiter < 1) stop("maxiter must be a positive integer.")
  if (!is.finite(tol) || tol <= 0) stop("tol must be positive and finite.")
  
  maxiter <- as.integer(maxiter)
  
  ## eta for S = (eta - n)I_n - sigma X X^T >= 0

  if (is.null(eta)) {
    svals <- svd(X, nu = 0, nv = 0)$d
    max_sv2 <- ifelse(length(svals) > 0, max(svals)^2, 0)
    eta <- n + sigma * max_sv2 + 1
  }
  
  if (!is.finite(eta) || eta <= 0) stop("eta must be positive and finite.")
  

  ## Initialization
  u  <- if (!is.null(init$u))  as.numeric(init$u)  else rep(0, n)
  v  <- if (!is.null(init$v))  as.numeric(init$v)  else rep(0, p)
  q  <- if (!is.null(init$q))  as.numeric(init$q)  else rep(0, r)
  xi <- if (!is.null(init$xi)) as.numeric(init$xi) else rep(0, p)
  
  if (length(u)  != n) stop("init$u has wrong length.")
  if (length(v)  != p) stop("init$v has wrong length.")
  if (length(q)  != r) stop("init$q has wrong length.")
  if (length(xi) != p) stop("init$xi has wrong length.")
  
  CtC <- if (r > 0) crossprod(C) else matrix(0, 0, 0)
  
  ## S u = (eta - n)u - sigma X X^T u
  Su_op <- function(uc) {
    as.numeric((eta - n) * uc - sigma * (X %*% (crossprod(X, uc))))
  }
  
  ## Safe linear solver
  solve_q_system <- function(A, b) {
    out <- tryCatch(
      as.numeric(solve(A, b)),
      error = function(e) as.numeric(qr.solve(A, b))
    )
    out
  }
  
  ## Dual KKT residual
  dual_kkt_residual <- function(u, v, q, xi) {
    Xt_u <- as.numeric(crossprod(X, u))
    Cq   <- if (r > 0) as.numeric(C %*% q) else rep(0, p)
    
    ## eta1:
    proj_u <- proj_linf_ball(
      (y - as.numeric(X %*% xi)) / n,
      tau / n
    )
    eta1 <- safe_ratio(
      vec_norm(u - proj_u),
      1 + vec_norm(u)
    )
    
    ## eta2:
    ## v = Proj_{B_inf^{lambda}}( v + xi )
    proj_v <- proj_linf_ball(v + xi, lambda1)
    eta2 <- safe_ratio(
      vec_norm(v - proj_v),
      1 + vec_norm(v)
    )
    
    Ct_xi <- if (r > 0) as.numeric(crossprod(C, xi)) else numeric(0)
    eta3 <- if (r > 0) {
      safe_ratio(vec_norm(Ct_xi), 1 + vec_norm(xi))
    } else {
      0
    }
    
    ## eta4:
    ## v - X^T u - C q = 0
    feas <- v - Xt_u - Cq
    eta4 <- safe_ratio(
      vec_norm(feas),
      1 + vec_norm(v) + vec_norm(Xt_u) + vec_norm(Cq)
    )
    
    eta_vec <- c(eta1, eta2, eta3, eta4)
    eta_d <- if (all(!is.finite(eta_vec))) {
      Inf
    } else {
      max(eta_vec, na.rm = TRUE)
    }
    
    list(
      eta1 = eta1,
      eta2 = eta2,
      eta3 = eta3,
      eta4 = eta4,
      eta_d = eta_d,
      feas = feas
    )
  }
  
  ## History
  history <- list(
    eta1    = rep(NA_real_, maxiter),
    eta2    = rep(NA_real_, maxiter),
    eta3    = rep(NA_real_, maxiter),
    eta4    = rep(NA_real_, maxiter),
    eta_d   = rep(NA_real_, maxiter),
    eta_min = rep(NA_real_, maxiter)
  )
  
  eta_best <- Inf
  stop_reason <- "maxiter"
  k_end <- 0
  
  for (k in 1:maxiter) {
  
    u_k  <- u
    v_k  <- v
    q_k  <- q
    xi_k <- xi
    
    Xt_u_k <- as.numeric(crossprod(X, u_k))
    
    if (r > 0) {
      A_q <- sigma * CtC + gamma * diag(r)
      
      rhs_q_half <- sigma * as.numeric(crossprod(C, v_k - Xt_u_k)) -
        as.numeric(crossprod(C, xi_k)) +
        gamma * q_k
      
      q_half <- solve_q_system(A_q, rhs_q_half)
    } else {
      q_half <- numeric(0)
    }
    
    Cq_half <- if (r > 0) as.numeric(C %*% q_half) else rep(0, p)
    
    a_half <- Xt_u_k + Cq_half
    v_new <- proj_linf_ball(a_half + xi_k / sigma, lambda1)
    
    if (r > 0) {
      rhs_q_new <- sigma * as.numeric(crossprod(C, v_new - Xt_u_k)) -
        as.numeric(crossprod(C, xi_k)) +
        gamma * q_k
      
      q_new <- solve_q_system(A_q, rhs_q_new)
    } else {
      q_new <- numeric(0)
    }
    
    Cq_new <- if (r > 0) as.numeric(C %*% q_new) else rep(0, p)
    
    Su_k   <- Su_op(u_k)
    X_xi_k <- as.numeric(X %*% xi_k)
    X_vnew <- as.numeric(X %*% v_new)
    X_Cq   <- if (r > 0) as.numeric(X %*% Cq_new) else rep(0, n)
    
    rhs_u <- y - X_xi_k +
      sigma * X_vnew +
      Su_k -
      sigma * X_Cq
    
    u_new <- proj_linf_ball(rhs_u / eta, tau / n)
    
    Xt_u_new <- as.numeric(crossprod(X, u_new))
    resid <- v_new - Xt_u_new - Cq_new
    xi_new <- xi_k - sigma * rho * resid
  
    if (any(!is.finite(c(u_new, v_new, q_new, xi_new)))) {
      stop_reason <- paste0("nonfinite_iter_at_", k)
      if (verbose) {
        cat(sprintf("Non-finite iterate encountered at iter %d. Algorithm stopped.\n", k))
      }
      k_end <- k
      break
    }
    
    ## Update variables
    u  <- u_new
    v  <- v_new
    q  <- q_new
    xi <- xi_new
    
    kkt <- dual_kkt_residual(u, v, q, xi)
    
    history$eta1[k]  <- kkt$eta1
    history$eta2[k]  <- kkt$eta2
    history$eta3[k]  <- kkt$eta3
    history$eta4[k]  <- kkt$eta4
    history$eta_d[k] <- kkt$eta_d
    
    eta_best <- min(eta_best, kkt$eta_d, na.rm = TRUE)
    history$eta_min[k] <- eta_best
    
    if (verbose && (k == 1 || k %% 50 == 0)) {
      cat(sprintf(
        "iter %4d: eta_d=%.3e, eta1=%.3e, eta2=%.3e, eta3=%.3e, eta4=%.3e, eta_min=%.3e\n",
        k,
        kkt$eta_d,
        kkt$eta1,
        kkt$eta2,
        kkt$eta3,
        kkt$eta4,
        eta_best
      ))
    }
    
    if (!is.finite(kkt$eta_d)) {
      stop_reason <- paste0("nonfinite_kkt_at_", k)
      if (verbose) {
        cat(sprintf("Non-finite KKT residual encountered at iter %d. Algorithm stopped.\n", k))
      }
      k_end <- k
      break
    }
    
    if (kkt$eta_d < tol) {
      stop_reason <- paste0("kkt_tol_at_", k)
      if (verbose) {
        cat("Converged by dual KKT residual at iter", k, "\n")
      }
      k_end <- k
      break
    }
    
    k_end <- k
  }
  
  if (k_end == 0) k_end <- maxiter
  
  Xt_u_final <- as.numeric(crossprod(X, u))
  Cq_final <- if (r > 0) as.numeric(C %*% q) else rep(0, p)
  
  beta_hat <- Xt_u_final + Cq_final + xi / sigma - v
  
  yhat <- as.numeric(X %*% beta_hat)
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)
  R2 <- 1 - ss_res / ss_tot
  
  history$eta1    <- history$eta1[1:k_end]
  history$eta2    <- history$eta2[1:k_end]
  history$eta3    <- history$eta3[1:k_end]
  history$eta4    <- history$eta4[1:k_end]
  history$eta_d   <- history$eta_d[1:k_end]
  history$eta_min <- history$eta_min[1:k_end]
  
  final_kkt <- dual_kkt_residual(u, v, q, xi)
  
  list(
    beta        = beta_hat,
    u           = u,
    v           = v,
    q           = q,
    xi          = xi,
    R2          = R2,
    history     = history,
    kkt         = final_kkt,
    iter        = k_end,
    stop_reason = stop_reason
  )
}

spadmm_lasso_C_with_stop <- spadmm_lasso_C_dual_kkt
