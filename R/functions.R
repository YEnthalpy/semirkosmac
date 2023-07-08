eres <- function(e, delta, pi, ind_km) {
  e_sub <- e[ind_km]
  ord_sub <- order(e_sub)
  km_sub <- km(e_sub[ord_sub], delta[ind_km][ord_sub], pi[ind_km][ord_sub])
  es_sub <- km_sub[[1]]
  s_sub <- km_sub[[2]]
  edif <- c(diff(es_sub), 0)
  int <- rev(cumsum(rev(edif * s_sub)))
  es_int <- int + s_sub * es_sub
  int <- int / s_sub + es_sub
  ehat <- approx(es_sub, int, e, method = "constant", ties = "ordered")$y
  ehat[is.na(ehat)] <- e[is.na(ehat)]
  return(list(ehat, es_int))
}

semi_rk_est <- function(x, y, delta, pi, n,
                        control = list(
                          init = lsfit(x, y, intercept = FALSE)$coefficient,
                          tol = 1e-5, maxit = 1000
                        )) {
  out <- nleqslv::nleqslv(
    x = control$init, fn = function(b) {
      colSums(gehan_smth(x, y, delta, pi, b, n))
    }, jac = function(b) {
      gehan_s_jaco(x, y, delta, pi, b, n)
    }, method = "Broyden", jacobian = FALSE,
    control = list(ftol = control$tol, xtol = 1e-20, maxit = control$maxit)
  )
  conv <- out$termcd
  coe <- out$x
  if (conv == 1) {
    conv <- 0
  } else if (conv %in% c(2, 4)) {
    conv <- 2
    coe <- rep(NA, ncol(x))
  } else {
    conv <- 1
    coe <- rep(NA, ncol(x))
  }
  names(coe) <- paste0("beta", seq_len(ncol(x)))
  return(list(
    coefficient = coe, converge = conv,
    iter = c(out$iter, out$njcnt, out$nfcnt)
  ))
}



## 5. get the optimal ssps
semi_rk_ssp <- function(x, y, delta, r0, ssp_type, alpha) {
  n <- nrow(x)
  if (ssp_type == "uniform") {
    ssp <- rep(1 / n, n)
    return(list(ssp = ssp, index.pilot = sample(n, r0, TRUE), 
                converge = 0))
  } else {
    ssp_pt <- rep(1 / n, r0)
    ind_pt <- sample(n, r0, TRUE) # index of pilot sample
    x_pt <- x[ind_pt, ] 
    y_pt <- y[ind_pt]
    dt_pt <- delta[ind_pt]
    mle_pt <- semi_rk_est(x_pt, y_pt, dt_pt, ssp_pt, n) # pilot estimator
    if (mle_pt$converge %in% c(1, 2)) {
      return(list(ssp = NA, index.pilot = NA, converge = mle_pt$converge))
    }
    bt_pt <- mle_pt$coefficient # pilot estimator
    # estimating function of all observations estimated by the pilot sample
    g <- gehan_s_mtg(x, y, delta, rep(1 / n, n), bt_pt, ind_pt - 1, n)
    # hessian matrix estimated by the pilot sample
    m <- gehan_s_jaco(x_pt, y_pt, dt_pt, ssp_pt, bt_pt, n)
    if (ssp_type == "optL") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - alpha) + alpha / n # optL SSPs
    } else if (ssp_type == "optA") {
      m_inv <- solve(m)
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - alpha) + alpha / n # optA SSPs
    }
    return(list(ssp = ssp, est.pilot = bt_pt, 
                index.pilot = ind_pt, M = m, converge = 0))
  }
}


## 6. Get estiamted coefficient and standard error
semi_rk_fit <- function(x, y, delta, r0, r, ssp_type, se = TRUE, alpha = 0.2) {
  n <- nrow(x)
  intercept <- (sum(x[, 1]) == n)
  # get rid of the intercept
  if (intercept) {
    x <- x[, -1]
  }
  # get optimal SSPs, the pilot sample, the pilot estimator
  # and the M estimated by the pilot sample
  t_ssp <- system.time(ssps <- semi_rk_ssp(x, y, delta, r0, ssp_type, alpha))
  if (ssps$converge != 0) {
    stop(paste0("Fail to get a converging pilot
    estimator. The converging code is ", ssps$converge))
  }
  pi <- ssps$ssp # optimal SSPs
  ind_r <- sample(n, r, prob = pi, replace = TRUE) # index of second-step subsample
  if (ssp_type == "uniform") {
    ind_st <- c(ind_r, ssps$index.pilot)
    # using uniform SSPs as the control group
    # combine pilot and second-step subsample to derive the final estimator
    t_est <- system.time(est <- semi_rk_est(x[ind_st, ], y[ind_st],
                                            delta[ind_st], rep(1/n, (r+r0)), n))
    coe_out <- as.vector(est$coefficient)
  }else {
    # second-step estimator
    t_est <- system.time(est <- semi_rk_est(x[ind_r, ], y[ind_r], 
                                            delta[ind_r], pi[ind_r], n))
    coe_sec <- as.vector(est$coefficient)
    # hessian matrix estimated by the second step subsample
    m_sec <- r * gehan_s_jaco(x[ind_r, ], y[ind_r], delta[ind_r],
                          pi[ind_r], coe_sec, n)
    m_pt <- r0 * ssps$M
    # aggregate the pilot and second-step estimator
    coe_out <- drop(solve(m_pt + m_sec) %*% (m_pt %*% ssps$est.pilot + m_sec %*% coe_sec))
  }
  iter <- est$iter
  if (est$converge %in% c(1, 2)) {
    stop(paste0("Fail to get a converging second-step
    estimator. The converging code is ", est$converge))
  }
  t_se <- system.time(
  {if (se) {
    ind_st <- c(ind_r, ssps$index.pilot) # combined subsample index
    pi_st <- c(pi[ind_r], rep(1/n, r0)) # combined SSPs
    # estimating function estimated by the combined subsample
    g <- gehan_s_mtg(
      x[ind_st, ], y[ind_st], delta[ind_st], pi_st,
      coe_out, seq_along(ind_st) - 1, n
    )
    # estimate vc
    vc <- crossprod(g) * r
    vc_add <- crossprod(pi[ind_st] * g, g) * r * n
    vc_amend <- vc + (r / n) * vc_add
    # inverse of the hessian matrix estimated by the second step subsample
    m_inv <- solve(gehan_s_jaco(
      x[ind_st, ], y[ind_st], delta[ind_st], pi_st, coe_out, n))
    vx <- m_inv %*% vc %*% m_inv / r # sandwich estimator for full estimator
    std <- c(sqrt(diag(vx)))
    vx_amend <- m_inv %*% vc_amend %*% m_inv / r # sandwich estimator for true coe
    std_amend <- c(sqrt(diag(vx_amend)))
  } else {
    std <- std_amend <- NA
  }})
  if (is.null(colnames(x))) {
    names(coe_out) <- paste0("Beta", seq_len(ncol(x)))
    names(std) <- names(std_amend) <- paste0("Beta", seq_len(ncol(x)))
  } else {
    names(coe_out) <- colnames(x)
    names(std) <- names(std_amend) <- colnames(x)
  }
  # add an intercept term if elements in the first column in x are 1's
  if (intercept) {
    coe_out <- c(max(eres((y[ind_st] - x[ind_st, ] %*% coe_out), 
                          delta[ind_st], pi_st, (r+r0))[[2]]), coe_out)
    names(coe_out) <- c("intercept", names(coe_out))
  }
  time <- cbind(t_ssp, t_est, t_se)
  colnames(time) <- c("SSPs", "Est", "SE")
  return(list(coe = coe_out, std = std, std_amend = std_amend, iter = iter, 
              converge = 0, time = time))
}
