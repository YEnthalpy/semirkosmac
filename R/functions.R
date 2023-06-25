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
                        control = list(init = lsfit(x[, -1], y, intercept = FALSE)$coefficient)) 
{
  out <- nleqslv::nleqslv(x = control$init, fn = function(b) {
    colSums(gehan_smth(x[, -1], y, delta, pi, b, n))
  }, jac = function(b) {
    gehan_s_jaco(x[, -1], y, delta, pi, b, n)
  }, method = "Broyden", jacobian = FALSE, 
  control = list(ftol = 1e-5, xtol = 1e-20, maxit = 1000))
  conv <- out$termcd
  coe <- out$x
  if (conv == 1) {
    conv <- 0
    coe <- c(max(eres((y - x[, -1] %*% coe), delta, pi, seq_along(y))[[2]]), coe)
  }
  else if (conv %in% c(2, 4)) {
    conv <- 2
    coe <- rep(NA, ncol(x))
  }
  else {
    conv <- 1
    coe <- rep(NA, ncol(x))
  }
  names(coe) <- c("intercept", paste0("beta", seq_len(ncol(x)-1)))
  return(list(coefficient = coe, converge = conv, 
              iter = c(out$iter, out$njcnt, out$nfcnt)))
}



## 5. get the optimal ssps
semi_rk_ssp <- function(x, y, delta, r0, ssp_type, alpha) {
  n <- nrow(x)
  if (ssp_type == "uniform") {
    ssp <- rep(1 / n, n)
    return(list(ssp = ssp, index.pilot = sample(n, r0, TRUE)))
  } else {
    ssp_pt <- rep(1 / n, r0)
    ind_pt <- sample(n, r0, TRUE)
    x_pt <- x[ind_pt, ]
    y_pt <- y[ind_pt]
    dt_pt <- delta[ind_pt]
    mle_pt <- semi_rk_est(x_pt, y_pt, dt_pt, ssp_pt, n)
    if (mle_pt$converge %in% c(1, 2)) {
      return(list(ssp = NA, index.pilot = NA, 
                  converge = c(3, mle_pt$converge)))
    }
    bt_pt <- mle_pt$coefficient
    e_pt <- y - x %*% bt_pt
    g <- gehan_s_mtg(x, y, delta, rep(1/n, n), bt_pt, ind_pt-1, n)
    g <- g[, -1]
    if (ssp_type == "mVc") {
      g_nm <- sqrt(rowSums(g^2))
      ssp <- g_nm / sum(g_nm) * (1 - alpha) + alpha / n
    } else if (ssp_type == "mMSE") {
      m_inv <- solve(gehan_s_jaco(x_pt[, -1], y_pt, dt_pt, ssp_pt, bt_pt[-1], n))
      m_mse <- sqrt(colSums((tcrossprod(m_inv, g))^2))
      ssp <- m_mse / sum(m_mse) * (1 - alpha) + alpha / n
    }
    return(list(ssp = ssp, index.pilot = ind_pt, converge = 0))
  }
}


## 6. Get estiamted coefficient and standard error
semi_rk_fit <- function(x, y, delta, r0, r, ssp_type, se = TRUE, alpha = 0.2) {
  n <- nrow(x)
  ssps <- semi_rk_ssp(x, y, delta, r0, ssp_type, alpha)
  if (is.na(ssps[[1]][1])) {
    return(list(coe = NA, std = NA, converge = ssps$converge))
  }
  pi <- ssps$ssp
  ind_pt <- ssps$index.pilot
  ind_r <- sample(n, r, prob = pi, replace = TRUE)
  ind <- c(ind_r, ind_pt)
  sec_ssp <- c(pi[ind_r], rep(1 / n, r0))
  est <- semi_rk_est(x[ind, ], y[ind], delta[ind], sec_ssp, n)
  if (est$converge %in% c(1, 2)) {
    return(list(coe = NA, std = NA, converge = est$converge))
  }
  coe <- as.vector(est$coefficient)
  iter <- est$iter
  if (se) {
    r <- r + r0
    g <- gehan_s_mtg(x[ind, ], y[ind], delta[ind], sec_ssp, 
                     coe, seq_along(ind)-1, n)
    vc <- crossprod(g) * r
    vc_amend <- crossprod(sec_ssp * g, g) * r * n
    vc <- vc[-1, -1] + (r / n) * vc_amend[-1, -1]
    m_inv <- solve(gehan_s_jaco(x[ind, -1], y[ind], delta[ind],
                                sec_ssp, coe[-1], n))
    vx <- m_inv %*% vc %*% m_inv / r
    std <- c(sqrt(diag(vx)))
  } else {
    std <- NA
  }
  return(list(coe = coe, std = std, iter = iter, converge = 0))
}
