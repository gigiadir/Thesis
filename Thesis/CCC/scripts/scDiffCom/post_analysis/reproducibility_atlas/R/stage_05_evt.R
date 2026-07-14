# Stage 5 — EVT / GPD tail calibration (Knijnenburg et al. 2009).
#
# The cross-gene null of R_g is not Gaussian, so a z-score p-value is wrong in
# the tail. By the Pickands-Balkema-de Haan theorem the exceedances above a high
# threshold of (almost) any distribution converge to a Generalized Pareto
# Distribution (GPD). We therefore fit a GPD to the null tail and read the
# p-value off the fitted tail — exactly the "fewer permutations, more accurate
# P-values" recipe of Knijnenburg et al. 2009, Bioinformatics 25(12):i161-i168.
#
# Self-contained: GPD MLE + Anderson-Darling goodness-of-fit (via parametric
# bootstrap, Choulakian & Stephens 2001 statistic) using base R only — no
# external EVT packages required.
#
# GPD parameterization (excess over threshold, x >= 0, beta > 0):
#   xi != 0 : F(x) = 1 - (1 + xi * x / beta) ^ (-1 / xi)
#   xi == 0 : F(x) = 1 - exp(-x / beta)
# For xi < 0 the support is bounded: x in [0, -beta/xi].

pgpd_local <- function(q, xi, beta) {
  q <- pmax(q, 0)
  if (abs(xi) < 1e-8) {
    p <- 1 - exp(-q / beta)
  } else {
    z <- 1 + xi * q / beta
    z <- pmax(z, 0)              # beyond the (xi<0) upper support -> F = 1
    p <- 1 - z^(-1 / xi)
  }
  pmin(pmax(p, 0), 1)
}

rgpd_local <- function(n, xi, beta) {
  u <- stats::runif(n)
  if (abs(xi) < 1e-8) {
    -beta * log(1 - u)
  } else {
    (beta / xi) * ((1 - u)^(-xi) - 1)
  }
}

gpd_negloglik <- function(par, x) {
  xi <- par[1]
  beta <- exp(par[2])
  n <- length(x)
  if (abs(xi) < 1e-8) {
    return(n * log(beta) + sum(x) / beta)
  }
  z <- 1 + xi * x / beta
  if (any(z <= 0)) return(1e10)
  n * log(beta) + (1 + 1 / xi) * sum(log(z))
}

gpd_fit_mle <- function(x) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 10L) return(list(xi = NA_real_, beta = NA_real_, conv = FALSE, n = n))
  m <- mean(x)
  v <- stats::var(x)
  # Method-of-moments starting values.
  xi0 <- 0.5 * (1 - m^2 / v)
  beta0 <- 0.5 * m * (m^2 / v + 1)
  if (!is.finite(beta0) || beta0 <= 0) beta0 <- max(m, 1e-6)
  if (!is.finite(xi0)) xi0 <- 0.1
  fit <- tryCatch(
    stats::optim(
      c(xi0, log(beta0)), gpd_negloglik, x = x,
      method = "Nelder-Mead", control = list(maxit = 1000L)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0L) {
    return(list(xi = NA_real_, beta = NA_real_, conv = FALSE, n = n))
  }
  list(xi = fit$par[1], beta = exp(fit$par[2]), conv = TRUE, n = n)
}

# Anderson-Darling statistic for a fitted GPD (Choulakian & Stephens 2001).
gpd_ad_statistic <- function(x, xi, beta) {
  x <- sort(x[is.finite(x) & x > 0])
  n <- length(x)
  if (n < 5L) return(NA_real_)
  z <- pgpd_local(x, xi, beta)
  eps <- 1e-12
  z <- pmin(pmax(z, eps), 1 - eps)
  i <- seq_len(n)
  -n - sum((2 * i - 1) * (log(z) + log(1 - rev(z)))) / n
}

# Parametric-bootstrap p-value for the AD GOF test (H0: excesses ~ GPD).
# Large p => cannot reject => GPD is an adequate fit.
gpd_ad_gof_pvalue <- function(x, xi, beta, B = 199L, seed = NULL) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 10L || is.na(xi) || is.na(beta)) return(NA_real_)
  A2_obs <- gpd_ad_statistic(x, xi, beta)
  if (!is.finite(A2_obs)) return(NA_real_)
  if (!is.null(seed)) set.seed(seed)
  count <- 0L
  valid <- 0L
  for (b in seq_len(B)) {
    xb <- rgpd_local(n, xi, beta)
    fb <- gpd_fit_mle(xb)
    if (!fb$conv) next
    A2b <- gpd_ad_statistic(xb, fb$xi, fb$beta)
    if (!is.finite(A2b)) next
    valid <- valid + 1L
    if (A2b >= A2_obs) count <- count + 1L
  }
  if (valid < 20L) return(NA_real_)
  (count + 1) / (valid + 1)
}

# Hybrid empirical / GPD p-value for a single statistic vs its null vector.
#   - bulk  (>= exceedance_min null values >= obs): empirical p (reliable)
#   - tail  (fewer exceedances): GPD extrapolation on the top n_exc null values,
#           shrinking n_exc by 10 whenever the AD GOF rejects the fit.
calibrate_evt_gpd <- function(obs, null_vec,
                              n_exc = 250L,
                              gof_alpha = 0.05,
                              exceedance_min = 10L,
                              gof_boot = 199L,
                              seed = NULL) {
  null_vec <- null_vec[is.finite(null_vec)]
  N <- length(null_vec)
  out <- list(
    obs = obs, n_null = N, empirical_p = NA_real_, p = NA_real_,
    method = "none", n_exc = NA_integer_, xi = NA_real_, beta = NA_real_,
    threshold = NA_real_, gof_p = NA_real_
  )

  if (!is.finite(obs) || N < 1L) return(out)

  M <- sum(null_vec >= obs)
  emp_p_floor <- (M + 1) / (N + 1)
  out$empirical_p <- emp_p_floor

  # Too few draws to attempt EVT — trust the (floored) empirical estimate.
  if (N < 50L) {
    out$p <- emp_p_floor
    out$method <- "empirical_smallN"
    return(out)
  }

  # Bulk: empirical p is reliable.
  if (M >= exceedance_min) {
    out$p <- emp_p_floor
    out$method <- "empirical"
    return(out)
  }

  # Tail: GPD extrapolation, shrinking n_exc on GOF failure.
  sorted_desc <- sort(null_vec, decreasing = TRUE)
  top <- min(as.integer(n_exc), N - 1L)
  nexc_seq <- seq(top, 10L, by = -10L)
  for (nexc in nexc_seq) {
    thr <- sorted_desc[nexc + 1L]              # threshold just below top nexc
    exceed <- null_vec[null_vec > thr] - thr
    if (length(exceed) < 10L) next
    fit <- gpd_fit_mle(exceed)
    if (!fit$conv) next
    gof_p <- gpd_ad_gof_pvalue(
      exceed, fit$xi, fit$beta, B = gof_boot,
      seed = if (!is.null(seed)) seed + nexc else NULL
    )
    # Accept the fit when GOF cannot reject it (or cannot be evaluated).
    if (is.na(gof_p) || gof_p > gof_alpha) {
      z <- obs - thr
      if (z <= 0) next
      tail_prob <- 1 - pgpd_local(z, fit$xi, fit$beta)
      p_gpd <- (length(exceed) / N) * tail_prob
      out$p <- max(p_gpd, .Machine$double.xmin)
      out$method <- "gpd"
      out$n_exc <- length(exceed)
      out$xi <- fit$xi
      out$beta <- fit$beta
      out$threshold <- thr
      out$gof_p <- gof_p
      return(out)
    }
  }

  # No acceptable GPD fit anywhere on the tail — fall back to empirical.
  out$p <- emp_p_floor
  out$method <- "empirical_fallback"
  out
}

# Apply calibrate_evt_gpd to each gene (row of null_mat vs obs_vec[g]).
calibrate_evt_genes <- function(obs_vec, null_mat,
                                n_exc = 250L, gof_alpha = 0.05,
                                exceedance_min = 10L, gof_boot = 199L,
                                seed = NULL, n_cores = 1L) {
  G <- length(obs_vec)
  worker <- function(g) {
    calibrate_evt_gpd(
      obs = obs_vec[g], null_vec = null_mat[g, ],
      n_exc = n_exc, gof_alpha = gof_alpha,
      exceedance_min = exceedance_min, gof_boot = gof_boot,
      seed = if (!is.null(seed)) seed + g else NULL
    )
  }
  res <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(seq_len(G), worker, mc.cores = n_cores)
  } else {
    lapply(seq_len(G), worker)
  }
  data.frame(
    empirical_p = vapply(res, function(x) x$empirical_p, numeric(1)),
    evt_p       = vapply(res, function(x) x$p, numeric(1)),
    evt_method  = vapply(res, function(x) x$method, character(1)),
    evt_n_exc   = vapply(res, function(x) as.integer(x$n_exc), integer(1)),
    evt_xi      = vapply(res, function(x) x$xi, numeric(1)),
    evt_beta    = vapply(res, function(x) x$beta, numeric(1)),
    evt_gof_p   = vapply(res, function(x) x$gof_p, numeric(1)),
    stringsAsFactors = FALSE
  )
}
