
library(ggplot2)
library(gridExtra)

setwd("C:/Users/christophe.sauser/Documents/Paper3/R")

load("OutColonies.RData")

# PART 1 — COLONY METADATA
meta <- data.frame(
  colony = c("Hornoya", "Bjornoya", "Runde", "Hjelmsoya", "Anda",
             "Skomer", "May", "Rost", "Grumant", "SorGjeslingan", "Brittany"),
  obj    = c("out_Hornoya", "out_Bjornoya", "out_Runde", "out_Hjelmsoya", "out_Anda",
             "out_Skomer", "out_May", "out_Rost", "out_Grumant", "out_SorGjeslingan", "out_Brittany"),
  latitude = c(70.39, 74.45, 62.40, 71.09, 69.05, 51.74, 56.19, 67.52, 78.18, 64.66, 48.06),
  stringsAsFactors = FALSE
)
meta <- meta[order(meta$latitude), ]

N_DRAWS <- 2000

# PART 2 — MATRIX FUNCTIONS
build_matrix <- function(phi, fe, sj = 0.6, b = 0.25) {
  A <- matrix(0, 4, 4)
  A[1, 4] <- (fe / 2) * sj
  A[2, 1] <- sj
  A[3, 2] <- phi * (1 - b)
  A[3, 3] <- phi * b
  A[4, 2] <- phi * b
  A[4, 3] <- phi * (1 - b)
  A[4, 4] <- phi
  A
}

matrix_analysis <- function(A) {
  ev  <- eigen(A)
  lam <- Re(ev$values[1])
  w   <- Re(ev$vectors[, 1]); w <- w / sum(w)
  evl <- eigen(t(A))
  v   <- Re(evl$vectors[, 1]); v <- v / sum(v * w)
  S_mat <- outer(v, w)
  E_mat <- (A / lam) * S_mat
  e_surv <- E_mat[3,2] + E_mat[3,3] + E_mat[4,2] + E_mat[4,3] + E_mat[4,4]
  e_recr <- E_mat[1,4] + E_mat[2,1]
  s_phi  <- S_mat[3,2] * 0.75 + S_mat[3,3] * 0.25 + S_mat[4,2] * 0.25 +
            S_mat[4,3] * 0.75 + S_mat[4,4]
  s_fe   <- S_mat[1,4] * 0.3
  list(lambda = lam, e_surv = e_surv, e_recr = e_recr,
       s_phi = s_phi, s_fe = s_fe)
}

summarize_vec <- function(x) {
  c(mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    lo   = quantile(x, 0.025, na.rm = TRUE),
    hi   = quantile(x, 0.975, na.rm = TRUE))
}

# PART 3 — PER-COLONY POSTERIOR EXTRACTION

post <- list()
for (i in seq_len(nrow(meta))) {
  col <- meta$colony[i]
  cat("Processing", col, "...\n")

  out <- tryCatch(get(meta$obj[i]), error = function(e) NULL)
  if (is.null(out)) { cat("  -> not found, skipping\n"); next }

  sl <- out$sims.list
  n_total <- length(sl$mean.phi)
  idx <- sample(seq_len(n_total), min(N_DRAWS, n_total))

  phi_draws <- sl$mean.phi[idx]
  fe_draws  <- sl$mean.fe[idx]
  sp_draws  <- sl$sigma.phi[idx]
  sf_draws  <- sl$sigma.fe[idx]

  # -- Posterior of elasticities and sensitivities --
  e_surv <- e_recr <- s_phi_post <- s_fe_post <- lam_det <- numeric(length(idx))
  for (k in seq_along(idx)) {
    res <- matrix_analysis(build_matrix(phi_draws[k], fe_draws[k]))
    e_surv[k]     <- res$e_surv
    e_recr[k]     <- res$e_recr
    s_phi_post[k] <- res$s_phi
    s_fe_post[k]  <- res$s_fe
    lam_det[k]    <- res$lambda
  }

  # -- Stochastic lambda: geometric mean over time series ---
  # sl$lambda is a matrix [n_iter x T_years]
  lam_mat <- sl$lambda[idx, , drop = FALSE]
  lam_s <- apply(log(lam_mat), 1, mean, na.rm = TRUE)

  # - Process variance in lambda (from posterior of sigma.phi + sigma.fe) -
  # Combined contribution: var(log lambda) approximated per draw via delta method
  # var(logit phi) = sigma.phi^2 -> var(phi) ~ (phi*(1-phi))^2 * sigma.phi^2
  phi_mean <- phi_draws
  var_phi_approx <- (phi_mean * (1 - phi_mean))^2 * sp_draws^2
  var_fe_approx  <- fe_draws^2 * sf_draws^2

  post[[col]] <- list(
    colony   = col,
    latitude = meta$latitude[i],
    phi_draws = phi_draws,
    fe_draws  = fe_draws,
    sp_draws  = sp_draws,
    sf_draws  = sf_draws,
    e_surv    = e_surv,
    e_recr    = e_recr,
    s_phi     = s_phi_post,
    s_fe      = s_fe_post,
    lam_det   = lam_det,
    lam_s     = lam_s,
    var_phi   = var_phi_approx,
    var_fe    = var_fe_approx,
    lam_mat   = lam_mat
  )
}

colonies_ok <- names(post)

# PART 4 — SUMMARY TABLE (one row per colony)

S <- do.call(rbind, lapply(post, function(p) {
  data.frame(
    colony   = p$colony,
    latitude = p$latitude,
    mean_phi    = mean(p$phi_draws),
    mean_phi_lo = quantile(p$phi_draws, 0.025),
    mean_phi_hi = quantile(p$phi_draws, 0.975),
    sd_phi      = sd(p$phi_draws),
    mean_fe    = mean(p$fe_draws),
    mean_fe_lo = quantile(p$fe_draws, 0.025),
    mean_fe_hi = quantile(p$fe_draws, 0.975),
    sd_fe      = sd(p$fe_draws),
    sigma_phi = mean(p$sp_draws),
    sigma_fe  = mean(p$sf_draws),
    e_surv    = mean(p$e_surv),
    e_surv_lo = quantile(p$e_surv, 0.025),
    e_surv_hi = quantile(p$e_surv, 0.975),
    e_recr    = mean(p$e_recr),
    e_recr_lo = quantile(p$e_recr, 0.025),
    e_recr_hi = quantile(p$e_recr, 0.975),
    s_phi    = mean(p$s_phi),
    s_phi_lo = quantile(p$s_phi, 0.025),
    s_phi_hi = quantile(p$s_phi, 0.975),
    s_fe    = mean(p$s_fe),
    s_fe_lo = quantile(p$s_fe, 0.025),
    s_fe_hi = quantile(p$s_fe, 0.975),
    lam_det    = mean(p$lam_det),
    lam_s      = mean(exp(p$lam_s)),
    lam_s_lo   = quantile(exp(p$lam_s), 0.025),
    lam_s_hi   = quantile(exp(p$lam_s), 0.975),
    log_lam_s  = mean(p$lam_s),
    var_phi_mean = mean(p$var_phi),
    var_fe_mean  = mean(p$var_fe),
    stringsAsFactors = FALSE
  )
}))
S <- S[order(S$latitude), ]
S$colony <- factor(S$colony, levels = S$colony)


# PART 5 — EIV REGRESSIONS (Monte Carlo uncertainty propagation)

eiv_regression <- function(post_list, x_var, y_var, n_draws = N_DRAWS) {
  cols <- names(post_list)
  n_col <- length(cols)
  slopes <- r2s <- numeric(n_draws)
  for (k in seq_len(n_draws)) {
    x <- sapply(cols, function(col) post_list[[col]][[x_var]][k])
    y <- sapply(cols, function(col) post_list[[col]][[y_var]][k])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4) { slopes[k] <- r2s[k] <- NA; next }
    m <- lm(y[ok] ~ x[ok])
    slopes[k] <- coef(m)[2]
    r2s[k]    <- summary(m)$r.squared
  }
  list(
    slope_mean = mean(slopes, na.rm = TRUE),
    slope_lo   = quantile(slopes, 0.025, na.rm = TRUE),
    slope_hi   = quantile(slopes, 0.975, na.rm = TRUE),
    slope_p    = mean(slopes < 0, na.rm = TRUE),  # P(slope < 0)
    r2_mean    = mean(r2s, na.rm = TRUE),
    slopes     = slopes,
    r2s        = r2s
  )
}

# Prediction 1: trade-off phi vs fe
eiv_tradeoff  <- eiv_regression(post, "phi_draws", "fe_draws")

# Prediction 2a: phi ~ latitude
eiv_phi_lat   <- eiv_regression(post, "phi_draws",
                   setNames(lapply(post, function(p) rep(p$latitude, N_DRAWS)), colonies_ok))

# Better: build a list with latitude repeated as draws for EIV
lat_draws <- lapply(post, function(p) rep(p$latitude, length(p$phi_draws)))

eiv_phi_lat_reg  <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         phi_draws = post[[col]]$phi_draws[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "phi_draws")

eiv_fe_lat_reg   <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         fe_draws = post[[col]]$fe_draws[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "fe_draws")

# Prediction 2b: sensitivity of lambda to phi ~ latitude
eiv_sphi_lat <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         s_phi = post[[col]]$s_phi[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "s_phi")

# Prediction 3: sensitivity of lambda to fe ~ latitude
eiv_sfe_lat  <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         s_fe = post[[col]]$s_fe[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "s_fe")

# Elasticity ~ latitude
eiv_esurv_lat <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         e_surv = post[[col]]$e_surv[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "e_surv")

eiv_erecr_lat <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         e_recr = post[[col]]$e_recr[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "e_recr")



write.csv(S, "colony_summary_v5.csv", row.names = FALSE)