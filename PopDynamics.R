library(ggplot2)
library(gridExtra)

setwd("C:/Users/christophe.sauser/Documents/Paper3/R")

load("OutColonies.RData")

# PART 1 — COLONY METADATA

meta <- data.frame(
  colony   = c("Hornoya", "Bjornoya", "Runde", "Hjelmsoya", "Anda",
               "Skomer", "May", "Rost", "Grumant", "SorGjeslingan", "Brittany"),
  obj      = c("out_Hornoya", "out_Bjornoya", "out_Runde", "out_Hjelmsoya", "out_Anda",
               "out_Skomer", "out_May", "out_Rost", "out_Grumant", "out_SorGjeslingan", "out_Brittany"),
  latitude = c(70.39, 74.45, 62.40, 71.09, 69.05, 51.74, 56.19, 67.52, 78.18, 64.66, 48.06),
  yr_start = c(1990,  2004,  2011,  2001,  2005,  1978,  1986,  2003,  2008,  2011,  1979),
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

# PART 3 — PER-COLONY POSTERIOR EXTRACTION

post <- list()
for (i in seq_len(nrow(meta))) {
  col <- meta$colony[i]
  cat("Processing", col, "...\n")

  out <- tryCatch(get(meta$obj[i]), error = function(e) NULL)
  if (is.null(out)) { cat("  -> not found, skipping\n"); next }

  sl      <- out$sims.list
  n_total <- length(sl$mean.phi)
  idx     <- sample(seq_len(n_total), min(N_DRAWS, n_total))

  phi_draws <- sl$mean.phi[idx]
  fe_draws  <- sl$mean.fe[idx]
  sp_draws  <- sl$sigma.phi[idx]
  sf_draws  <- sl$sigma.fe[idx]

  e_surv <- e_recr <- s_phi_post <- s_fe_post <- lam_det <- numeric(length(idx))
  for (k in seq_along(idx)) {
    res           <- matrix_analysis(build_matrix(phi_draws[k], fe_draws[k]))
    e_surv[k]     <- res$e_surv
    e_recr[k]     <- res$e_recr
    s_phi_post[k] <- res$s_phi
    s_fe_post[k]  <- res$s_fe
    lam_det[k]    <- res$lambda
  }

  lam_mat    <- sl$lambda[idx, , drop = FALSE]
  lam_s_obs  <- apply(log(lam_mat), 1, mean, na.rm = TRUE)
  omega_mat  <- sl$omega[idx, , drop = FALSE]
  omega_mean <- apply(omega_mat, 1, mean, na.rm = TRUE)

  N_YEARS <- 500
  phi_sim <- matrix(plogis(rnorm(length(idx) * N_YEARS, qlogis(phi_draws), sp_draws)),
                    nrow = length(idx), ncol = N_YEARS)
  fe_sim  <- matrix(exp(rnorm(length(idx) * N_YEARS, log(fe_draws), sf_draws)),
                    nrow = length(idx), ncol = N_YEARS)

  lam_s_sim <- sapply(seq_along(idx), function(k) {
    mean(sapply(seq_len(N_YEARS), function(yr) {
      log(Re(eigen(build_matrix(phi_sim[k, yr], fe_sim[k, yr]))$values[1]))
    }))
  })

  var_phi_approx <- (phi_draws * (1 - phi_draws))^2 * sp_draws^2
  var_fe_approx  <- fe_draws^2 * sf_draws^2

  # Annual time series posterior summaries
  yr_start <- meta$yr_start[meta$colony == col]
  n_occ    <- ncol(sl$phi.a) + 1
  yr_cmr_v <- yr_start:(yr_start + n_occ - 1)

  phi_ts    <- apply(sl$phi.a,  2, mean,            na.rm = TRUE)
  phi_ts_lo <- apply(sl$phi.a,  2, quantile, 0.025, na.rm = TRUE)
  phi_ts_hi <- apply(sl$phi.a,  2, quantile, 0.975, na.rm = TRUE)
  yr_phi    <- yr_cmr_v[seq_len(ncol(sl$phi.a))]

  fe_ts    <- apply(sl$fe,  2, mean,            na.rm = TRUE)
  fe_ts_lo <- apply(sl$fe,  2, quantile, 0.025, na.rm = TRUE)
  fe_ts_hi <- apply(sl$fe,  2, quantile, 0.975, na.rm = TRUE)
  yr_fe    <- yr_cmr_v[seq_len(ncol(sl$fe))]

  N4_ts    <- apply(sl$N4,  2, mean,            na.rm = TRUE)
  N4_ts_lo <- apply(sl$N4,  2, quantile, 0.025, na.rm = TRUE)
  N4_ts_hi <- apply(sl$N4,  2, quantile, 0.975, na.rm = TRUE)
  yr_N4    <- yr_cmr_v[seq_len(ncol(sl$N4))]

  lam_ts    <- apply(sl$lambda, 2, mean,            na.rm = TRUE)
  lam_ts_lo <- apply(sl$lambda, 2, quantile, 0.025, na.rm = TRUE)
  lam_ts_hi <- apply(sl$lambda, 2, quantile, 0.975, na.rm = TRUE)
  yr_lam    <- yr_cmr_v[seq_len(ncol(sl$lambda))]

  # Temporal correlations on posterior mean time series
  T_common <- min(length(phi_ts), length(lam_ts), length(fe_ts))
  ok_t     <- is.finite(phi_ts[1:T_common]) &
               is.finite(fe_ts[1:T_common])  &
               is.finite(lam_ts[1:T_common])

  r_lam_phi <- if (sum(ok_t) >= 4) cor(lam_ts[1:T_common][ok_t], phi_ts[1:T_common][ok_t]) else NA
  r_lam_fe  <- if (sum(ok_t) >= 4) cor(lam_ts[1:T_common][ok_t], fe_ts[1:T_common][ok_t])  else NA
  r_phi_fe  <- if (sum(ok_t) >= 4) cor(phi_ts[1:T_common][ok_t], fe_ts[1:T_common][ok_t])  else NA

  # LTRE variance decomposition on posterior mean time series
  s_phi_m <- mean(s_phi_post)
  s_fe_m  <- mean(s_fe_post)
  phi_ok  <- phi_ts[1:T_common][ok_t]
  fe_ok   <- fe_ts[1:T_common][ok_t]

  C_phi  <- s_phi_m^2 * var(phi_ok)
  C_fe   <- s_fe_m^2  * var(fe_ok)
  C_cov  <- 2 * s_phi_m * s_fe_m * cov(phi_ok, fe_ok)
  C_tot  <- C_phi + C_fe + C_cov

  post[[col]] <- list(
    colony     = col,
    latitude   = meta$latitude[i],
    yr_start   = yr_start,
    phi_draws  = phi_draws, fe_draws  = fe_draws,
    sp_draws   = sp_draws,  sf_draws  = sf_draws,
    e_surv     = e_surv,    e_recr    = e_recr,
    s_phi      = s_phi_post, s_fe     = s_fe_post,
    lam_det    = lam_det,
    lam_s_obs  = lam_s_obs, lam_s_sim = lam_s_sim,
    var_phi    = var_phi_approx, var_fe = var_fe_approx,
    omega_mean = omega_mean,
    phi_ts = phi_ts, phi_ts_lo = phi_ts_lo, phi_ts_hi = phi_ts_hi, yr_phi = yr_phi,
    fe_ts  = fe_ts,  fe_ts_lo  = fe_ts_lo,  fe_ts_hi  = fe_ts_hi,  yr_fe  = yr_fe,
    N4_ts  = N4_ts,  N4_ts_lo  = N4_ts_lo,  N4_ts_hi  = N4_ts_hi,  yr_N4  = yr_N4,
    lam_ts = lam_ts, lam_ts_lo = lam_ts_lo, lam_ts_hi = lam_ts_hi, yr_lam = yr_lam,
    r_lam_phi = r_lam_phi, r_lam_fe = r_lam_fe, r_phi_fe = r_phi_fe,
    ltre_phi = C_phi / C_tot, ltre_fe = C_fe / C_tot, ltre_cov = C_cov / C_tot
  )
}

colonies_ok <- names(post)

# PART 4 — SUMMARY TABLE

S <- do.call(rbind, lapply(post, function(p) {
  data.frame(
    colony   = p$colony, latitude = p$latitude,
    mean_phi    = mean(p$phi_draws),
    mean_phi_lo = quantile(p$phi_draws, 0.025),
    mean_phi_hi = quantile(p$phi_draws, 0.975),
    sd_phi      = sd(p$phi_draws),
    mean_fe    = mean(p$fe_draws),
    mean_fe_lo = quantile(p$fe_draws, 0.025),
    mean_fe_hi = quantile(p$fe_draws, 0.975),
    sd_fe      = sd(p$fe_draws),
    sigma_phi   = mean(p$sp_draws),
    sigma_fe    = mean(p$sf_draws),
    e_surv = mean(p$e_surv), e_surv_lo = quantile(p$e_surv, 0.025), e_surv_hi = quantile(p$e_surv, 0.975),
    e_recr = mean(p$e_recr), e_recr_lo = quantile(p$e_recr, 0.025), e_recr_hi = quantile(p$e_recr, 0.975),
    s_phi  = mean(p$s_phi),  s_phi_lo  = quantile(p$s_phi,  0.025), s_phi_hi  = quantile(p$s_phi,  0.975),
    s_fe   = mean(p$s_fe),   s_fe_lo   = quantile(p$s_fe,   0.025), s_fe_hi   = quantile(p$s_fe,   0.975),
    lam_det      = mean(p$lam_det),
    lam_s_obs    = mean(exp(p$lam_s_obs)),
    lam_s_obs_lo = quantile(exp(p$lam_s_obs), 0.025),
    lam_s_obs_hi = quantile(exp(p$lam_s_obs), 0.975),
    lam_s_sim    = mean(exp(p$lam_s_sim)),
    lam_s_sim_lo = quantile(exp(p$lam_s_sim), 0.025),
    lam_s_sim_hi = quantile(exp(p$lam_s_sim), 0.975),
    delta_lam_s  = mean(p$lam_s_obs - p$lam_s_sim),
    var_phi_mean = mean(p$var_phi), var_fe_mean = mean(p$var_fe),
    omega    = mean(p$omega_mean),
    omega_lo = quantile(p$omega_mean, 0.025),
    omega_hi = quantile(p$omega_mean, 0.975),
    r_lam_phi = p$r_lam_phi, r_lam_fe = p$r_lam_fe, r_phi_fe = p$r_phi_fe,
    ltre_phi  = p$ltre_phi,  ltre_fe  = p$ltre_fe,  ltre_cov = p$ltre_cov,
    stringsAsFactors = FALSE
  )
}))
S <- S[order(S$latitude), ]
S$colony <- factor(S$colony, levels = S$colony)

# PART 5 — EIV REGRESSIONS

eiv_regression <- function(post_list, x_var, y_var, n_draws = N_DRAWS) {
  cols   <- names(post_list)
  slopes <- r2s <- numeric(n_draws)
  for (k in seq_len(n_draws)) {
    x  <- sapply(cols, function(col) post_list[[col]][[x_var]][k])
    y  <- sapply(cols, function(col) post_list[[col]][[y_var]][k])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4) { slopes[k] <- r2s[k] <- NA; next }
    m         <- lm(y[ok] ~ x[ok])
    slopes[k] <- coef(m)[2]
    r2s[k]    <- summary(m)$r.squared
  }
  list(slope_mean = mean(slopes, na.rm = TRUE),
       slope_lo   = quantile(slopes, 0.025, na.rm = TRUE),
       slope_hi   = quantile(slopes, 0.975, na.rm = TRUE),
       slope_p    = mean(slopes < 0, na.rm = TRUE),
       r2_mean    = mean(r2s, na.rm = TRUE),
       slopes = slopes, r2s = r2s)
}

eiv_tradeoff    <- eiv_regression(post, "phi_draws", "fe_draws")

eiv_phi_lat_reg <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         phi_draws = post[[col]]$phi_draws[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "phi_draws")

eiv_fe_lat_reg  <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         fe_draws = post[[col]]$fe_draws[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "fe_draws")

eiv_sphi_lat <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         s_phi = post[[col]]$s_phi[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "s_phi")

eiv_sfe_lat  <- eiv_regression(
  setNames(lapply(colonies_ok, function(col)
    list(lat = rep(post[[col]]$latitude, N_DRAWS),
         s_fe = post[[col]]$s_fe[seq_len(N_DRAWS)])), colonies_ok),
  "lat", "s_fe")

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

write.csv(S, "colony_summary.csv", row.names = FALSE)
