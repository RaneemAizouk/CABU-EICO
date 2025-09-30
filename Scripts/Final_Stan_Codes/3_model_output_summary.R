###############################################################################
# THIS CODE PRODUCED SUMMARY PLOTS FOR THE MODEL OUTPUT
###############################################################################

# Date created: 12 August 2025
# Date last updated: 12 August 2025
# Author: Esther van Kleef

rm(list=ls())

#------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------
pacman::p_load(posterior, bayesplot, dplyr, shinystan, ggplot2, tibble, tidyr, lubridate, gridExtra,grid, 
               rstan, loo, flextable, purrr, cowplot, forcats, future.apply, wesanderson)

#------------------------------------------------------------------------------
# Load data & fit
#------------------------------------------------------------------------------
data_source <- "observed"  # "simulated" or "observed"
intervention <- "Two_step_" #"Two_step_" or "One_step"
scen    <- "seasonality" # "seasonality" or "noseasonality"
seasonality <- "sine_" # "sine_" or "spline_"
intervention_effect <- "NonAdd_NonCol_" # "" or "NonAdd_" or "NonAdd_NonCol_" or ""

scenario = paste0(intervention, seasonality, scen)
  
if (data_source == "simulated") {
  scen_data <- paste0("Simulated_data_", ifelse(scenario == paste0("Two_step_", seasonality,intervention_effect, "seasonality"), "seasonality", "noseasonality"))
  stan_data_fit <- readRDS(paste0("./Data/Simulated_data/", scen_data, "_stan_data.rds"))
  sim_df        <- readRDS(paste0("./Data/Simulated_data/", scen_data, ".rds"))
  data_fit      <- sim_df
} else {
  stan_data_fit <- readRDS("./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")
  data_fit      <- read.csv("./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")
}

str(stan_data_fit)

# Used for setting up spline basis
#------------------------------------------------------------------------------
num_knots <- 5
knots <- seq(0, 1, length.out = num_knots)

knots <- as.numeric(knots)
spline_degree <- 3
num_basis <- num_knots + spline_degree 


dir = ifelse(data_source=="simulated", "Simulated_data/","Observed_data/")
fit <- readRDS(paste0("./Output/Model_results/", dir,intervention, seasonality, intervention_effect, scen, data_source,".rds"))


#------------------------------------------------------------------------------
# Traceplots & ShinyStan
#------------------------------------------------------------------------------
# For spline
if (trimws(seasonality) == "spline_") {
  trace_pars <- c("a_1_2","q_1_2_raw","q_2_1_raw",
                  "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                  "beta_int1_1","beta_int1_2",#"beta_int2_1","beta_int2_2",
                  "sigma_u","sigma_q_1_2","sigma_q_2_1")
} else { # for sine
 trace_pars <- c("a1","phi","q_1_2_raw","q_2_1_raw",#"q_1_2_base","q_2_1_base",
                 "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                 "beta_int1_1","beta_int1_2",#"beta_int2_1","beta_int2_2",
                 "sigma_u","sigma_q_1_2","sigma_q_2_1")
}

trace_pars

p = rstan::traceplot(fit, pars = trace_pars, inc_warmup = FALSE, nrow = 3)
p

# Function to run diagnostics
# NUTS stats table (like Shiny) 
nuts_stats_table_stanfit <- function(fit, stat = c("Mean","SD","Max","Min"),
                                     inc_warmup = FALSE, digits = 4,
                                     log_lik_param = "log_lik") {
  stopifnot(inherits(fit, "stanfit"))
  stat <- match.arg(stat)
  
  # --- sampler diagnostics (like Shiny) ---
  L <- rstan::get_sampler_params(fit, inc_warmup = inc_warmup)  # list per chain
  cols_base <- c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "energy__")
  has_div <- "divergent__" %in% colnames(L[[1]])
  
  aggfun <- switch(stat,
                   Mean = function(v) mean(v, na.rm = TRUE),
                   SD   = function(v)  sd(v,  na.rm = TRUE),
                   Max  = function(v)  suppressWarnings(max(v, na.rm = TRUE)),
                   Min  = function(v)  suppressWarnings(min(v, na.rm = TRUE))
  )
  
  # Per-chain table (base stats)
  per_chain_base <- lapply(seq_along(L), function(i) {
    x <- L[[i]][, intersect(cols_base, colnames(L[[i]])), drop = FALSE]
    out <- as.data.frame(lapply(as.data.frame(x), aggfun))
    names(out) <- sub("__$", "", names(out))
    cbind(chain = paste0("chain", i), out, row.names = NULL)
  })
  tab <- do.call(rbind, per_chain_base)
  
  # Per-chain divergences (count & %)
  if (has_div) {
    n_div_chain   <- sapply(L, function(m) sum(m[, "divergent__"] > 0))
    iters_chain   <- sapply(L, nrow)
    pct_div_chain <- 100 * n_div_chain / iters_chain
    tab$n_divergent   <- n_div_chain
    tab$pct_divergent <- pct_div_chain
  } else {
    tab$n_divergent   <- NA_integer_
    tab$pct_divergent <- NA_real_
  }
  
  # "All chains" row (base stats)
  S <- do.call(rbind, L)
  all_mat <- S[, intersect(cols_base, colnames(S)), drop = FALSE]
  all_row <- as.data.frame(lapply(as.data.frame(all_mat), aggfun))
  names(all_row) <- sub("__$", "", names(all_row))
  if (has_div) {
    n_div_total   <- sum(S[, "divergent__"] > 0)
    iters_total   <- nrow(S)
    pct_div_total <- 100 * n_div_total / iters_total
  } else {
    n_div_total <- NA_integer_
    pct_div_total <- NA_real_
  }
  
  # R-hat min/max (exclude lp__), guard missing column
  sm_sum <- rstan::summary(fit)$summary
  if (!is.null(sm_sum) && nrow(sm_sum)) {
    if ("lp__" %in% rownames(sm_sum)) sm_sum <- sm_sum[setdiff(rownames(sm_sum), "lp__"), , drop = FALSE]
    if ("Rhat" %in% colnames(sm_sum)) {
      rhat_min <- suppressWarnings(min(sm_sum[, "Rhat"], na.rm = TRUE))
      rhat_max <- suppressWarnings(max(sm_sum[, "Rhat"], na.rm = TRUE))
    } else {
      rhat_min <- NA_real_; rhat_max <- NA_real_
    }
  } else {
    rhat_min <- NA_real_; rhat_max <- NA_real_
  }
  
  # --- LOO block with robust guards ---
  loo_obj <- NULL  # will remain NULL if unavailable
  
  ll <- tryCatch(
    loo::extract_log_lik(fit, parameter_name = log_lik_param, merge_chains = FALSE),
    error = function(e) NULL
  )
  
  # Identify obs entirely NA across all draws/chains (e.g., first obs per individual)
  if (!is.null(ll)) {
    drop_idx <- which(vapply(seq_len(dim(ll)[3]),
                             function(n) all(is.na(ll[,,n])),
                             logical(1)))
    ll_for_loo <- if (length(drop_idx)) ll[, , -drop_idx, drop = FALSE] else ll
    
    if (length(dim(ll_for_loo)) == 3 && dim(ll_for_loo)[3] > 0) {
      loo_obj <- tryCatch(
        withCallingHandlers(
          loo::loo.array(ll_for_loo),
          warning = function(w) {
            # Quiet PSIS k>0.7 etc., but keep the returned object
            tryInvokeRestart("muffleWarning")
          }
        ),
        error = function(e) NULL
      )
    }
  }
  
  # Prepare LOO/Rhat fields (default NA)
  loo_fields <- list(
    Rhat_min = rhat_min, Rhat_max = rhat_max,
    elpd_loo = NA_real_, se_elpd_loo = NA_real_,
    LOOIC = NA_real_,   se_LOOIC = NA_real_,
    p_loo = NA_real_,   se_p_loo = NA_real_,
    pareto_k_max = NA_real_, pareto_k_gt_0.5 = NA_integer_,
    pareto_k_gt_0.7 = NA_integer_, n_points = NA_integer_
  )
  
  if (!is.null(loo_obj)) {
    est <- loo_obj$estimates
    k   <- loo_obj$diagnostics$pareto_k
    if (!is.null(est) && all(c("elpd_loo","p_loo") %in% rownames(est))) {
      elpd <- est["elpd_loo","Estimate"]; se_e <- est["elpd_loo","SE"]
      ploo <- est["p_loo","Estimate"];    se_p <- est["p_loo","SE"]
      
      loo_fields$elpd_loo   <- elpd
      loo_fields$se_elpd_loo <- se_e
      loo_fields$LOOIC      <- -2 * elpd
      loo_fields$se_LOOIC   <-  2 * se_e
      loo_fields$p_loo      <- ploo
      loo_fields$se_p_loo   <- se_p
    }
    if (!is.null(k)) {
      loo_fields$pareto_k_max   <- suppressWarnings(max(k, na.rm = TRUE))
      loo_fields$pareto_k_gt_0.5 <- sum(k > 0.5, na.rm = TRUE)
      loo_fields$pareto_k_gt_0.7 <- sum(k > 0.7, na.rm = TRUE)
      loo_fields$n_points        <- length(k)
    }
  }
  
  # Bind: LOO/Rhat + total divergence only on "All chains" row
  all_out <- cbind(
    chain = "All chains",
    all_row,
    n_divergent   = n_div_total,
    pct_divergent = pct_div_total,
    as.data.frame(loo_fields)
  )
  
  per_chain_out <- cbind(tab, as.data.frame(lapply(loo_fields, function(x) NA)))
  out <- rbind(all_out, per_chain_out)
  
  # Round numeric columns
  is_num <- vapply(out, is.numeric, TRUE)
  out[is_num] <- lapply(out[is_num], function(x) ifelse(is.finite(x), round(x, digits), x))
  
  rownames(out) <- NULL
  out
}

diagnostics <- nuts_stats_table_stanfit(fit, stat = "Mean", log_lik_param = "log_lik")

ft <- flextable(diagnostics[,c(1:8)])
ft <- autofit(ft)
ft <- set_caption(ft, "NUTS sampler diagnostics")
ft

ft_f <- flextable(diagnostics[1,c(9:20)])
ft_f <- autofit(ft_f)
ft_f <- set_caption(ft_f, "LOO statistics")
ft_f

# Ran these subsequently, changing the scenario in the first lines of code
ft_f_sin1 = ft_f # One_Step_sine
ft_f_sin2 = ft_f # Two_Step_sine
ft_f_spline1 = ft_f # One_Step_spline
ft_f_spline2 = ft_f # Two_Step_spline

# Compare LOOIC
ft_f_sin1; ft_f_sin2; ft_f_spline1; ft_f_spline2

D <- posterior::as_draws_df(fit)

#------------------------------------------------------------
# PLOT SEASONAL PATTERN
#------------------------------------------------------------

global_interval_start <- ifelse(data_source=="simulated", as.Date("2022-10-01"), as.Date("2022-10-03"))
global_interval_end   <- as.Date("2024-02-19")
interval_length <- 28

# Convert to numeric
global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric   <- as.numeric(global_interval_end)

# Number of intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)

# Interval starts/ends
interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)
interval_ends <- interval_starts + (interval_length - 1)

# Midpoints (special rule for short last interval)
X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2
if ((global_interval_end_numeric - interval_starts[num_intervals]) + 1 < interval_length) {
  X_midpoints[num_intervals] <- interval_starts[num_intervals]
}

phase_grid <- seq(0, 1, length.out = 400)  # smooth one-year grid

# ---- One-year calendar helper (anchors plotting year to Jan 1 of year_anchor) ----
make_calendar <- function(start_date, end_date, interval_len) {
  gs <- as.numeric(start_date); ge <- as.numeric(end_date)
  n  <- ceiling((ge - gs) / interval_len)
  starts <- seq(gs, gs + (n - 1) * interval_len, by = interval_len)
  mids <- (starts + pmin(starts + interval_len - 1, ge)) / 2
  if ((ge - starts[n]) + 1 < interval_len) mids[n] <- starts[n]
  year_anchor <- 2023
  date0 <- as.Date(paste0(year_anchor, "-01-01"))
  list(midpoints = mids, date0 = date0)
}
cal <- make_calendar(global_interval_start, global_interval_end, interval_length)

summarise_draws <- function(mat, probs = c(.05, .5, .95)) {
  qs <- apply(mat, 2, quantile, probs = probs)
  tibble(q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,])
}

#------------------------------------------------------------
# Plot season (HR vs yearly median = 1)
#------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Phase/date helpers
year_len <- 365  # must match Stan
phase_to_date <- function(phase, date0, year_len = year_len) date0 + round(phase * year_len)
date_to_phase <- function(d, origin, year_len = year_len) ((as.numeric(d) - as.numeric(origin)) / year_len) %% 1
plot_phase_curve <- function(
    df,                       # columns: Date, q_low, q50, q_hi (on LOG scale)
    title,
    subtitle,
    ylab = "Hazard ratio (reference = yearly median)",
    rainy = NULL,             # data.frame with columns start, end (Date)
    ref = c("mean","median"),
    show_peak_trough = TRUE,
    show_hr1 = TRUE,          # draw HR=1 reference
    year_median_qs = NULL,    # optional LOG-scale c(0.025,0.5,0.975)
    peak_qs = NULL,           # optional LOG-scale c(0.025,0.5,0.975)
    trough_qs = NULL,         # optional LOG-scale c(0.025,0.5,0.975)
    base_size = 18
) {
  ref <- match.arg(ref)
  fmt <- function(x) formatC(x, format = "f", digits = 2, drop0trailing = TRUE)
  
  # --- LOG -> HR for plotting ---
  dfp <- df
  dfp$q_low <- exp(dfp$q_low); dfp$q50 <- exp(dfp$q50); dfp$q_hi <- exp(dfp$q_hi)
  if (!is.null(year_median_qs)) year_median_qs <- exp(year_median_qs)
  if (!is.null(peak_qs))        peak_qs        <- exp(peak_qs)
  if (!is.null(trough_qs))      trough_qs      <- exp(trough_qs)
  
  # Base layers
  p <- ggplot(dfp, aes(Date, q50)) +
    { if (!is.null(rainy))
      geom_rect(
        data = rainy,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = "forestgreen", alpha = 0.10
      )
    } +
    geom_ribbon(aes(ymin = q_low, ymax = q_hi), alpha = 0.25, fill = "cyan4") +
    geom_line(linewidth = 2, color = "cyan4") +
    { if (show_hr1) geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.9) }
  
  # ------- Dynamic axis limits with extra space BELOW 0 -------
  ymax <- max(dfp$q_hi, na.rm = TRUE)
  ymin <- min(dfp$q_low, na.rm = TRUE)
  
  # how much space to show below zero (both relative and absolute cushions)
  bottom_frac <- 0.10                 # 10% of top scale as room below 0
  bottom_abs  <- 0.08                 # or at least 0.08 HR units
  bottom_pad  <- max(bottom_frac * ymax, bottom_abs)
  
  # lower limit goes BELOW 0 to create visual breathing room
  lower_limit <- min(0 - bottom_pad, ymin * 0.92)
  upper_limit <- ymax * 1.12
  
  # --- rainy label: bottom-left of (visible) band; never clipped ---
  if (!is.null(rainy) && nrow(rainy) > 0) {
    xmin <- min(dfp$Date, na.rm = TRUE)
    xmax <- max(dfp$Date, na.rm = TRUE)
    rainy2 <- subset(rainy, end >= xmin & start <= xmax)
    if (nrow(rainy2) > 0) {
      y_lab <- lower_limit + 0.02 * (upper_limit - lower_limit)  # 6% above bottom
      rainy_labs <- transform(rainy2, x_lab = pmax(start, xmin))
      p <- p + geom_text(
        data = rainy_labs,
        aes(x = x_lab, y = y_lab, label = "Rainy season"),
        inherit.aes = FALSE,
        hjust = 0, vjust = 0,
        color = "darkgreen", fontface = "bold", size = base_size/3.2
      )
    }
  }
  
  # --- peak & trough labels (never overlap the curve) ---
  if (show_peak_trough) {
    # Peak ABOVE ribbon
    pk_ix   <- which.max(dfp$q50)
    pk_date <- dfp$Date[pk_ix]
    pk_top  <- dfp$q_hi[pk_ix] * 1.04
    pk_txt  <- if (!is.null(peak_qs)) {
      sprintf("Peak month: %s\nHR = %s [%s–%s]",
              format(pk_date, "%b"), fmt(peak_qs[2]), fmt(peak_qs[1]), fmt(peak_qs[3]))
    } else {
      sprintf("Peak month: %s\nHR = %s", format(pk_date, "%b"), fmt(dfp$q50[pk_ix]))
    }
    p <- p +
      geom_vline(xintercept = pk_date, linetype = "dotdash",
                 col = "darkred", linewidth = 0.9, alpha = 0.7) +
      annotate("label", x = pk_date, y = pk_top,
               label = pk_txt, vjust = 0, hjust = 0.5,
               label.size = 0, alpha = 0.95, color = "black", size = base_size/3.2)
    
    # Trough BELOW ribbon but safely above axis
    tr_ix   <- which.min(dfp$q50)
    tr_date <- dfp$Date[tr_ix]
    
    # safe floor a bit above the bottom axis
    safe_floor <- lower_limit + 0.08 * (upper_limit - lower_limit)  # 8% above bottom
    tr_bot     <- max(dfp$q_low[tr_ix] * 0.96, safe_floor)
    
    tr_txt  <- if (!is.null(trough_qs)) {
      sprintf("Trough month: %s\nHR = %s [%s–%s]",
              format(tr_date, "%b"), fmt(trough_qs[2]), fmt(trough_qs[1]), fmt(trough_qs[3]))
    } else {
      sprintf("Trough month: %s\nHR = %s", format(tr_date, "%b"), fmt(dfp$q50[tr_ix]))
    }
    p <- p +
      geom_vline(xintercept = tr_date, linetype = "dotdash", alpha = 0.7,
                 col = "darkred", linewidth = 0.9) +
      annotate("label", x = tr_date, y = tr_bot,
               label = tr_txt, vjust = 1, hjust = 0.5,   # TOP of label at tr_bot
               label.size = 0, alpha = 0.95, color = "black", size = base_size/3.2)
  }
  
  # Axis & theme (allow negative bottom for space; prevent clipping)
  p +
    scale_x_date(date_labels = "%b", date_breaks = "1 month") +
    scale_y_continuous(limits = c(lower_limit, upper_limit),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Calendar month",
      y = ylab,
      title = title,
      subtitle = subtitle,
      caption = "Dotted = HR=1 reference"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 2),
      plot.caption  = element_text(size = base_size - 3),
      panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
      panel.grid.minor = element_line(linetype = "dotted", color = "grey85"),
      plot.margin = grid::unit(c(10, 10, 18, 10), "pt")
    )
}

# Rainy season bands (adjust as needed)
rainy <- data.frame(start = as.Date("2023-06-01"), end = as.Date("2023-10-31"))

#-----------------------------------------------------------------
# Model branches (spline / sine)
#-----------------------------------------------------------------

if (trimws(seasonality) == "spline_") {
  
  # Spline 
  K_fit   <- num_knots
  deg_fit <- spline_degree
  knots_vec <- knots
  stopifnot(length(knots_vec) == K_fit)
  nbasis_fit <- K_fit + deg_fit
  
  # Stan-matching periodic B-spline basis
  build_b_spline_vec <- function(t, ext_knots, ind, order) {
    if (order == 1) {
      as.numeric(ext_knots[ind] <= t & t < ext_knots[ind + 1])
    } else {
      d1 <- ext_knots[ind + order - 1] - ext_knots[ind]
      d2 <- ext_knots[ind + order]     - ext_knots[ind + 1]
      w1 <- if (d1 > 0) (t - ext_knots[ind]) / d1 else 0
      w2 <- if (d2 > 0) 1 - (t - ext_knots[ind + 1]) / d2 else 0
      w1 * build_b_spline_vec(t, ext_knots, ind, order - 1) +
        w2 * build_b_spline_vec(t, ext_knots, ind + 1, order - 1)
    }
  }
  
  make_B_matrix <- function(X_in, knots, spline_degree) {
    X_mod <- pmin(X_in %% 1, 1 - 1e-12)
    K  <- length(knots)
    nb <- K + spline_degree
    ext_knots <- c(knots[(K - spline_degree + 1):K] - 1,
                   knots,
                   knots[1:(spline_degree + 1)] + 1)
    order <- spline_degree + 1
    B <- sapply(1:nb, function(ind) build_b_spline_vec(X_mod, ext_knots, ind, order))
    as.matrix(B)
  }
  
  # Build basis on the phase grid
  B_phase <- make_B_matrix(phase_grid, knots_vec, deg_fit)
  
  vars_a  <- sprintf("a_1_2[%d]", seq_len(nbasis_fit))
  A_mat   <- as.matrix(D[, vars_a])  # draws × nbasis
  stopifnot(ncol(A_mat) == ncol(B_phase))
  
  # draws × n_phase (LOG scale)
  Y_log <- A_mat %*% t(B_phase)
  
  # Center per draw → yearly median HR = 1 after exp()
  row_med <- apply(Y_log, 1, median)
  Y_log   <- sweep(Y_log, 1, row_med, "-")
  
  # Timewise 95% CrIs (LOG)
  qs <- apply(Y_log, 2, quantile, probs = c(.05, .5, .95))
  spline_df <- tibble(
    phase = phase_grid,
    q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,]
  ) |>
    mutate(Date = cal$date0 + round(phase * year_len))
  
  # Peak & Trough across draws (LOG) for annotations
  peak_draws   <- apply(Y_log, 1, max)
  trough_draws <- apply(Y_log, 1, min)
  peak_qs      <- quantile(peak_draws,   c(.025, .5, .975), na.rm = TRUE)
  trough_qs    <- quantile(trough_draws, c(.025, .5, .975), na.rm = TRUE)
  
  p_season <- plot_phase_curve(
    df = spline_df,
    title = "Estimated seasonal variation in ESBL-E acquisition risk",
    subtitle = "Cubic spline seasonal function (median hazard ratio with 95% CrI)",
    ylab = NULL,
    rainy = rainy,
    ref = "median",
    show_peak_trough = TRUE,
    show_hr1 = TRUE,
    peak_qs = peak_qs,
    trough_qs = trough_qs
  )
  
} else {
  
  # Sine (evaluate on calendar; auto-align φ convention) 
  jan1 <- cal$date0
  date_grid <- seq(jan1, by = "1 day", length.out = year_len)
  
  # Model origin used to define phase in Stan (change to X_midpoints[1] if that’s what Stan used)
  model_origin <- as.Date(global_interval_start_numeric, origin = "1970-01-01")
  
  # Phase per calendar date
  phase_date <- ((as.numeric(date_grid) - as.numeric(model_origin)) / year_len) %% 1  # length = n_dates
  
  # Draws: amplitude (log scale) and φ
  s <- D[, c("a1","phi")]
  n_draws <- nrow(s); n_dates <- length(date_grid)
  
  # Try common φ conventions; pick the one peaking nearest August
  build_Y <- function(convention = c("sin_plus","sin_minus","cos_plus","sin_plus_cycles")) {
    convention <- match.arg(convention)
    angle_base <- matrix(2*pi*phase_date, nrow = n_dates, ncol = n_draws)
    phi_mat <- matrix(rep(s$phi, each = n_dates), nrow = n_dates, ncol = n_draws)
    a1_mat  <- matrix(rep(s$a1,  each = n_dates), nrow = n_dates, ncol = n_draws)
    Y <- switch(convention,
                sin_plus         = sin(angle_base +  phi_mat) * a1_mat,
                sin_minus        = sin(angle_base -  phi_mat) * a1_mat,
                cos_plus         = cos(angle_base +  phi_mat) * a1_mat,
                sin_plus_cycles  = { phi_rad <- 2*pi*phi_mat; sin(angle_base + phi_rad) * a1_mat }
    )
    t(Y)  # draws × dates
  }
  
  conventions <- c("sin_plus","sin_minus","cos_plus","sin_plus_cycles")
  Y_list <- lapply(conventions, build_Y)
  
  # Center per draw
  Y_list <- lapply(Y_list, function(Y) {
    rm <- apply(Y, 1, median)
    sweep(Y, 1, rm, "-")
  })
  
  # Choose convention whose median curve peaks closest to August
  target_month <- 8L
  med_curves <- lapply(Y_list, function(Y) apply(Y, 2, median))
  peak_dist <- function(dates, y_med) {
    pk <- dates[which.max(y_med)]
    pm <- as.integer(format(pk, "%m"))
    abs(((pm - target_month + 6) %% 12) - 6)
  }
  dists <- vapply(med_curves, function(m) peak_dist(date_grid, m), numeric(1))
  Y_log <- Y_list[[which.min(dists)]]
  
  # Timewise 95% CrIs (LOG)
  qs <- apply(Y_log, 2, quantile, probs = c(.05, .5, .95), na.rm = TRUE)
  sine_df <- tibble(Date = date_grid, q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,])
  
  # Peak & Trough across draws (LOG) for annotations
  peak_draws   <- apply(Y_log, 1, max)
  trough_draws <- apply(Y_log, 1, min)
  peak_qs      <- quantile(peak_draws,   c(.025, .5, .975), na.rm = TRUE)
  trough_qs    <- quantile(trough_draws, c(.025, .5, .975), na.rm = TRUE)
  
  p_season <- plot_phase_curve(
    df = sine_df,
    title = "Estimated seasonal variation in ESBL-E acquisition risk",
    subtitle = "Sinusoidal seasonal function (median hazard ratio with 95% CrI)",
    ylab = NULL,
    rainy = rainy,
    ref = "median",
    show_peak_trough = TRUE,
    show_hr1 = TRUE,
    peak_qs = peak_qs,
    trough_qs = trough_qs
  )
}

p_season

#----------------------------------------------------------------------------------------------
# PLOT RATES OVER TIME
#----------------------------------------------------------------------------------------------

# Ensure row index matches Stan input order
data_fit <- data_fit %>% mutate(stan_index = row_number())

# Summarise a vector-valued Stan var "[n]" to median/CrI per n
summarise_stan_var <- function(D, varname, probs = c(0.05, 0.50, 0.95)) {
  # Pull all columns matching varname[<index>]
  M <- as.matrix(D[, grep(paste0("^", varname, "\\["), names(D)), drop = FALSE])
  #M <- get_mat(D, regex = paste0("^", trace_pars, "\\["))
  idx <- as.integer(gsub(sprintf("^%s\\[(\\d+)\\]$", varname), "\\1", colnames(M)))
  
  # Fast column-wise quantiles
  qs <- matrixStats::colQuantiles(M, probs = probs, na.rm = TRUE)
  tibble::tibble(
    index = idx,
    !!paste0(varname, "_q05")    := qs[, 1],
    !!paste0(varname, "_median") := qs[, 2],
    !!paste0(varname, "_q95")    := qs[, 3]
  ) |> arrange(index)
}

# Summarise acquisition, decolonisation, and y_hat 
s12 <- summarise_stan_var(D, "log_lambda_1_2_out")
s21 <- summarise_stan_var(D, "log_lambda_2_1_out")
yhat <- summarise_stan_var(D, "Y_hat_1_2_out")  # assuming y_hat_out is on probability scale already

if(data_source == "simulated"){
data_fit$date.use = data_fit$Date
data_fit$intervention_date = data_fit$Intervention_Date1
data_fit$intervention_date2 = data_fit$Intervention_Date2
data_fit$intervention_village = data_fit$Intervention
data_fit$round = data_fit$Round
}

# Merge with data (index should align with N rows in data_fit) 
df_all <- data_fit %>%
  mutate(index = row_number()) %>%
  left_join(s12, by = "index") %>%
  left_join(s21, by = "index") %>%
  #left_join(yhat, by = "index") %>%
  mutate(
    lambda_1_2       = exp(log_lambda_1_2_out_median),
    lambda_1_2_lower = exp(log_lambda_1_2_out_q05),
    lambda_1_2_upper = exp(log_lambda_1_2_out_q95),
    lambda_2_1       = exp(log_lambda_2_1_out_median),
    lambda_2_1_lower = exp(log_lambda_2_1_out_q05),
    lambda_2_1_upper = exp(log_lambda_2_1_out_q95),
    data.use = as.Date(date.use, format = "%Y-%m-%d"),
    Date = as.Date(date.use, format = "%Y-%m-%d"),
    Intervention = factor(intervention_village, levels=c(0,1), labels=c("No intervention", "Intervention")),
    Week = floor_date(Date, unit = "week"),
    Intervention_Date1 = as.Date(intervention_date,  origin = "1970-01-01"),
    Intervention_Date2 = as.Date(intervention_date2, origin = "1970-01-01"),
    Time_Period = case_when(
      intervention_village == 1 & date.use >= Intervention_Date1 & date.use < Intervention_Date2 ~ "Post-Intervention 1",
      intervention_village == 1 & date.use >= Intervention_Date2 ~ "Post-Intervention 2",
      intervention_village == 1 & date.use < Intervention_Date1 ~ "Pre-Intervention",
      intervention_village == 0 ~ "Control"
    ),
    Time_Period = factor(Time_Period, levels = c("Control","Pre-Intervention","Post-Intervention 1","Post-Intervention 2"))
  )

lambda_week <- df_all %>%
  group_by(Time_Period, Intervention, round, Week) %>%
  summarise(
    Mean_lambda_1_2 = mean(lambda_1_2, na.rm = TRUE),
    Mean_lambda_2_1 = mean(lambda_2_1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("Mean_lambda_"),
    names_to = "Rate",
    values_to = "Mean"
  ) %>%
  mutate(
    Rate = factor(Rate,
                  levels = c("Mean_lambda_1_2", "Mean_lambda_2_1"),
                  labels = c("Acquisition rate (λ12)", "Decolonisation rate (λ21)")),
    round = factor(round, levels = c(1,2,3,4), labels = c("3-months\npre-intervention", "start intervention",
                                                          "3-months\npost-intervention", "9-months\npost-intervention"))
  )

lambda_summary_stats <- df_all %>%
  group_by(Time_Period, round) %>%
  summarise(
    Acquisition    = sprintf("%.3f (%.3f–%.3f)",
                             median(lambda_1_2, na.rm = TRUE),
                             quantile(lambda_1_2, 0.25, na.rm = TRUE),
                             quantile(lambda_1_2, 0.75, na.rm = TRUE)),
    Decolonisation = sprintf("%.3f (%.3f–%.3f)",
                             median(lambda_2_1, na.rm = TRUE),
                             quantile(lambda_2_1, 0.25, na.rm = TRUE),
                             quantile(lambda_2_1, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Time_Period, round) #%>%
  #mutate(
  #  round = factor(round, levels = c(1,2,3,4), labels = c("3-months\npre-intervention", "start intervention",
  #                                                        "3-months\npost-intervention", "9-months\npre-intervention"))
  #)

lambda_summary_stats2 <- df_all %>%
  group_by(Time_Period) %>%
  summarise(
    Acquisition    = sprintf("%.3f (%.3f–%.3f)",
                             median(lambda_1_2, na.rm = TRUE),
                             quantile(lambda_1_2, 0.25, na.rm = TRUE),
                             quantile(lambda_1_2, 0.75, na.rm = TRUE)),
    Decolonisation = sprintf("%.3f (%.3f–%.3f)",
                             median(lambda_2_1, na.rm = TRUE),
                             quantile(lambda_2_1, 0.25, na.rm = TRUE),
                             quantile(lambda_2_1, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(Time_Period)

ft2 = flextable(lambda_summary_stats)
ft2 = autofit(ft2)
ft2 = set_caption(ft2, "Model estimates for acquisition and decolonisation - per round")
ft2

ft3 = flextable(lambda_summary_stats2)
ft3 = autofit(ft3)
ft3 = set_caption(ft3, "Model estimates for acquisition and decolonisation - per intervention period")
ft3

# plot acquisition & decolonisation by time period 
p_period <- ggplot(lambda_week %>% filter(round!="3-months\npre-intervention"), aes(x = Time_Period, y = Mean, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap(~ Rate, scales = "free_y") +
  labs(title = "Weekly mean rates by intervention period",
       x = "Time Period", y = expression(lambda ~ "(per day)")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_brewer(palette="Set1")
p_period

p_round <- ggplot(lambda_week %>% filter(round!="3-months\npre-intervention"), aes(x = factor(round), y = Mean, fill = Intervention)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap(~ Rate, scales = "free_y") +
  labs(
    title = "Intervention effect (intention to treat)",
    subtitle = expression("Daily transition rate - weekly average (" * lambda * ")"),
    y =  expression(lambda ~ "(per day)"),
    x = "")+
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)  # centered
  ) +
  scale_fill_brewer(palette = "Set1")
p_round

# Yhat
interval_df <- tibble(
  index      = 1:num_intervals,
  start_num  = interval_starts,
  end_num    = pmin(interval_ends, as.numeric(global_interval_end)),  # cap last partial interval
  mid_num    = X_midpoints
) %>%
  mutate(
    start_date = as.Date(round(start_num), origin = "1970-01-01"),
    end_date   = as.Date(round(end_num),   origin = "1970-01-01"),
    mid_date   = as.Date(round(mid_num),   origin = "1970-01-01"),
    year       = year(mid_date),
    month      = month(mid_date),
    month_lab  = format(mid_date, "%b %Y")  # e.g., "Jan 2023"
  ) %>%
  select(index, start_date, end_date, mid_date, year, month, month_lab)

#---------------------------
# Link yhat to dates
#---------------------------
yhat_linked <- yhat %>%
  left_join(interval_df, by = "index")

# Quick check
head(yhat_linked)

#-------------------------------------------------
# Plot yhat by interval midpoints (ribbon)
#-------------------------------------------------
p_mid <- ggplot(yhat_linked, aes(mid_date, Y_hat_1_2_out_median)) +
  geom_ribbon(aes(ymin = Y_hat_1_2_out_q05, ymax = Y_hat_1_2_out_q95), alpha = 0.25) +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(x = "Interval midpoint", y = expression(hat(Y)[1%->%2]),
       title = "Seasonal effect by 28-day intervals (median & 95% CrI)",
       subtitle = "Based on yhat") +
  theme_minimal() +
  #ylim(-5,0.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_mid

df_step <- yhat_linked %>%
  transmute(
    start = start_date,
    end   = end_date,
    lo    = Y_hat_1_2_out_q05,
    med   = Y_hat_1_2_out_median,
    hi    = Y_hat_1_2_out_q95
  )

p_step = ggplot(df_step) +
  # shaded 90% CrI per interval
  geom_rect(aes(xmin = start, xmax = end, ymin = lo, ymax = hi),
            alpha = 0.25) +
  # median as a step (flat within each interval)
  geom_segment(aes(x = start, xend = end, y = med, yend = med), linewidth = 0.7) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", expand = c(0, 2)) +
  labs(x = "Interval", y = expression(hat(Y)[1 %->% 2]),
       title = "Seasonal effect by 28-day interval (median with shaded 95% CrI)",
       subtitle = "Based on yhat") +
  theme_minimal() +
  ylim(-2,1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p_step

#-----------------------------------------------------------------------
# Plot caterpillar
#-----------------------------------------------------------------------
# choose which parameters to show

# Spline
#--------------------------------------------------------------------
if (trimws(seasonality) == "spline_") {
  keep <- c(
    "q_1_2_base","q_2_1_base",
    "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
    "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
    "sigma","sigma_u","a_1_2"
  )
} else { 
# Sine
#--------------------------------------------------------------------
  keep <- c("q_1_2_base","q_2_1_base",
                "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                "beta_int1_1","beta_int1_2","a1","phi")
}

keep = trace_pars
columns = ifelse(seasonality =="spline_", 34, 12)
dd_sel = D[, which(names(D)%in%keep)] # Select variables of interest

# Turn into matrix
#Xmat <- as.matrix(dd[, c(1:36), drop = FALSE])  # draws × params (numeric)

# summarize (med, 80% & 95% CrI) 
qs <- apply(dd_sel, 2, function(x) quantile(x, c(0.025, 0.10, 0.50, 0.90, 0.975)))
sumdf <- tibble(
  variable = colnames(dd_sel),
  q025 = qs[1,], q10 = qs[2,], q50 = qs[3,], q90 = qs[4,], q975 = qs[5,]
) %>%
  mutate(
    param = sub("\\[.*$", "", variable),
    index = suppressWarnings(as.integer(sub("^.*\\[(\\d+)\\].*$", "\\1", variable))),
    label = ifelse(is.na(index), param, paste0(param, "[", index, "]"))
  )

# overlay true values: supply a named numeric vector 'truth_named' matching Stan names
# True parameters
true_params <- list(
  beta_1_2_age = 0.01,         # Age effect on transition from state 1 to 2 (acquisition)
  beta_2_1_age = -0.03,        # Age effect on transition from state 2 to 1 (clearance)
  beta_1_2_sexe = 0.03,        # Sex effect on transition from 1 to 2
  beta_2_1_sexe = -0.02,       # Sex effect on transition from 2 to 1
  #u_raw = seq(-1, 1, length.out = H),  # Raw household-level random effects (length H)
  #sigma_u = 0.3,               # Standard deviation of household random effects
  beta_int1_1 = -1.5,          # Effect of the 1st intervention on transition rate λ₁₂
  beta_int1_2 = -1.6,          # Effect of the 2nd intervention on transition rate λ₁₂
  beta_int2_1 = 1.0,           # Effect of the 1st intervention on transition rate λ₂₁
  beta_int2_2 = 0.15,          # Effect of the 2nd intervention on transition rate λ₂₁
  #sigma_q_1_2 = 0.3,           # Standard deviation of baseline λ₁₂
  #sigma_q_2_1 = 0.3,           # Standard deviation of baseline λ₂₁
  q_sum = 0.02,                # Sum of baseline transition rates (λ₁₂ + λ₂₁)
  alpha = 0.5,                 # Mixing parameter between global and individual hazard contributions
  a1 = ifelse(scen =="seasonality", 1.2,0),                    # Amplitude of seasonal effect on λ₁₂
  phi = ifelse(scen =="seasonality", pi/2,0)                   # Phase shift of seasonal effect
)
true_params

# Convert truths; keep only names present in the 'base' column
truth_named <- unlist(true_params)
truth_named <- truth_named[names(truth_named) %in% unique(sumdf$variable)]

# Attach truth by base (recycled to all indices of that base)
sumdf <- sumdf %>%
  dplyr::mutate(truth = dplyr::if_else(variable %in% names(truth_named),
                                       truth_named[variable], NA_real_))

# which params considered "seasonal" parameters (for plotting them seperately from other parameters)
season_params <- c("a_1_2[1]","a_1_2[2]","a_1_2[3]","a_1_2[4]","a_1_2[5]","a_1_2[6]","a_1_2[7]","a_1_2[8]",
                   "a_raw_1_2[1]","a_raw_1_2[2]","a_raw_1_2[3]","a_raw_1_2[4]","a_raw_1_2[5]","a_raw_1_2[6]","a_raw_1_2[7]","a_raw_1_2[8]",
                   "a_raw_1_2_free[1]","a_raw_1_2_free[2]","a_raw_1_2_free[3]","a_raw_1_2_free[4]","a_raw_1_2_free[5]","a_raw_1_2_free[6]","a_raw_1_2_free[7]","a_raw_1_2_free[8]",
                  "log_tau_raw_1_2", "phi", "a1")

sumdf_plot <- sumdf %>%
  mutate(
    group = ifelse(variable %in% season_params, "Seasonal (λ12)", "Non-seasonal")
  ) %>%
  group_by(group, variable) %>%
  arrange(q50, .by_group = TRUE) %>%
  mutate(label_ord = factor(label, levels = label)) %>%
  ungroup()

p_cat <- ggplot(sumdf_plot, aes(x = q50, y = label_ord)) +
  geom_linerange(aes(xmin = q025, xmax = q975), linewidth = 0.6, alpha = 0.75) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, col="red", linetype=2) +
  #geom_point(aes(x = truth), colour="red", shape = 4, stroke = 1, size = 2.2, na.rm = TRUE) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    x = "Posterior (median with 95% CrI) log-scale",
    y = NULL,
    title = "Caterpillar plot of posterior summaries"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

p_cat

p_cat_exp <- ggplot(sumdf_plot%>%filter(group!="Seasonal (λ12)"), aes(x = exp(q50), y = label_ord)) +
  geom_linerange(aes(xmin = exp(q025), xmax = exp(q975)), linewidth = 0.6, alpha = 0.75) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 1, col="red", linetype=2) +
#  geom_point(aes(x = truth), colour="red", shape = 4, stroke = 1, size = 2.2, na.rm = TRUE) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    x = "Posterior (median with 95% CrI) normal-scale",
    y = NULL,
    title = "Caterpillar plot of posterior summaries"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

p_cat_exp

#----------------------------------------------------------------
# FOR PAPER: Acquisitions per 1000 population days
#----------------------------------------------------------------

#-----------------------------------
# Set denominator
#-----------------------------------

use_uncolonised_pt <- FALSE   # TRUE = person-days-at-risk, FALSE = all person-days

# ===============================
# Calendar-aligned rates + DiD
# ===============================

# This approach does ignore those intervals which have pre- and post-intervention both falling
# in the same interval, as based on weighted mean lambda12 averaged over the full intervention length
# Did this as storing seperate loglambdas for each period adds a lot to memory of the stan output file

#-------------------------------------------------------
# Extract draws for log lambda12 per observation (N rows)
#-------------------------------------------------------
draws <- rstan::extract(fit, pars = "log_lambda_1_2_out")$log_lambda_1_2_out  # [draw, N]
stopifnot(length(dim(draws)) == 2)
n_draws <- dim(draws)[1]
N        <- dim(draws)[2]

#-------------------------------------------------------
# Build row-level dataframe from stan_data_fit
#    (one row per observation interval)
#-------------------------------------------------------
df_rows <- tibble(
  row_id         = seq_len(N),
  date_use       = stan_data_fit$date_use,            # numeric day index
  observed_state = stan_data_fit$observed_state,      # 1 = uncolonised, 2 = colonised
  person         = stan_data_fit$menage_id_member,    # unit id
  group_flag     = stan_data_fit$intervention,        # 1 = treated, 0 = control
  t_star1        = stan_data_fit$intervention_date,   # first intervention date (if treated)
  t_star2        = stan_data_fit$intervention_date2   # second intervention date (if treated)
) %>%
  arrange(person, row_id) %>%
  group_by(person) %>%
  mutate(
    # interval length (days) and midpoint
    dt_days = date_use - dplyr::lag(date_use),
    mid_date = date_use - dt_days/2,
    
    # simple "at-risk" indicator at start (optional; not used by default below)
    p_U_start = case_when(
      is.na(lag(observed_state)) ~ NA_real_,
      lag(observed_state) == 1   ~ 1.0,
      lag(observed_state) == 2   ~ 0.0,
      TRUE                       ~ NA_real_
    )
  ) %>%
  ungroup()

# Keep only valid intervals (positive length & defined at-risk indicator if needed)
df_valid <- df_rows %>%
  filter(!is.na(dt_days), dt_days > 0)

#-------------------------------------------------------
# Define global calendar windows from treated units
#    (so controls can be aligned in calendar time)
#-------------------------------------------------------
treated <- df_valid %>% filter(group_flag == 1)

rng_or_na <- function(x) {
  if (length(x) == 0 || all(!is.finite(x))) c(NA_real_, NA_real_) else range(x, na.rm = TRUE)
}

win_pre <- rng_or_na(treated$mid_date[treated$mid_date < treated$t_star1])
win_post1 <- rng_or_na(treated$mid_date[
  treated$mid_date >= treated$t_star1 &
    (is.na(treated$t_star2) | treated$mid_date < treated$t_star2)
])
win_post2 <- rng_or_na(treated$mid_date[
  !is.na(treated$t_star2) & treated$mid_date >= treated$t_star2
])

in_win <- function(x, rng) {
  is.finite(rng[1]) & is.finite(rng[2]) & x >= rng[1] & x <= rng[2]
}

df_valid <- df_valid %>%
  mutate(
    Phase_cal = case_when(
      in_win(mid_date, win_post2) ~ "Post-Intervention 2",
      in_win(mid_date, win_post1) ~ "Post-Intervention 1",
      in_win(mid_date, win_pre)   ~ "Pre-Intervention",
      TRUE                        ~ NA_character_
    ),
    Group = ifelse(group_flag == 1, "Intervention", "Control"),
    Phase_cal = factor(Phase_cal, levels = c("Pre-Intervention","Post-Intervention 1","Post-Intervention 2"))
  )

# drop rows outside any window
df_valid <- df_valid %>% filter(!is.na(Phase_cal))

#-------------------------------------------------------
# Rates per 1,000 person-days (per draw -> CrIs)
#    Using ALL person-time (not only uncolonised)
#-------------------------------------------------------
use_uncolonised_pt <- FALSE  # set TRUE to switch denominator to uncolonised time approximation

rate_one_draw <- function(d) {
  lambda12_d <- exp(draws[d, ])  # per-day hazard for draw d (length N)
  
  tmp <- df_valid %>%
    mutate(
      lambda12 = lambda12_d[row_id],
      pt       = if (use_uncolonised_pt) p_U_start * dt_days else dt_days,
      exp_acq  = pt * lambda12
    )
  
  tmp %>%
    group_by(Group, Phase_cal) %>%
    summarise(
      events  = sum(exp_acq, na.rm = TRUE),
      pt_days = sum(pt,      na.rm = TRUE),
      ID12_per1000 = 1000 * events / pt_days,
      .groups = "drop"
    ) %>%
    mutate(draw = d)
}

# (can thin draws for speed/memory, e.g. draw_ids <- sample(seq_len(n_draws), 1000))
draw_ids <- seq_len(n_draws)

rates_draws <- bind_rows(lapply(draw_ids, rate_one_draw))

rates_cri <- rates_draws %>%
  group_by(Group, Phase_cal) %>%
  summarise(
    Median = median(ID12_per1000, na.rm = TRUE),
    Lower  = quantile(ID12_per1000, 0.025, na.rm = TRUE),
    Upper  = quantile(ID12_per1000, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(Group, levels = c("Control","Intervention")), Phase_cal)

rates_cri$phase_label = factor(rates_cri$Phase_cal, labels=c("Pre-intervention",
                                                             "3-months\npost-intervention",
                                                             "9-months\npost-intervention"))

rates_draws$phase_label = factor(rates_draws$Phase_cal, labels=c("Pre-intervention",
                                                             "3-months\npost-intervention",
                                                             "9-months\npost-intervention"))

#------------------------------------------------------------------------
# Proper DiD on log scale (per draw -> CrIs)
#    DiD = (T_post - T_pre) - (C_post - C_pre) on log-rate
#------------------------------------------------------------------------
did_one_draw <- function(d) {
  log_lambda_d <- draws[d, ]  # average log hazard per row for draw d
  
  tmp <- df_valid %>%
    mutate(
      log_rate = log_lambda_d[row_id],  # already a log-rate per day
      w        = dt_days                # person-time weight
    ) %>%
    group_by(Group, Phase_cal) %>%
    summarise(
      log_rate = sum(log_rate * w, na.rm = TRUE) / sum(w, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = Group, values_from = log_rate)
  
  getI <- function(ph) tmp$Intervention[tmp$Phase_cal == ph]
  getC <- function(ph) tmp$Control[tmp$Phase_cal == ph]
  
  pre <- "Pre-Intervention"; p1 <- "Post-Intervention 1"; p2 <- "Post-Intervention 2"
  
  did_p1 <- (getI(p1) - getI(pre)) - (getC(p1) - getC(pre))
  did_p2 <- (getI(p2) - getI(pre)) - (getC(p2) - getC(pre))
  
  tibble(draw = d, DiD_post1_log = did_p1, DiD_post2_log = did_p2)
}

# DRAWS
#---------------------------------------------------------------------------------------
did_draws <- bind_rows(lapply(draw_ids, did_one_draw))

# geometric mean on HR scale
#---------------------------------------------------------------------------------------
did_draws <- did_draws %>%
  mutate(
    DiD_overall_log = 0.5 * (DiD_post1_log + DiD_post2_log) 
  )

DiD_summary <- did_draws %>%
  summarise(
    post1_log_med = median(DiD_post1_log, na.rm=TRUE),
    post1_log_lo  = quantile(DiD_post1_log, 0.025, na.rm=TRUE),
    post1_log_hi  = quantile(DiD_post1_log, 0.975, na.rm=TRUE),
    
    post2_log_med = median(DiD_post2_log, na.rm=TRUE),
    post2_log_lo  = quantile(DiD_post2_log, 0.025, na.rm=TRUE),
    post2_log_hi  = quantile(DiD_post2_log, 0.975, na.rm=TRUE),
    
    overall_log_med = median(DiD_overall_log, na.rm=TRUE),
    overall_log_lo  = quantile(DiD_overall_log, 0.025, na.rm=TRUE),
    overall_log_hi  = quantile(DiD_overall_log, 0.975, na.rm=TRUE)
  ) %>%
  mutate(
    post1_HR_med = exp(post1_log_med),
    post1_HR_lo  = exp(post1_log_lo),
    post1_HR_hi  = exp(post1_log_hi),
    
    post2_HR_med = exp(post2_log_med),
    post2_HR_lo  = exp(post2_log_lo),
    post2_HR_hi  = exp(post2_log_hi),
    
    overall_HR_med = exp(overall_log_med),
    overall_HR_lo  = exp(overall_log_lo),
    overall_HR_hi  = exp(overall_log_hi)
  )


# PRE-INTERVENTION HR between groups 
# This assumes `rates_draws` has one row per posterior draw with:
#   - a draw identifier column named `draw` (if not, replace below with whatever identifies draws)
#   - columns: phase_label == "Pre-intervention", Group (e.g., "Intervention","Control"),
#              and a rate metric like ID12_per1000 (any positive rate works; scale cancels in ratio)
pre_gap <- rates_draws %>%
   filter(phase_label == "Pre-intervention") %>%
   select(draw, Group, ID12_per1000) %>%
   tidyr::pivot_wider(names_from = Group, values_from = ID12_per1000) %>%
   # rename to match your actual group labels if different:
   # e.g., rename(Intervention = TreatmentArmName, Control = ControlArmName)
   mutate(
     pre_HR  = Intervention / Control,
     pre_log = log(pre_HR)
   )
 
 pre_summary <- pre_gap %>%
   summarise(
     pre_HR_med = median(pre_HR, na.rm=TRUE),
     pre_HR_lo  = quantile(pre_HR, 0.025, na.rm=TRUE),
     pre_HR_hi  = quantile(pre_HR, 0.975, na.rm=TRUE)
   )

pre_gap <- rates_draws %>%
  dplyr::filter(phase_label == "Pre-intervention") %>%
  dplyr::select(draw, Group, ID12_per1000) %>%
  tidyr::pivot_wider(names_from = Group, values_from = ID12_per1000)

pre_gap <- pre_gap %>%
  dplyr::mutate(pre_HR = .data$Intervention / .data$Control,
                pre_log = log(pre_HR))

# Helper to get med + central 50/80/95% intervals from a vector on the HR scale
summ_ci <- function(x) {
  tibble::tibble(
    HR      = stats::median(x, na.rm = TRUE),
    HR_lo50 = stats::quantile(x, 0.25,  na.rm = TRUE),
    HR_hi50 = stats::quantile(x, 0.75,  na.rm = TRUE),
    HR_lo80 = stats::quantile(x, 0.10,  na.rm = TRUE),
    HR_hi80 = stats::quantile(x, 0.90,  na.rm = TRUE),
    HR_lo   = stats::quantile(x, 0.025, na.rm = TRUE),
    HR_hi   = stats::quantile(x, 0.975, na.rm = TRUE)
  )
}

# Compute HR draws for each phase
post1_hr <- exp(did_draws$DiD_post1_log)
post2_hr <- exp(did_draws$DiD_post2_log)
overall_hr <- exp(did_draws$DiD_overall_log)
pre_hr <- pre_gap$pre_HR

# Assemble the summary table in desired order & labels
summ <- dplyr::bind_rows(
  dplyr::bind_cols(label = "Pre-intervention HR (between groups)", summ_ci(pre_hr)),
  dplyr::bind_cols(label = "Post-Intervention 1 (DiD)",           summ_ci(post1_hr)),
  dplyr::bind_cols(label = "Post-Intervention 2 (DiD)",           summ_ci(post2_hr)),
  dplyr::bind_cols(label = "Overall (Post1 & Post2)",             summ_ci(overall_hr))
) %>%
  dplyr::mutate(
    label = factor(
      label,
      levels = c(
        "Pre-intervention HR (between groups)",
        "Post-Intervention 1 (DiD)",
        "Post-Intervention 2 (DiD)",
        "Overall (Post1 & Post2)"
      )
    )
  )

# Build long data for plotting Post1, Post2, OVERALL, and PRE
#---------------------------------------------------------------------------
did_long <- did_draws %>%
  summarise(
    Post1   = median(exp(DiD_post1_log), na.rm=TRUE),
    Post2   = median(exp(DiD_post2_log), na.rm=TRUE),
    Overall = median(exp(DiD_overall_log), na.rm=TRUE),
    
    Post1_lo   = quantile(exp(DiD_post1_log), 0.025, na.rm=TRUE),
    Post1_hi   = quantile(exp(DiD_post1_log), 0.975, na.rm=TRUE),
    Post2_lo   = quantile(exp(DiD_post2_log), 0.025, na.rm=TRUE),
    Post2_hi   = quantile(exp(DiD_post2_log), 0.975, na.rm=TRUE),
    Overall_lo = quantile(exp(DiD_overall_log), 0.025, na.rm=TRUE),
    Overall_hi = quantile(exp(DiD_overall_log), 0.975, na.rm=TRUE)
  ) %>%
  tidyr::pivot_longer(everything()) %>%
  mutate(
    phase = dplyr::case_when(
      grepl("^Post1", name)   ~ "Post-Intervention 1",
      grepl("^Post2", name)   ~ "Post-Intervention 2",
      grepl("^Overall", name) ~ "Overall (Post1 & Post2)",
      TRUE ~ NA_character_
    ),
    stat  = dplyr::case_when(
      grepl("_lo$", name) ~ "lo",
      grepl("_hi$", name) ~ "hi",
      TRUE                ~ "med"
    )
  ) %>%
  select(phase, stat, value) %>%
  tidyr::pivot_wider(names_from = stat, values_from = value) %>%
  filter(!is.na(phase))

# Bind PRE row
did_long <- bind_rows(
  did_long,
  tibble::tibble(
    phase = "Pre-intervention (gap)",
    med = pre_summary$pre_HR_med,
    lo  = pre_summary$pre_HR_lo,
    hi  = pre_summary$pre_HR_hi
  )
) %>%
  mutate(
    phase = factor(
      phase,
      levels = c(
        "Pre-intervention (gap)",
        "Post-Intervention 1",
        "Post-Intervention 2",
        "Overall (Post1 & Post2)"
      )
    ),
    phase_labels = factor(phase,
                          levels = c(
                            "Pre-intervention (gap)",
                            "Post-Intervention 1",
                            "Post-Intervention 2",
                            "Overall (Post1 & Post2)"
                          ),
                          labels=c("Pre-intervention",
                                   "3-months\npost-intervention",
                                   "9-months\npost-intervention",
                                   "Overall"))
                          )                     


#-------------------------------------------------------
# Plots
#-------------------------------------------------------

# Incidence per 1,000 PD by Group and intervention Phase (CrIs)
pd <- position_dodge(width = 0.5)
p_rates <- ggplot(rates_cri, aes(x = phase_label, y = Median, color = Group, group = Group)) +
  coord_flip() +
  geom_point(position = pd, size = 6, shape = 18) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = pd, width = 0.1, size = 0.8) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Acquisition incidence by intervention period",
    subtitle = if (use_uncolonised_pt)
      "Per 1,000 uncolonised person-days"
    else
      "Per 1,000 person-days",
    y = "Incidence density (per 1,000 person-days)",
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
p_rates

#-----------------------------------------------------------------
# CURRENTLY USED IN OTN PRESENTION
#-----------------------------------------------------------------
# as boxplot
col_set = wes_palette("Darjeeling1", n = 5)
col_set

col_set2 = wes_palette("GrandBudapest1", n = 4)
col_set2

rates_draws$Group=factor(rates_draws$Group, levels=c("Intervention", "Control"))

# Add an "Overall" phase (pooled across all phases), keeping Group separate
rates_draws_overall <- bind_rows(
  rates_draws,
  rates_draws %>% mutate(phase_label = "Overall")
) %>%
  mutate(
    phase_label = fct_relevel(phase_label, "Overall", after = 0), # put Overall at the top (with coord_flip)
    Group = factor(Group, levels = c("Intervention", "Control"))
  )

rates_draws_overall$phase_label=factor(rates_draws_overall$phase_label, levels=c("Pre-intervention",
         "3-months\npost-intervention",
         "9-months\npost-intervention",
         "Overall"))

# Switch to rates_draws_overall if want to plot overall too
base_size <- 16
p_rates_boxplot <- ggplot(rates_draws, aes(x = phase_label, y = ID12_per1000, fill = Group)) +
  coord_flip()+
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  labs(
    title = "ESBL-E acquisition incidence by intervention period",
    subtitle = if (use_uncolonised_pt) "Per 1,000 uncolonised person-days" else "Per 1,000 person-days",
    y = "Incidence density (per 1,000 person-days)",
    x = ""
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 150, 25), limits = c(0, 150)) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
  ) +
  scale_fill_manual(values = c(col_set[3], col_set[2])) +
  theme(
    plot.title    = element_text(size = base_size + 4, face = "bold"),
    plot.subtitle = element_text(size = base_size + 2),
    axis.title    = element_text(size = base_size + 2),
    axis.text     = element_text(size = base_size + 2),
    plot.caption  = element_text(size = base_size - 3),
    legend.text   = element_text(size = 16),
    legend.title  = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor = element_line(linetype = "dotted", color = "grey85")
  )

p_rates_boxplot
#-----------------------------------------------------------------
# CURRENTLY USED IN OTN PRESENTATION
#-----------------------------------------------------------------

did_long <- did_long %>%
  mutate(hr_label = sprintf("HR = %.2f [%.2f–%.2f]", med, lo, hi))

p_did <-  ggplot(did_long, aes(x = phase_labels, y = med)) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.1, size = 1.2,
                color = c("black")) +
  geom_point(size = 10, shape = c(18,18,17,18),
             color = c("black", "black","darkblue", "black")) +
  # fixed text position at x-axis value 2.7
  geom_text(aes(y = 2.3, label = hr_label),
            hjust = 0, size = 6) +
  scale_y_log10() +
  scale_y_continuous(breaks = seq(0, 2.25, 0.25),
                     limits = c(0, 3.0)) + # allow space for labels
  labs(
    title = "Hazard ratio of ESBL-acquisition incidence",
    subtitle = "Intervention/Control (median with 95% CrI)",
    x = "",
    y = "Hazard ratio (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title    = element_text(size = 18),
    axis.text     = element_text(size = 18),
    panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor = element_line(linetype = "dotted", color = "grey85"),
    axis.title.y = element_text(size = 18, vjust = 0.5) 
  )

p_did

# DiD HRs (log scale), shown as rectangles
#------------------------------------------------------------------------
if (!exists("title_text")) title_text <- "Difference-in-Differences (HR) with Pre-intervention Check"

p <- ggplot2::ggplot(summ, ggplot2::aes(x = HR, y = label)) +
  ggplot2::geom_rect(
    ggplot2::aes(xmin = HR_lo, xmax = HR_hi,
                 ymin = as.numeric(label) - 0.3,
                 ymax = as.numeric(label) + 0.3,
                 fill = "95% CI"),
    alpha = 0.4, inherit.aes = FALSE
  ) +
  ggplot2::geom_rect(
    ggplot2::aes(xmin = HR_lo80, xmax = HR_hi80,
                 ymin = as.numeric(label) - 0.2,
                 ymax = as.numeric(label) + 0.2,
                 fill = "80% CI"),
    alpha = 0.6, inherit.aes = FALSE
  ) +
  ggplot2::geom_rect(
    ggplot2::aes(xmin = HR_lo50, xmax = HR_hi50,
                 ymin = as.numeric(label) - 0.1,
                 ymax = as.numeric(label) + 0.1,
                 fill = "50% CI"),
    alpha = 0.8, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi), height = 0.2) +
  scale_x_continuous(breaks = seq(0, 3, 0.25), limits=c(0,3), expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_manual(
    name = "Shading Level",
    values = c("50% CI" = "#3182BD", "80% CI" = "#9ECAE1", "95% CI" = "#DEEBF7")
  ) +
  ggplot2::labs(
    title = title_text,
    x = "Hazard ratio",
    y = NULL,
    caption = "Points = posterior median; bars = 95% CrI"
  ) +
  theme_cowplot() +
  theme(
    axis.text.y = ggplot2::element_text(size = 10, hjust = 0),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 9),
    plot.title = ggplot2::element_text(face = "bold")
  )

list(plot = p, table = summ)

# Print plots
p_rates
p_rates_boxplot
p_did

#-------------------------------------------------------
# What to report/use downstream
#-------------------------------------------------------
# - 'rates_cri'  : phase- & group-specific incidence per 1,000 PD (median & 95% CrI)
# - 'DiD_summary': log-scale DiD and HR-scale DiD (median & 95% CrI)
# - 'p_rates' and 'p_did' plots


#--------------------------------------------------------
# NOW USING MODEL PARAMETERS (MORE ACCURATE THAN LAMBDAOUT)
#-----------------------------------------------------------



#-------------------------------------------------------------------------------------------------------------
# STORE OUTPUT
#-------------------------------------------------------------------------------------------------------------

pdf(file = paste0("./Output/Model_results/Model_summaries/", dir,"m_output_", intervention, seasonality, intervention_effect, scen, "_",data_source,".pdf"), width = 15, height = 8)

grid.newpage()
grid.table(diagnostics[,c(1:8)])

grid.newpage()
grid.table(diagnostics[1,c(9:20)])

print(p)

grid.newpage()
grid.table(lambda_summary_stats)

grid.newpage()
grid.table(lambda_summary_stats2)

# Loop over a list of plots
plots <- list(p_period,p_round,p_season, p_mid,p_step, p_cat, p_cat_exp, p_rates,p_rates_boxplot, p_did)
for (plt in plots) {
  print(plt)
}

dev.off()


# Store plots
#-------------------------------------------------------------------------------------

ggsave(
  filename = paste0("./Output/Model_results/Figures/", dir, 
                    "Figure4_Seasonal_term - ", intervention, seasonality, intervention_effect, 
                    scen, "_", data_source, ".png"),
  plot     = p_season,
  width    = 12,   # in inches
  height   = 9,    # in inches
  dpi      = 300,   
  bg = "white"
)

ggsave(
  filename = paste0("./Output/Model_results/Figures/", dir, 
                    "Figure3_p_acq_rates - ", intervention, seasonality, intervention_effect, 
                    scen, "_", data_source, ".png"),
  plot     = p_rates_boxplot,
  width    = 12,   # in inches
  height   = 9,    # in inches
  dpi      = 300,   
  bg = "white"
)


ggsave(
  filename = paste0("./Output/Model_results/Figures/", dir, 
                    "Figure3_p_acq_hr_did - ", intervention, seasonality, intervention_effect, 
                    scen, "_", data_source, ".png"),
  plot     = p_did,
  width    = 12,   # in inches
  height   = 9,    # in inches
  dpi      = 300,   
  bg = "white"
)










