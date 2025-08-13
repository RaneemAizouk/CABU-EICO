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
pacman::p_load(posterior, bayesplot, dplyr, shinystan, ggplot2, tibble, tidyr, lubridate, gridExtra)

#------------------------------------------------------------------------------
# Load data & fit
#------------------------------------------------------------------------------
data_source <- "observed"  # "simulated" or "observed"
intervention <- "Two_step_" #"Two_step_" or "One_step"
scen    <- "seasonality" # "seasonality" or "noseasonality"
seasonality <- "spline_" # "sine_" or "spline_"

scenario = paste0(intervention, seasonality, scen)
  
if (data_source == "simulated") {
  scen_data <- paste0("Simulated_data_", ifelse(scenario == paste0("Two_step_", seasonality, "seasonality"), "seasonality", "noseasonality"))
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
fit <- readRDS(paste0("./Output/Model_results/", dir,scenario, data_source,".rds"))

#------------------------------------------------------------------------------
# Traceplots & ShinyStan
#------------------------------------------------------------------------------
# For spline
if (trimws(seasonality) == "spline_") {
  trace_pars <- c("a_1_2","q_1_2_raw","q_2_1_raw",
                  "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                  "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
                  "sigma","sigma_u","sigma_q_1_2","sigma_q_2_1")
} else { # for sine
trace_pars <- c("a1","phi","q_1_2_raw","q_2_1_raw",
                "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
                "sigma_u","sigma_q_1_2","sigma_q_2_1")
}

trace_pars

p = rstan::traceplot(fit, pars = trace_pars, inc_warmup = FALSE, nrow = 3)
p

#y_rep = fit@par_dims$y_rep

#launch_shinystan(fit)

D <- posterior::as_draws_df(fit)

# PLOT SEASONAL PATTERN
#-----------------------------------------------------------------------

# Inputs
global_interval_start <- ifelse(data_source=="simulated", as.Date("2022-10-01"), as.Date("2022-10-03"))  # observed data start date of the study
global_interval_end   <- as.Date("2024-02-19")  # end date of the study
interval_length <- 28   

# # Convert to numeric
global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric   <- as.numeric(global_interval_end)

# # Number of intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)

# Interval starts and ends
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

# Develop function that sets up the calendar for one year
make_calendar <- function(start_date, end_date, interval_len) {
  gs <- as.numeric(start_date) # Global interval start
  ge <- as.numeric(end_date)  # Global interval end
  n  <- ceiling((ge - gs) / interval_len)
  starts <- seq(gs, gs + (n - 1) * interval_len, by = interval_len)
  mids <- (starts + pmin(starts + interval_len - 1, ge)) / 2
  if ((ge - starts[n]) + 1 < interval_len) mids[n] <- starts[n]
  M1 <- mids[1]; medA <- median(mids - M1)
  #date0 <- as.Date(round(medA + M1), origin = "1970-01-01")
  year_anchor <- 2023  # choose the reference year
  date0 <- as.Date(paste0(year_anchor, "-01-01"))
  list(midpoints = mids, date0 = date0)
}
cal <- make_calendar(global_interval_start, global_interval_end, interval_length)
cal

summarise_draws <- function(mat, probs = c(.05, .5, .95)) {
  qs <- apply(mat, 2, quantile, probs = probs)
  tibble(q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,])
}

plot_phase_curve <- function(df, ylab, title) {
  ggplot(df, aes(Date, q50)) +
    geom_ribbon(aes(ymin = q_low, ymax = q_hi), alpha = 0.25) +
    geom_line(linewidth = 0.9) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
    labs(x = "Calendar time (one-year seasonal cycle from January)",
         y = ylab, title = title) +
    theme_minimal()
}

phase_to_date <- function(phase, date0) date0 + round(phase * 365)

# For spline seasonality
# --------------------------------------------------------------------------
if (trimws(seasonality) == "spline_") {
  
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
    X_mod <- pmin(X_in %% 1, 1 - 1e-12)  # avoid right-open boundary at 1
    K  <- length(knots)
    nb <- K + spline_degree
    ext_knots <- c(knots[(K - spline_degree + 1):K] - 1,
                   knots,
                   knots[1:(spline_degree + 1)] + 1)
    order <- spline_degree + 1
    B <- sapply(1:nb, function(ind) build_b_spline_vec(X_mod, ext_knots, ind, order))
    as.matrix(B)  # n_points × nb
  }
  
  # Build basis once on the phase grid
  B_phase <- make_B_matrix(phase_grid, knots_vec, deg_fit)
  
  vars_a  <- sprintf("a_1_2[%d]", seq_len(nbasis_fit))
  #A_df    <- posterior::as_draws_df(fit, variables = vars_a) # Extracts draws for a12 parameter
  #A_mat   <- as.matrix(A_df[, vars_a, drop = FALSE])   # draws × nbasis_fit
  A_df <- D[,vars_a]
  A_mat   <- as.matrix(A_df)
  stopifnot(ncol(A_mat) == ncol(B_phase))
  
  Y_spline <- A_mat %*% t(B_phase)  # draws × n_phase
  
  qs <- apply(Y_spline, 2, quantile, probs = c(.05, .5, .95))
  spline_df <- tibble(
    phase = phase_grid,
    q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,]
  ) |>
    mutate(Date = cal$date0 + round(phase * 365))
  
  
  
  p_season <- plot_phase_curve(spline_df, "Spline seasonal effect (log scale)",
                               "Periodic B-spline for λ12 (median & CrI)")
  
} else {
# For sine seasonality
#-------------------------------------------------------------------------------
  # Phase grid (avoid t = 1 right-open boundary) 
  phase_grid <- seq(0, 1, length.out = 400)
  
  # expects a1, phi in `fit
  s <- D[, c("a1","phi")]
  #s <- as_draws_df(fit, variables = c("a1","phi"))
  Y_sine <- outer(2*pi*phase_grid, rep(1, nrow(s)))  # angle per phase
  
  # Each draw: a1*sin(2*pi*x + phi)
  sine_mat <- sweep(sin(Y_sine + matrix(rep(s$phi, each = length(phase_grid)),
                                      nrow = length(phase_grid))), 2, s$a1, `*`)
  sine_df <- bind_cols(tibble(phase = phase_grid),
                     summarise_draws(t(sine_mat))) |>
  mutate(Date = phase_to_date(phase, cal$date0))
  
  p_season <- plot_phase_curve(sine_df, "Seasonal effect (λ12)", 
                           "Sinusoidal seasonal pattern (median & CrI)")
}

p_season

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
                  labels = c("Acquisition rate (λ12)", "Decolonisation rate (λ21)"))
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
  arrange(Time_Period, round)

grid.table(lambda_summary_stats)

# plot acquisition & decolonisation by time period 
p_period <- ggplot(lambda_week, aes(x = Time_Period, y = Mean, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap(~ Rate, scales = "free_y") +
  labs(title = "Weekly mean rates by intervention period",
       x = "Time Period", y = "Rate (per day)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_brewer(palette="Set1")
p_period

p_round <- ggplot(lambda_week, aes(x = factor(round), y = Mean, fill =Intervention)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap( ~ Rate, scales = "free_y") +
  labs(title = "Weekly mean rates by sample period",
       x = "Round", y = "Rate (per day)") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 1))+
  scale_fill_brewer(palette="Set1")
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
    "q_1_2_raw","q_2_1_raw",
    "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
    "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
    "sigma","sigma_u","sigma_q_1_2","sigma_q_2_1","a_1_2"
  )
} else { 
# Sine
#--------------------------------------------------------------------
  keep <- c("q_1_2_raw","q_2_1_raw",
                "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
                "sigma_u","sigma_q_1_2","sigma_q_2_1","a1","phi")
}

columns = ifelse(seasonality =="spline_", 36, 15)
dd_sel = D[, c(1:columns)] # Select variables of interest

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
  sigma_u = 0.3,               # Standard deviation of household random effects
  beta_int1_1 = -1.5,          # Effect of the 1st intervention on transition rate λ₁₂
  beta_int1_2 = -1.6,          # Effect of the 2nd intervention on transition rate λ₁₂
  beta_int2_1 = 1.0,           # Effect of the 1st intervention on transition rate λ₂₁
  beta_int2_2 = 0.15,          # Effect of the 2nd intervention on transition rate λ₂₁
  sigma_q_1_2 = 0.3,           # Standard deviation of baseline λ₁₂
  sigma_q_2_1 = 0.3,           # Standard deviation of baseline λ₂₁
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
  geom_point(aes(x = truth), colour="red", shape = 4, stroke = 1, size = 2.2, na.rm = TRUE) +
  facet_wrap(~ group, scales = "free_y") +
  labs(
    x = "Posterior (median with 95% CrI) log-scale",
    y = NULL,
    title = "Caterpillar plot of posterior summaries"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

p_cat

#-------------------------------------------------------------------------------------------------------------
# STORE OUTPUT
#-------------------------------------------------------------------------------------------------------------


pdf(file = paste0("./Output/Figures/", dir,"m_output_", scenario,"_", data_source,".pdf"), width = 15, height = 8)

# Loop over a list of plots
grid.table(lambda_summary_stats)
plots <- list(p, p_period,p_round,p_season, p_mid,p_step, p_cat)
for (plt in plots) {
  print(plt)
}

dev.off()
