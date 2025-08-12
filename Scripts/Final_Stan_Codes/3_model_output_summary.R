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
pacman::pload(posterior, bayesplot, dplyr, shinystan, ggplot2, tibble, tidyr, lubridate, gridExtra)

#------------------------------------------------------------------------------
# Load data & fit
#------------------------------------------------------------------------------
data_source <- "observed"
scenario    <- "Two_step_sine_seasonality"

if (data_source == "simulated") {
  scen_data <- paste0("Simulated_data_", ifelse(scenario == "Two_step_sine_seasonality", "seasonality", "noseasonality"))
  sim_stan_data <- readRDS(paste0("./Data/Simulated_data/", scen_data, "_stan_data.rds"))
  sim_df        <- readRDS(paste0("./Data/Simulated_data/", scen_data, ".rds"))
  data_fit      <- sim_df
} else {
  stan_data_fit <- readRDS("./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")
  data_fit      <- read.csv("./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")
}

fit <- readRDS("./Output/Model_results/Observed_data/Two_step_sine_seasonalityobserved.rds")

#------------------------------------------------------------------------------
# Traceplots & ShinyStan
#------------------------------------------------------------------------------
trace_pars <- c("a1","phi","q_1_2_raw","q_2_1_raw",
                "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
                "sigma_u","sigma_q_1_2","sigma_q_2_1")

rstan::traceplot(fit, pars = trace_pars, inc_warmup = FALSE, nrow = 3)

launch_shinystan(fit)


d <- as_draws_df(fit, variables = c("a1","phi","q_1_2_raw","q_2_1_raw",
                                     "beta_1_2_age","beta_2_1_age","beta_1_2_sexe","beta_2_1_sexe",
                                     "beta_int1_1","beta_int1_2","beta_int2_1","beta_int2_2",
                                     "sigma_u","sigma_q_1_2","sigma_q_2_1"))



# PLOT SEASONAL PATTERN
#-----------------------------------------------------------------------

# Inputs
global_interval_start <- as.Date("2022-10-03")  # start date of the study
global_interval_end   <- as.Date("2024-02-19")  # end date of the study
interval_length <- 28                           # length in days

# Convert to numeric
global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric   <- as.numeric(global_interval_end)

# Number of intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)
max_middle    <- num_intervals - 1

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

# Create scaled X
X <- X_midpoints
X <- X - X[1]                  # start from 0
X <- (X - median(X)) / 365     # center & scale to years
num_data <- length(X)

# Function to convert scaled X back to dates
x_to_date <- function(X_scaled, X_midpoints) {
  M1    <- X_midpoints[1]
  medA  <- median(X_midpoints - M1)
  M_num <- X_scaled * 365 + medA + M1
  as.Date(round(M_num), origin = "1970-01-01")
}

# Example — convert X back to dates
dates_from_X <- x_to_date(X, X_midpoints)

# View results
data.frame(X_scaled = X, Date = dates_from_X)


#------------------------------------------------------------------------------
# Seasonal pattern plot (sinusoidal)
#------------------------------------------------------------------------------

# posterior draws
s <- as_draws_df(fit, variables = c("a1","phi"))

# smooth grid over one seasonal cycle (one year in your X-units)
x <- seq(0, 1, length.out = 200)

# posterior summary of a1*sin(2*pi*x + phi)
qfun <- function(xi) {
  y <- s$a1 * sin(2*pi*xi + s$phi)
  quantile(y, probs = c(0.05, 0.5, 0.95))
}
summ <- sapply(x, qfun)

df <- tibble(
  X = x,
  q05 = summ[1,],
  median = summ[2,],
  q95 = summ[3,]
)

# Map phase X to dates
M1   <- X_midpoints[1]
medA <- median(X_midpoints - M1)
date0 <- as.Date(round(medA + M1), origin = "1970-01-01")  # date where X=0 (median study midpoint)

df <- df %>%
  mutate(Date = date0 + round(X * 365))  # one seasonal year from the median date

# Plot with calendar months on x-axis
p_season = ggplot(df, aes(Date, median)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25) +
  geom_line() +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(x = "Calendar time (one-year seasonal cycle from study median)",
       y = "Seasonal effect (λ12)",
       title = "Sinusoidal seasonal pattern (posterior median & 95% CrI)",
       subtitle = "Based on estimates of phi and a1") +
  theme_minimal()
p_season



# PLOT RATES OVER TIME
#----------------------------------------------------------------------------------------------

# Ensure row index matches Stan input order
data_fit <- data_fit %>% mutate(stan_index = row_number())

# Summarise a vector-valued Stan var "[n]" to median/CrI per n
summarise_stan_var <- function(fit, varname) {
  posterior::as_draws_df(fit) %>%
    posterior::subset_draws(regex = paste0("^", varname, "(\\[|$)")) %>%
    as.data.frame() %>%
    select(.draw, matches(paste0("^", varname, "\\["))) %>%
    pivot_longer(
      cols = - .draw,
      names_to = c("param","index"),
      names_pattern = "^(.*)\\[(\\d+)\\]$",
      values_to  = "value"
    ) %>%
    mutate(index = as.integer(index)) %>%
    group_by(index) %>%
    summarise(
      q05    = quantile(value, 0.05),
      median = quantile(value, 0.50),
      q95    = quantile(value, 0.95),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(varname, "_", .), c("q05","median","q95"))
}

# Summarise acquisition, decolonisation, and y_hat 
s12 <- summarise_stan_var(fit, "log_lambda_1_2_out")
s21 <- summarise_stan_var(fit, "log_lambda_2_1_out")
yhat <- summarise_stan_var(fit, "Y_hat_1_2_out")  # assuming y_hat_out is on probability scale already


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
    Date = as.Date(date.use, format="%Y-%m-%d"),
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
  group_by(Time_Period, round, Week) %>%
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
p <- ggplot(lambda_week, aes(x = Time_Period, y = Mean, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap(~ Rate, scales = "free_y") +
  labs(title = "Weekly mean rates by intervention period",
       x = "Time Period", y = "Rate (per day)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_brewer(palette="Set1")
p

p_round <- ggplot(lambda_week, aes(x = factor(round), y = Mean, fill = factor(round))) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  facet_wrap(~ Rate, scales = "free_y") +
  labs(title = "Weekly mean rates by intervention period",
       x = "Round", y = "Rate (per day)") +
  theme_minimal() +
  theme(legend.position = "none",
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

#---------------------------
# Plot yhat by interval midpoints (ribbon)
#---------------------------
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

#-------------------------------------------------------------------------------------------------------------
# STORE OUTPUT
#-------------------------------------------------------------------------------------------------------------

pdf(file = "./Output/Figures/Observed_data/m_output_two_step_sine_seasonality.pdf", width = 15, height = 8)

# Loop over a list of plots
grid.table(lambda_summary_stats)
plots <- list(p_season,p_mid, p_step, p,p_round)
for (plt in plots) {
  print(plt)
}

dev.off()
