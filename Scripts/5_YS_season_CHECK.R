# Load required libraries
library(dplyr)
library(lubridate)
library(rstan)
set.seed(123)

# ------------------------------
# Define global interval structure
# ------------------------------
global_interval_start <- as.Date("2010-01-01")
global_interval_end <- as.Date("2014-12-31")

global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric <- as.numeric(global_interval_end)

interval_length <- 28

num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)

interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)
interval_ends <- interval_starts + (interval_length - 1)

X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2

if ((global_interval_end_numeric - interval_starts[num_intervals]) + 1 < interval_length) {
  X_midpoints[num_intervals] <- interval_starts[num_intervals]
}

X <- X_midpoints
X <- X - X[1]
X <- (X - median(X)) / 365

interval_starts_date <- as.Date(interval_starts, origin = "1970-01-01")
interval_ends_date <- as.Date(interval_ends, origin = "1970-01-01")

# ------------------------------
# Define individuals and households
# ------------------------------
num_individuals <- 200
num_households <- 40
obs_per_year <- 4
num_years <- 5
obs_per_individual <- obs_per_year * num_years

individuals <- data.frame(
  individual_number = 1:num_individuals,
  HouseID = rep(1:num_households, length.out = num_individuals)
) %>%
  mutate(
    menage_id_member = HouseID * 100 + individual_number,
    household_id = paste0("HH", sprintf("%02d", HouseID))
  )

years <- 2010:2014
base_month_days <- c("01-01", "04-01", "07-01", "10-01")
visit_dates <- as.Date(paste0(rep(years, each = obs_per_year), "-", rep(base_month_days, times = num_years)))

observation_df <- individuals %>%
  slice(rep(1:n(), each = obs_per_individual)) %>%
  mutate(
    visit_date = rep(visit_dates, times = num_individuals),
    day_of_year = yday(visit_date),
    season = case_when(
      day_of_year <= 60 | day_of_year >= 335 ~ "winter",
      day_of_year >= 150 & day_of_year <= 250 ~ "summer",
      TRUE ~ "shoulder"
    ),
    observed_state = case_when(
      season == "winter" ~ rbinom(n(), 1, prob = 0.75),
      season == "summer" ~ rbinom(n(), 1, prob = 0.40),
      season == "shoulder" ~ rbinom(n(), 1, prob = 0.50)
    )
  )

# ------------------------------
# Map visit_date to intervals and spline X
# ------------------------------
get_interval_index <- function(date) {
  idx <- which(date >= interval_starts_date & date <= interval_ends_date)
  if (length(idx) > 0) return(idx[1]) else return(NA)
}

interval_indices <- sapply(observation_df$visit_date, get_interval_index)
observation_df$interval_start <- interval_starts_date[interval_indices]
observation_df$interval_end <- interval_ends_date[interval_indices]
observation_df$X <- X[interval_indices]

# ------------------------------
# Final formatted data
# ------------------------------
observation_df <- observation_df %>%
  arrange(household_id, individual_number, visit_date)

# ------------------------------
# Prepare Stan data
# ------------------------------
N <- nrow(observation_df)
H <- num_households
num_data <- length(X)

X_index <- sapply(observation_df$X, function(x_val) {
  which.min(abs(X - x_val))
})

num_knots <- 4
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1
knots <- seq(min(X), max(X), length.out = num_knots)

date_use <- as.numeric(observation_df$visit_date)
observed_state <- observation_df$observed_state
max_middle <- num_data - 1
intervention_village1 <- rep(0, N)
intervention_village2 <- rep(0, N)

stan_data_small <- list(
  N = N,
  H = H,
  observed_state = as.integer(observed_state),
  menage_id_member = observation_df$menage_id_member,
  HouseID = observation_df$HouseID,
  #age = sample(20:70, N, replace = TRUE),
  round = sample(1:3, N, replace = TRUE),
  #sexe = sample(0:1, N, replace = TRUE),
  global_interval_start = global_interval_start_numeric,
  global_interval_end = global_interval_end_numeric,
  interval_length = interval_length,
  date_use = date_use,
  num_knots = num_knots,
  knots = as.vector(knots),
  spline_degree = spline_degree,
  num_basis = num_basis,
  num_data = num_data,
  X = as.vector(X),
  max_middle = max_middle,
  num_intervals = num_intervals,
  X_index = as.integer(X_index)
)



range(X)  # Should span approximately 5 years (2010–2014) in year units

# transformed parameters {
#   real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
#   real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
#   vector[H] u = u_raw * sigma_u;
#   
#   real tau_1_2 = exp(log_tau_raw_1_2);
#   
#   row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;
#   vector[num_data] Y_hat_1_2 =   to_vector(a_1_2 * B);
#   // vector[num_data] Y_hat_1_2_centered = (Y_hat_1_2 - mean(Y_hat_1_2)) / sd(Y_hat_1_2);
#   
#   real tau_2_1 = exp(log_tau_raw_2_1);
#   
#   row_vector[num_basis] a_2_1 = a_raw_2_1 * tau_2_1;
#   vector[num_data] Y_hat_2_1 = to_vector(a_2_1 * B);
#   
# }

# Prepare data for Stan
stan_model_seasonality_only <- "
functions {
  matrix transition_matrix(real t, real lambda_1_2, real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;
    real exp_term = exp(-total_lambda * t);
    matrix[2,2] P;
    P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;
    P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;
    P[1,2] = (lambda_2_1 / total_lambda) * (1 - exp_term);
    P[2,1] = (lambda_1_2 / total_lambda) * (1 - exp_term);
    return P;
  }

  int wrap_index(int i, int N) {
    return ((i - 1) % N) + 1;
  }

  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order, real period, int num_basis) {
    int M = size(t);
    vector[M] b_spline = rep_vector(0, M);
    vector[M] t_wrapped;

    for (i in 1:M) {
      t_wrapped[i] = t[i] - period * floor(t[i] / period);
    }

    if (order == 1) {
      for (i in 1:M) {
        b_spline[i] = (ext_knots[ind] <= t_wrapped[i]) && (t_wrapped[i] < ext_knots[ind + 1]);
      }
    } else {
      int ind_next = wrap_index(ind + 1, num_basis);
      real denom1 = ext_knots[ind + order - 1] - ext_knots[ind];
      real denom2 = ext_knots[ind + order] - ext_knots[ind + 1];

      vector[M] w1 = (denom1 > 0) ? (t_wrapped - rep_vector(ext_knots[ind], M)) / denom1 : rep_vector(0, M);
      vector[M] w2 = (denom2 > 0) ? 1 - (t_wrapped - rep_vector(ext_knots[ind + 1], M)) / denom2 : rep_vector(0, M);

      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order - 1, period, num_basis)
               + w2 .* build_b_spline(t, ext_knots, ind_next, order - 1, period, num_basis);
    }

    return b_spline;
  }
}

data {
  int<lower=1> N;
  int<lower=1> H;
  array[N] int<lower=0, upper=1> observed_state; // Changed to 0/1
  array[N] int menage_id_member;
  array[N] int HouseID;
  //array[N] int age;
  array[N] int round;
  //array[N] int sexe;
  real global_interval_start;
  real global_interval_end;
  int<lower=1> interval_length;
  array[N] real date_use;
  int<lower=1> num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int<lower=1> num_basis;
  int<lower=1> num_data;
  array[num_data] real X;
  int<lower=1> max_middle;
  int<lower=1> num_intervals;
  array[N] int X_index;
}

transformed data {
  matrix[num_basis, num_data] B;
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2 * spline_degree + num_knots] ext_knots;
  real period = 1;

  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));

  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(
      build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1, period, num_basis)
    );
  }
}

parameters {
  real q_1_2_raw;
  real q_2_1_raw;
  //real beta_1_2_age;
  //real beta_2_1_age;
  //real beta_1_2_sexe;
  //real beta_2_1_sexe;
  //vector[H] u_raw;
  //real<lower=0> sigma_u;
  real<lower=0> sigma_q_1_2;
  real<lower=0> sigma_q_2_1;
  row_vector[num_basis] a_raw_1_2;
  real beta_season_1_2;
  real beta_season_2_1;


  real log_tau_raw_1_2;
  row_vector[num_basis] a_raw_2_1;
 

  real log_tau_raw_2_1;
}


transformed parameters {
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
 // vector[H] u = u_raw * sigma_u;

  real tau_1_2 = exp(log_tau_raw_1_2);
  row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;
  vector[num_data] Y_hat_raw_1_2 = to_vector(a_1_2 * B);
  vector[num_data] Y_hat_1_2 = (Y_hat_raw_1_2 - mean(Y_hat_raw_1_2)) / sd(Y_hat_raw_1_2);  // normalize

  real tau_2_1 = exp(log_tau_raw_2_1);
  row_vector[num_basis] a_2_1 = a_raw_2_1 * tau_2_1;
  vector[num_data] Y_hat_raw_2_1 = to_vector(a_2_1 * B);
  vector[num_data] Y_hat_2_1 = (Y_hat_raw_2_1 - mean(Y_hat_raw_2_1)) / sd(Y_hat_raw_2_1);  // normalize
}

model {
  q_1_2_raw ~ normal(0, 0.3);
  q_2_1_raw ~ normal(0, 0.38);
  //beta_1_2_age ~ normal(0, 1);
  //beta_2_1_age ~ normal(0, 1);
  //beta_1_2_sexe ~ normal(0, 1);
  //beta_2_1_sexe ~ normal(0, 1);
  //u_raw ~ normal(0, 0.5);
  //sigma_u ~ normal(0, 0.5);
  sigma_q_1_2 ~ normal(0, 0.5);
  sigma_q_2_1 ~ normal(0, 0.5);
  a_raw_1_2 ~ normal(0, 1);
  a_raw_2_1 ~ normal(0, 1);
   beta_season_1_2 ~ normal(0, 1);
  beta_season_2_1 ~ normal(0, 1);
  log_tau_raw_1_2 ~ normal(0, 1);
  log_tau_raw_2_1 ~ normal(0, 2);

  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n - 1]) {
      real lambda12 = q_1_2_base  ;
      real lambda21 = q_2_1_base  ;
// + u[HouseID[n]]+ beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n], + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n]+ u[HouseID[n]]
      int i = X_index[n];
      real s12 = Y_hat_1_2[i];
      //real s12 = Y_hat_1_2_centered[i])

      real s21 = Y_hat_2_1[i];

      matrix[2,2] P_total = transition_matrix(date_use[n] - date_use[n - 1],
                                              exp(lambda12 + beta_season_1_2 * s12),
                                              exp(lambda21 + beta_season_2_1 * s21));

      // Adjust indices for 0/1 states
      target += log(P_total[observed_state[n - 1] + 1, observed_state[n] + 1] + 1e-9);
    }
  }
}
generated quantities {
  vector[num_data] Y_hat_1_2_out = Y_hat_1_2;
  vector[num_data] Y_hat_2_1_out = Y_hat_2_1;
  vector[num_data] seasonal_effect_1_2 = beta_season_1_2 * Y_hat_1_2_out;
  vector[num_data] seasonal_effect_2_1 = beta_season_2_1 * Y_hat_2_1_out;
}



"
stan_model <- stan_model(model_code = stan_model_seasonality_only)

fit_test <- sampling(
  stan_model,
  data = stan_data_small,
  iter = 500,
  chains = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 123
)
print(fit_test )

# stan_fit <- sampling(stan_model, data = stan_data, 
#                      iter = 500, warmup = 250, chains = 1, cores = 1,
#                      control = list(adapt_delta = 0.99, max_treedepth = 15),
#                      verbose = TRUE)
# 
df <- data.frame(
  date = as.Date(global_interval_start_numeric + X * 365, origin = "1970-01-01"),
  Y_hat_mean = apply(as.matrix(fit_test, pars = "Y_hat_1_2_out"), 2, mean)
)

ggplot(df, aes(date, Y_hat_mean)) +
  geom_line(color = "steelblue") +
  theme_minimal() +
  labs(title = "Estimated Seasonality Spline (Y_hat_1_2)")




#####
library(shinystan)
shinystan_obj <- as.shinystan(fit_test)
launch_shinystan(shinystan_obj)

#####################


###############
# Assuming you have Stan output already extracted
posterior_samples <- rstan::extract(fit_test)

# Example: spline estimates from Stan model output
Y_hat_1_2_df <- as.data.frame(posterior_samples$Y_hat_1_2_out)
Y_hat_2_1_df <- as.data.frame(posterior_samples$Y_hat_2_1_out)

# Summarize posterior means and 95% intervals
spline_summary <- data.frame(
  time = stan_data_small$X,
  mean_1_2 = apply(Y_hat_1_2_df, 2, mean),
  lwr_1_2 = apply(Y_hat_1_2_df, 2, function(x) quantile(x, 0.025)),
  upr_1_2 = apply(Y_hat_1_2_df, 2, function(x) quantile(x, 0.975)),
  
  mean_2_1 = apply(Y_hat_2_1_df, 2, mean),
  lwr_2_1 = apply(Y_hat_2_1_df, 2, function(x) quantile(x, 0.025)),
  upr_2_1 = apply(Y_hat_2_1_df, 2, function(x) quantile(x, 0.975))
)

# Plot for infection rate spline (Y_hat_1_2)
ggplot(spline_summary, aes(x = time, y = mean_1_2)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = lwr_1_2, ymax = upr_1_2), fill = "red", alpha = 0.2) +
  labs(title = "Estimated Seasonal Effect on Infection Rate (λ₁₂)",
       x = "Time (scaled years)", y = "Log Rate") +
  theme_minimal()
################################
posterior_samples <- rstan::extract(fit_test)

# For a_raw_1_2
a_raw_1_2_samples <- posterior_samples$a_raw_1_2  # matrix: iterations × num_basis

# For tau_1_2
tau_1_2_samples <- exp(posterior_samples$log_tau_raw_1_2)  # since log_tau was sampled
# Mean, SD, 95% Credible Interval for tau
summary_tau <- c(
  mean = mean(tau_1_2_samples),
  sd = sd(tau_1_2_samples),
  quantile(tau_1_2_samples, probs = c(0.025, 0.975))
)

print(summary_tau)

# Plot histogram of tau
hist(tau_1_2_samples, breaks = 50, main = "Posterior of tau_1_2", xlab = "tau_1_2")

# Optional: summary of a_raw_1_2 coefficients
apply(a_raw_1_2_samples, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))

library(dplyr)
library(lubridate)
library(ggplot2)

observation_df %>%
  mutate(
    year_month = floor_date(visit_date, "month")
  ) %>%
  group_by(year_month) %>%
  summarise(
    positive = sum(observed_state == 1),
    total = n(),
    prevalence = positive / total
  ) %>%
  ggplot(aes(x = year_month, y = prevalence)) +
  geom_line(color = "darkred", size = 1.2) +
  geom_point(color = "darkred", size = 2) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Monthly Prevalence of Positive Carriage (Observed State = 1)",
    x = "Date",
    y = "Prevalence"
  )
##################################
library(splines2)
library(ggplot2)
library(reshape2)

# Define parameters
X_plot <- seq(-0.5, 1.5, length.out = 400)  # Covers more than one cycle
degree <- 3
period <- 1

# Internal knots (exclude boundaries)
internal_knots <- seq(0.2, 0.8, length.out = 4)  # 4 internal knots

# Define boundary
boundary <- c(0, 1)

# Generate periodic B-spline basis matrix
B <- bSpline(
  x = X_plot %% period,
  knots = internal_knots,
  Boundary.knots = boundary,
  degree = degree,
  intercept = TRUE,
  periodic = TRUE
)

# Convert for plotting
basis_df <- as.data.frame(B)
basis_df$X <- X_plot
long_df <- melt(basis_df, id.vars = "X")

# Plot the basis functions
ggplot(long_df, aes(x = X, y = value, color = variable)) +
  geom_line(size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Periodic B-Spline Basis Functions (1-Year Cycle)",
    x = "Years",
    y = "Basis Value"
  ) +
  theme(legend.position = "none")
######################################
# Extract seasonal effect samples from Stan fit
seasonal_effect_samples <- rstan::extract(fit_test, pars = "seasonal_effect_1_2")$seasonal_effect_1_2

# Compute mean and 95% CI across posterior draws
seasonal_mean <- apply(seasonal_effect_samples, 2, mean)
seasonal_lower <- apply(seasonal_effect_samples, 2, quantile, probs = 0.025)
seasonal_upper <- apply(seasonal_effect_samples, 2, quantile, probs = 0.975)

# Prepare plotting data
df_seasonal <- data.frame(
  X = X,
  time = seq(as.Date("2010-01-01"), as.Date("2014-12-31"), length.out = length(X)),
  mean = seasonal_mean,
  lower = seasonal_lower,
  upper = seasonal_upper
)
library(lubridate)
df_seasonal <- df_seasonal %>%
  mutate(
    month = month(time),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "Winter",
      month %in% c(6, 7, 8)   ~ "Summer",
      TRUE                   ~ "Other"
    )
  )
library(ggplot2)

ggplot(df_seasonal, aes(x = time, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.4) +
  geom_line(color = "darkblue", size = 1) +
  geom_rect(data = subset(df_seasonal, season == "Winter"),
            aes(xmin = time, xmax = time + 28, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "skyblue", alpha = 0.2) +
  geom_rect(data = subset(df_seasonal, season == "Summer"),
            aes(xmin = time, xmax = time + 28, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  labs(
    title = "Posterior Seasonal Effect (Infection Hazard)",
    subtitle = "Shaded areas: Winter (blue), Summer (orange)",
    x = "Date",
    y = "Effect size (log scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################################
library(dplyr)
library(lubridate)
library(ggplot2)
# Extract posterior samples from your Stan model
posterior_samples <- rstan::extract(fit_test)

# Confirm seasonal effect exists
str(posterior_samples$seasonal_effect_1_2)
# Compute summaries over MCMC draws
seasonal_mean <- apply(posterior_samples$seasonal_effect_1_2, 2, mean)
seasonal_lwr <- apply(posterior_samples$seasonal_effect_1_2, 2, quantile, probs = 0.025)
seasonal_upr <- apply(posterior_samples$seasonal_effect_1_2, 2, quantile, probs = 0.975)

# Assume `df` is your full simulation data with `observed_state` and `date_use`
# And `seasonal_effect_df` contains posterior seasonal effect over time (aligned with `X_midpoints`)

# Step 1: Compute monthly observed prevalence
monthly_prevalence <- observation_df %>%
  mutate(year_month = floor_date(as.Date(date_use, origin = "1970-01-01"), "month")) %>%
  group_by(year_month) %>%
  summarise(
    prevalence = mean(observed_state),
    .groups = "drop"
  )

X_dates <- as.Date(global_interval_start_numeric + X * 365, origin = "1970-01-01")

posterior_df <- data.frame(
  date = X_dates,
  seasonal_mean = seasonal_mean,
  seasonal_lwr = seasonal_lwr,
  seasonal_upr = seasonal_upr
)
monthly_prevalence <- observation_df %>%
  mutate(year_month = floor_date(as.Date(date_use, origin = "1970-01-01"), "month")) %>%
  group_by(year_month) %>%
  summarise(
    prevalence = mean(observed_state),
    .groups = "drop"
  )


ggplot() +
  geom_ribbon(data = posterior_df, aes(x = date, ymin = seasonal_lwr, ymax = seasonal_upr),
              fill = "steelblue", alpha = 0.2) +
  geom_line(data = posterior_df, aes(x = date, y = seasonal_mean), color = "blue", size = 1.2) +
  geom_point(data = monthly_prevalence, aes(x = year_month, y = prevalence), color = "darkred", size = 2) +
  geom_line(data = monthly_prevalence, aes(x = year_month, y = prevalence), color = "darkred", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior Seasonal Effect vs. Observed Prevalence",
    x = "Date",
    y = "Log Hazard (Effect) / Prevalence"
  )
#################
hist(posterior_samples$beta_season_1_2, breaks = 30, main = "Posterior of β_season_1_2")
##########
