# ------------------------------------------------------------------------------
# Seasonality Setup – Simulation vs. Fitting 
# SIMULATED DATA: No seasonality
# MODEL FITTING: Seasonality with spline
# ------------------------------------------------------------------------------
# Date created: 10 July 2025
# Date last updated: 17 July 2025
# Author: Raneem Aizouk

rm(list=ls())

# Load libraries
library(rstan); library(V8); library(ggplot2); library(dplyr); library(wesanderson);library(patchwork); 
library(lubridate); library(tidyr)

# Simulation:
# Seasonality Disabled (i.e. no seasonality, this to check if a model with seasonality will then also fit no seasonal trend):
#  λ₁₂ (infection/acquisition) and λ₂₁ (recovery) rates do not vary with time.
#  Seasonality terms (sin/cos or splines) removed from simulation.
#  Rates depend only on age, sex, household random effects, and baseline values.
# ------------------------------------------------------------------------------
# Fitting:
#  Acquisition seasonality modeled using B-splines (flexible).
#  Recovery modeled as constant (no seasonal component).
# Goal:
#  Evaluate if spline-based fitting can recover sinusoidal acquisition pattern.

#---------------------------------------------
# Simulate Individual-Level Inputs
# ---------------------------------------------
set.seed(123)
V <- 50; houses_per_village <- 15; individuals_per_house <- 5
H <- V * houses_per_village
N_individuals <- H * individuals_per_house
n_rounds <- 4
N <- N_individuals * n_rounds  

# ---------------------------------------------
# Assign VillageID and Intervention Assignment
# ---------------------------------------------
Village_base <- rep(1:V, each = houses_per_village)
village_ids_per_individual <- rep(Village_base, each = individuals_per_house)
VillageID <- rep(village_ids_per_individual, each = n_rounds)

# Randomly assign 25 intervention villages
set.seed(123)
intervention_villages <- sample(1:V, 25)
village_intervention_map <- ifelse(1:V %in% intervention_villages, 1, 0)
Intervention <- village_intervention_map[VillageID]

##########################
HouseID_base <- rep(1:H, each = individuals_per_house)
menage_id_member_base <- HouseID_base * 10 + rep(1:individuals_per_house, times = H)
HouseID <- rep(HouseID_base, each = n_rounds)
menage_id_member <- rep(menage_id_member_base, each = n_rounds)
rounds <- rep(1:n_rounds, times = N_individuals)

#####################
# 1. Each group of 15 houses belongs to one village
house_to_village <- rep(1:V, each = houses_per_village)  # length = H

# 2. Assign each individual (via their HouseID) to the correct village
VillageID <- house_to_village[HouseID]  # HouseID is of length N (repeats per round)

# Individual characteristics
age_indiv <- sample(1:90, N_individuals, replace = TRUE)
sexe_indiv <- sample(0:1, N_individuals, replace = TRUE)
age <- rep(age_indiv, each = n_rounds)
sexe <- rep(sexe_indiv, each = n_rounds)

# Time intervals in days between rounds
d1 <- 90  # R1 -> R2
d2 <- 90 # R2 -> R3
d3 <- 180  # R3 -> R4

# Simulate start dates
baseline_date <- as.Date("2022-10-01")
start_dates <- baseline_date + sample(0:d1, N_individuals, replace = TRUE)

# Construct observation dates with constrained noise
set.seed(123)
date_use <- rep(as.Date(NA), N_individuals * n_rounds)
for (i in 1:N_individuals) {
  base <- start_dates[i]
  idx <- ((i - 1) * n_rounds + 1):(i * n_rounds)
  noise <- sample(-20:20, 3, replace = TRUE)  # Reduced noise to ensure monotonicity
  dates <- base + c(0, d1 + noise[1], d1 + d2 + noise[2], d1 + d2 + d3 + noise[3])
  date_use[idx] <- sort(dates)  # Ensure dates are in order
}

# Time Grid & Seasonality
# ─────────────────────────────────────────────────────────────────────────────
global_interval_start <- as.Date("2022-10-01")
global_interval_end <- as.Date("2024-02-19")
interval_length <- 28
interval_starts <- seq(as.numeric(global_interval_start),
                       as.numeric(global_interval_end),
                       by = interval_length)
interval_ends <- interval_starts + interval_length - 1
interval_starts_dates <- as.Date(interval_starts, origin = "1970-01-01")
interval_ends_dates <- as.Date(interval_ends, origin = "1970-01-01")

X_midpoints <- (interval_starts + pmin(interval_ends, as.numeric(global_interval_end))) / 2
if ((as.numeric(global_interval_end) - tail(interval_starts, 1)) < interval_length) {
  X_midpoints[length(X_midpoints)] <- tail(interval_starts, 1)
}

X <- (X_midpoints - X_midpoints[1])
X <- (X - median(X)) / 365
num_data <- length(X)

# True Parameters (for seasonal effect)
true_params <- list(
  q_1_2_raw = 0.01,
  q_2_1_raw = 0.01,
  sigma_q_1_2 = 0.3,
  sigma_q_2_1 = 0.3,
  beta_1_2_age = 0.01,
  beta_2_1_age = -0.03,
  beta_1_2_sexe = 0.03,
  beta_2_1_sexe = -0.02,
  u_raw = seq(-1, 1, length.out = H),
  sigma_u = 0.3,
  beta_int1_1 = -0.5,
  beta_int1_2 = -0.6,
  beta_int2_1 = 0.01,
  beta_int2_2 = 0.1
  
)
true_params$u <- true_params$u_raw * true_params$sigma_u


#####################

# Validate dates
if (any(is.na(date_use)) || any(date_use < global_interval_start) || any(date_use > global_interval_end)) {
  stop("Invalid observation dates detected.")
}
# Global start for intervention phase (e.g., after most R2 dates)
global_intervention_date1 <- as.Date("2023-03-01")
global_intervention_date2 <- as.Date("2023-06-01")

# Add small noise (e.g., ±5 days)
set.seed(123)
noise1 <- sample(-5:5, N_individuals, replace = TRUE)
noise2 <- sample(-5:5, N_individuals, replace = TRUE)

# Generate intervention dates per individual
simulated_intervention_date <- rep(global_intervention_date1, N_individuals) + noise1
simulated_intervention_date2 <- rep(global_intervention_date2, N_individuals) + noise2


simulated_intervention_date <- rep(simulated_intervention_date, each = n_rounds)
simulated_intervention_date2 <- rep(simulated_intervention_date2, each = n_rounds)

# # Simulate intervention dates
# simulated_intervention_date <- rep(
#   sample(data_complete$intervention_date, N_individuals, replace = TRUE),
#   each = n_rounds
# )
# simulated_intervention_date2 <- rep(
#   sample(data_complete$intervention_date2, N_individuals, replace = TRUE),
#   each = n_rounds
# )
simulated_intervention_date <- as.Date(simulated_intervention_date)
simulated_intervention_date2 <- as.Date(simulated_intervention_date2)
if (any(is.na(simulated_intervention_date)) || any(is.na(simulated_intervention_date2))) {
  stop("Invalid intervention dates detected.")
}
simulated_intervention_date <- as.numeric(simulated_intervention_date - global_interval_start)
simulated_intervention_date2 <- as.numeric(simulated_intervention_date2 - global_interval_start)

# Compute time grid indices
total_days <- as.numeric(global_interval_end - global_interval_start + 1)
num_intervals <- ceiling(total_days / interval_length)
max_middle <- num_intervals - 2
global_interval_start_numeric <- as.numeric(global_interval_start)
raw_idx_first <- floor((as.numeric(date_use) - global_interval_start_numeric) / interval_length) + 1
idx_first <- pmax(1, pmin(num_data, raw_idx_first))
raw_idx_last <- floor((as.numeric(date_use) - global_interval_start_numeric) / interval_length) + 1
idx_last <- pmax(1, pmin(num_data, raw_idx_last))

# Correct assignment of first/last/middle subintervals
# Initialize
first_subinterval <- rep(0, N)
last_subinterval <- rep(0, N)
num_middle_subintervals <- rep(0, N)


#######################

# Loop over individuals
for (i in 1:N_individuals) {
  idx <- ((i - 1) * n_rounds + 1):(i * n_rounds)
  
  d1 <- as.numeric(date_use[idx[1]])
  d2 <- as.numeric(date_use[idx[2]])
  d3 <- as.numeric(date_use[idx[3]])
  d4 <- as.numeric(date_use[idx[4]])
  
  # ---  First Subinterval (R1 → R2) → store at round 1
  days12 <- d2 - d1 + 1
  idx_first <- floor((d1 - global_interval_start_numeric) / interval_length) + 1
  interval_end_first <- global_interval_start_numeric + idx_first * interval_length - 1
  first_subinterval[idx[1]] <- pmin(interval_end_first - d1 + 1, days12,interval_length - 1)
  
  # --- Middle Subinterval (R2 → R3) ---
  total_days_23 <- d3 - d2 + 1
  num_middle_subintervals[idx[3]] <- pmax(0, ceiling(total_days_23 / interval_length) - 2)  # ✅ at Round 3
  
  # --- Last Subinterval (R3 → R4) ---
  total_days_34 <- d4 - d3 + 1
  idx_last <- floor((d4 - global_interval_start_numeric) / interval_length) + 1
  interval_start_last <- global_interval_start_numeric + (idx_last - 1) * interval_length
  last_subinterval[idx[4]] <- pmin(d4 - interval_start_last + 1,interval_length - 1)  # at Round 4
}


# Compute middle interval indices
idx_middle <- matrix(1L, nrow = N, ncol = max_middle)
raw_start_idx <- floor((as.numeric(date_use) - global_interval_start_numeric) / interval_length)
for (i in 1:N) {
  if (num_middle_subintervals[i] > 0) {
    for (m in 1:num_middle_subintervals[i]) {
      interval_start_m <- global_interval_start_numeric + (raw_start_idx[i] + m - 1) * interval_length
      interval_end_m <- interval_start_m + interval_length - 1
      interval_mid_m <- (interval_start_m + interval_end_m) / 2
      raw_midm <- floor((interval_mid_m - global_interval_start_numeric) / interval_length) + 1
      idx_middle[i, m] <- max(1, min(num_data, raw_midm))
    }
  }
}

# Final time index calculations
raw_idx_first <- floor((as.numeric(date_use) - global_interval_start_numeric) / interval_length) + 1
idx_first <- pmax(1, pmin(num_data, raw_idx_first))

raw_idx_last <- floor((as.numeric(date_use) - global_interval_start_numeric) / interval_length) + 1
idx_last <- pmax(1, pmin(num_data, raw_idx_last))

# Transition rates (final)
lambda12_vec <- true_params$q_1_2_raw +
  true_params$u[HouseID] +
  true_params$beta_1_2_age * age +
  true_params$beta_1_2_sexe * sexe# +
#Y_hat_1_2[idx_last]

lambda21_vec <- true_params$q_2_1_raw +
  true_params$u[HouseID] +
  true_params$beta_2_1_age * age +
  true_params$beta_2_1_sexe * sexe


########################
# Number of individuals
#N_individuals <- length(date_use) / 4

# Indices for Round 2 and Round 3
idx_round2 <- seq(2, by = 4, length.out = N_individuals)
idx_round3 <- seq(3, by = 4, length.out = N_individuals)

# Get dates
dates_round2 <- as.Date(date_use[idx_round2])
dates_round3 <- as.Date(date_use[idx_round3])

# Compute duration between R2 and R3
interval_r2_r3 <- as.numeric(dates_round3 - dates_round2)

print(interval_r2_r3)
print(min(interval_r2_r3))
print(max(interval_r2_r3))

# Assemble Stan Data
stan_data <- list(
  N = N,
  H = H,
  N_persons = length(unique(menage_id_member)),
  
  menage_id_member = menage_id_member,
  
  round = rounds,
  HouseID = HouseID,
  age = age,
  sexe = sexe,
  date_use = as.numeric(date_use),
  intervention_date = simulated_intervention_date,
  intervention_date2 = simulated_intervention_date2,
  global_interval_start = global_interval_start_numeric,
  global_interval_end = as.numeric(global_interval_end),
  interval_length = interval_length,
  
  num_data = num_data,
  X = X,
  
  q_1_2_base = true_params$q_1_2_raw,
  q_2_1_base = true_params$q_2_1_raw,
  beta_1_2_age = true_params$beta_1_2_age,
  beta_2_1_age = true_params$beta_2_1_age,
  beta_1_2_sexe = true_params$beta_1_2_sexe,
  beta_2_1_sexe = true_params$beta_2_1_sexe,
  beta_int1_1 = true_params$beta_int1_1,
  beta_int1_2 = true_params$beta_int1_2,
  beta_int2_1 = true_params$beta_int2_1,
  beta_int2_2 = true_params$beta_int2_2,
  u = true_params$u,
  
  max_middle = max_middle,
  idx_first = idx_first,
  idx_last = idx_last,
  first_subinterval = pmax(0, first_subinterval),
  last_subinterval = pmax(0, last_subinterval),
  num_middle_subintervals = num_middle_subintervals,
  idx_middle = idx_middle,
  
  lambda12_vec = as.numeric(lambda12_vec),
  lambda21_vec = as.numeric(lambda21_vec)
)


# Stan Model
stan_code <- "
functions {
  matrix transition_matrix(real delta, real rate12, real rate21) {
    matrix[2, 2] P;
    real total_rate = rate12 + rate21;
    real p12 = total_rate > 0 ? rate12 / total_rate : 0.5;
    real p21 = total_rate > 0 ? rate21 / total_rate : 0.5;
    real exp_neg = exp(-delta * total_rate);

    P[1,1] = p21 + p12 * exp_neg;
    P[1,2] = p12 - p12 * exp_neg;
    P[2,1] = p21 - p21 * exp_neg;
    P[2,2] = p12 + p21 * exp_neg;

    return P;
  }
}

data {
  int<lower=1> N;
  int<lower=1> H;
  int<lower=1> N_persons;
  array[N] int<lower=1> menage_id_member;
  array[N] int<lower=1> HouseID;
  array[N] int<lower=0> age;
  array[N] int<lower=0, upper=1> sexe;
  array[N] real date_use;
  array[N] real intervention_date;
  array[N] real intervention_date2;
  real global_interval_start;
  int interval_length;
  int<lower=1> num_data;
  array[num_data] real X;

  real global_interval_end;
  int<lower=0> max_middle;
  real q_1_2_base;
  real q_2_1_base;
  real beta_1_2_age;
  real beta_2_1_age;
  real beta_1_2_sexe;
  real beta_2_1_sexe;
  real beta_int1_1;
  real beta_int1_2;
  real beta_int2_1;
  real beta_int2_2;
  array[H] real u;
  array[N] real lambda12_vec;
  array[N] real lambda21_vec;
  array[N] int<lower=1> idx_first;
  array[N] int<lower=1> idx_last;
  array[N] int<lower=0> first_subinterval;
  array[N] int<lower=0> last_subinterval;
  array[N] int<lower=0, upper=max_middle> num_middle_subintervals;
  array[N, max_middle] int<lower=1> idx_middle;
}

generated quantities {
  vector[N]   log_lambda_1_2_out;
  vector[N]   log_lambda_2_1_out;

       int    y_rep[N];

   for (n in 1:N) {
    real log12 = q_1_2_base + u[HouseID[n]] + beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n];
    real log21 = q_2_1_base + u[HouseID[n]] + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n];

    if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
      y_rep[n] = categorical_rng([0.5, 0.5]');
      log_lambda_1_2_out[n] = log12 ;
      log_lambda_2_1_out[n] = log21;
      continue;
    }

    matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];

    // --- First Subinterval ---
    if (first_subinterval[n] > 0) {
      int i1 = idx_first[n];
     // real s12 = Y_hat_1_2_out[i1];
      real t0 = date_use[n - 1];
      real t1 = t0 + first_subinterval[n] - 1;

      if (t0 < t_star1 && t_star1 < t1) {
        real d1a = t_star1 - t0;
        real d1b = t1 - t_star1 + 1;
        P_total *= transition_matrix(d1a, exp(log12 ), exp(log21));
        int flag1 = 1;
        int flag2 = ((t_star1 + t1) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(d1b, exp(log12 + flag1 * beta_int1_1 + flag2 * beta_int1_2 ),
                                          exp(log21 + flag1 * beta_int2_1 + flag2 * beta_int2_2));
      } else {
        int flag1 = ((t0 + t1) / 2 >= t_star1) ? 1 : 0;
        int flag2 = ((t0 + t1) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(first_subinterval[n], exp(log12 + flag1 * beta_int1_1 + flag2 * beta_int1_2 ),
                                                   exp(log21 + flag1 * beta_int2_1 + flag2 * beta_int2_2));
      }
    }

    // --- Middle Subintervals ---
    int raw_start_idx = to_int(floor((date_use[n - 1] - global_interval_start) / interval_length));
    for (m in 1:num_middle_subintervals[n]) {
      int im = idx_middle[n, m];
      //real s12m = Y_hat_1_2_out[im];
      real t0m = global_interval_start + (raw_start_idx + m - 1) * interval_length;
      real t1m = t0m + interval_length - 1;

      if (t0m < t_star1 && t_star1 < t1m) {
        real d2a = t_star1 - t0m;
        real d2b = t1m - t_star1 + 1;
        P_total *= transition_matrix(d2a, exp(log12 ), exp(log21));
        int flag1 = 1;
        int flag2 = ((t_star1 + t1m) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(d2b, exp(log12 + flag1 * beta_int1_1 + flag2 * beta_int1_2 ),
                                          exp(log21 + flag1 * beta_int2_1 + flag2 * beta_int2_2));
      } else {
        int flag1 = ((t0m + t1m) / 2 >= t_star1) ? 1 : 0;
        int flag2 = ((t0m + t1m) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(interval_length, exp(log12 + flag1 * beta_int1_1 + flag2 * beta_int1_2),
                                                    exp(log21 + flag1 * beta_int2_1 + flag2 * beta_int2_2));
      }
    }

    // --- Last Subinterval ---
    if (last_subinterval[n] > 0) {
      int il = idx_last[n];
      //real s12l = Y_hat_1_2_out[il];
      real t1l = date_use[n];
      real t0l = t1l - last_subinterval[n] + 1;

      if (t0l < t_star2 && t_star2 < t1l) {
        real d3a = t_star2 - t0l;
        real d3b = t1l - t_star2 + 1;
        P_total *= transition_matrix(d3a, exp(log12 + beta_int1_1 ), exp(log21 + beta_int2_1));
        P_total *= transition_matrix(d3b, exp(log12 + beta_int1_1 + beta_int1_2 ),
                                          exp(log21 + beta_int2_1 + beta_int2_2));
      } else {
        int flag2 = ((t0l + t1l) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(last_subinterval[n], exp(log12 + beta_int1_1 + flag2 * beta_int1_2 ),
                                                     exp(log21 + beta_int2_1 + flag2 * beta_int2_2));
      }
    }

    // Simulate next state
    y_rep[n] = categorical_rng(P_total[y_rep[n - 1]]' / sum(P_total[y_rep[n - 1]]));

    // Output transition rates
    log_lambda_1_2_out[n] = log12 ;
    log_lambda_2_1_out[n] = log21;
    }
}
"

# Compile and Run Stan Model
stan_model <- stan_model(model_code = stan_code)

set.seed(100)

stan_fit <- sampling(
  object = stan_model,
  data = stan_data,
  algorithm = "Fixed_param",
  iter = 1,
  chains = 1,
  warmup = 0,
  refresh = 0
)

#print(stan_fit)


# 1. Extract simulated data from stan_fit
sim_data <- rstan::extract(stan_fit)

# Extract generated quantities
observed_state_sim <- sim_data$y_rep[1, ]  # Simulated states (y_rep is the correct variable)
log_lambda_1_2_out <- sim_data$log_lambda_1_2_out[1, ]  # Log acquisition rates
log_lambda_2_1_out <- sim_data$log_lambda_2_1_out[1, ]  # Log recovery rates
#Y_hat_1_2_out <- sim_data$Y_hat_1_2_out[1, ]  # Seasonal effect
observed_state_sim <- sim_data$y_rep[1, ]

# 2. Create data frame combining simulated and input data
simulated_dataset <- data.frame(
  Observation = 1:N,
  MenageID = menage_id_member,
  HouseID = HouseID,
  VillageID = VillageID,
  Intervention = Intervention,
  Round = rounds,
  Date = as.Date(date_use, origin = "1970-01-01"),
  Observed_State_Sim = observed_state_sim,  # Changed from Simulated_State
  Age = age,
  Sexe = sexe,
  Intervention_Date1 = as.Date(simulated_intervention_date + as.numeric(global_interval_start), origin = "1970-01-01"),
  Intervention_Date2 = as.Date(simulated_intervention_date2 + as.numeric(global_interval_start), origin = "1970-01-01"),
  Idx_First = stan_data$idx_first,
  Idx_Last = stan_data$idx_last,
  First_Subinterval = stan_data$first_subinterval,
  Last_Subinterval = stan_data$last_subinterval,
  Num_Middle_Subintervals = stan_data$num_middle_subintervals,
  Log_Lambda_1_2 = log_lambda_1_2_out,
  Log_Lambda_2_1 = log_lambda_2_1_out
)


head(simulated_dataset)
colnames(simulated_dataset)
print(min(simulated_dataset$Intervention_Date1))
print(max(simulated_dataset$Intervention_Date1))
print(min(simulated_dataset$Intervention_Date2))
print(max(simulated_dataset$Intervention_Date2))

# PLOT SIMULATED DATASET
# Colours
palette <- wes_palette("Darjeeling1", n = 5)
palette2 <- wes_palette("BottleRocket2", n = 1)
palette3 <- wes_palette("GrandBudapest1", n = 2)[2]
palette4 <- wes_palette("BottleRocket2", n = 2)[2]
palette5 = c(palette3, palette[2],palette2,palette[5], palette[4],palette4)

# Define colors for each pair of shaded areas
shaded_colors <- c("red", "blue", "black")

p = ggplot(simulated_dataset%>%filter(Intervention==1), aes(x = Date, fill =as.character(Round))) +
  geom_bar() +
  facet_grid(VillageID~., scales = "free_x") +  # Use facet_grid with free x-scales
  labs(x = "Date", y = "Count", title = "Village sampling times",
       subtitle="Intervention group", fill="Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "4 week", date_labels = "%Y-%m-%d") +  # Set weekly breaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = "white")) +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = min(simulated_dataset$Intervention_Date1), 
                              xmax = max(simulated_dataset$Intervention_Date1), 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors[1]),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2) +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = min(simulated_dataset$Intervention_Date2), 
                              xmax = max(simulated_dataset$Intervention_Date2), 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors[2]),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2)+
  scale_fill_manual(values = c("1" = palette5[1], 
                               "2" = palette5[2], 
                               "3" =palette5[3], 
                               "4" = palette5[4])) 
p   

p2 = ggplot(simulated_dataset%>%filter(Intervention==0), aes(x = Date, fill =as.character(Round))) +
  geom_bar() +
  facet_grid(VillageID~., scales = "free_x") +  # Use facet_grid with free x-scales
  labs(x = "Date", y = "Count", title = "Village sampling times",
       subtitle="Control group", fill="Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "4 week", date_labels = "%Y-%m-%d") +  # Set weekly breaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = "white")) +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = min(simulated_dataset$Intervention_Date1), 
                              xmax = max(simulated_dataset$Intervention_Date1), 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors[1]),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2) +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = min(simulated_dataset$Intervention_Date2), 
                              xmax = max(simulated_dataset$Intervention_Date2), 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors[2]),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2)+
  scale_fill_manual(values = c("1" = palette5[1], 
                               "2" = palette5[2], 
                               "3" =palette5[3], 
                               "4" = palette5[4])) 
p2   

combined = p+p2 +plot_layout(ncol=2)  

#ggsave(file = "./Output/Figures/Simulated_data/samplingtimes_villages_sim.png", 
#       plot = combined, width = 20, height = 13) 

# Summarize observed prevalence over time (control group)
prev_df <- simulated_dataset %>%
  mutate(month_year = floor_date(Date, "month")) %>%
  group_by(month_year, Intervention) %>%
  summarise(
    n_colonized = sum(Observed_State_Sim == 2),
    n_total = n(),
    prevalence = n_colonized / n_total,
    .groups = "drop"
  )

# Plot without rescaling
p3 = ggplot(prev_df, aes(x = month_year, y = prevalence, color=factor(Intervention))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Observed Prevalence",
    subtitle = "Control vs Intervention group",
    x = "Time", y = "Prevalence") +
  ylim(0,1.5) +
  #facet_wrap(~ Intervention, ncol=2) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = "white"))+
  scale_color_manual(values = c("black", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed"))
p3

ggsave(file = "./Output/Figures/Simulated_data/Cases_no_seasonalpattern_sim_seasonality.png", 
       plot = p3, width = 8, height = 7) 

###################
num_knots <- 5
knots <- quantile(X, probs = seq(0, 1, length.out = num_knots))
knots <- as.numeric(knots)
spline_degree <- 3
num_basis <- num_knots + spline_degree - 1

# ─── Ensure all required columns exist in simulated_dataset ─────────────────────
simulated_dataset$age <- if (!"age" %in% names(simulated_dataset)) age else simulated_dataset$age
simulated_dataset$sexe <- if (!"sexe" %in% names(simulated_dataset)) sexe else simulated_dataset$sexe
simulated_dataset$round <- if (!"Round" %in% names(simulated_dataset)) rounds else simulated_dataset$Round

# Dates need to be numeric for Stan
simulated_dataset$date_use <- as.numeric(simulated_dataset$Date)
simulated_dataset$intervention_date <- as.numeric(simulated_dataset$Intervention_Date1)
simulated_dataset$intervention_date2 <- as.numeric(simulated_dataset$Intervention_Date2)
stan_data$intervention_village <- as.integer(simulated_dataset$Intervention)

# Use these simulated variables from your simulation
stan_data <- list(
  # Basic counts
  N = nrow(simulated_dataset),
  H = length(unique(simulated_dataset$HouseID)),
  
  # Observed state sequence (must be 1 or 2)
  observed_state = as.integer(simulated_dataset$Observed_State_Sim),
  
  # IDs and individual-level covariates
  menage_id_member = simulated_dataset$MenageID,
  HouseID = simulated_dataset$HouseID,
  age = simulated_dataset$age,
  sexe = simulated_dataset$sexe,
  round = simulated_dataset$round,
  
  # Time settings
  global_interval_start = as.numeric(global_interval_start),
  global_interval_end = as.numeric(global_interval_end),
  interval_length = interval_length,
  date_use = simulated_dataset$date_use,
  
  # Interventions
  intervention_date = simulated_dataset$intervention_date,
  intervention_date2 = simulated_dataset$intervention_date2,
  
  # Spline / seasonal components
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  num_basis = num_basis,
  num_data = length(X),
  X = X,
  
  # Binning / indexing
  max_middle = max_middle,
  num_intervals = num_intervals,
  
  # Optional: Pre-indexed subinterval mapping
  idx_first = simulated_dataset$Idx_First,
  idx_last = simulated_dataset$Idx_Last,
  first_subinterval = simulated_dataset$First_Subinterval,
  last_subinterval = simulated_dataset$Last_Subinterval,
  num_middle_subintervals = simulated_dataset$Num_Middle_Subintervals
)
stan_data$intervention_village <- simulated_dataset$Intervention

str(stan_data)

# Prepare data for Stan
# Prepare data for Stan
stan_model_full <- "
functions {
  // 2×2 CTMC transition matrix for a 2-state chain
  matrix transition_matrix(real t,
                           real lambda_1_2,
                           real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;
    real exp_term     = exp(-total_lambda * t);

    matrix[2,2] P;
    P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;
    P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;
    P[1,2] = (lambda_2_1 / total_lambda) * (1 - exp_term);
    P[2,1] = (lambda_1_2 / total_lambda) * (1 - exp_term);
    return P;
  }

  // Recursive function to compute B-spline basis functions
  vector build_b_spline(real[] t,
                        real[] ext_knots,
                        int   ind,
                        int   order) {
    int M = size(t);
    vector[M] b_spline = rep_vector(0, M);
    if (order == 1) {
      for (i in 1:M)
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    } else {
      real denom1 = ext_knots[ind+order-1] - ext_knots[ind];
      real denom2 = ext_knots[ind+order]   - ext_knots[ind+1];
      vector[M] w1 = denom1 > 0
        ? (to_vector(t) - rep_vector(ext_knots[ind], M)) / denom1
        : rep_vector(0, M);
      vector[M] w2 = denom2 > 0
        ? 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], M)) / denom2
        : rep_vector(0, M);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind,     order-1)
               + w2 .* build_b_spline(t, ext_knots, ind + 1, order-1);
    }
    return b_spline;
  }


 }

data {
  int<lower=1> N;
  int<lower=1> H;
  array[N] int<lower=1, upper=2> observed_state;
  array[N] int menage_id_member;
  array[N] int HouseID;
  array[N] int age;
  array[N] int round;
  array[N] int sexe;
  real global_interval_start;
  real global_interval_end;
  int<lower=1> interval_length;
  array[N] real date_use;
  array[N] int intervention_village;
  int<lower=1> num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int<lower=1> num_basis;
  int<lower=1> num_data;
  array[num_data] real X;

  //array[N] real Intervention_start_date;
  array[N] real intervention_date;  // or int if stored as day numbers
  array[N] real intervention_date2; 
  int<lower=1> max_middle;
  int<lower=1> num_intervals;


}

transformed data {
  matrix[num_basis, num_data] B;

  // Extend knots
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2 * spline_degree + num_knots] ext_knots;
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots      = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));

  // Build B-spline basis
  for (ind in 1:num_basis) {
    B[ind, :] = to_row_vector(
      build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1)
    );
  }
  // Boundary constraint
  B[num_knots + spline_degree - 1, num_data] = 1;

  array[N] real first_subinterval;
  array[N] real middle_subinterval;
  array[N] real last_subinterval;
  array[N] int  num_middle_subintervals;
  array[N] int  global_interval_index_start;
  array[N] int  global_interval_index_end;
  array[N] real interval_start_first;
  array[N] real interval_end_first;

  array[N, max_middle] real interval_start_middle;
  array[N, max_middle] real interval_end_middle;

  array[N] real interval_start_last;
  array[N] real interval_end_last;
  // Declare and initialize the index arrays to 1
  array[N] int idx_first;
  array[N, max_middle] int idx_middle;
  array[N] int idx_last;
  // Fill every position with “1” so that unused slots remain in range
  idx_first  = rep_array(1, N);
  idx_middle = rep_array(1, N, max_middle);
  idx_last   = rep_array(1, N);

  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n-1]) {
      real date1 = date_use[n-1];
      real date2 = date_use[n];

      // 1) map date1/date2 into global-interval bins
      int raw_start = to_int(floor((date1 - global_interval_start) / interval_length));
      int raw_end   = to_int(ceil((date2 - global_interval_start) / interval_length));

      // 1.1) clamp into [1, num_intervals]
      global_interval_index_start[n] = max(1, min(num_intervals, raw_start + 1));
      global_interval_index_end  [n] = max(1, min(num_intervals, raw_end));

      // 2) number of full intervals strictly between
      num_middle_subintervals[n] = global_interval_index_end[n]
                                - global_interval_index_start[n]
                                - 1;

      // 3) first subinterval calendar start/end
      interval_start_first[n] = global_interval_start
                              + (global_interval_index_start[n] - 1) * interval_length;
      interval_end_first  [n] = interval_start_first[n] + interval_length - 1;

      // 4) middle subintervals calendar start/end
      for (m in 1:num_middle_subintervals[n]) {
        int global_idx = global_interval_index_start[n] + m;
        interval_start_middle[n,m] = global_interval_start
                                   + (global_idx - 1) * interval_length;
        interval_end_middle  [n,m] = interval_start_middle[n,m] + interval_length - 1;
      }

      // 5) last subinterval calendar start/end
      interval_start_last[n] = global_interval_start + (raw_end - 1) * interval_length;
      interval_end_last  [n] = date2;

      // 6) durations
      first_subinterval[n]  = interval_end_first[n] - date1 + 1;
      last_subinterval [n]  = date2
                            - (global_interval_start
                               + (global_interval_index_end[n] - 1) * interval_length)
                            + 1;
      middle_subinterval[n] = num_middle_subintervals[n] * interval_length;
    // 7) For each subinterval, compute the midpoint and map it to the corresponding global‐interval index (1…num_data):
    // 7.1) first midpoint → idx_first
      {
        real mid1 = (date1 + interval_end_first[n]) / 2;
        int raw_mid1 = to_int(floor((mid1 - global_interval_start) / interval_length)) + 1;
        idx_first[n] = max(1, min(num_data, raw_mid1));
      }

      // 7.2) middle midpoints → idx_middle
      for (m in 1:num_middle_subintervals[n]) {
        real ms   = interval_start_middle[n,m];
        real me   = interval_end_middle  [n,m];
        real midm = (ms + me) / 2;
        int raw_midm = to_int(floor((midm - global_interval_start) / interval_length)) + 1;
        idx_middle[n,m] = max(1, min(num_data, raw_midm));
      }

      // 7.3) last midpoint → idx_last
      {
        real ml = (interval_start_last[n] + date2) / 2;
        int raw_ml = to_int(floor((ml - global_interval_start) / interval_length)) + 1;
        idx_last[n] = max(1, min(num_data, raw_ml));
      }
    }
  }
}


 
parameters {
  // Transition rate base parameters (non-centered)
  real q_1_2_raw;
  real q_2_1_raw;

  real<lower=0> sigma_q_1_2;
  real<lower=0> sigma_q_2_1;

  // Random effects
  vector[H] u_raw;
  real<lower=0> sigma_u;

  // Covariate effects
  real beta_1_2_age;
  real beta_2_1_age;
  real beta_1_2_sexe;
  real beta_2_1_sexe;

  // Intervention and seasonal covariate effects
  real beta_int1;             // first intervention
  real beta_int2;             // second intervention
  real beta_season_1_2;
  real beta_season_2_1;

  // Infection spline (periodic)
  row_vector[num_basis - 3] a_raw_1_2_free;
  real log_tau_raw_1_2;

  // Recovery spline (periodic)
  row_vector[num_basis - 3] a_raw_2_1_free;
  real log_tau_raw_2_1;

  // Linear trend components added to seasonal splines
  real a0_raw_1_2;
  real<lower=0> sigma_a0_1_2;

  real a0_raw_2_1;
  real<lower=0> sigma_a0_2_1;

  // Noise parameter for observation model
  real<lower=0> sigma;
}



transformed parameters {
  // ─── Base transition log-rates ───────────────────────────────────────────────
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;

  vector[H] u = u_raw * sigma_u;

  // ─── Infection spline coefficients ───────────────────────────────────────────
  row_vector[num_basis] a_raw_1_2;

  // Complete the periodic spline by repeating the first 3 basis weights at the end
  for (i in 1:(num_basis - 3))
    a_raw_1_2[i] = a_raw_1_2_free[i];

  for (j in 1:3)
    a_raw_1_2[num_basis - 3 + j] = a_raw_1_2_free[j];

  // Scale by global variance
  real tau_1_2 = exp(log_tau_raw_1_2);
  row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;

  // Linear trend term
  real a0_1_2 = a0_raw_1_2 * sigma_a0_1_2;

  // Final seasonal effect for infection (with linear + periodic terms)
  vector[num_data] Y_hat_1_2 = a0_1_2 * to_vector(X) + to_vector(a_1_2 * B);

  // ─── Recovery spline coefficients ────────────────────────────────────────────
  row_vector[num_basis] a_raw_2_1;

  for (i in 1:(num_basis - 3))
    a_raw_2_1[i] = a_raw_2_1_free[i];

  for (j in 1:3)
    a_raw_2_1[num_basis - 3 + j] = a_raw_2_1_free[j];

  real tau_2_1 = exp(log_tau_raw_2_1);
  row_vector[num_basis] a_2_1 = a_raw_2_1 * tau_2_1;

  real a0_2_1 = a0_raw_2_1 * sigma_a0_2_1;

  // Final seasonal effect for recovery (with linear + periodic terms)
  vector[num_data] Y_hat_2_1 = a0_2_1 * to_vector(X) + to_vector(a_2_1 * B);
}





model {
  // ─────────────────────────────────────────────────────────────────────────────
  // PRIORS
  // ─────────────────────────────────────────────────────────────────────────────
  // Base CTMC transition coefficients (raw scale → centered in transformed parameters)
  q_1_2_raw ~ normal(0, 0.3);
  q_2_1_raw ~ normal(0, 0.38);
  
  // Age and sex effects on infection/recovery hazards
  beta_1_2_age  ~ normal(0, 1);
  beta_2_1_age  ~ normal(0, 1);
  beta_1_2_sexe ~ normal(0, 1);
  beta_2_1_sexe ~ normal(0, 1);

  // Household‐level random effects
  u_raw   ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);

  // Spline parameters for seasonal variation (infection and recovery)
  a_raw_1_2     ~ normal(0, 1);
  a_raw_2_1     ~ normal(0, 1);
 
  a0_raw_1_2    ~ normal(0, 0.5);
  a0_raw_2_1    ~ normal(0, 0.5);
  sigma_a0_1_2  ~ normal(0, 0.1);
  sigma_a0_2_1  ~ normal(0, 0.05);
 
  log_tau_raw_1_2 ~ normal(0, 1);
  log_tau_raw_2_1 ~ normal(0, 2);

  // Noise on continuous‐time transitions
  sigma ~ normal(0, 1);

  // Intervention‐specific effects (first and second)
  beta_int1 ~ normal(0, 1);
  beta_int2 ~ normal(0, 1);
  beta_season_1_2~ normal(0, 1);
  beta_season_2_1~ normal(0, 1);
  // ─────────────────────────────────────────────────────────────────────────────
  // LIKELIHOOD: FOR EACH PAIR OF SUCCESSIVE OBSERVATIONS (n−1 → n) WITHIN A PERSON
  // ─────────────────────────────────────────────────────────────────────────────
  for (n in 2:N) {
    // Only compute transitions if this row (n) belongs to the same individual
    // as the previous row (n−1).  Otherwise, skip (no baseline → no transition).
    if (menage_id_member[n] == menage_id_member[n-1]) {
      // 1) Compute baseline log‐hazards for infection (1→2) and recovery (2→1):
      //    q_1_2_base and q_2_1_base are the “population‐average” log‐rates,
      //    u[HouseID[n]] is the household‐specific random intercept,
      //    beta_*_age and beta_*_sexe are individual‐level covariates.
      real lambda12 = q_1_2_base
                     + u[HouseID[n]]
                     + beta_1_2_age  * age[n]
                     + beta_1_2_sexe * sexe[n];
      real lambda21 = q_2_1_base
                     + u[HouseID[n]]
                     + beta_2_1_age  * age[n]
                     + beta_2_1_sexe * sexe[n];

      // Start with an “identity” 2×2 matrix P_total; we will multiply in pieces.
      matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));

      // Intervention dates (first and second) in numeric form (days since epoch)
      real t_star1 = intervention_date[n];
      real t_star2 = intervention_date2[n];

      // ─────────────────────────────────────────────────────────────────────────
      // (A) FIRST SUBINTERVAL: from (date_use[n−1]) up to FIRST interval boundary
      // ─────────────────────────────────────────────────────────────────────────
      if (first_subinterval[n] > 0) {
        // idx_first[n] gives the spline index (1…num_data) for this piece’s midpoint.
        int    i1  = idx_first[n];
        real   s12 = Y_hat_1_2[i1];  // seasonal adjustment at that midpoint for infection
       // real   s21 = Y_hat_2_1[i1];  // seasonal adjustment at that midpoint for recovery

        // t0 = the actual calendar date of observation (n−1)
        real t0    = date_use[n-1];
        // t1 = end of “first subinterval,” measured in days.
        // (Because first_subinterval[n] = length of that piece in days.)
        real t1    = t0 + first_subinterval[n] - 1;

        // (A.1) If t_star1 lies strictly inside [t0, t1], we split this subinterval
        if (t0 < t_star1 && t_star1 < t1) {
          real d1a = t_star1 - t0;          // # days before first‐intervention
          real d1b = t1      - t_star1 + 1; // # days after first‐intervention (but before second)

          // — Pre‐first piece (no intervention active):
          P_total *= transition_matrix(
            d1a,
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 +  s12),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2 )
          );

          // — Post‐first piece (first intervention active, but second still not yet active):
          // Compute midpoint of [t_star1, t1]
          real midpoint_post1 = (t_star1 + t1) / 2;
          // Because we are after t_star1, flag1_post1 = 1
          int  flag1_post1    = 1;
          // flag2_post1 = 1 only if that midpoint ≥ t_star2 (rare, since t_star2 >> t_star1)
          int  flag2_post1    = (midpoint_post1 >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            d1b,
            exp(lambda12 + flag1_post1*beta_int1 + flag2_post1*beta_int2 +  s12),
            exp(lambda21 + flag1_post1*beta_int1 + flag2_post1*beta_int2 )
          );

        } else {
          // (A.2) No split by t_star1 inside this “first” piece
          real midpoint1 = (t0 + t1) / 2;
          int  flag1_1   = (midpoint1 >= t_star1) ? 1 : 0;
          int  flag2_1   = (midpoint1 >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            first_subinterval[n],
            exp(lambda12 + flag1_1*beta_int1 + flag2_1*beta_int2 +  s12),
            exp(lambda21 + flag1_1*beta_int1 + flag2_1*beta_int2 )
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (B) MIDDLE SUBINTERVALS: intervals fully between 1st and last global bins
      // ─────────────────────────────────────────────────────────────────────────
      for (m in 1:num_middle_subintervals[n]) {
        int    im = idx_middle[n,m];
        real   s12_m = Y_hat_1_2[im];
      //  real   s21_m = Y_hat_2_1[im];

        // t0m = start of the m-th middle bin, t1m = end of that bin
        real t0m = global_interval_start
                   + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;

        // (B.1) If the first intervention date t_star1 falls inside this middle bin
        if (t0m < t_star1 && t_star1 < t1m) {
          real d2a = t_star1 - t0m;          // time before first‐intervention in this bin
          real d2b = t1m     - t_star1 + 1;  // time after first‐intervention (but before second)

          // — Pre‐first in this middle bin (no interventions yet):
          P_total *= transition_matrix(
            d2a,
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 +  s12_m),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2 )
          );

          // — Post‐first in this same bin (first active, second not yet active):
          real midpoint_post1m = (t_star1 + t1m) / 2;
          int  flag1_post1m    = 1;
          int  flag2_post1m    = (midpoint_post1m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            d2b,
            exp(lambda12 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 +  s12_m),
            exp(lambda21 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 )
          );

        } else {
          // (B.2) No split by t_star1 inside this middle bin
          real midpoint2m = (t0m + t1m) / 2;
          int  flag1_2m   = (midpoint2m >= t_star1) ? 1 : 0;
          int  flag2_2m   = (midpoint2m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            interval_length,
            exp(lambda12 + flag1_2m*beta_int1 + flag2_2m*beta_int2 +  s12_m),
            exp(lambda21 + flag1_2m*beta_int1 + flag2_2m*beta_int2 )
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (C) LAST SUBINTERVAL: from the last global‐bin boundary up to date_use[n]
      // ─────────────────────────────────────────────────────────────────────────
      if (last_subinterval[n] > 0) {
        int    il = idx_last[n];
        real   s12_l = Y_hat_1_2[il];
       // real   s21_l = Y_hat_2_1[il];

        real t1l = date_use[n];
        real t0l = t1l - last_subinterval[n] + 1;

        // By design, t0l > t_star1, since the “last subinterval” always starts
        // after the end of the first‐intervention window.  Therefore, we skip
        // any check for t_star1 here.  The only possible split in this last piece
        // is around t_star2 (September 2023), which may or may not lie inside [t0l, t1l].

        // (C.1) If second intervention falls inside this last subinterval:
        if (t0l < t_star2 && t_star2 < t1l) {
          real d3a = t_star2 - t0l;          // time from t0l until just before t_star2
          real d3b = t1l     - t_star2 + 1;  // time from t_star2 through t1l

          // — Pre‐second piece:  first intervention is already active on entire last piece —
          //   (Because t0l > t_star1 implies we are already post‐first for all days in [t0l, t_star2).)
          real flag1_pre2 = 1;
          real flag2_pre2 = 0;
          P_total *= transition_matrix(
            d3a,
            exp(lambda12 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 +  s12_l),
            exp(lambda21 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 )
          );

          // — Post‐second piece: both interventions active —
          real flag1_post2 = 1;
          real flag2_post2 = 1;
          P_total *= transition_matrix(
            d3b,
            exp(lambda12 + flag1_post2*beta_int1 + flag2_post2*beta_int2 +  s12_l),
            exp(lambda21 + flag1_post2*beta_int1 + flag2_post2*beta_int2 )
          );

        } else {
          // (C.2) No split by t_star2 inside this last piece:
          //    first intervention is already active (flag1 = 1).
          //    second intervention may or may not yet be active, depending on midpoint.
          real midpoint_l = (t0l + t1l) / 2;
          int  flag1_l    = 1;  // always post‐first
          int  flag2_l    = (midpoint_l >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            last_subinterval[n],
            exp(lambda12 + flag1_l*beta_int1 + flag2_l*beta_int2 +  s12_l),
            exp(lambda21 + flag1_l*beta_int1 + flag2_l*beta_int2 )
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (D) ADD LOG‐LIKELIHOOD CONTRIBUTION FOR THIS TRANSITION
      // ─────────────────────────────────────────────────────────────────────────
      // P_total is now the 2×2 transition probability matrix for going from
      // state at (n−1) to state at n, after multiplying through every “piece”
      // (first/middle/last) with the appropriate interventions flags.
      target += log(
        P_total[ observed_state[n-1], observed_state[n] ] + 1e-9
     );
    }
  }
}

generated quantities {
  vector[N]        log_lambda_1_2_out;
  vector[N]        log_lambda_2_1_out;
  vector[num_data] Y_hat_1_2_out;
  vector[num_data] Y_hat_2_1_out;
  int              y_rep[N];

  // Expose spline weights
  row_vector[num_basis] a_1_2_out = a_1_2;
  row_vector[num_basis] a_2_1_out = a_2_1;

  // Recompute full seasonal effect (spline + linear trend)
  Y_hat_1_2_out = a0_1_2 * to_vector(X) + to_vector(a_1_2 * B);
  Y_hat_2_1_out = a0_2_1 * to_vector(X) + to_vector(a_2_1 * B);

  for (n in 1:N) {
    // (1) Compute baseline log-hazards
    real log12 = q_1_2_base
                 + u[HouseID[n]]
                 + beta_1_2_age * age[n]
                 + beta_1_2_sexe * sexe[n];
    real log21 = q_2_1_base
                 + u[HouseID[n]]
                 + beta_2_1_age * age[n]
                 + beta_2_1_sexe * sexe[n];

    // (2) First observation: skip transition
    if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
      y_rep[n] = 1;  // arbitrary initial state
      int i0 = idx_first[n];
      log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[i0];
      log_lambda_2_1_out[n] = log21 + Y_hat_2_1_out[i0];
      continue;
    }

    // (3) Init transition matrix
    matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];

    // (A) First subinterval
    if (first_subinterval[n] > 0) {
      int i1 = idx_first[n];
      real s12 = Y_hat_1_2_out[i1];
      real s21 = Y_hat_2_1_out[i1];
      real t0 = date_use[n - 1];
      real t1 = t0 + first_subinterval[n] - 1;

      if (t0 < t_star1 && t_star1 < t1) {
        real d1a = t_star1 - t0;
        real d1b = t1 - t_star1 + 1;
        P_total *= transition_matrix(d1a,
          exp(log12 + beta_season_1_2 * s12),
          exp(log21 + beta_season_2_1 * s21));
        int flag1 = 1;
        int flag2 = ((t_star1 + t1) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(d1b,
          exp(log12 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12),
          exp(log21 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21));
      } else {
        real midpoint = (t0 + t1) / 2;
        int flag1 = (midpoint >= t_star1) ? 1 : 0;
        int flag2 = (midpoint >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(first_subinterval[n],
          exp(log12 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12),
          exp(log21 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21));
      }
    }

    // (B) Middle subintervals
    for (m in 1:num_middle_subintervals[n]) {
      int im = idx_middle[n, m];
      real s12m = Y_hat_1_2_out[im];
      real s21m = Y_hat_2_1_out[im];
      real t0m = global_interval_start + (global_interval_index_start[n] + m - 1) * interval_length;
      real t1m = t0m + interval_length - 1;

      if (t0m < t_star1 && t_star1 < t1m) {
        real d2a = t_star1 - t0m;
        real d2b = t1m - t_star1 + 1;
        P_total *= transition_matrix(d2a,
          exp(log12 + beta_season_1_2 * s12m),
          exp(log21 + beta_season_2_1 * s21m));
        int flag1 = 1;
        int flag2 = ((t_star1 + t1m) / 2 >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(d2b,
          exp(log12 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12m),
          exp(log21 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21m));
      } else {
        real midpoint = (t0m + t1m) / 2;
        int flag1 = (midpoint >= t_star1) ? 1 : 0;
        int flag2 = (midpoint >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(interval_length,
          exp(log12 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12m),
          exp(log21 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21m));
      }
    }

    // (C) Last subinterval
    if (last_subinterval[n] > 0) {
      int il = idx_last[n];
      real s12l = Y_hat_1_2_out[il];
      real s21l = Y_hat_2_1_out[il];
      real t1l = date_use[n];
      real t0l = t1l - last_subinterval[n] + 1;

      if (t0l < t_star2 && t_star2 < t1l) {
        real d3a = t_star2 - t0l;
        real d3b = t1l - t_star2 + 1;
        P_total *= transition_matrix(d3a,
          exp(log12 + beta_int1 + beta_season_1_2 * s12l),
          exp(log21 + beta_int1 + beta_season_2_1 * s21l));
        P_total *= transition_matrix(d3b,
          exp(log12 + beta_int1 + beta_int2 + beta_season_1_2 * s12l),
          exp(log21 + beta_int1 + beta_int2 + beta_season_2_1 * s21l));
      } else {
        real midpoint = (t0l + t1l) / 2;
        int flag2 = (midpoint >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(last_subinterval[n],
          exp(log12 + beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12l),
          exp(log21 + beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21l));
      }
    }

    // (4) Sample state
    vector[2] probs = to_vector(P_total[observed_state[n - 1]]);
    real total = sum(probs);
    if (is_nan(total) || total <= 0 || is_nan(probs[1]) || is_nan(probs[2])) {
      y_rep[n] = 0;
    } else {
      y_rep[n] = categorical_rng(probs / total);
    }

    // (5) Save log-hazards
    log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[idx_last[n]];
    log_lambda_2_1_out[n] = log21 + Y_hat_2_1_out[idx_last[n]];
  }
}
"

compiled_model <- stan_model(model_code = stan_model_full)

# OR USE THIS (MODEL STORED IN .STAN object)
#compiled_model <- rstan::stan_model(file = "two_step_model_spline.stan")

fit <- sampling(
  compiled_model,
  data = stan_data,
  iter = 1000,
  warmup = 200,
  chains = 4,
  cores = 4,
  seed = 123
)


saveRDS(file = "./Output/Validation_with_simulated_data/two_step_spline_noseason_simulated.rds")




