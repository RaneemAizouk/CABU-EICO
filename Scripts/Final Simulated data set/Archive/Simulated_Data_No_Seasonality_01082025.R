# Clear environment
rm(list = ls())

# Libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(wesanderson)
library(patchwork)

# ------------------- Setup -------------------
set.seed(122)
V <- 50  # Number of villages
houses_per_village <- 15
individuals_per_house <- 5
H <- V * houses_per_village
n_rounds <- 4
N_individuals <- H * individuals_per_house
N <- N_individuals * n_rounds

# Assign IDs
Village_base <- rep(1:V, each = houses_per_village)
village_ids_per_individual <- rep(Village_base, each = individuals_per_house)
VillageID <- rep(village_ids_per_individual, each = n_rounds)
str(VillageID)


# Assign interventions (25 intervention villages, rest control)
intervention_villages <- sample(1:V, 25)
village_intervention_map <- ifelse(1:V %in% intervention_villages, 1, 0)
Intervention <- village_intervention_map[VillageID]

# Verify VillageID and Intervention assignment
length(VillageID)       # Should be 15000 (N)
length(Intervention)    # Should be 15000 (N)
intervention_check <- data.frame(VillageID = VillageID, Intervention = Intervention) %>%
  group_by(VillageID) %>%
  summarise(intervention_status = unique(Intervention), .groups = "drop")
print(intervention_check)

# Household and individual IDs
HouseID_base <- rep(1:H, each = individuals_per_house)
menage_id_member_base <- HouseID_base * 10 + rep(1:individuals_per_house, times = H)
HouseID <- rep(HouseID_base, each = n_rounds)
menage_id_member <- rep(menage_id_member_base, each = n_rounds)
rounds <- rep(1:n_rounds, times = N_individuals)

# Individual-level attributes
age_indiv <- sample(1:90, N_individuals, replace = TRUE)
sexe_indiv <- sample(0:1, N_individuals, replace = TRUE)
age <- rep(age_indiv, each = n_rounds)
sexe <- rep(sexe_indiv, each = n_rounds)

# Time intervals between rounds (in days) - Adjusted to capture period between interventions
d1 <- 90
d2 <- 150  # Increased to cover Feb 2023 to Aug 2023
d3 <- 180
baseline_date <- as.Date("2022-10-01")
start_dates <- baseline_date + sample(0:30, N_individuals, replace = TRUE)

# Observation dates for each round
date_use <- rep(as.Date(NA), N)
for (i in 1:N_individuals) {
  base <- start_dates[i]
  idx <- ((i - 1) * n_rounds + 1):(i * n_rounds)
  noise <- sample(-20:20, 3, replace = TRUE)
  dates <- base + c(0, 
                    d1 + noise[1], 
                    d1 + d2 + noise[2], 
                    d1 + d2 + d3 + noise[3])
  date_use[idx] <- sort(dates)
}

# Validate dates
global_interval_start <- as.Date("2022-10-01")
global_interval_end <- as.Date("2024-04-01")
if (any(is.na(date_use)) || any(date_use < global_interval_start) || any(date_use > global_interval_end)) {
  stop("Invalid observation dates detected.")
}

# Convert dates to numeric
date_use_num <- as.numeric(date_use)
global_interval_start_num <- as.numeric(global_interval_start)
global_interval_end_num <- as.numeric(global_interval_end)

# Time grid for piecewise-constant intervals
interval_length <- 28
interval_starts_num <- seq(global_interval_start_num, global_interval_end_num, by = interval_length)
interval_ends_num <- pmin(interval_starts_num + interval_length - 1, global_interval_end_num)
X_midpoints <- (interval_starts_num + interval_ends_num) / 2
X_midpoints[length(X_midpoints)] <- tail(interval_starts_num, 1)
X <- (X_midpoints - X_midpoints[1]) / 365
num_data <- length(X)

# True parameters
true_params <- list(
  beta_1_2_age = 0.01,
  beta_2_1_age = -0.03,
  beta_1_2_sexe = 0.03,
  beta_2_1_sexe = -0.02,
  u_raw = seq(-1, 1, length.out = H),
  sigma_u = 0.3,
  beta_int1_1 = -1.5,   # effect of *1st* intervention on λ₁₂
  beta_int1_2 = -1.6,   # effect of *2nd* intervention on λ₁₂
  beta_int2_1 = 1.0,    # effect of *1st* intervention on λ₂₁
  beta_int2_2 = 0.15,   # effect of *2nd* intervention on λ₂₁
  sigma_q_1_2 = 0.3,
  sigma_q_2_1 = 0.3,
  q_sum = 0.02,
  alpha = 0.5,
  a1 = 0.0,
  phi=0
)

# Compute base transition rates
q_sum <- true_params$q_sum
alpha <- true_params$alpha
sigma_q_1_2 <- true_params$sigma_q_1_2
sigma_q_2_1 <- true_params$sigma_q_2_1
a1 <- true_params$a1
q_1_2_raw <- q_sum * alpha
q_2_1_raw <- q_sum * (1 - alpha)
q_1_2_base <- -3.5 + q_1_2_raw * sigma_q_1_2
q_2_1_base <- 4.71 + q_2_1_raw * sigma_q_2_1
true_params$q_1_2_raw <- q_1_2_raw
true_params$q_2_1_raw <- q_2_1_raw
true_params$q_1_2_base <- q_1_2_base
true_params$q_2_1_base <- q_2_1_base
true_params$u <- true_params$u_raw * true_params$sigma_u

# Intervention dates
global_intervention_date1 <- as.Date("2023-02-01")
global_intervention_date2 <- as.Date("2023-08-01")
noise1 <- sample(-10:10, N_individuals, replace = TRUE)
noise2 <- sample(-10:10, N_individuals, replace = TRUE)
simulated_intervention_date_orig <- rep(as.Date(NA), N_individuals)
simulated_intervention_date2_orig <- rep(as.Date(NA), N_individuals)

treated_indices <- which(village_intervention_map[village_ids_per_individual] == 1)



# Expand intervention indicator to match rounds
Intervention_expanded <- Intervention


# Assign intervention dates only for treated (intervention) villages
simulated_intervention_date_orig[treated_indices]  <- global_intervention_date1 + noise1[treated_indices]
simulated_intervention_date2_orig[treated_indices] <- global_intervention_date2 + noise2[treated_indices]

# Expand to rounds
simulated_intervention_date  <- rep(simulated_intervention_date_orig,  each = n_rounds)
simulated_intervention_date2 <- rep(simulated_intervention_date2_orig, each = n_rounds)

# Convert to numeric relative to global start
simulated_intervention_date_num  <- as.numeric(simulated_intervention_date  - global_interval_start)
simulated_intervention_date2_num <- as.numeric(simulated_intervention_date2 - global_interval_start)

# For control villages only: replace with -999
simulated_intervention_date_num[Intervention_expanded == 0]  <- -999
simulated_intervention_date2_num[Intervention_expanded == 0] <- -999

simulated_intervention_date[Intervention_expanded == 0]  <- -999
simulated_intervention_date2[Intervention_expanded == 0] <- -999

# Original intervention flags (corrected to account for -999 and control villages)
Intervention1_active <- ifelse(Intervention_expanded == 1 & simulated_intervention_date != -999 & date_use >= simulated_intervention_date, 1, 0)
Intervention2_active <- ifelse(Intervention_expanded == 1 & simulated_intervention_date2 != -999 & date_use >= simulated_intervention_date2, 1, 0)

# Interval calculations
total_days <- global_interval_end_num - global_interval_start_num + 1
num_intervals <- ceiling(total_days / interval_length)
max_middle <- num_intervals - 2
first_subinterval_sim <- rep(0, N)
last_subinterval_sim <- rep(0, N)
num_middle_subintervals_sim <- rep(0, N)
idx_first_sim <- rep(1L, N)
idx_last_sim <- rep(1L, N)
idx_middle_sim <- matrix(1L, nrow = N, ncol = max_middle)
global_interval_index_start <- rep(1L, N)

for (i in 1:N_individuals) {
  idx <- ((i - 1) * n_rounds + 1):(i * n_rounds)
  d1 <- date_use_num[idx[1]]
  d2 <- date_use_num[idx[2]]
  d3 <- date_use_num[idx[3]]
  d4 <- date_use_num[idx[4]]
  days_12 <- d2 - d1 + 1
  idx_first_temp <- floor((d1 - global_interval_start_num) / interval_length) + 1
  interval_end_first <- global_interval_start_num + idx_first_temp * interval_length - 1
  first_subinterval_sim[idx[2]] <- pmin(interval_end_first - d1 + 1, days_12, interval_length)
  total_days_23 <- d3 - d2 + 1
  num_middle_subintervals_sim[idx[3]] <- pmax(0, ceiling(total_days_23 / interval_length) - 2)
  total_days_34 <- d4 - d3 + 1
  idx_last_temp <- floor((d4 - global_interval_start_num) / interval_length) + 1
  interval_start_last <- global_interval_start_num + (idx_last_temp - 1) * interval_length
  last_subinterval_sim[idx[4]] <- pmin(d4 - interval_start_last + 1, interval_length)
}

for (n in 2:N) {
  if (menage_id_member[n] == menage_id_member[n - 1] && num_middle_subintervals_sim[n] > 0) {
    t_prev <- date_use_num[n - 1]
    raw_start_idx <- floor((t_prev - global_interval_start_num) / interval_length)
    for (m in 1:num_middle_subintervals_sim[n]) {
      interval_start_m <- global_interval_start_num + (raw_start_idx + m) * interval_length
      interval_end_m <- min(interval_start_m + interval_length - 1, global_interval_end_num)
      interval_mid_m <- (interval_start_m + interval_end_m) / 2
      raw_mid_idx <- floor((interval_mid_m - global_interval_start_num) / interval_length) + 1
      idx_middle_sim[n, m] <- max(1L, min(num_intervals, raw_mid_idx))
    }
  }
}

for (n in 1:N) {
  if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
    idx_first_sim[n] <- floor((date_use_num[n] - global_interval_start_num) / interval_length) + 1
    idx_last_sim[n] <- idx_first_sim[n]
  } else {
    idx_first_sim[n] <- floor((date_use_num[n - 1] - global_interval_start_num) / interval_length) + 1
    idx_last_sim[n] <- floor((date_use_num[n] - global_interval_start_num) / interval_length) + 1
  }
  idx_first_sim[n] <- max(1L, min(num_intervals, idx_first_sim[n]))
  idx_last_sim[n] <- max(1L, min(num_intervals, idx_last_sim[n]))
}

global_interval_index_start[1] <- floor((date_use_num[1] - global_interval_start_num) / interval_length) + 1
global_interval_index_start[1] <- max(1L, min(num_intervals, global_interval_index_start[1]))
for (n in 2:N) {
  if (menage_id_member[n] != menage_id_member[n - 1]) {
    global_interval_index_start[n] <- floor((date_use_num[n] - global_interval_start_num) / interval_length) + 1
  } else {
    t_prev <- date_use_num[n - 1]
    global_interval_index_start[n] <- floor((t_prev - global_interval_start_num) / interval_length) + 1
  }
  global_interval_index_start[n] <- max(1L, min(num_intervals, global_interval_index_start[n]))
}

# Diagnostic table
n_check <- 5
get_round_indices <- function(i) {
  ((i - 1) * n_rounds + 1):(i * n_rounds)
}
diagnostic_table <- data.frame()
for (i in 1:n_check) {
  idx <- get_round_indices(i)
  df <- data.frame(
    Individual = i,
    Round = 1:n_rounds,
    Date = as.Date(date_use_num[idx], origin = "1970-01-01"),
    Intervention_Group = Intervention[idx],
    Interv1_active = Intervention1_active[idx],
    Interv2_active = Intervention2_active[idx],
    Idx_First = idx_first_sim[idx],
    Idx_Last = idx_last_sim[idx],
    First_Sub = first_subinterval_sim[idx],
    Middle_Sub = num_middle_subintervals_sim[idx] * interval_length,
    Last_Sub = last_subinterval_sim[idx],
    Num_Middle = num_middle_subintervals_sim[idx],
    Start_Idx = global_interval_index_start[idx],
    # correct here
    # End_Idx = idx_last_sim[idx],
    stringsAsFactors = FALSE
  )
  round3_idx <- idx[3]
  if (num_middle_subintervals_sim[round3_idx] > 0) {
    mid_idx <- idx_middle_sim[round3_idx, 1:num_middle_subintervals_sim[round3_idx]]
    df$Middle_Idx <- c(NA, NA, paste(mid_idx, collapse = ","), NA)
  } else {
    df$Middle_Idx <- NA
  }
  diagnostic_table <- rbind(diagnostic_table, df)
}
print(diagnostic_table)

if (any(diagnostic_table$First_Sub < 0 | diagnostic_table$Last_Sub < 0)) {
  stop("Negative subinterval durations detected.")
}
if (any(diagnostic_table$Idx_First > num_intervals | diagnostic_table$Idx_Last > num_intervals)) {
  stop("Interval indices out of bounds.")
}

# Stan Data
stan_data <- list(
  N = N,
  H = H,
  N_persons = length(unique(menage_id_member)),
  menage_id_member = menage_id_member,
  round = rounds,
  HouseID = HouseID,
  age = age,
  sexe = sexe,
  date_use = date_use_num,
  intervention_date = as.numeric(simulated_intervention_date),
  intervention_date2 = as.numeric(simulated_intervention_date2),
  
  global_interval_start = global_interval_start_num,
  global_interval_end = global_interval_end_num,
  interval_length = interval_length,
  num_data = num_data,
  X = X,
  q_1_2_base = true_params$q_1_2_base,
  q_2_1_base = true_params$q_2_1_base,
  beta_1_2_age = true_params$beta_1_2_age,
  beta_2_1_age = true_params$beta_2_1_age,
  beta_1_2_sexe = true_params$beta_1_2_sexe,
  beta_2_1_sexe = true_params$beta_2_1_sexe,
  beta_int1_1 = true_params$beta_int1_1,
  beta_int1_2 = true_params$beta_int1_2,
  beta_int2_1 = true_params$beta_int2_1,
  beta_int2_2 = true_params$beta_int2_2,
  u = true_params$u,
  a1 = true_params$a1,
  max_middle = max_middle,
  idx_first_sim = idx_first_sim,
  idx_last_sim = idx_last_sim,
  first_subinterval_sim = pmax(0, first_subinterval_sim),
  last_subinterval_sim = pmax(0, last_subinterval_sim),
  num_middle_subintervals_sim = num_middle_subintervals_sim,
  idx_middle_sim = idx_middle_sim,
  global_interval_index_start = global_interval_index_start,
  intervention = Intervention,
  Date = as.Date(date_use, origin = "1970-01-01"),
  date_use_num = as.numeric(date_use),
  phi = true_params$phi
)

# Stan Model
stan_code <- "
functions {
  matrix transition_matrix(real t, real lambda_1_2, real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;
    real exp_term = total_lambda > 0 ? exp(-total_lambda * t) : 1.0;
    matrix[2,2] P;
    if (total_lambda == 0) {
      P[1,1] = 1.0; P[1,2] = 0.0;
      P[2,1] = 0.0; P[2,2] = 1.0;
    } else {
      P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;
      P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;
      P[1,2] = (lambda_1_2 / total_lambda) * (1 - exp_term);
      P[2,1] = (lambda_2_1 / total_lambda) * (1 - exp_term);
    }
    return P;
  }
}
data {
  int<lower=1> N;
  int<lower=1> H;
  int<lower=1> N_persons;
  int menage_id_member[N];
  int round[N];
  int HouseID[N];
  int age[N];
  int sexe[N];
  real date_use[N];
  real intervention_date[N];
  real intervention_date2[N];
  real global_interval_start;
  real global_interval_end;
  real interval_length;
  int<lower=1> num_data;
  real X[num_data];
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
  real u[H];
  real a1;
  real phi;
  int<lower=1> max_middle;
  int<lower=1> idx_first_sim[N];
  int<lower=1> idx_last_sim[N];
  real<lower=0> first_subinterval_sim[N];
  real<lower=0> last_subinterval_sim[N];
  int<lower=0> num_middle_subintervals_sim[N];
  int<lower=1> idx_middle_sim[N, max_middle];
  int<lower=1> global_interval_index_start[N];
  int<lower=0, upper=1> intervention[N];
}
generated quantities {
  vector[N] log_lambda_1_2_out;
  vector[N] log_lambda_2_1_out;
  int observed_state[N];
  int y_rep[N];

  vector[num_data] Y_hat_1_2_out;

  real log_lambda_1_2;
  real log_lambda_2_1;
 

  // Generate seasonal covariate
  for (i in 1:num_data) {
     Y_hat_1_2_out[i] = a1 * sin(2 * pi() * X[i] + phi);  
  }

  for (n in 1:N) {
    real log12_base = q_1_2_base + u[HouseID[n]] + beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n];
    real log21_base = q_2_1_base + u[HouseID[n]] + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n];
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];
   

    // Handle first observation in household
    if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
      observed_state[n] = 1;
      y_rep[n] = 1;
     

      int i0 = idx_first_sim[n];
      real s12 = Y_hat_1_2_out[i0];

      log_lambda_1_2 = log12_base;
      log_lambda_2_1 = log21_base;

      matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
      // continue;  // prevents illegal access to observed_state[n - 1]
    } else {

      // Cumulative transition matrix
      matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
         
      // --- First subinterval
if (first_subinterval_sim[n] > 0) {
  int i1 = idx_first_sim[n];
  real s12 = Y_hat_1_2_out[i1];
  real t0 = date_use[n - 1];
  real t1 = t0 + first_subinterval_sim[n] - 1;

  // In the first subinterval:
  // If the first intervention date (t_star1) lies in this interval,
  // and intervention 1 is actually applied, we split into:
  //    1. Pre-intervention (baseline only)
  //    2. Post-intervention 1 (effect of the first intervention)
  // If intervention 1 does NOT belong to this subinterval, then
  // intervention 2 also cannot happen (due to 3 months separation),
  // so we only use baseline + seasonality for the whole subinterval.
  if (t0 < t_star1 && t_star1 < t1 && intervention[n] == 1) {
    real d1a = t_star1 - t0;         // duration before first intervention
    real d1b = t1 - t_star1 + 1;     // duration after first intervention

    // --- Pre-intervention 1: baseline (no interventions)
    log_lambda_1_2 = log12_base + s12;
    log_lambda_2_1 = log21_base;
    P_total *= transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));

    // --- Post-intervention 1: add first intervention effect
    log_lambda_1_2 = log12_base + s12 + intervention[n] * beta_int1_1;
    log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1;
    P_total *= transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));

  } else {
    // --- No first intervention effect in this subinterval:
    // apply only baseline + seasonality
    log_lambda_1_2 = log12_base + s12;
    log_lambda_2_1 = log21_base;
    P_total *= transition_matrix(first_subinterval_sim[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}

      // --- Middle subintervals
      for (m in 1:num_middle_subintervals_sim[n]) {
        int im = idx_middle_sim[n, m];
        real s12m = Y_hat_1_2_out[im];
        real t0m = global_interval_start + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;
      //If the first intervention lies in the subinterval :
     //Split the interval into two parts: pre-intervention and post-intervention.
    //In post-intervention part, apply only the first intervention effect (no need to check the second intervention because 3 months separation ensures it can't happen in the same subinterval).
        if (t0m < t_star1 && t_star1 < t1m) {
        //// CASE 1: First intervention lies inside the subinterval -> split
          real d2a = t_star1 - t0m;
          real d2b = t1m - t_star1 + 1;
  // pre intervention 1(baseline )
          log_lambda_1_2 = log12_base + s12m;
          log_lambda_2_1 = log21_base;
          
          P_total *= transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
  // Post-intervention part (only first intervention effect)
    log_lambda_1_2 = log12_base + s12m + intervention[n] * beta_int1_1;
    log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1;
     P_total *= transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
     
       } else {
       // CASE 2: First intervention does NOT lie inside -> use midpoint logic
    real midpoint = (t0m + t1m) / 2;

    if (midpoint >= t_star2) {
      // Both interventions cumulative
      log_lambda_1_2 = log12_base + s12m + intervention[n] * beta_int1_1 + intervention[n] * beta_int1_2;
      log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1 +intervention[n] *  beta_int2_2;
    } else if (midpoint >= t_star1) {
      // Only first intervention
      log_lambda_1_2 = log12_base + s12m + intervention[n] * beta_int1_1;
      log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1;
    } else {
      // Baseline
      log_lambda_1_2 = log12_base + s12m;
      log_lambda_2_1 = log21_base;
    }

    P_total *= transition_matrix(interval_length, exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}
      // --- Last subinterval
      // If t_star2 (second intervention) lies within the subinterval:
// Split into two parts:
// Pre-second intervention: Apply the first intervention only .
// Post-second intervention: Apply both interventions.

// --- Last subinterval
if (last_subinterval_sim[n] > 0) {
  int il = idx_last_sim[n];
  real s12l = Y_hat_1_2_out[il];
  real t1l = date_use[n];
  real t0l = t1l - last_subinterval_sim[n] + 1;

  if (t0l < t_star2 && t_star2 < t1l) {
    // Split at second intervention
    real d3a = t_star2 - t0l;       // before second intervention
    real d3b = t1l - t_star2 + 1;   // after second intervention

    // --- Pre-second intervention: only first intervention effect
    log_lambda_1_2 = log12_base + s12l + intervention[n] * beta_int1_1;
    log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1;
    P_total *= transition_matrix(d3a, exp(log_lambda_1_2), exp(log_lambda_2_1));

    // --- Post-second intervention: both interventions
    log_lambda_1_2 = log12_base + s12l + intervention[n] * beta_int1_1 + intervention[n] * beta_int1_2;
    log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1 +intervention[n] *  beta_int2_2;
    P_total *= transition_matrix(d3b, exp(log_lambda_1_2), exp(log_lambda_2_1));

  } else {
    // No split: use midpoint check
    real midpoint = (t0l + t1l) / 2;

    if (midpoint >= t_star2) {
      // Both interventions cumulative
      log_lambda_1_2 = log12_base + s12l + intervention[n] * beta_int1_1 + intervention[n] * beta_int1_2;
      log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1 +intervention[n] *  beta_int2_2;
    } else if (midpoint >= t_star1) {
      // Only first intervention
      log_lambda_1_2 = log12_base + s12l + intervention[n] * beta_int1_1;
      log_lambda_2_1 = log21_base + intervention[n] *  beta_int2_1;
    } else {
      // Baseline (no interventions)
      log_lambda_1_2 = log12_base + s12l;
      log_lambda_2_1 = log21_base;
    }

    P_total *= transition_matrix(last_subinterval_sim[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}


    
      // --- Final update (safe after continue)
      //real midpoint_final = (date_use[n - 1] + date_use[n]) / 2;
      //int flag1_final = (intervention[n] == 1 && midpoint_final >= t_star1) ? 1 : 0;
      //int flag2_final = (intervention[n] == 1 && midpoint_final >= t_star2) ? 1 : 0;
    
    
      vector[2] probs = to_vector(P_total[observed_state[n - 1]]);
      observed_state[n] = categorical_rng(probs / sum(probs));
      y_rep[n] = observed_state[n];
    
      //int idx_log = idx_last_sim[n];
      //real s12_final = Y_hat_1_2_out[idx_log];
    }
    print(\"n: \", n, 
          \", beta_int1_1: \", beta_int1_1, \", beta_int1_2: \", beta_int1_2, 
          \", log_lambda_1_2: \", log_lambda_1_2);
    

    log_lambda_1_2_out[n] = log_lambda_1_2; // 
    log_lambda_2_1_out[n] = log_lambda_2_1; // 
  }
}


"

# Compile and run Stan model
compiled_model <- stan_model(model_code = stan_code, verbose = TRUE)

stan_fit <- sampling(
  object = compiled_model,
  data = stan_data,
  algorithm = "Fixed_param",  # REQUIRED
  iter = 1,                   # Only 1 iteration needed
  chains = 1,                 # One simulation run
  seed = 123,
  verbose = TRUE
)

# Extract simulated data
sim_data <- rstan::extract(stan_fit)
observed_state_sim <- as.vector(sim_data$y_rep)
log_lambda_1_2_out <- as.vector(sim_data$log_lambda_1_2_out)
log_lambda_2_1_out <- as.vector(sim_data$log_lambda_2_1_out)

# Create simulated_dataset
simulated_dataset <- data.frame(
  Observation = 1:N,
  MenageID = menage_id_member,
  HouseID = HouseID,
  VillageID = VillageID,
  Intervention = Intervention,
  Round = rounds,
  Date = date_use,
  date_use_num = as.numeric(date_use),
  Age = age,
  Sexe = sexe,
  Intervention_Date1 = simulated_intervention_date,
  Intervention_Date2 = simulated_intervention_date2,
  Idx_First_Sim = stan_data$idx_first_sim,
  Idx_Last_Sim = stan_data$idx_last_sim,
  First_Subinterval_Sim = stan_data$first_subinterval_sim,
  Last_Subinterval_Sim = stan_data$last_subinterval_sim,
  Num_Middle_Subintervals_Sim = stan_data$num_middle_subintervals_sim,
  Log_Lambda_1_2 = log_lambda_1_2_out,
  Log_Lambda_2_1 = log_lambda_2_1_out,
  Observed_State_Sim = observed_state_sim,
  
  stringsAsFactors = FALSE
)


# Ensure date columns are in Date format
simulated_dataset$Date <- as.Date(simulated_dataset$Date)
simulated_dataset$Intervention_Date1 <- as.Date(simulated_dataset$Intervention_Date1, origin = "1970-01-01")
simulated_dataset$Intervention_Date2 <- as.Date(simulated_dataset$Intervention_Date2, origin = "1970-01-01")

simulated_dataset$Intervention1_active <- ifelse(
  simulated_dataset$Intervention == 1 &
    !is.na(simulated_dataset$Intervention_Date1) &
    simulated_dataset$Date >= simulated_dataset$Intervention_Date1,
  1, 0
)

simulated_dataset$Intervention2_active <- ifelse(
  simulated_dataset$Intervention == 1 &
    !is.na(simulated_dataset$Intervention_Date2) &
    simulated_dataset$Date >= simulated_dataset$Intervention_Date2,
  1, 0
)

any(simulated_dataset$Intervention2_active == 1 & simulated_dataset$Intervention1_active == 0)
# Should return FALSE



# Extract log hazard outputs
log_lambda <- rstan::extract(stan_fit, pars = "log_lambda_1_2_out")$log_lambda_1_2_out
log_lambda2 <- rstan::extract(stan_fit, pars = "log_lambda_2_1_out")$log_lambda_2_1_out
log_lambda_mean <- apply(log_lambda, 2, mean)
log_lambda2_mean <- apply(log_lambda2, 2, mean)
simulated_dataset$Log_Lambda_1_2 <- log_lambda_mean
simulated_dataset$Log_Lambda_2_1 <- log_lambda2_mean
simulated_dataset$lambda_1_2 <- exp(simulated_dataset$Log_Lambda_1_2)
simulated_dataset$lambda_2_1 <- exp(simulated_dataset$Log_Lambda_2_1)

# Define group membership
treated1_only <- simulated_dataset$Intervention == 1 & simulated_dataset$Intervention2_active == 0 & simulated_dataset$Intervention1_active == 1
treated2 <- simulated_dataset$Intervention == 1 & simulated_dataset$Intervention2_active == 1 & simulated_dataset$Intervention1_active == 1

control1_only <- simulated_dataset$Intervention == 0 & simulated_dataset$Intervention2_active == 0 & simulated_dataset$Intervention1_active == 1
control2 <- simulated_dataset$Intervention == 0 & simulated_dataset$Intervention2_active == 1 & simulated_dataset$Intervention1_active == 1
table(treated1_only, treated2) # Check overlaps
summary(simulated_dataset$lambda_1_2[treated1_only])
summary(simulated_dataset$lambda_1_2[treated2])

#simulated_dataset$Intervention2_active == 1
# With this more accurate logic:
control <- simulated_dataset$Intervention == 0 
# Compute group means
control_mean <- mean(simulated_dataset$lambda_1_2[control], na.rm = TRUE)
treated1_mean <- mean(simulated_dataset$lambda_1_2[treated1_only], na.rm = TRUE)
treated2_mean <- mean(simulated_dataset$lambda_1_2[treated2], na.rm = TRUE)

# Output results
cat(" λ₁₂ Transition Rates:\n")
cat(" - Control group:              ", round(control_mean, 4), "\n")
cat(" - After 1st intervention:     ", round(treated1_mean, 4), "\n")
cat(" - After 2nd intervention:     ", round(treated2_mean, 4), "\n")
cat("\n Observation Counts:\n")

lambda_df <- data.frame(
  Type = c("lambda_1_2 Control", "lambda_1_2 Intervention 1", "lambda_1_2 Intervention 2"),
  Value = c(control_mean, treated1_mean, treated2_mean)
)

ggplot(lambda_df, aes(x = Type, y = Value, fill = Type)) +
  geom_col() +
  labs(title = "λ Transition Rates", x = "Transition", y = "Rate")

# Verify observation timing
treated_data <- simulated_dataset[simulated_dataset$Intervention == 1, ]
treated_data$in_between <- treated_data$Date >= treated_data$Intervention_Date1 &
  treated_data$Date < treated_data$Intervention_Date2
cat("Observations between 1st and 2nd intervention: ", sum(treated_data$in_between, na.rm = TRUE), "\n")

#############################
####Villages and Intervention Allocation
library(ggplot2)
village_summary <- data.frame(
  VillageID = 1:V,
  Intervention = village_intervention_map
)

ggplot(village_summary, aes(x = factor(Intervention), fill = factor(Intervention))) +
  geom_bar() +
  scale_fill_manual(values = wes_palette("GrandBudapest1", 2)) +
  labs(title = "Intervention vs Control Villages", x = "Intervention (1 = Treated)", y = "Count")

#################################Observation Dates Per Round
obs_dates <- data.frame(Date = as.Date(date_use, origin = "1970-01-01"), Round = rounds)
ggplot(obs_dates, aes(x = Date, fill = factor(Round))) +
  geom_histogram(binwidth = 30, position = "identity", alpha = 0.7) +
  labs(title = "Distribution of Observation Dates by Round", x = "Date", y = "Count")

##########Intervention Activation Over Time
library(tidyr)
library(ggplot2)

# Prepare intervention data in long format
intervention_status <- data.frame(
  Date = as.Date(date_use, origin = "1970-01-01"),
  Interv1_active = Intervention1_active,
  Interv2_active = Intervention2_active
) %>%
  pivot_longer(cols = c(Interv1_active, Interv2_active),
               names_to = "Intervention",
               values_to = "Active")

# Aggregate by date (sum of active individuals per day)
intervention_summary <- intervention_status %>%
  group_by(Date, Intervention) %>%
  summarise(Active_Count = sum(Active), .groups = "drop")

# Plot correctly
ggplot(intervention_summary, aes(x = Date, y = Active_Count, color = Intervention)) +
  geom_line() +
  labs(
    title = "Activation of Interventions Over Time",
    x = "Date",
    y = "Active Individuals"
  ) +
  scale_color_manual(values = c("Interv1_active" = "red", "Interv2_active" = "blue"))

###############
age_df <- data.frame(Age = age)
ggplot(age_df, aes(x = Age)) +
  geom_histogram(binwidth = 5, fill = "darkgreen") +
  labs(title = "Age Distribution of Individuals", x = "Age", y = "Count")
###############
# subintervals <- data.frame(
#   First = first_subinterval_sim,
#   Middle = num_middle_subintervals_sim * interval_length,
#   Last = last_subinterval_sim
# )
# 
# subintervals_long <- tidyr::pivot_longer(subintervals, cols = everything(), names_to = "Type", values_to = "Duration")
# subintervals_long <- subintervals_long[subintervals_long$Duration != 0, ]
# ggplot(subintervals_long, aes(x = Duration, fill = Type)) +
#   geom_histogram(position = "dodge", bins = 30) +
#   labs(title = "Distribution of Subinterval Durations", x = "Duration (days)", y = "Count")
#############
q_df <- data.frame(
  Type = c("q_1_2_base", "q_2_1_base"),
  Value = c(true_params$q_1_2_base, true_params$q_2_1_base)
)

ggplot(q_df, aes(x = Type, y = Value, fill = Type)) +
  geom_col() +
  labs(title = "Base Transition Rates", x = "Transition", y = "Rate")

