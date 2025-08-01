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
  if (t0 < t_star1 && t_star1 < t1 ) {
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
#########################
#################### Plot for Plot 1: Raw count of positives.
# Plot 2: Positivity rate (positives / total).
library(dplyr)
library(ggplot2)
library(scales)

# 1. Create positivity indicator (based on Observed_State)
df <- simulated_dataset %>%
  mutate(
    positive = ifelse(Observed_State_Sim == 2, 1, 0),
    year_month = format(Date, "%Y-%m"),
    month = as.numeric(format(Date, "%m")),
    Intervention_Group = Intervention
  )

# 2. Summarize data per month & group
df_monthly <- df %>%
  group_by(year_month, month, Intervention_Group) %>%
  summarise(
    positives = sum(positive, na.rm = TRUE),
    total_samples = n(),
    positivity_rate = positives / total_samples,
    .groups = "drop"
  )

# 3. Ensure chronological order
df_monthly$year_month <- factor(df_monthly$year_month,
                                levels = unique(df_monthly$year_month[order(df_monthly$year_month)]))

# ---------------------------
# Plot 1: Raw positives
# ---------------------------
plot_positives <- ggplot(df_monthly, aes(x = year_month, y = positives, fill = factor(Intervention_Group))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = positives),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkorange"),
                    labels = c("Control", "Intervention"),
                    name = "Group") +
  labs(
    title = "Monthly Positive Counts: Control vs Intervention",
    x = "Month",
    y = "Number of Positives"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------
# Plot 2: Positivity rate
# ---------------------------
plot_positivity_rate <- ggplot(df_monthly, aes(x = year_month, y = positivity_rate, fill = factor(Intervention_Group))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = scales::percent(positivity_rate, accuracy = 0.1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "darkorange"),
                    labels = c("Control", "Intervention"),
                    name = "Group") +
  labs(
    title = "Monthly Positivity Rate: Control vs Intervention",
    x = "Month",
    y = "Positivity Rate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the two plots
plot_positives
plot_positivity_rate
####################################


###################
stan_data <- list(
  menage_id_member = as.integer(simulated_dataset$MenageID),
  HouseID = as.integer(simulated_dataset$HouseID),
  H = length(unique(simulated_dataset$HouseID)),
  VillageID = as.integer(simulated_dataset$VillageID),
  intervention = simulated_dataset$Intervention,   # Corrected spacing
  round = simulated_dataset$Round,
  age = as.integer(simulated_dataset$Age),
  sexe = simulated_dataset$Sexe,
  intervention_date = as.numeric(simulated_dataset$Intervention_Date1),  # Corrected spacing
  intervention_date2 = as.numeric(simulated_dataset$Intervention_Date2),
  
  N = N,
  
  observed_state = as.integer(observed_state_sim),
  global_interval_start = as.numeric(global_interval_start),
  global_interval_end = as.numeric(global_interval_end),
  interval_length = interval_length,
  date_use = as.numeric(simulated_dataset$date_use),
  
  num_data = length(X),
  X = X,  
  max_middle = max_middle,
  num_intervals = num_intervals  # Removed extra comma
)



stan_data$observed_state <- as.integer(observed_state_sim)
stan_data$menage_id_member   <- as.integer(simulated_dataset$MenageID)
stan_data$age                <- as.integer(simulated_dataset$Age)
stan_data$round              <- as.integer(simulated_dataset$Round)
stan_data$sexe               <- as.integer(simulated_dataset$Sexe)
stan_data$date_use           <- as.numeric(simulated_dataset$Date)
stan_data$intervention_date  <- as.numeric(simulated_dataset$Intervention_Date1)
stan_data$intervention_date2 <- as.numeric(simulated_dataset$Intervention_Date2)


str(stan_data)
str(simulated_dataset)

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
  
 }

data {
  int<lower=1> N;
  array[N] int HouseID;
  int<lower=1> H;
  array[N] int<lower=1, upper=2> observed_state;
  array[N] int menage_id_member;

  array[N] int age;
  array[N] int round;
  array[N] int sexe;
  real global_interval_start;
  real global_interval_end;
  int<lower=1> interval_length;
  array[N] real date_use;
  //array[N] int intervention_village;
  vector[N] intervention;        // intervention variable
  int<lower=1> num_data;
  array[num_data] real X;

  //array[N] real Intervention_start_date;
  array[N] real intervention_date;  // or int if stored as day numbers
  array[N] real intervention_date2; 
  int<lower=1> max_middle;
  int<lower=1> num_intervals;


}

transformed data {
  
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
      // First subinterval duration
      first_subinterval[n] = fmin(fmin(interval_end_first[n] - date1 + 1,
                                 date2 - date1 + 1),
                            interval_length * 1.0);

       // Last subinterval duration
       last_subinterval[n]  = fmin(date2 
                            - (global_interval_start 
                               + (global_interval_index_end[n] - 1) * interval_length) 
                                + 1, interval_length);

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
  real q_1_2_raw;
  real q_2_1_raw;

 
  real beta_1_2_age;
  real beta_2_1_age;
  real beta_1_2_sexe;
  real beta_2_1_sexe;

  vector[H] u_raw;

  real<lower=0> sigma_u;
  real<lower=0> sigma_q_1_2;
  real<lower=0> sigma_q_2_1;

  real<lower=0> sigma;
 

  // Replace spline parameters with sinusoid coefficients:
  real a1;  // Amplitude for sin(2piX)
  real phi;
 
  real beta_int1_1;
  real beta_int1_2;
  real beta_int2_1;
  real beta_int2_2;


}

transformed parameters {
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
  vector[H] u = u_raw * sigma_u;
  
  // Define seasonality as sinusoidal:
  vector[num_data] Y_hat_1_2;
  for (i in 1:num_data) {
    Y_hat_1_2[i] = a1 * sin(2 * pi() * X[i] + phi);  
  }


 
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
   beta_int1_1 ~ normal(0, 1);
  beta_int1_2 ~ normal(0, 1);
  beta_int2_1 ~ normal(0, 1);
  beta_int2_2 ~ normal(0, 1);
  sigma_q_1_2 ~ normal(0, 1);  
  sigma_q_2_1 ~ normal(0, 1);  
  // Household‐level random effects
  u_raw   ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);

  
  // Noise on continuous‐time transitions
  sigma ~ normal(0, 1);


  // New priors for sinusoidal coefficients:
   a1 ~ normal(0, 0.5);
  

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
      real log12_base = q_1_2_base
                     + u[HouseID[n]]
                     + beta_1_2_age  * age[n]
                     + beta_1_2_sexe * sexe[n];
      real log21_base = q_2_1_base
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
      

        // t0 = the actual calendar date of observation (n−1)
        
          real t0  = date_use[n-1];
          
          
        // t1 = end of “first subinterval,” measured in days.
        // (Because first_subinterval[n] = length of that piece in days.)
        
        
        real t1  = t0 + first_subinterval[n] - 1;

        // (A.1) If t_star1 lies strictly inside [t0, t1], we split this subinterval
        
        real log_lambda_1_2;
        real log_lambda_2_1;
        
        if (t0 < t_star1 && t_star1 < t1) {
        real d1a = t_star1 - t0;        // duration before first intervention
        real d1b = t1 - t_star1 + 1;    // duration after first intervention

         // --- Pre-intervention 1: baseline (no interventions)
           log_lambda_1_2 = log12_base + s12;
            log_lambda_2_1 = log21_base;
         P_total *= transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));


           // Post-intervention 1: add first intervention effect
          log_lambda_1_2 = log12_base + s12 + intervention[n] * beta_int1_1;
          log_lambda_2_1 = log21_base +intervention[n] *  beta_int2_1;
           P_total *= transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));

        } else {
      // No first intervention effect in this subinterval:
    // apply only baseline + seasonality as if the first intervention doesnt belong to the first subinterval the 
          log_lambda_1_2 = log12_base + s12;
           log_lambda_2_1 = log21_base;
    P_total *= transition_matrix(first_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}
      // ─────────────────────────────────────────────────────────────────────────
      // (B) MIDDLE SUBINTERVALS: intervals fully between 1st and last global bins
      // ─────────────────────────────────────────────────────────────────────────
      for (m in 1:num_middle_subintervals[n]) {
        int    im = idx_middle[n,m];
        real   s12m = Y_hat_1_2[im];

        // t0m = start of the m-th middle bin, t1m = end of that bin
        real t0m = global_interval_start
                   + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;

        // (B.1) If the first intervention date t_star1 falls inside this middle bin
        // Declare log_lambda variables at the start of the loop
         real log_lambda_1_2;
         real log_lambda_2_1;
        
        if (t0m < t_star1 && t_star1 < t1m) {
          real d2a = t_star1 - t0m;          // time before first‐intervention in this bin
          real d2b = t1m - t_star1 + 1;  // time after first‐intervention (but before second)
          
       
          // Pre‐first in this middle bin (no interventions yet/ baseline):
          
           log_lambda_1_2 = log12_base + s12m;
              log_lambda_2_1 = log21_base;
          
          P_total *= transition_matrix(
            d2a,
            exp( log_lambda_1_2 ),
            exp(log_lambda_2_1 )
          );

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
      // ─────────────────────────────────────────────────────────────────────────
      // (C) LAST SUBINTERVAL: from the last global‐bin boundary up to date_use[n]
      // ─────────────────────────────────────────────────────────────────────────
      // --- Last subinterval
      // If t_star2 (second intervention) lies within the subinterval:
// Split into two parts:
// Pre-second intervention: Apply the first intervention only .
// Post-second intervention: Apply both interventions.

// --- Last subinterval
if (last_subinterval[n] > 0) {
  int il = idx_last[n];
  real s12l = Y_hat_1_2[il];
  real t1l = date_use[n];
  real t0l = t1l - last_subinterval[n] + 1;
  real log_lambda_1_2;
  real log_lambda_2_1;

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

    P_total *= transition_matrix(last_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
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
  vector[num_data] seasonality;
for (i in 1:num_data) {
  seasonality[i] = a1 * sin(2 * pi() * X[i] + phi);
}
  }


"

compiled_model <- stan_model(model_code = stan_model_full)

# OR USE THIS (MODEL STORED IN .STAN object)
#compiled_model <- rstan::stan_model(file = "two_step_model_sine.stan")

fit <- sampling(
  compiled_model,
  data = stan_data,
  iter = 500,
  warmup = 100,
  chains = 1,
  cores = 1,
  seed = 123
)

# Extract posterior samples
posterior <- rstan::extract(fit, pars = c(
  "q_1_2_raw", "q_2_1_raw",
  "beta_1_2_age", "beta_2_1_age",
  "beta_1_2_sexe", "beta_2_1_sexe",
  "beta_int1_1", "beta_int1_2", "beta_int2_1", "beta_int2_2",
  "sigma_u", "sigma_q_1_2", "sigma_q_2_1", "sigma",
  "a1", "phi",
  "u_raw"
))
# Summary of posterior distributions
posterior_summary <- summary(fit, pars = names(posterior))$summary
print(posterior_summary)

# Extract mean and 95% credible intervals
param_means <- posterior_summary[, "mean"]
param_lower <- posterior_summary[, "2.5%"]
param_upper <- posterior_summary[, "97.5%"]

# Compute posterior for transformed parameters
q_1_2_base_posterior <- -3.5 + posterior$q_1_2_raw * posterior$sigma_q_1_2
q_2_1_base_posterior <- 4.71 + posterior$q_2_1_raw * posterior$sigma_q_2_1
u_posterior <- sweep(posterior$u_raw, 1, posterior$sigma_u, "*")

# Summarize transformed parameters
q_1_2_base_mean <- mean(q_1_2_base_posterior)
q_1_2_base_lower <- quantile(q_1_2_base_posterior, 0.025)
q_1_2_base_upper <- quantile(q_1_2_base_posterior, 0.975)

q_2_1_base_mean <- mean(q_2_1_base_posterior)
q_2_1_base_lower <- quantile(q_2_1_base_posterior, 0.025)
q_2_1_base_upper <- quantile(q_2_1_base_posterior, 0.975)

# Create comparison data frame
comparison_df <- data.frame(
  Parameter = c("q_1_2_raw", "q_2_1_raw", "q_1_2_base", "q_2_1_base",
                "beta_1_2_age", "beta_2_1_age", "beta_1_2_sexe", "beta_2_1_sexe",
                "beta_int1_1", "beta_int1_2", "beta_int2_1", "beta_int2_2",
                "sigma_u", "sigma_q_1_2", "sigma_q_2_1", "sigma", "a1", "phi"),
  True_Value = c(true_params$q_1_2_raw, true_params$q_2_1_raw, true_params$q_1_2_base, true_params$q_2_1_base,
                 true_params$beta_1_2_age, true_params$beta_2_1_age, true_params$beta_1_2_sexe, true_params$beta_2_1_sexe,
                 true_params$beta_int1_1, true_params$beta_int1_2, true_params$beta_int2_1, true_params$beta_int2_2,
                 true_params$sigma_u, true_params$sigma_q_1_2, true_params$sigma_q_2_1, 1.0, true_params$a1, true_params$phi),
  Posterior_Mean = c(param_means["q_1_2_raw"], param_means["q_2_1_raw"], q_1_2_base_mean, q_2_1_base_mean,
                     param_means["beta_1_2_age"], param_means["beta_2_1_age"], param_means["beta_1_2_sexe"], param_means["beta_2_1_sexe"],
                     param_means["beta_int1_1"], param_means["beta_int1_2"], param_means["beta_int2_1"], param_means["beta_int2_2"],
                     param_means["sigma_u"], param_means["sigma_q_1_2"], param_means["sigma_q_2_1"], param_means["sigma"],
                     param_means["a1"], param_means["phi"]),
  Lower_95CI = c(param_lower["q_1_2_raw"], param_lower["q_2_1_raw"], q_1_2_base_lower, q_2_1_base_lower,
                 param_lower["beta_1_2_age"], param_lower["beta_2_1_age"], param_lower["beta_1_2_sexe"], param_lower["beta_2_1_sexe"],
                 param_lower["beta_int1_1"], param_lower["beta_int1_2"], param_lower["beta_int2_1"], param_lower["beta_int2_2"],
                 param_lower["sigma_u"], param_lower["sigma_q_1_2"], param_lower["sigma_q_2_1"], param_lower["sigma"],
                 param_lower["a1"], param_lower["phi"]),
  Upper_95CI = c(param_upper["q_1_2_raw"], param_upper["q_2_1_raw"], q_1_2_base_upper, q_2_1_base_upper,
                 param_upper["beta_1_2_age"], param_upper["beta_2_1_age"], param_upper["beta_1_2_sexe"], param_upper["beta_2_1_sexe"],
                 param_upper["beta_int1_1"], param_upper["beta_int1_2"], param_upper["beta_int2_1"], param_upper["beta_int2_2"],
                 param_upper["sigma_u"], param_upper["sigma_q_1_2"], param_upper["sigma_q_2_1"], param_upper["sigma"],
                 param_upper["a1"], param_upper["phi"])
)

# Print comparison
print(comparison_df)

# Plot
library(ggplot2)
ggplot(comparison_df, aes(x = Parameter)) +
  geom_point(aes(y = Posterior_Mean), color = "blue", size = 3) +
  geom_errorbar(aes(ymin = Lower_95CI, ymax = Upper_95CI), width = 0.2, color = "blue") +
  geom_point(aes(y = True_Value), color = "red", shape = 17, size = 3) +
  labs(title = "Posterior Estimates vs True Parameters",
       x = "Parameter",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



summary(q_1_2_base_posterior)
summary(q_2_1_base_posterior)
summary(posterior$sigma_q_1_2)
summary(posterior$sigma_q_2_1)

