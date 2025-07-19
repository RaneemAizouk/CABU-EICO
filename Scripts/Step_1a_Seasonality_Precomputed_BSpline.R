
rm(list=ls())
#####Library loading 

library(splines2)

library(ggplot2)
library(rstan)
library(tidyr)  # For pivot_longer function
library(bayesplot)
library(gridExtra)  # For arranging plots together
library(dplyr)
library(GGally)  # Load the package
library(lubridate)
# load dataset----------------------------
load("/Users/raizouk/Desktop/bf_esbl0123_long_all.rda")
#load("./Data/latest_data_version/bf_esbl0123_long_all.rda")
data <- dfls0
#print(ls())
library(readr)
setwd("/Users/raizouk/Desktop/")
timing_interventions <- read.csv("timing_interventions.csv")
# Check
head(timing_interventions)
#############################

# Data pre-processing

data$observed_state <- data$esble
data <- data %>%
  relocate(observed_state, .after = esble)

# HouseID is the house number in which The factor() function assigns a unique integer to each distinct value of menage_id in the order it appears
data$HouseID <- as.numeric(factor(data$menage_id))
data <- data %>%
  relocate(HouseID, .after = menage_id)

### individual_number, Individual number :Create a new identifier for menage_id_member based on HouseID and individual_number.

data$individual_number <- as.numeric(sub(".*-", "", data$menage_id_member))

##Scales the HouseID by 100 to create a base value, so if the HouseID =1 and individual number is 1 the menage_id_member=101 in which last two digits are saved for the individual number as some individual have 2 digits number.
multiplier <- 100
data$menage_id_member <- data$HouseID * multiplier + data$individual_number

##### Extract the round number from the redcap_event_name
data$round <- as.integer(gsub("round_(\\d+)_arm_1", "\\1", data$redcap_event_name)) +1
data <- data %>%
  relocate(round, .after = redcap_event_name)


#########################
# Handle missing data-----------------------------------------
# Removes individuals with all missing age or sexe values, then fills remaining missing values within individuals using fill(.direction = "downup").
missing_data_individuals <- data %>%
  group_by(menage_id_member) %>%
  summarise(
    missing_age = all(is.na(age)),
    missing_sexe = all(is.na(sexe))
  ) %>%
  filter(missing_age | missing_sexe)

data_filtered <- data %>%
  filter(!menage_id_member %in% missing_data_individuals$menage_id_member)

data_filled <- data_filtered %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()

missing_after_fill <- data_filled %>%
  filter(is.na(age) | is.na(sexe))

if (nrow(missing_after_fill) > 0) {
  print("Rows with remaining missing values:")
  print(missing_after_fill)
} else {
  print("All missing values for age and sexe have been filled.")
}

num_individuals_after_filling <- n_distinct(data_filled$menage_id_member)
print(paste("Number of individuals after filling missing age and sexe:", num_individuals_after_filling))

###################################################
###If we want to use the full data set just filter out this section.
# Filter individuals who have participated in all rounds

# 1) Count how many rounds each person actually appears in
round_counts <- data_filled %>%
  group_by(menage_id_member) %>%
  summarise(rounds_participated = n_distinct(round))

# 2) Find the true maximum of that count
max_participation <- max(round_counts$rounds_participated)

# 3) Keep only individuals who hit that maximum
individuals_in_all_rounds <- round_counts %>%
  filter(rounds_participated == max_participation) %>%
  pull(menage_id_member)

# 4) Subset your main data
data_complete <- data_filled %>%
  filter(menage_id_member %in% individuals_in_all_rounds)
# 5) Check that each remaining person appears exactly max_participation times
table(data_complete$menage_id_member) %>% 
  table()   # counts how many individuals have 1 obs, 2 obs, etc

########################################

# Create a binary variable for intervention villages with date condition
data_complete <- data_complete %>%
  mutate(
    intervention_village = case_when(
      intervention.text == "intervention"  ~ 1,
      TRUE ~ 0
    )
  )


##########################

### Convert sexe and esble =observed_state to binary 0,1 for male and female respectively.

# First, ensure 'sexe' is a character type
data_complete$sexe <- as.character(data_complete$sexe)

# Now, correctly convert 'sexe' to binary (1 = Female, 0 = Male)
data_complete$sexe <- ifelse(data_complete$sexe == "Female", 1, 0)

# Confirm conversion worked
unique(data_complete$sexe)  # Should print only 0 and 1

# Recode esble to 1 or 2 only (assuming it's 0/1)
data_complete$observed_state <- ifelse(data_complete$esble == 0, 1,
                                       ifelse(data_complete$esble == 1, 2, NA))


##########################################################
# Global Interval Setup
# Defines a global time interval (Oct 2022 to Feb 2024) divided into 28-day intervals.
# Computes midpoints for spline and step function evaluation, adjusting for short final intervals.
#  Scales time (X) to years, centered at the median.
#  Sets up 5 knots for cubic B-splines (spline_degree = 3).

# Step 1: Define global interval structure
global_interval_start <- as.Date("2022-10-01")
global_interval_end <- as.Date("2024-02-19")

global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric <- as.numeric(global_interval_end)

interval_length <- 28

# Step 2: Compute number of full intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)
# Pass only the maximum possible number of full middle‐intervals to Stan
max_middle <- num_intervals - 2

# Step 3: Compute global interval starts,seq(from, to, by)
interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)
#step 4: compute global interval end, we dont need to use seq() as it will create a seqence of values for each start point while this way is vectorised and will create a vector for all the points the start intervention points sequence.
interval_ends <- interval_starts + (interval_length - 1)

# Step 5: Compute seasonal time points (X) — use midpoint of each full global interval,
# but if the last interval is short, use the start of that interval instead
X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2

# Step 6: For short last interval (e.g., <28 days), use start instead of midpoint
if ((global_interval_end_numeric - interval_starts[num_intervals]) + 1 < interval_length) {
  X_midpoints[num_intervals] <- interval_starts[num_intervals]  # Use start for last short interval
}

# 2. Normalize midpoints to [0, 1]
X <- (X_midpoints - min(X_midpoints)) / (max(X_midpoints) - min(X_midpoints))

# this table for debugging 
# Convert numeric dates back to Date format
interval_table <- tibble(
  `Interval Index` = 1:num_intervals,
  `Start Date` = as.Date(interval_starts, origin = "1970-01-01"),
  `End Date` = as.Date(interval_ends, origin = "1970-01-01"),
  `Midpoint (Date)` = as.Date(X_midpoints, origin = "1970-01-01"),
  `interval_start (numeric)` = interval_starts,
  `interval_end (numeric)` = interval_ends,
  `X_midpoints (numeric)` = X_midpoints
)

print(interval_table)


# Step 6: 
#num_knots <- 5   # Desired number of knots.

# Ensure the knots are selected from global intervals
#knots <- seq(min(X), max(X), length.out = 5)  # Define knots over global intervals
#knots <- quantile(X, probs = seq(0, 1, length.out = num_knots))

# Convert knots to numeric
#knots <- as.numeric(knots)

# Spline degree should be defined explicitly
#spline_degree <- 3  # Example: cubic splines

# Number of basis functions
#num_basis <- num_knots + spline_degree - 1  # Correctly computed basis functions

# Print to verify
#print(knots)
#print(num_basis)

# Normalize X into [0, 1] range for periodic splines
# Normalize X into [0, 1]
# X <- X_midpoints
# X <- (X - min(X)) / (max(X) - min(X))
# #For full time plot
# # 1. Rescale X back to original numeric scale (Unix date format)
# X_rescaled <- X * (max(X_midpoints) - min(X_midpoints)) + min(X_midpoints)
# 
# # 2. Convert numeric to actual Date format
# X_dates <- as.Date(X_rescaled, origin = "1970-01-01")
# 
# # 3. Add this to your season_df for plotting
# season_df$Date <- X_dates
# 
# 
# # Choose number of internal knots and spline degree
# num_knots <- 10
# spline_degree <- 3
# 
# # Define internal knots strictly inside (0,1), not including boundaries
# internal_knots <- seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)]
# 
# # Construct periodic B-spline basis
# library(splines2)
# B_periodic <- bSpline(
#   x = X,
#   knots = internal_knots,           # internal knots only
#   degree = spline_degree,
#   intercept = TRUE,
#   Boundary.knots = c(0, 1),
#   periodic = TRUE
# )
# 
# dim(B_periodic)         # [num_data, num_basis]
# head(B_periodic)        # preview basis matrix
# 
# num_data <- nrow(B_periodic)
# num_basis <- ncol(B_periodic)
# print(num_data)
# print(X)
# 1. Compute midpoints


# 3. Define internal knots (strictly inside)
num_knots <- 10
spline_degree <- 3
internal_knots <- seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)]

# 4. Create periodic B-spline basis
library(splines2)
B_periodic <- bSpline(
  x = X,
  knots = internal_knots,
  degree = spline_degree,
  intercept = TRUE,
  Boundary.knots = c(0, 1),
  periodic = TRUE
)

num_data <- nrow(B_periodic)
num_basis <- ncol(B_periodic)

# 5. Rescale X to dates (for plotting)
X_dates <- as.Date(X * (max(X_midpoints) - min(X_midpoints)) + min(X_midpoints), origin = "1970-01-01")


# X <- X - X[1]                  # Start from 0
# X <- (X - median(X)) / 365     # Center and scale to years
#num_data <- length(X)  # matches number of spline evaluation points
# Define 5 equally spaced knots from global intervals

##########
#Intervention Date Integration

# Step 1: Clean and prepare timing_interventions
timing_interventions <- timing_interventions %>%
  mutate(
    Village_name = tolower(trimws(Village_name)),
    Round = as.integer(Round),
    Intervention_start_date = dmy(Intervention_start_date)  # handles dd/mm/yyyy safely
  ) %>%
  mutate(data_round = Round + 1) %>%  # shift rounds by +1 to match data_complete
  select(Village_name, data_round, Intervention_start_date)

# Step 2: Clean and prepare data_complete
data_complete <- data_complete %>%
  mutate(
    village_name = tolower(trimws(village_name)),
    round = as.integer(round)
  )

# Step 3: Join intervention dates to data_complete by village and adjusted round
data_complete <- data_complete %>%
  left_join(
    timing_interventions,
    by = c("village_name" = "Village_name", "round" = "data_round")
  )

# Step 4: Fill round 1 (baseline) with fixed pre-intervention date
data_complete <- data_complete %>%
  mutate(
    Intervention_start_date = case_when(
      round == 1 ~ as.Date("2022-10-01"),
      TRUE ~ Intervention_start_date
    )
  )

# Step 5: Relocate the column after 'round'
data_complete <- data_complete %>%
  relocate(Intervention_start_date, .after = round)

# step6: Compute one intervention date per person: the first non-baseline date
intervention_dates <- data_complete %>%
  filter(round > 1) %>%                # ignore the artificial round-1 baseline
  group_by(menage_id_member) %>%
  summarize(
    intervention_date = first(Intervention_start_date)
  )

# Join it back in
data_complete <- data_complete %>%
  left_join(intervention_dates, by = "menage_id_member") %>%
  relocate(intervention_date, .after = round)

# Check
data_complete %>%
  select(menage_id_member, round, Intervention_start_date, intervention_date) %>%
  distinct() %>%
  head(8)

##########
#### fill Na values in the intervention date  for control group 
# 1) Compute the last sampling date across everyone:
max_obs_date <- max(data_complete$date.use, na.rm = TRUE)

# 2) For control villages (where intervention_date is NA), 
#    set intervention_date = max_obs_date + 1 interval_length
data_complete <- data_complete %>%
  mutate(
    intervention_date = as.Date(intervention_date),
    intervention_date = if_else(
      is.na(intervention_date),
      max_obs_date + interval_length,  # e.g. 28 days past last date
      intervention_date
    )
  )

# 3) Check that there are no more NAs
stopifnot(!any(is.na(data_complete$intervention_date)))

# 4) Re‐convert to numeric days since epoch for Stan
data_complete$intervention_date <- as.numeric(data_complete$intervention_date)


max_date <- max(data_complete$date.use, na.rm = TRUE)
# fill the dates of the first intervention for control villages with big number to get zero
data_complete <- data_complete %>%
  mutate(
    Intervention_start_date = case_when(
      is.na(Intervention_start_date) ~ max_date + interval_length,  # instead of crazy 9999
      TRUE ~ Intervention_start_date
    )
  )



sum(is.na(data_complete$Intervention_start_date))


# Make sure Intervention_start_date is in Date format first
data_complete$Intervention_start_date <- as.Date(data_complete$Intervention_start_date)




##########################
###Conversion to stan units 


# Ensure HouseID is correctly converted to an integer sequence
data_complete$HouseID <- as.integer(factor(data_complete$HouseID))
# Scale age before creating the stan_data list
data_complete$age_scaled <- as.integer(data_complete$age - mean(data_complete$age)) / sd(data_complete$age)

data_complete$month <- as.integer(data_complete$month)

# Ensure all necessary conversions are done

data_complete$date.use <- as.Date(data_complete$date.use)  # Ensure it's in Date format
date_use <- as.numeric(data_complete$date.use)  # Convert to numeric

# Convert categorical variables
data_complete$round <- as.integer(data_complete$round)
data_complete$intervention_village <- as.integer(data_complete$intervention_village)

##################################
####Check why im using x_midpoints for global intervals and then later in stan we are using mid point for subintervals one of them is enough.
#### Do I need the intervention definition to be 0 or 1?
######
# Define the final stan_data list
stan_data <- list(
  N = nrow(data_complete),
  observed_state = as.integer(data_complete$observed_state),
  menage_id_member = data_complete$menage_id_member,
  HouseID = data_complete$HouseID,
  H = length(unique(data_complete$HouseID)),
  age = as.integer(data_complete$age),
  round = as.integer(data_complete$round),
  sexe = as.integer(data_complete$sexe),
  intervention_village = as.integer(data_complete$intervention_village),
  
  # Interval length and date range
  num_intervals = num_intervals,  # Pass precomputed value
  global_interval_start = global_interval_start_numeric,
  global_interval_end = global_interval_end_numeric,
  interval_length = interval_length,
  
  # Date use (converted to numeric)
  date_use = as.numeric(as.Date(data_complete$date.use)),
  
  # Spline settings (now based on global intervals)
  #num_knots = num_knots,
  #knots = knots,
  spline_degree = spline_degree,
  
  X = X,# Time variable based on global intervals
  num_data = length(X), # Number of global intervals
  max_middle = max_middle,
  Intervention_start_date = as.numeric(data_complete$Intervention_start_date),
  intervention_date  = as.numeric(data_complete$intervention_date),


  num_basis = ncol(B_periodic),
  B = t(B_periodic)  # Stan expects [num_basis, num_data]
  
  
)
stan_data$num_knots <- length(internal_knots)
stan_data$knots <- internal_knots  # must be a numeric vector

# Convert Intervention_start_date to numeric (days since 1970-01-01)

# Check structure
str(stan_data)



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
int<lower=1> num_data;
  int<lower=1> num_basis;
  vector[num_data] X;
  matrix[num_basis, num_data] B;
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
 
  array[N] real Intervention_start_date;
  array[N] real intervention_date;  // or int if stored as day numbers
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
  // 3) **new**: one global‐bin index per transition
  array[N] int    bin_global;

  // Fill every position with “1” so that unused slots remain in range
  idx_first  = rep_array(1, N);
  idx_middle = rep_array(1, N, max_middle);
  idx_last   = rep_array(1, N);
  //new**
  bin_global = rep_array(1, N);

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
      // — **new**: compute **one** global‐bin index at the midpoint of the **whole** transition
      {
        real mid_all = 0.5 * (date1 + date2);
        int raw_bin = to_int( ceil((mid_all - global_interval_start)
                                   / interval_length) );
        bin_global[n] = max(1, min(num_data, raw_bin));
      }  // close inner bin_global block
    }    // close if (menage_id_member … )
 }      // close for (n in 2:N)
}        // close transformed data block

 

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

  // Infection spline
  row_vector[num_basis] a_raw_1_2;
  real a0_raw_1_2;

  real log_tau_raw_1_2;

  // Recovery spline
  row_vector[num_basis] a_raw_2_1;


  real log_tau_raw_2_1;
  // Linear trend components added to seasonal splines

  real<lower=0> sigma_a0_1_2;

  real a0_raw_2_1;
  real<lower=0> sigma_a0_2_1;

  real<lower=0> sigma;
  real beta_int1;  // intervention effect on acquetion rate
  real beta_int2;//intervention effect on decolonisation rate


}

transformed parameters {
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
  vector[H] u = u_raw * sigma_u;

  real tau_1_2 = exp(log_tau_raw_1_2);
  row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;
  real a0_1_2 = a0_raw_1_2 * sigma_a0_1_2;

  vector[num_data] Y_hat_1_2;
  Y_hat_1_2 = a0_1_2 * to_vector(X) + to_vector(a_1_2 * B);

}

  

model {
  // Priors for model parameters
  q_1_2_raw ~ normal(0, 0.3);
  q_2_1_raw ~ normal(0, 0.38);
 
  
  beta_1_2_age ~ normal(0, 1);
  beta_2_1_age ~ normal(0, 1);
  beta_1_2_sexe ~ normal(0, 1);
  beta_2_1_sexe ~ normal(0, 1);

  u_raw ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);
  sigma_q_1_2 ~ normal(0, 0.5); // Weakly informative prior
  sigma_q_2_1 ~ normal(0, 0.5);
  
 // added spline parameters
 
  a_raw_1_2 ~ normal(0, 1);
  a_raw_2_1 ~ normal(0, 1);

  a0_raw_1_2    ~ normal(0, 0.5);
  a0_raw_2_1    ~ normal(0, 0.5);
  sigma_a0_1_2  ~ normal(0, 0.1);
  sigma_a0_2_1  ~ normal(0, 0.05);




log_tau_raw_1_2 ~ normal(0, 1);
log_tau_raw_2_1 ~ normal(0, 2);


  sigma ~ normal(0, 1);


  beta_int1 ~ normal(0, 1);  // Prior for acquisition effect
  beta_int2 ~ normal(0, 1);  // Prior for decolonization effect
  
 

// ——————————————————————————————————————————————————————
// Model block: likelihood contribution for all observations
// ——————————————————————————————————————————————————————

  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n-1]) {
      // 1) Baseline log‐rates
      real lambda12 = q_1_2_base
                     + u[HouseID[n]]
                     + beta_1_2_age  * age[n]
                     + beta_1_2_sexe * sexe[n];
      real lambda21 = q_2_1_base
                     + u[HouseID[n]]
                     + beta_2_1_age  * age[n]
                     + beta_2_1_sexe * sexe[n];

      matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));
      real t_star = intervention_date[n];
      // One seasonal effect for the acquestion transition:
       real s12 = Y_hat_1_2[bin_global[n]];


      // — First subinterval —
      if (first_subinterval[n] > 0) {
        // use the precomputed midpoint‐index
        int  i1 = idx_first[n];
       //   Single‐bin seasonal effect

        real dt1, dt2;
        real t0 = date_use[n-1];
        real t1 = t0 + first_subinterval[n] - 1;

        if (t0 < t_star && t_star < t1) {
          dt1 = t_star - t0;
          dt2 = t1 - t_star + 1;
          // pre
          P_total *= transition_matrix(
            dt1,
            exp(lambda12 + /*pre*/0 * beta_int1 +  s12),
            exp(lambda21 + /*pre*/0 * beta_int2 )
          );
          // post
          P_total *= transition_matrix(
            dt2,
            exp(lambda12 + /*post*/1 * beta_int1 + s12),
            exp(lambda21 + /*post*/1 * beta_int2 )
          );
        } else {
          int ie = ((t0 + t1) / 2 >= t_star) ? 1 : 0;
          P_total *= transition_matrix(
            first_subinterval[n],
            exp(lambda12 + ie * beta_int1 +  s12),
            exp(lambda21 + ie * beta_int2 )
          );
        }
      }

      // — Middle subintervals —
      for (m in 1:num_middle_subintervals[n]) {
        int im = idx_middle[n, m];
       

        real t0 = global_interval_start
                  + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1 = t0 + interval_length - 1;

        if (t0 < t_star && t_star < t1) {
          real dt1 = t_star - t0;
          real dt2 = t1     - t_star + 1;
          // split
          P_total *= transition_matrix(
            dt1,
            exp(lambda12 + 0 * beta_int1 + s12),
            exp(lambda21 + 0 * beta_int2 )
          );
          P_total *= transition_matrix(
            dt2,
            exp(lambda12 + 1 * beta_int1 +  s12),
            exp(lambda21 + 1 * beta_int2 )
          );
        } else {
          int ie = ((t0 + t1) / 2 >= t_star) ? 1 : 0;
          P_total *= transition_matrix(
            interval_length,
            exp(lambda12 + ie * beta_int1 + s12),
            exp(lambda21 + ie * beta_int2 )
          );
        }
      }

      // — Last subinterval —
      if (last_subinterval[n] > 0) {
        int il = idx_last[n];
    
        real t1 = date_use[n];
        real t0 = t1 - last_subinterval[n] + 1;

        if (t0 < t_star && t_star < t1) {
          real dt1 = t_star - t0;
          real dt2 = t1 - t_star + 1;
          // use the two pieces as above, but each with the same s12_l/s21_l
          P_total *= transition_matrix(
            dt1,
            exp(lambda12 + /*pre*/0 * beta_int1 +  s12),
            exp(lambda21 + /*pre*/0 * beta_int2 )
          );
          P_total *= transition_matrix(
            dt2,
            exp(lambda12 + /*post*/1 * beta_int1 + s12),
            exp(lambda21 + /*post*/1 * beta_int2 )
          );
        } else {
          int ie = ((t0 + t1) / 2 >= t_star) ? 1 : 0;
          P_total *= transition_matrix(
            last_subinterval[n],
            exp(lambda12 + ie * beta_int1 +  s12),
            exp(lambda21 + ie * beta_int2 )
          );
        }
      }

      // 4) add to log‐likelihood
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

  int              y_rep[N];

  // Expose spline parameters

  row_vector[num_basis] a_1_2_out  = a_1_2;


  // Recompute seasonal splines
  Y_hat_1_2_out = a0_1_2 * to_vector(X)+ to_vector(a_1_2 * B);
 
  for (n in 1:N) {
    real log12 = q_1_2_base + u[HouseID[n]] 
               + beta_1_2_age  * age[n] 
               + beta_1_2_sexe * sexe[n];
    real log21 = q_2_1_base + u[HouseID[n]] 
               + beta_2_1_age  * age[n] 
               + beta_2_1_sexe * sexe[n];

    // First ever obs: no transition, just record baseline rate
    if (n == 1 || menage_id_member[n] != menage_id_member[n-1]) {
      y_rep[n] = -1;
      // pick the first‐interval seasonal index
      int i0 = idx_first[n];
      log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[i0];
      log_lambda_2_1_out[n] = log21 ;
      continue;
    }

    // Otherwise build the same P_total as in the model
    matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));
    real t_star  = intervention_date[n];

    //  SINGLE‐BIN seasonal effect for the entire transition 
    real s12 = Y_hat_1_2_out[bin_global[n]];


    // 1) FIRST SUBINTERVAL
    if (first_subinterval[n] > 0) {
      real t0 = date_use[n-1];
      real t1 = t0 + first_subinterval[n] - 1;
      if (t0 < t_star && t_star < t1) {
        real dt1 = t_star - t0;
        real dt2 = t1     - t_star + 1;
        P_total *= transition_matrix(
          dt1,
          exp(log12 + 0 * beta_int1 +  s12),
          exp(log21 + 0 * beta_int2 )
        );
        P_total *= transition_matrix(
          dt2,
          exp(log12 + 1 * beta_int1+  s12),
          exp(log21 + 1 * beta_int2 )
        );
      } else {
        int ie = ((t0 + t1)/2 >= t_star) ? 1 : 0;
        P_total *= transition_matrix(
          first_subinterval[n],
          exp(log12 + ie * beta_int1 +  s12),
          exp(log21 + ie * beta_int2 )
        );
      }
    }

    // 2) MIDDLE SUBINTERVALS
    for (m in 1:num_middle_subintervals[n]) {
      real t0 = global_interval_start
                + (global_interval_index_start[n] + m - 1) * interval_length;
      real t1 = t0 + interval_length - 1;
      if (t0 < t_star && t_star < t1) {
        real dt1 = t_star - t0;
        real dt2 = t1     - t_star + 1;
        P_total *= transition_matrix(
          dt1,
          exp(log12 + 0 * beta_int1 +  s12),
          exp(log21 + 0 * beta_int2 )
        );
        P_total *= transition_matrix(
          dt2,
          exp(log12 + 1 * beta_int1 +  s12),
          exp(log21 + 1 * beta_int2 )
        );
      } else {
        int ie = ((t0 + t1)/2 >= t_star) ? 1 : 0;
        P_total *= transition_matrix(
          interval_length,
          exp(log12 + ie * beta_int1 +  s12),
          exp(log21 + ie * beta_int2 )
        );
      }
    }

    // 3) LAST SUBINTERVAL
    if (last_subinterval[n] > 0) {
      real t1 = date_use[n];
      real t0 = t1 - last_subinterval[n] + 1;
      if (t0 < t_star && t_star < t1) {
        real dt1 = t_star - t0;
        real dt2 = t1     - t_star + 1;
        P_total *= transition_matrix(
          dt1,
          exp(log12 + 0 * beta_int1 + s12),
          exp(log21 + 0 * beta_int2 )
        );
        P_total *= transition_matrix(
          dt2,
          exp(log12 + 1 * beta_int1 +   s12),
          exp(log21 + 1 * beta_int2 )
        );
      } else {
        int ie = ((t0 + t1)/2 >= t_star) ? 1 : 0;
        P_total *= transition_matrix(
          last_subinterval[n],
          exp(log12 + ie * beta_int1 +  s12),
          exp(log21 + ie * beta_int2 )
        );
      }
    }

    // Draw the next state
    y_rep[n] = categorical_rng(
      to_vector(P_total[ observed_state[n-1] ]) 
      / sum(P_total[ observed_state[n-1] ])
    );

    // Record the log‐rates at the very end of this transition:
    log_lambda_1_2_out[n] = log12 +  s12;
    log_lambda_2_1_out[n] = log21 ;
  }
}



"



# 1. Compile
stan_model <- stan_model(model_code = stan_model_full)

# 2. Sample (4 chains, 4 cores)
stan_fit <- sampling(
  object      = stan_model,
  data        = stan_data,
  iter        = 500,
  warmup      = 250,
  chains      = 1,
  #cores       = 4,
  control     = list(adapt_delta = 0.9995, max_treedepth = 20)
)

# 3. Quick summary
print(stan_fit)
print(summary(stan_fit)$summary)

# 4. Extract parameter names
params <- names(rstan::extract(stan_fit, permuted = TRUE))
print(params)

# 5. Launch Shinystan
if (!requireNamespace("shinystan", quietly = TRUE)) {
  install.packages("shinystan")
}
library(shinystan)

# Convert and launch
shinystan_obj <- as.shinystan(stan_fit)
launch_shinystan(shinystan_obj)


###############
#### full study period 
post <- rstan::extract(stan_fit)
a_draws <- post$a_1_2_out       # [iterations, num_basis]
a0_draws <- post$a0_1_2         # [iterations]
B <- t(stan_data$B)             # [num_data, num_basis]

# Compute Y_hat across all posterior samples
Y_hat_draws <- a_draws %*% t(B)  # [250, 19]
Y_hat_draws <- Y_hat_draws + matrix(a0_draws, nrow = length(a0_draws), ncol = ncol(Y_hat_draws))


# Posterior summary: mean, 95% CI
season_stats <- t(apply(Y_hat_draws, 2, function(x) {
  c(
    Mean  = mean(x, na.rm = TRUE),
    Lower = quantile(x, 0.025, na.rm = TRUE),
    Upper = quantile(x, 0.975, na.rm = TRUE)
  )
}))
colnames(season_stats) <- c("Mean", "Lower", "Upper")

# X_dates already defined from X_midpoints
season_df <- data.frame(
  Date  = X_dates,
  Mean  = season_stats[, "Mean"],
  Lower = season_stats[, "Lower"],
  Upper = season_stats[, "Upper"]
)

# Plot
library(ggplot2)
ggplot(season_df, aes(x = Date, y = Mean)) +
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "skyblue", alpha = 0.3) +
  labs(
    title = expression("Posterior Seasonality of Infection Rate " ~ lambda[12]),
    x = "Date",
    y = expression(hat(Y)[12] ~ "(Seasonal Component)")
  ) +
  theme_minimal()



########one year 
library(rstan)
class(stan_fit)
# For example, extract Y_hat_1_2_out
post <- rstan::extract(stan_fit)

# 1. Extract posterior matrix
Y_hat_1_2_mat <- post$Y_hat_1_2_out  # [iterations, num_data]

# 2. Compute posterior mean and 95% CI across iterations
season_stats <- t(apply(Y_hat_1_2_mat, 2, function(x) {
  c(mean = mean(x), lwr = quantile(x, 0.025), upr = quantile(x, 0.975))
}))
colnames(season_stats) <- c("mean", "lwr", "upr")

season_df <- data.frame(
  X = stan_data$X,
  Mean = season_stats[, "mean"],
  Lower = season_stats[, "lwr"],
  Upper = season_stats[, "upr"]
)



ggplot(season_df, aes(x = X, y = Mean)) +
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "skyblue") +
  labs(
    title = expression("Posterior Seasonality of Infection Rate " ~ lambda[12]),
    x = "Time Midpoint (X)",
    y = expression(hat(Y)[12] ~ "(Seasonal Component)")
  ) +
  theme_minimal()





















#######
# Extract parameter names from the Stan fit
params <- rstan::extract(stan_fit, permuted = TRUE)
param_names <- names(params)

# Show unique top-level parameters (remove indexing like [1], [2], etc.)
unique_param_names <- unique(gsub("\\[.*\\]", "", param_names))

# Print the unique parameter names
print(unique_param_names)
#Extract the values of the parameters from S1_onestep.rds
S1_onestep <- readRDS("/Users/raizouk/Desktop/S1_onestep.rds")
class(S1_onestep)  # Should return "stanfit" (for rstan) or similar
library(rstan)  # Ensure rstan is loaded
post_mat <- as.matrix(S1_onestep)
# Mean for each parameter
means <- colMeans(post_mat)

# Standard deviation for each parameter
sds <- apply(post_mat, 2, sd)

# 95% credible intervals (2.5th and 97.5th percentiles)
ci_95 <- t(apply(post_mat, 2, quantile, probs = c(0.025, 0.975)))
summary_df <- data.frame(
  Parameter = colnames(post_mat),
  Mean = means,
  SD = sds,
  CI_2.5 = ci_95[, 1],
  CI_97.5 = ci_95[, 2]
)
print(summary_df)
############################
####Save the results
# Save the compiled stan_model object
save(stan_model, file = "stan_model.rda")

# Save the fitted model object
save(fit_test, file = "fit_test.rda")

# Save multiple objects (e.g., data and summaries)
save(stan_data_small, spline_summary, file = "stan_data_summary.rda")



################
##checking code
library(dplyr)
library(lubridate)

# ---- Step 1–6: Define global interval structure ----
global_interval_start <- as.Date("2022-10-01")
global_interval_end <- as.Date("2024-02-19")

global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric <- as.numeric(global_interval_end)

interval_length <- 28
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)
max_middle <- num_intervals - 2

interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)

interval_ends <- interval_starts + (interval_length - 1)

X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2

if ((global_interval_end_numeric - interval_starts[num_intervals]) + 1 < interval_length) {
  X_midpoints[num_intervals] <- interval_starts[num_intervals]
}

X <- (X_midpoints - min(X_midpoints)) / (max(X_midpoints) - min(X_midpoints))
num_data <- length(X)

# ---- Simulated individual observation data ----
date_use <- as.numeric(as.Date(c("2022-10-01", "2022-10-20", "2022-11-10", "2022-12-01", "2022-12-25")))
menage_id_member <- c(1, 1, 1, 2, 2)
N <- length(date_use)

# ---- Stan-style logic replicated in R ----
results <- list()

for (n in 2:N) {
  if (menage_id_member[n] == menage_id_member[n - 1]) {
    date1 <- date_use[n - 1]
    date2 <- date_use[n]
    
    raw_start <- floor((date1 - global_interval_start_numeric) / interval_length)
    raw_end   <- ceiling((date2 - global_interval_start_numeric) / interval_length)
    
    global_interval_index_start <- max(1, min(num_intervals, raw_start + 1))
    global_interval_index_end   <- max(1, min(num_intervals, raw_end))
    num_middle_subintervals     <- max(0, global_interval_index_end - global_interval_index_start - 1)
    
    interval_start_first <- global_interval_start_numeric + (global_interval_index_start - 1) * interval_length
    interval_end_first   <- interval_start_first + interval_length - 1
    
    interval_start_last  <- global_interval_start_numeric + (raw_end - 1) * interval_length
    interval_end_last    <- date2
    
    first_subinterval  <- interval_end_first - date1 + 1
    last_subinterval   <- date2 - (global_interval_start_numeric + (global_interval_index_end - 1) * interval_length) + 1
    middle_subinterval <- num_middle_subintervals * interval_length
    
    # Midpoint-based indexing
    mid1 <- (date1 + interval_end_first) / 2
    raw_mid1 <- floor((mid1 - global_interval_start_numeric) / interval_length) + 1
    idx_first <- max(1, min(num_data, raw_mid1))
    
    idx_middle <- c()
    for (m in 1:num_middle_subintervals) {
      global_idx <- global_interval_index_start + m
      ms <- global_interval_start_numeric + (global_idx - 1) * interval_length
      me <- ms + interval_length - 1
      midm <- (ms + me) / 2
      raw_midm <- floor((midm - global_interval_start_numeric) / interval_length) + 1
      idx_middle[m] <- max(1, min(num_data, raw_midm))
    }
    
    ml <- (interval_start_last + date2) / 2
    raw_ml <- floor((ml - global_interval_start_numeric) / interval_length) + 1
    idx_last <- max(1, min(num_data, raw_ml))
    
    mid_all <- (date1 + date2) / 2
    raw_bin <- ceiling((mid_all - global_interval_start_numeric) / interval_length)
    bin_global <- max(1, min(num_data, raw_bin))
    
    results[[length(results) + 1]] <- tibble(
      n = n,
      global_interval_index_start = global_interval_index_start,
      global_interval_index_end = global_interval_index_end,
      num_middle_subintervals = num_middle_subintervals,
      first_subinterval = first_subinterval,
      middle_subinterval = middle_subinterval,
      last_subinterval = last_subinterval,
      idx_first = idx_first,
      idx_middle = paste(idx_middle, collapse = ","),
      idx_last = idx_last,
      bin_global = bin_global,
      global_interval_start=global_interval_start,
      global_interval_end =global_interval_end
    )
  }
}

# Combine and view results
result_df <- bind_rows(results)
print(result_df)
