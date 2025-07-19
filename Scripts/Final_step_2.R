
#####Library loading 

library(ggplot2)
library(rstan)
library(tidyr)  # For pivot_longer function
library(bayesplot)
library(gridExtra)  # For arranging plots together
library(dplyr)
library(GGally)  # Load the package
library(lubridate)
# load dataset----------------------------
#load("/Users/raizouk/Desktop/Oxford/Projects/data_set/bf_esbl0123_long_all.rda")
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

# Data pre-processing --------------------------------------------------------

# 1. Create an explicit "observed_state" variable from the original carriage indicator
#    This makes it easier to recode or rename later without touching the raw "esble" column.
data$observed_state <- data$esble
data <- data %>%
  relocate(observed_state, .after = esble)

# 2. Assign each household a unique integer ID
#    factor() converts the string menage_id into a factor with levels in the order encountered,
#    and as.numeric() then maps those levels 1,2,3,… to the HouseID column.
data$HouseID <- as.numeric(factor(data$menage_id))
data <- data %>%
  relocate(HouseID, .after = menage_id)

# 3. Extract each person’s within-household number from menage_id_member
#    The pattern ".*-" strips off everything up to the dash, leaving the numeric suffix.
data$individual_number <- as.numeric(sub(".*-", "", data$menage_id_member))

# 4. Combine HouseID and individual_number into a single integer identifier:
#    multiply HouseID by 100 (reserving two digits) then add the individual number.
#    For example, HouseID=1 & individual_number=5 → 100*1 + 5 = 105
multiplier <- 100
data$menage_id_member <- data$HouseID * multiplier + data$individual_number

# 5. Parse the survey round from the RedCap event name
#    Extract the number after "round_" and before "_arm_1", then add 1 so that "round_1" → 2, etc.
data$round <- as.integer(
  gsub("round_(\\d+)_arm_1", "\\1", data$redcap_event_name)
) + 1
data <- data %>%
  relocate(round, .after = redcap_event_name)

#########################
# Handle missing data-----------------------------------------
# Removes individuals with all missing age or sexe values, then fills remaining missing values within individuals using fill(.direction = "downup").
# 1. Identify individuals with all‐missing age or sex
missing_data_individuals <- data %>%
  group_by(menage_id_member) %>%  # for each person…
  summarise(
    missing_age = all(is.na(age)),    #   is age always missing?
    missing_sexe = all(is.na(sexe))  #   is sex always missing?
  ) %>%
  filter(missing_age | missing_sexe)  # keep only those with a fully missing field

# 2. Drop those completely uninformative individuals(we dont have many of them around 4)
data_filtered <- data %>%
  filter(!menage_id_member %in% missing_data_individuals$menage_id_member)
# 3. Carry forward/backward any remaining sporadic NAs
#    (we assume age/sex are constant within a person)
data_filled <- data_filtered %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()
# 4. Check that no NAs remain for age or sex
missing_after_fill <- data_filled %>%
  filter(is.na(age) | is.na(sexe))

if (nrow(missing_after_fill) > 0) {
  print("Rows with remaining missing values:")
  print(missing_after_fill)
} else {
  print("All missing values for age and sexe have been filled.")
}
# 5. Report how many individuals survived this filtering
num_individuals_after_filling <- n_distinct(data_filled$menage_id_member)
print(paste("Number of individuals after filling missing age and sexe:", num_individuals_after_filling))
# s<- unique(length(villages))
# print(s)

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
max_middle <- num_intervals - 1

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

# Step 6: Scale X for spline stability (centered at median, in years)
X <- X_midpoints
X <- X - X[1]                  # Start from 0
X <- (X - median(X)) / 365     # Center and scale to years

num_data <- length(X)  # matches number of spline evaluation points


# Define 5 equally spaced knots from global intervals
num_knots <- 5   # Desired number of knots.

# Ensure the knots are selected from global intervals
#knots <- seq(min(X), max(X), length.out = 5)  # Define knots over global intervals
knots <- quantile(X, probs = seq(0, 1, length.out = num_knots))

# Convert knots to numeric
knots <- as.numeric(knots)

# Spline degree should be defined explicitly
spline_degree <- 3  # Example: cubic splines

# Number of basis functions
num_basis <- num_knots + spline_degree - 1  # Correctly computed basis functions

# Print to verify
print(knots)
print(num_basis)

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
# Step 7: Compute second intervention date per person (from round 3, applied to rounds > 3)
# Extract the Intervention_start_date for round 3 as the second intervention
intervention_dates2 <- data_complete %>%
  filter(round >= 3) %>%  # Focus on round 3 to get the second intervention date
  group_by(menage_id_member) %>%
  summarize(
    intervention_date2 = first(Intervention_start_date)
  )

# Step 8: Join both intervention dates back into data_complete
data_complete <- data_complete %>%
  left_join(intervention_dates, by = "menage_id_member") %>%
  left_join(intervention_dates2, by = "menage_id_member") %>%
  relocate(intervention_date, intervention_date2, .after = Intervention_start_date)

# Step 9: Check the result
data_complete %>%
  select(menage_id_member, round, Intervention_start_date, intervention_date, intervention_date2) %>%
  distinct() %>%
  head(10)
sum(is.na(data_complete$intervention_date2))
sum(is.na(data_complete$Intervention_start_date))
class(data_complete$intervention_date)
class(data_complete$intervention_date2)

##########
###########################
# (Assuming you have already computed:)
#   data_complete$intervention_date    # may contain NAs for controls
#   data_complete$intervention_date2   # may contain NAs for controls
#   data_complete$date.use             # each observation’s date
#   interval_length                     # e.g. 28

#### fill Na values in the intervention date  for control group 
# 1) Compute the last sampling date across everyone:
# 1) Compute the last sampling date across all observations:
max_obs_date <- max(data_complete$date.use, na.rm = TRUE)

# 2) Fill NA in intervention_date and intervention_date2 (both are Date already):
data_complete <- data_complete %>%
  mutate(
    intervention_date  = if_else(
      is.na(intervention_date),
      max_obs_date + interval_length,
      intervention_date
    ),
    intervention_date2 = if_else(
      is.na(intervention_date2),
      max_obs_date + interval_length,
      intervention_date2
    )
  ) %>%
  # 3) Sanity‐check: no NAs should remain
  { stopifnot(!any(is.na(.$intervention_date))); . } %>%
  { stopifnot(!any(is.na(.$intervention_date2))); . } %>%
  # 4) Convert back to numeric days‐since‐1970 for Stan
  mutate(
    intervention_date  = as.numeric(intervention_date),
    intervention_date2 = as.numeric(intervention_date2)
  )

# 5) If you also need to fill Intervention_start_date (the first‐intervention column), do likewise:
data_complete <- data_complete %>%
  mutate(
    Intervention_start_date = if_else(
      is.na(Intervention_start_date),
      max_obs_date + interval_length,
      Intervention_start_date
    )
  ) %>%
  { stopifnot(!any(is.na(.$Intervention_start_date))); . } %>%
  mutate(
    Intervention_start_date = as.Date(Intervention_start_date)
  )

# Verify:
sum(is.na(data_complete$intervention_date))    # should print 0
sum(is.na(data_complete$intervention_date2))   # should print 0
sum(is.na(data_complete$Intervention_start_date))  # should print 0


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

summary(data_complete$ date.use)



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
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  num_basis = num_basis,
  X = X,# Time variable based on global intervals
  num_data = length(X), # Number of global intervals
  max_middle = max_middle,
  Intervention_start_date = as.numeric(data_complete$Intervention_start_date),
  intervention_date  = as.numeric(data_complete$intervention_date),
  intervention_date2 = data_complete$intervention_date2
 )
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

  array[N] real Intervention_start_date;
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

  // Noise parameter for observation model
  real<lower=0> sigma;
}


transformed parameters {
  // -- Base parameters --
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
  vector[H] u = u_raw * sigma_u;

  // -- Periodic B-spline coefficients: a_raw_1_2 --
  row_vector[num_basis] a_raw_1_2;

  for (i in 1:(num_basis - 3))
    a_raw_1_2[i] = a_raw_1_2_free[i];

  for (j in 1:3)
    a_raw_1_2[num_basis - 3 + j] = a_raw_1_2_free[j];

  real tau_1_2 = exp(log_tau_raw_1_2);
  row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;
  vector[num_data] Y_hat_1_2 = to_vector(a_1_2 * B);

  // -- Periodic B-spline coefficients: a_raw_2_1 --
  row_vector[num_basis] a_raw_2_1;

  for (i in 1:(num_basis - 3))
    a_raw_2_1[i] = a_raw_2_1_free[i];

  for (j in 1:3)
    a_raw_2_1[num_basis - 3 + j] = a_raw_2_1_free[j];

  real tau_2_1 = exp(log_tau_raw_2_1);
  row_vector[num_basis] a_2_1 = a_raw_2_1 * tau_2_1;
  vector[num_data] Y_hat_2_1 = to_vector(a_2_1 * B);
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
        real   s21 = Y_hat_2_1[i1];  // seasonal adjustment at that midpoint for recovery

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
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 + beta_season_1_2 * s12),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2 + beta_season_2_1* s21)
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
            exp(lambda12 + flag1_post1*beta_int1 + flag2_post1*beta_int2 + beta_season_1_2 * s12),
            exp(lambda21 + flag1_post1*beta_int1 + flag2_post1*beta_int2 + beta_season_2_1* s21)
          );

        } else {
          // (A.2) No split by t_star1 inside this “first” piece
          real midpoint1 = (t0 + t1) / 2;
          int  flag1_1   = (midpoint1 >= t_star1) ? 1 : 0;
          int  flag2_1   = (midpoint1 >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            first_subinterval[n],
            exp(lambda12 + flag1_1*beta_int1 + flag2_1*beta_int2 + beta_season_1_2 * s12),
            exp(lambda21 + flag1_1*beta_int1 + flag2_1*beta_int2 + beta_season_2_1* s21)
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (B) MIDDLE SUBINTERVALS: intervals fully between 1st and last global bins
      // ─────────────────────────────────────────────────────────────────────────
      for (m in 1:num_middle_subintervals[n]) {
        int    im = idx_middle[n,m];
        real   s12_m = Y_hat_1_2[im];
        real   s21_m = Y_hat_2_1[im];

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
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 + beta_season_1_2 * s12_m),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2 + beta_season_2_1 * s21_m)
          );

          // — Post‐first in this same bin (first active, second not yet active):
          real midpoint_post1m = (t_star1 + t1m) / 2;
          int  flag1_post1m    = 1;
          int  flag2_post1m    = (midpoint_post1m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            d2b,
            exp(lambda12 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 + beta_season_1_2 * s12_m),
            exp(lambda21 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 + beta_season_2_1* s21_m)
          );

        } else {
          // (B.2) No split by t_star1 inside this middle bin
          real midpoint2m = (t0m + t1m) / 2;
          int  flag1_2m   = (midpoint2m >= t_star1) ? 1 : 0;
          int  flag2_2m   = (midpoint2m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            interval_length,
            exp(lambda12 + flag1_2m*beta_int1 + flag2_2m*beta_int2 + beta_season_1_2 * s12_m),
            exp(lambda21 + flag1_2m*beta_int1 + flag2_2m*beta_int2 + beta_season_2_1*s21_m)
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (C) LAST SUBINTERVAL: from the last global‐bin boundary up to date_use[n]
      // ─────────────────────────────────────────────────────────────────────────
      if (last_subinterval[n] > 0) {
        int    il = idx_last[n];
        real   s12_l = Y_hat_1_2[il];
        real   s21_l = Y_hat_2_1[il];

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
            exp(lambda12 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 + beta_season_1_2 * s12_l),
            exp(lambda21 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 +beta_season_2_1*  s21_l)
          );

          // — Post‐second piece: both interventions active —
          real flag1_post2 = 1;
          real flag2_post2 = 1;
          P_total *= transition_matrix(
            d3b,
            exp(lambda12 + flag1_post2*beta_int1 + flag2_post2*beta_int2 + beta_season_1_2 * s12_l),
            exp(lambda21 + flag1_post2*beta_int1 + flag2_post2*beta_int2 + beta_season_2_1* s21_l)
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
            exp(lambda12 + flag1_l*beta_int1 + flag2_l*beta_int2 + beta_season_1_2 * s12_l),
            exp(lambda21 + flag1_l*beta_int1 + flag2_l*beta_int2 + beta_season_2_1*s21_l)
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

  // Expose spline parameters
  row_vector[num_basis] a_1_2_out = a_1_2;
  row_vector[num_basis] a_2_1_out = a_2_1;

  // Recompute seasonal splines
  Y_hat_1_2_out = to_vector(a_1_2 * B);
  Y_hat_2_1_out = to_vector(a_2_1 * B);

  for (n in 1:N) {
    // (1) Compute baseline log-rates
    real log12 = q_1_2_base
                 + u[HouseID[n]]
                 + beta_1_2_age * age[n]
                 + beta_1_2_sexe * sexe[n];
    real log21 = q_2_1_base
                 + u[HouseID[n]]
                 + beta_2_1_age * age[n]
                 + beta_2_1_sexe * sexe[n];

    // (2) First observation of each individual: no transition
    if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
      y_rep[n] = 1;  // placeholder state
      int i0 = idx_first[n];
      log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[i0];
      log_lambda_2_1_out[n] = log21 + Y_hat_2_1_out[i0];
      continue;
    }

    // (3) Initialize transition matrix
    matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];

    // ───── (A) FIRST SUBINTERVAL ─────
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

    // ───── (B) MIDDLE SUBINTERVALS ─────
    for (m in 1:num_middle_subintervals[n]) {
      int im = idx_middle[n, m];
      real s12m = Y_hat_1_2_out[im];
      real s21m = Y_hat_2_1_out[im];
      real t0m = global_interval_start
                 + (global_interval_index_start[n] + m - 1) * interval_length;
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

    // ───── (C) LAST SUBINTERVAL ─────
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
          exp(log12 + 1 * beta_int1 + 0 * beta_int2 + beta_season_1_2 * s12l),
          exp(log21 + 1 * beta_int1 + 0 * beta_int2 + beta_season_2_1 * s21l));
        P_total *= transition_matrix(d3b,
          exp(log12 + 1 * beta_int1 + 1 * beta_int2 + beta_season_1_2 * s12l),
          exp(log21 + 1 * beta_int1 + 1 * beta_int2 + beta_season_2_1 * s21l));
      } else {
        real midpoint = (t0l + t1l) / 2;
        int flag1 = 1;
        int flag2 = (midpoint >= t_star2) ? 1 : 0;
        P_total *= transition_matrix(last_subinterval[n],
          exp(log12 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_1_2 * s12l),
          exp(log21 + flag1 * beta_int1 + flag2 * beta_int2 + beta_season_2_1 * s21l));
      }
    }
   // (4) Draw the next state
     vector[2] probs = to_vector(P_total[observed_state[n - 1]]);
     real total = sum(probs);
if (is_nan(total) || total <= 0 || is_nan(probs[1]) || is_nan(probs[2])) {
  y_rep[n] = 0;  // fallback or placeholder
} else {
  y_rep[n] = categorical_rng(probs / total);
}


    // (5) Save log-hazards for output
    log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[idx_last[n]];
    log_lambda_2_1_out[n] = log21 + Y_hat_2_1_out[idx_last[n]];
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
# 5. Launch Shinystan

library(shinystan)

# Convert and launch
shinystan_obj <- as.shinystan(stan_fit)
launch_shinystan(shinystan_obj)


#####
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Extract posterior samples
posterior_samples <- rstan::extract(stan_fit)

# Convert both seasonal effects to data frames
Y_hat_1_2_df <- as.data.frame(posterior_samples$Y_hat_1_2_out)
Y_hat_2_1_df <- as.data.frame(posterior_samples$Y_hat_2_1_out)

# Summarize posterior draws for acquisition (1→2)
Y_hat_1_2_summary <- data.frame(
  mean = apply(Y_hat_1_2_df, 2, mean),
  lwr = apply(Y_hat_1_2_df, 2, function(x) quantile(x, 0.025)),
  upr = apply(Y_hat_1_2_df, 2, function(x) quantile(x, 0.975)),
  time = stan_data$X
)

# Summarize posterior draws for decolonization (2→1)
Y_hat_2_1_summary <- data.frame(
  mean = apply(Y_hat_2_1_df, 2, mean),
  lwr = apply(Y_hat_2_1_df, 2, function(x) quantile(x, 0.025)),
  upr = apply(Y_hat_2_1_df, 2, function(x) quantile(x, 0.975)),
  time = stan_data$X
)

# Plot: Seasonal effect on Acquisition (1→2)
p1 <- ggplot(Y_hat_1_2_summary, aes(x = time, y = mean)) +
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "red") +
  labs(
    title = "Seasonal Effect on Acquisition Rate (1→2)",
    x = "Time (scaled years)",
    y = "Log Effect"
  ) +
  theme_minimal()

# Plot: Seasonal effect on Recovery (2→1)
p2 <- ggplot(Y_hat_2_1_summary, aes(x = time, y = mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "blue") +
  labs(
    title = "Seasonal Effect on Recovery Rate (2→1)",
    x = "Time (scaled years)",
    y = "Log Effect"
  ) +
  theme_minimal()

# Display plots
print(p1)
print(p2)
