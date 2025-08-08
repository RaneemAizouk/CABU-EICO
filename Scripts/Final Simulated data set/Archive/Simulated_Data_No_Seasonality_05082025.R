# Clear environment

rm(list = ls())

# Libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(wesanderson)
library(patchwork)
library(tidyr)
library(lubridate)
library(tibble)
library(purrr)
library(gridExtra)
library(posterior)

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
table(intervention_check$intervention_status)

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
  beta_1_2_age = 0.01,         # Age effect on transition from state 1 to 2 (acquisition)
  beta_2_1_age = -0.03,        # Age effect on transition from state 2 to 1 (clearance)
  beta_1_2_sexe = 0.03,        # Sex effect on transition from 1 to 2
  beta_2_1_sexe = -0.02,       # Sex effect on transition from 2 to 1
  u_raw = seq(-1, 1, length.out = H),  # Raw household-level random effects (length H)
  sigma_u = 0.3,               # Standard deviation of household random effects
  beta_int1_1 = -1.5,          # Effect of the 1st intervention on transition rate λ₁₂
  beta_int1_2 = -1.6,          # Effect of the 2nd intervention on transition rate λ₁₂
  beta_int2_1 = 1.0,           # Effect of the 1st intervention on transition rate λ₂₁
  beta_int2_2 = 0.15,          # Effect of the 2nd intervention on transition rate λ₂₁
  sigma_q_1_2 = 0.3,           # Standard deviation of baseline λ₁₂
  sigma_q_2_1 = 0.3,           # Standard deviation of baseline λ₂₁
  q_sum = 0.02,                # Sum of baseline transition rates (λ₁₂ + λ₂₁)
  alpha = 0.5,                 # Mixing parameter between global and individual hazard contributions
  a1 = 0                    # Amplitude of seasonal effect on λ₁₂
  
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
q_2_1_base <- -4.71 + q_2_1_raw * sigma_q_2_1
true_params$q_1_2_raw <- q_1_2_raw
true_params$q_2_1_raw <- q_2_1_raw
true_params$q_1_2_base <- q_1_2_base
true_params$q_2_1_base <- q_2_1_base
true_params$u <- true_params$u_raw * true_params$sigma_u

# Intervention dates
global_intervention_date1 <- as.Date("2023-02-01")
global_intervention_date2 <- as.Date("2023-08-01")


# Intervention dates at village level
V_intervention <- sum(village_intervention_map == 1)

# Generate village-level noise for intervention villages only
set.seed(123)
noise_village1 <- sample(-10:10, V_intervention, replace = TRUE)
noise_village2 <- sample(-10:10, V_intervention, replace = TRUE)

# Create empty vectors to store village-level intervention dates
simulated_intervention_date_orig_village  <- rep(as.Date(NA), V)
simulated_intervention_date2_orig_village <- rep(as.Date(NA), V)

# Assign noisy intervention dates only to intervention villages
intervention_village_ids <- which(village_intervention_map == 1)

simulated_intervention_date_orig_village[intervention_village_ids]  <- global_intervention_date1 + noise_village1
simulated_intervention_date2_orig_village[intervention_village_ids] <- global_intervention_date2 + noise_village2

# Expand to individuals: use each person's village ID to assign their village's intervention date
simulated_intervention_date_orig  <- simulated_intervention_date_orig_village[village_ids_per_individual]
simulated_intervention_date2_orig <- simulated_intervention_date2_orig_village[village_ids_per_individual]

# Check that each individual in a treated village got the same date as others from their village
chv = data.frame(cbind(village_ids_per_individual, as.character(simulated_intervention_date_orig),as.character(simulated_intervention_date2_orig)))
table(chv$V2,chv$village_ids_per_individual) # 75 individuals per village
table(chv$V3,chv$village_ids_per_individual) # 75 individuals per village

# Expand to rounds
simulated_intervention_date  <- rep(simulated_intervention_date_orig,  each = n_rounds)
simulated_intervention_date2 <- rep(simulated_intervention_date2_orig, each = n_rounds)

# Convert to numeric relative to global start
simulated_intervention_date_num  <- as.numeric(simulated_intervention_date  - global_interval_start)
simulated_intervention_date2_num <- as.numeric(simulated_intervention_date2 - global_interval_start)

Intervention_expanded <- Intervention

# For control villages only: replace with -999
simulated_intervention_date_num[Intervention_expanded == 0]  <- -999
simulated_intervention_date2_num[Intervention_expanded == 0] <- -999

simulated_intervention_date[Intervention_expanded == 0]  <- -999
simulated_intervention_date2[Intervention_expanded == 0] <- -999

# Original intervention flags (corrected to account for -999 and control villages)
Intervention1_active <- ifelse(Intervention_expanded == 1 & simulated_intervention_date != -999 & date_use >= simulated_intervention_date, 1, 0)
Intervention2_active <- ifelse(Intervention_expanded == 1 & simulated_intervention_date2 != -999 & date_use >= simulated_intervention_date2, 1, 0) 

table(Intervention1_active)
table(Intervention2_active)# 15000/2 intervention groups/4 rounds = 1875 observations in last round with active dates



# Interval calculations
# we took two consecutive observation dates d_start and d_end in two consecutive rounds for the same individual in which "idx" isolates that individual’s observations.
# global_end :is the last day of the interval that contains d_start
# first_subinterval_sim: The number of days from the observation date in the first round (d_start) until the end of the interval ( 28-day block) in which that date lies.
# num_middle_subintervals_sim:The number of full intervals (e.g. 28-day blocks) that exist between the end of the first subinterval and the date of the second observation (d_end).
# last_subinterval_sim: The number of days from the start of the interval in which the second observation date (d_end) lies, up to d_end itself.
# idx_first_sim: The index of the interval (28-day block) that contains the start date of the period between two consecutive observations.
# idx_last_sim: The index of the interval (28-day block) that contains the second observation date (end of the period between two rounds).
total_days <- global_interval_end_num - global_interval_start_num + 1
num_intervals <- ceiling(total_days / interval_length)
max_middle <- num_intervals - 2

split_intervals <- function(N_individuals, 
                            n_rounds, 
                            date_use_num, 
                            global_interval_start_num, 
                            global_interval_end_num, 
                            interval_length) {
  
  # Derived values
  total_days <- global_interval_end_num - global_interval_start_num + 1
  num_intervals <- ceiling(total_days / interval_length)
  max_middle <- num_intervals - 2
  N <- N_individuals * n_rounds
  
  # Initialize storage
  first_subinterval_sim <- rep(0, N)
  last_subinterval_sim <- rep(0, N)
  num_middle_subintervals_sim <- rep(0, N)
  idx_first_sim <- rep(1L, N)
  idx_last_sim <- rep(1L, N)
  idx_middle_sim <- matrix(1L, nrow = N, ncol = max_middle)
  global_interval_index_start <- rep(1L, N)
  
  # Loop through individuals and rounds
  for (i in 1:N_individuals) {
    idx <- ((i - 1) * n_rounds + 1):(i * n_rounds)
    
    for (j in 1:(n_rounds - 1)) {
      n <- idx[j + 1]
      
      d_start <- date_use_num[idx[j]]
      d_end <- date_use_num[idx[j + 1]]
      
      # Identify the interval in which the first observation lies
      interval_idx <- floor((d_start - global_interval_start_num) / interval_length)
      global_end <- global_interval_start_num + (interval_idx + 1) * interval_length - 1
      
      # FIRST subinterval
      first_subinterval_sim[n] <- pmin(global_end - d_start + 1, interval_length)
      idx_first_sim[n] <- interval_idx + 1  # direct mapping
      
      # Remaining days
      remaining_days <- d_end - global_end
      
      # MIDDLE subintervals
      if (remaining_days > 0) {
        num_middle <- remaining_days %/% interval_length
        num_middle_subintervals_sim[n] <- num_middle
        
        if (num_middle > 0) {
          for (m in 1:num_middle) {
            idx_middle_sim[n, m] <- interval_idx + 1 + m
          }
        }
        
        # LAST subinterval
        last_subinterval_sim[n] <- remaining_days %% interval_length
        idx_last_sim[n] <- interval_idx + 1 + num_middle + ifelse(last_subinterval_sim[n] > 0, 1, 0)
      } else {
        last_subinterval_sim[n] <- 0
        idx_last_sim[n] <- idx_first_sim[n]
      }
      
      # Global interval start index
      global_interval_index_start[n] <- idx_first_sim[n]
    }
  }
  
  # Return results as a list
  return(list(
    first_subinterval_sim = first_subinterval_sim,
    last_subinterval_sim = last_subinterval_sim,
    num_middle_subintervals_sim = num_middle_subintervals_sim,
    idx_first_sim = idx_first_sim,
    idx_last_sim = idx_last_sim,
    idx_middle_sim = idx_middle_sim,
    global_interval_index_start = global_interval_index_start
  ))
}

subintervals <- split_intervals(
  N_individuals = N_individuals,
  n_rounds = n_rounds,
  date_use_num = date_use_num,
  global_interval_start_num = global_interval_start_num,
  global_interval_end_num = global_interval_end_num,
  interval_length = interval_length
)

# Extract results into variables
idx_first_sim <- subintervals$idx_first_sim
idx_last_sim <- subintervals$idx_last_sim
first_subinterval_sim <- subintervals$first_subinterval_sim
last_subinterval_sim <- subintervals$last_subinterval_sim
num_middle_subintervals_sim <- subintervals$num_middle_subintervals_sim
idx_middle_sim <- subintervals$idx_middle_sim
global_interval_index_start <- subintervals$global_interval_index_start




n_check <- 10  # Number of individuals to check
diagnostic_table <- data.frame()

get_round_indices <- function(i) {
  ((i - 1) * n_rounds + 1):(i * n_rounds)
}

for (i in 1:n_check) {
  idx <- get_round_indices(i)
  
  for (r in 2:n_rounds) {
    n <- idx[r]
    n_prev <- idx[r - 1]
    
    date_n <- as.Date(date_use_num[n], origin = "1970-01-01")
    date_prev <- as.Date(date_use_num[n_prev], origin = "1970-01-01")
    days_between <- as.numeric(date_n - date_prev)
    
    # Interval breakdown
    first <- first_subinterval_sim[n]
    middle <- num_middle_subintervals_sim[n] * interval_length
    last <- last_subinterval_sim[n]
    total_from_subs <- first + middle + last
    
    # Middle interval indices
    if (num_middle_subintervals_sim[n] > 0) {
      mids <- paste(idx_middle_sim[n, 1:num_middle_subintervals_sim[n]], collapse = ",")
    } else {
      mids <- NA
    }
    
    row <- data.frame(
      Individual = i,
      Round = r,
      Date_prev = date_prev,
      Date_n = date_n,
      Days_Between = days_between,
      First_Sub = first,
      Middle_Sub = middle,
      Last_Sub = last,
      Total_From_Subs = total_from_subs,
      Match = days_between == total_from_subs,
      Idx_First = idx_first_sim[n],
      Idx_Last = idx_last_sim[n],
      Num_Middle = num_middle_subintervals_sim[n],
      Middle_Idx = mids,
      stringsAsFactors = FALSE
    )
    
    diagnostic_table <- rbind(diagnostic_table, row)
  }
}

print(diagnostic_table)

if (any(diagnostic_table$First_Sub < 0 | diagnostic_table$Last_Sub < 0)) {
  stop("Negative subinterval durations detected.")
}
if (any(diagnostic_table$Idx_First > num_intervals | diagnostic_table$Idx_Last > num_intervals)) {
  stop("Interval indices out of bounds.")
}

# First observation is defined (as can not be modelled, i.e. no previous observation) based on prevalence of 60%
n_individuals <- length(unique(menage_id_member))

# One observation per individual (0 or 1)
first_obs <- sample(c(1, 2), size = n_individuals, replace = TRUE, prob = c(0.4, 0.6))

# Create vector of NAs of length N
obs <- rep(1, N)

# Assign the first observation per individual
first_indices <- match(unique(menage_id_member), menage_id_member)
obs[first_indices] <- first_obs

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
  obs = obs, # for first observation
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
  date_use_num = as.numeric(date_use)
  
)

# CHECK DATAFRAME
stan_df <- data.frame(
  VillageID = VillageID,
  menage_id_member = stan_data$menage_id_member,
  round = stan_data$round,
  HouseID = stan_data$HouseID,
  age = stan_data$age,
  sexe = stan_data$sexe,
  date_use = date_use,
  obs = obs,
  intervention = stan_data$intervention,
  intervention_date = as.Date(simulated_intervention_date, origin = "1970-01-01"),
  intervention_date2 = as.Date(simulated_intervention_date2,origin = "1970-01-01"),
  Intervention1_active = Intervention1_active,
  Intervention2_active = Intervention2_active,
  Date = stan_data$Date
)

# Create a unique intervention period per village
village_intervention_periods <- stan_df %>%
  filter(round==4, intervention==1 ) %>%
  select(VillageID, intervention_date, intervention_date2) %>%
  distinct() %>%
  filter(!is.na(intervention_date), !is.na(intervention_date2)) %>%
  mutate(ymin = -Inf, ymax = Inf)

village_intervention_periods2 <- stan_df %>%
  filter(round == 4, intervention==1 ) %>%
  filter(!is.na(intervention_date2), !is.na(Date)) %>%
  group_by(VillageID, intervention_date2) %>%
  summarise(Date = max(Date), .groups = "drop") %>%
  mutate(ymin = -Inf, ymax = Inf)

unique(village_intervention_periods$VillageID)
unique(village_intervention_periods2$VillageID)
unique(intervention_check$VillageID[intervention_check$intervention_status==1])

p <- stan_df %>%filter(intervention==1)%>%
  ggplot(aes(x = Date, fill = as.factor(round))) +
  # Add shaded areas specific to each village
  geom_rect(data = village_intervention_periods,
            aes(xmin = intervention_date, xmax = intervention_date2,
                ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "coral", alpha=0.2)+
  geom_rect(data = village_intervention_periods2,
            aes(xmin = intervention_date2, xmax =Date,
                ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "darkgreen", alpha=0.2)+
  geom_bar()  +
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~ VillageID, ncol=5) +  # Facet per village
  labs(x = "Date", y = "Count", title = "Village Sampling Times",
       subtitle = "Intervention group", fill = "Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = "white")
  ) 

p

p2 <- stan_df %>%filter(intervention==0)%>%
  ggplot(aes(x = Date, fill = as.factor(round))) +
  
  geom_bar()  +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~ VillageID, ncol=5) +  # Facet per village
  labs(x = "Date", y = "Count", title = "Village Sampling Times",
       subtitle = "Control group", fill = "Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = "white")
  ) 

p2


# CHECK INTERVENTION ACTIVE
unique(stan_df$Intervention1_active[stan_df$intervention==0]) # should be 0
unique(stan_df$Intervention2_active[stan_df$intervention==0]) # should be 0

unique(stan_df$Intervention1_active[stan_df$intervention==1&stan_df$round==1]) # should be 0
unique(stan_df$Intervention2_active[stan_df$intervention==1&stan_df$round%in%c(1,2)]) # should be 0


stan_code <- "
functions {
  matrix transition_matrix(real t, real lambda_1_2, real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;                           // Total transition rate (acquisitions and decolonisations)
    real exp_term = total_lambda > 0 ? exp(-total_lambda * t) : 1.0;      // Exponential decay term
    matrix[2,2] P;
    if (total_lambda == 0) {
      P[1,1] = 1.0; P[1,2] = 0.0;                                          // No transition possible
      P[2,1] = 0.0; P[2,2] = 1.0;
    } else {
      P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;  // Prob staying in state 1
      P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;  // Prob staying in state 2
      P[1,2] = (lambda_1_2 / total_lambda) * (1 - exp_term);                           // Prob moving from 1 to 2
      P[2,1] = (lambda_2_1 / total_lambda) * (1 - exp_term);                           // Prob moving from 2 to 1
    }
    return P;
  }
}
data {
  int<lower=1> N;                           // Number of observations
  int<lower=1> H;                           // Number of households
  int<lower=1> N_persons;                  // Number of individuals
  int menage_id_member[N];                // Household ID for each individual
  int round[N];                            // Round number
  int HouseID[N];                          // Household ID
  int age[N];                              // Age of individual
  int sexe[N];                             // Sex of individual
  int<lower=1, upper=2> obs[N];             // For first observation
  real date_use[N];                        // Observation date
  real intervention_date[N];              // First/second intervention date
  real intervention_date2[N];             // Third intervention date
  real global_interval_start;             // Study start date
  real global_interval_end;               // Study end date
  real interval_length;                   // Length of one interval (28 days)
  int<lower=1> num_data;                  // Number of intervals (e.g., time grid size)
  real X[num_data];                        // Scaled time points (in years) for seasonal effect
  real q_1_2_base;                         // Baseline acquisition rate
  real q_2_1_base;                         // Baseline decolonisation rate
  real beta_1_2_age;                      // Age effect on acquisition
  real beta_2_1_age;                      // Age effect on decolonisation
  real beta_1_2_sexe;                     // Sex effect on acquisition
  real beta_2_1_sexe;                     // Sex effect on decolonisation
  real beta_int1_1;                       // Effect of interventions 1+2 on acquisition
  real beta_int1_2;                       // Effect of intervention 3 on acquisition
  real beta_int2_1;                       // Effect of interventions 1+2 on decolonisation
  real beta_int2_2;                       // Effect of intervention 3 on decolonisation
  real u[H];                               // Household-level random effects
  real a1;                                 // Amplitude of seasonal effect
 
  int<lower=1> max_middle;                // Max number of middle intervals
  int<lower=1> idx_first_sim[N];         // Index of first time interval per obs
  int<lower=1> idx_last_sim[N];          // Index of last time interval per obs
  real<lower=0> first_subinterval_sim[N]; // Proportion of first interval contributing
  real<lower=0> last_subinterval_sim[N];  // Proportion of last interval contributing
  int<lower=0> num_middle_subintervals_sim[N]; // Number of full intervals per obs
  int<lower=1> idx_middle_sim[N, max_middle]; // Indices of middle intervals
  int<lower=1> global_interval_index_start[N]; // Global time interval for observation start
  int<lower=0, upper=1> intervention[N];        // Whether intervention was applied (binary, i.e. intervention or control group)
}
generated quantities {
  vector[N] log_lambda_1_2_out;          // Log acquisition rate per observation
  vector[N] log_lambda_2_1_out;          // Log decolonisation rate per observation
  int<lower=1, upper=2> observed_state[N]; // Simulated colonisation state (1 or 2)
  int y_rep[N];                          // Replicated outcome for posterior predictive checks
  vector[num_data] Y_hat_1_2_out;        // Seasonal intervention effect over time (i.e. what seasonal effect to expect at the time of the observation)
  int second_intervention_used[N];       // Whether second intervention was active
  real log_lambda_1_2;                   // Placeholder for a global acquisition rate
  real log_lambda_2_1;                   // Placeholder for a global decolonisation rate
  int acquisitions[N] = rep_array(0, N);
  int decolonisations[N] = rep_array(0, N);
  int at_risk_acquisition[N] = rep_array(0, N);
  int at_risk_decolonisation[N] = rep_array(0, N);


  // Generate seasonal covariate
  # Seasonal intervention effect modeled as a sine wave with:
  # - a1: amplitude of the seasonal effect (max effect size)
  # - X[i]: time (in years) since the start of the observation period, where 0 is the start of the study period and 1.4 the end of the study period (i.e. 1.4 years)
  # - phi: phase shift (in radians), determines when the seasonal peak occurs
  # The sine function completes one full cycle per year.
  # For example, if phi = 0, the effect peaks at the start of the study.

 

  for (n in 1:N) {
    real log12_base = q_1_2_base + u[HouseID[n]] + beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n];
    real log21_base = q_2_1_base + u[HouseID[n]] + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n];
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];
   

    // Handle first observation in household
    if (n == 1 || menage_id_member[n] != menage_id_member[n - 1]) {
      observed_state[n] = obs[n];
      y_rep[n] = observed_state[n];
     

      int i0 = idx_first_sim[n];
     
      log_lambda_1_2 = log12_base;
      log_lambda_2_1 = log21_base;

      matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
      continue;  // prevents illegal access to observed_state[n - 1]
    } else {

      // Cumulative transition matrix
      matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
      
      int current_state = observed_state[n - 1];
      int next_state;
         
      // --- First subinterval
if (first_subinterval_sim[n] > 0) {
  int i1 = idx_first_sim[n]; 

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
    log_lambda_1_2 = log12_base ;
    log_lambda_2_1 = log21_base;
    
{
      matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      current_state = next_state;
      P_total *= P;
    }
    
    // --- Post-intervention 1: add first intervention effect
    log_lambda_1_2 = log12_base  + beta_int1_1;
    log_lambda_2_1 = log21_base + beta_int2_1;
    
  {
      matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }
    
  } else {
    // --- No first intervention effect in this subinterval:
    // apply only baseline + seasonality
    log_lambda_1_2 = log12_base ;
    log_lambda_2_1 = log21_base;
    
  {
      matrix[2,2] P = transition_matrix(first_subinterval_sim[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }
  }
}
// --- Middle subintervals
for (m in 1:num_middle_subintervals_sim[n]) {
  int im = idx_middle_sim[n, m];

  real t0m = global_interval_start + (global_interval_index_start[n] + m - 1) * interval_length;
  real t1m = t0m + interval_length - 1;

  if (t0m < t_star1 && t_star1 < t1m && intervention[n] == 1) {
    // CASE 1: First intervention lies inside the subinterval → split
    real d2a = t_star1 - t0m;
    real d2b = t1m - t_star1 + 1;

    // --- Pre-intervention: baseline + seasonality
    log_lambda_1_2 = log12_base ;
    log_lambda_2_1 = log21_base;
    {
      matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }

    // --- Post-intervention 1
    log_lambda_1_2 = log12_base +  beta_int1_1;
    log_lambda_2_1 = log21_base + beta_int2_1;
    {
      matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }

  } else {
    // CASE 2: First intervention does NOT lie inside → use midpoint logic
    real midpoint = (t0m + t1m) / 2;

    if (midpoint >= t_star2 && intervention[n] == 1) {
      log_lambda_1_2 = log12_base + beta_int1_1 + beta_int1_2;
      log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
    } else if (midpoint >= t_star1 && intervention[n] == 1) {
      log_lambda_1_2 = log12_base  + beta_int1_1;
      log_lambda_2_1 = log21_base + beta_int2_1;
    } else {
      log_lambda_1_2 = log12_base ;
      log_lambda_2_1 = log21_base;
    }

    {
      matrix[2,2] P = transition_matrix(interval_length, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }
  }
}
// --- Last subinterval
if (last_subinterval_sim[n] > 0) {
  int il = idx_last_sim[n];

  real t1l = date_use[n];
  real t0l = t1l - last_subinterval_sim[n] + 1;

  if (t0l < t_star2 && t_star2 < t1l && intervention[n] == 1) {
    // Split at second intervention (intervention group only)
    real d3a = t_star2 - t0l;       // Before second intervention
    real d3b = t1l - t_star2 + 1;   // After second intervention

    // Pre-second intervention: first intervention only
    log_lambda_1_2 = log12_base  + beta_int1_1;
    log_lambda_2_1 = log21_base + beta_int2_1;
    {
      matrix[2,2] P = transition_matrix(d3a, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }

    // Post-second intervention: both interventions
    log_lambda_1_2 = log12_base  + beta_int1_1 + beta_int1_2;
    log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
    {
      matrix[2,2] P = transition_matrix(d3b, exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }

  } else {
    // No split: use midpoint check
    real midpoint = (t0l + t1l) / 2;

    if (midpoint >= t_star2 && intervention[n] == 1) {
      log_lambda_1_2 = log12_base  + beta_int1_1 + beta_int1_2;
      log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
    } else if (midpoint >= t_star1 && intervention[n] == 1) {
      log_lambda_1_2 = log12_base  + beta_int1_1;
      log_lambda_2_1 = log21_base + beta_int2_1;
    } else {
      log_lambda_1_2 = log12_base ;
      log_lambda_2_1 = log21_base;
    }

    {
      matrix[2,2] P = transition_matrix(last_subinterval_sim[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
      next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
      if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
      if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
      if (current_state == 1) at_risk_acquisition[n] += 1;
      if (current_state == 2) at_risk_decolonisation[n] += 1;
      
      current_state = next_state;
      P_total *= P;
    }
  }
}

       // --- Final update (safe after continue)
      //real midpoint_final = (date_use[n - 1] + date_use[n]) / 2;
      //int flag1_final = (intervention[n] == 1 && midpoint_final >= t_star1) ? 1 : 0;
      //int flag2_final = (intervention[n] == 1 && midpoint_final >= t_star2) ? 1 : 0;
      
      vector[2] probs = to_vector(P_total[observed_state[n - 1]]);
      int final_state = categorical_rng(probs / sum(probs));
      if (observed_state[n - 1] == 1 && final_state == 2) acquisitions[n] += 1;
      if (observed_state[n - 1] == 2 && final_state == 1) decolonisations[n] += 1;
      if (observed_state[n - 1] == 1) at_risk_acquisition[n] += 1;
      if (observed_state[n - 1] == 2) at_risk_decolonisation[n] += 1;

      observed_state[n] = final_state;
      y_rep[n] = observed_state[n];
      log_lambda_1_2_out[n] = log_lambda_1_2;
      log_lambda_2_1_out[n] = log_lambda_2_1;
      
      //int idx_log = idx_last_sim[n];
    
    }
    print(\"n: \", n, 
          \", beta_int1_1: \", beta_int1_1, \", beta_int1_2: \", beta_int1_2, 
          \", log_lambda_1_2: \", log_lambda_1_2);
    

    log_lambda_1_2_out[n] = log_lambda_1_2; // log12_base  + flag1_final * beta_int1_1 + flag2_final * beta_int1_2;
    log_lambda_2_1_out[n] = log_lambda_2_1; // log21_base + flag1_final * beta_int2_1 + flag2_final * beta_int2_2;
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




get_expected_rates <- function(age, sex, date, household_id, intervention,
                               first_subinterval, last_subinterval,
                               intervention_date1, intervention_date2,
                               true_params) {
  q_sum <- true_params$q_sum
  alpha <- true_params$alpha
  q_1_2_raw <- q_sum * alpha
  q_2_1_raw <- q_sum * (1 - alpha)
  q_1_2_base <- -3.5 + q_1_2_raw * true_params$sigma_q_1_2
  q_2_1_base <- -4.71 + q_2_1_raw * true_params$sigma_q_2_1
  
  t_obs <- as.numeric(date - as.Date(global_interval_start))
  t_start <- t_obs - first_subinterval + 1
  t_end <- t_obs
  t_star1 <- as.numeric(intervention_date1 - as.Date(global_interval_start))
  t_star2 <- as.numeric(intervention_date2 - as.Date(global_interval_start))
  
  u <- true_params$u[household_id]
  
  log_lambda_1_2 <- q_1_2_base + true_params$beta_1_2_age * age + true_params$beta_1_2_sexe * sex + u
  log_lambda_2_1 <- q_2_1_base + true_params$beta_2_1_age * age + true_params$beta_2_1_sexe * sex + u
  
  if (intervention == 1) {
    # First subinterval
    if (first_subinterval > 0 && t_start < t_star1 && t_star1 < t_end) {
      duration_pre <- t_star1 - t_start
      duration_post <- first_subinterval - duration_pre
      lambda_1_2_pre <- exp(log_lambda_1_2)
      lambda_2_1_pre <- exp(log_lambda_2_1)
      lambda_1_2_post <- exp(log_lambda_1_2 + true_params$beta_int1_1)
      lambda_2_1_post <- exp(log_lambda_2_1 + true_params$beta_int2_1)
      # Effective rate as weighted average of hazard rates
      log_lambda_1_2 <- log((lambda_1_2_pre * duration_pre + lambda_1_2_post * duration_post) / first_subinterval)
      log_lambda_2_1 <- log((lambda_2_1_pre * duration_pre + lambda_2_1_post * duration_post) / first_subinterval)
    } else if (t_obs >= t_star1) {
      log_lambda_1_2 <- log_lambda_1_2 + true_params$beta_int1_1
      log_lambda_2_1 <- log_lambda_2_1 + true_params$beta_int2_1
    }
    # Last subinterval
    if (last_subinterval > 0 && t_start < t_star2 && t_star2 < t_end) {
      duration_pre <- t_star2 - (t_end - last_subinterval + 1)
      duration_post <- last_subinterval - duration_pre
      lambda_1_2_pre <- exp(log_lambda_1_2)
      lambda_2_1_pre <- log_lambda_2_1
      lambda_1_2_post <- exp(log_lambda_1_2 + true_params$beta_int1_2)
      lambda_2_1_post <- exp(log_lambda_2_1 + true_params$beta_int2_2)
      log_lambda_1_2 <- log((lambda_1_2_pre * duration_pre + lambda_1_2_post * duration_post) / last_subinterval)
      log_lambda_2_1 <- log((lambda_2_1_pre * duration_pre + lambda_2_1_post * duration_post) / last_subinterval)
    } else if (t_obs >= t_star2) {
      log_lambda_1_2 <- log_lambda_1_2 + true_params$beta_int1_2
      log_lambda_2_1 <- log_lambda_2_1 + true_params$beta_int2_2
    }
  }
  
  return(list(
    lambda_1_2 = exp(log_lambda_1_2),
    lambda_2_1 = exp(log_lambda_2_1),
    log_lambda_1_2 = log_lambda_1_2,
    log_lambda_2_1 = log_lambda_2_1
  ))
}



# Extract simulated data
sim_data <- rstan::extract(stan_fit)
observed_state_sim <- as.vector(sim_data$y_rep)
log_lambda_1_2_out <- as.vector(sim_data$log_lambda_1_2_out)
log_lambda_2_1_out <- as.vector(sim_data$log_lambda_2_1_out)
acquisitions <- as.vector(sim_data$acquisitions)
at_risk_acquisitions <- as.vector(sim_data$at_risk_acquisition)
decolonisations <- as.vector(sim_data$decolonisations)
at_risk_decolonisations <- as.vector(sim_data$at_risk_decolonisation)

# Create simulated_dataset
simulated_dataset <- data.frame(
  Observation = 1:N,
  MenagememberID = menage_id_member,
  HouseID = HouseID,
  VillageID = VillageID,
  Intervention = Intervention,
  Round = rounds,
  Date = date_use,
  date_use_num = as.numeric(date_use),
  Age = age,
  Sexe = sexe,
  Intervention_Date1 = as.Date(simulated_intervention_date,  origin = "1970-01-01"),
  Intervention_Date2 = as.Date(simulated_intervention_date2, origin = "1970-01-01"),
  Idx_First_Sim = stan_data$idx_first_sim,
  Idx_Last_Sim = stan_data$idx_last_sim,
  #  Idx_Mid_Sim = stan_data$idx_middle_sim,
  First_Subinterval_Sim = stan_data$first_subinterval_sim,
  Last_Subinterval_Sim = stan_data$last_subinterval_sim,
  Num_Middle_Subintervals_Sim = stan_data$num_middle_subintervals_sim,
  Log_Lambda_1_2 = log_lambda_1_2_out,
  Log_Lambda_2_1 = log_lambda_2_1_out,
  Observed_State_Sim = observed_state_sim,
  Acquisitions = acquisitions,
  Decolonisations = decolonisations,
  AtRisk_Acquisitions = at_risk_acquisitions,
  AtRisk_Decolonisations = at_risk_decolonisations,
  stringsAsFactors = FALSE
)

# IDENTIFY ACQUISITIONS AND DECOLONISATIONS
simulated_dataset <- simulated_dataset %>%
  mutate(
    lambda_1_2 = exp(Log_Lambda_1_2),
    lambda_2_1 = exp(Log_Lambda_2_1),
    Intervention_Date1 = ifelse(Intervention==0, NA, as.character(Intervention_Date1)),
    Intervention_Date1 = as.Date(Intervention_Date1, "%Y-%m-%d"),
    Intervention_Date2 = ifelse(Intervention==0, NA, as.character(Intervention_Date2)),
    Intervention_Date2 = as.Date(Intervention_Date2, "%Y-%m-%d"),
    seasonal_t = as.numeric(Date - as.Date(global_interval_start)) / 365,
    month = month(Date)
  ) %>%
  arrange(MenagememberID, Date) %>%
  group_by(MenagememberID) %>%
  mutate(
    Prev_State = lag(Observed_State_Sim),
    Prev_Date = lag(Date),
    Acquisition = Acquisitions,
    Decolonisation = Decolonisations,
    AtRisk_Acquisitions = AtRisk_Acquisitions,
    AtRisk_Decolonisations = AtRisk_Decolonisations,
    #Acquisition = ifelse(Prev_State == 1 & Observed_State_Sim == 2, 1, 0),
    #Decolonisation = ifelse(Prev_State == 2 & Observed_State_Sim == 1, 1, 0),
    Time_Period = case_when(
      Intervention == 1 & Date >= Intervention_Date1 & Date < Intervention_Date2 ~ "Post-Intervention 1",
      Intervention == 1 & Date >= Intervention_Date2 ~ "Post-Intervention 2",
      Intervention == 1 & Date < Intervention_Date1 ~ "Pre-Intervention",
      Intervention == 0 ~ "Control"
    )
  ) %>%
  ungroup()

write.csv(simulated_dataset, "./Output/Model_results/Validation_with_simulated_data/no_season_simdata_raneem.csv")
write.csv(diagnostic_table, "./Output/Model_results/Validation_with_simulated_data/diagnostics_table_raneem.csv")

# Apply expected rates function row-by-row
sim_exp_dataset <- simulated_dataset %>%
  filter(Round != 1) %>%
  mutate(expected = pmap(
    list(age = Age,
         sex = Sexe,
         date = Date,
         household_id = HouseID,
         intervention = Intervention,
         first_subinterval = First_Subinterval_Sim,
         last_subinterval = Last_Subinterval_Sim,
         intervention_date1 = Intervention_Date1,
         intervention_date2 = Intervention_Date2),
    function(age, sex, date, household_id, intervention,
             first_subinterval, last_subinterval,
             intervention_date1, intervention_date2) {
      get_expected_rates(age, sex, date, household_id, intervention,
                         first_subinterval, last_subinterval,
                         intervention_date1, intervention_date2,
                         true_params)
    }
  )) %>%
  mutate(
    expected_lambda_1_2 = map_dbl(expected, "lambda_1_2"),
    expected_lambda_2_1 = map_dbl(expected, "lambda_2_1"),
    diff_lambda_1_2 = expected_lambda_1_2 - lambda_1_2,
    diff_lambda_2_1 = expected_lambda_2_1 - lambda_2_1
  ) %>%
  select(-expected) %>%
  mutate(Time_Period = factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))

round2_post1 <- sim_exp_dataset %>% filter(Round == 2 & Time_Period == "Post-Intervention 1")
head(round2_post1 %>% select(lambda_1_2, expected_lambda_1_2, diff_lambda_1_2))

# Summarise λ₁₂ (acquisition rate) and λ₂₁ (decolonisation rate)
lambda_summary <- simulated_dataset %>%
  mutate(Week = floor_date(Date, unit = "week")) %>%
  group_by(Time_Period, Week) %>%
  summarise(
    Mean_lambda_1_2 = mean(lambda_1_2, na.rm = TRUE),
    Mean_lambda_2_1 = mean(lambda_2_1, na.rm = TRUE),
    N = n()
  ) %>%
  arrange(factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))

lambda_summary <- lambda_summary %>%
  mutate(Time_Period = factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))


lambda_summary_stats <- simulated_dataset %>%
  mutate(Week = floor_date(Date, unit = "week")) %>%
  group_by(Time_Period, Week) %>%
  summarise(
    lambda_1_2 = mean(lambda_1_2, na.rm = TRUE),
    lambda_2_1 = mean(lambda_2_1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Time_Period) %>%
  summarise(
    lambda_1_2_summary = paste0(
      round(median(lambda_1_2, na.rm = TRUE), 3), 
      " (", 
      round(quantile(lambda_1_2, 0.25, na.rm = TRUE), 3), 
      "–", 
      round(quantile(lambda_1_2, 0.75, na.rm = TRUE), 3), 
      ")"
    ),
    lambda_2_1_summary = paste0(
      round(median(lambda_2_1, na.rm = TRUE), 3), 
      " (", 
      round(quantile(lambda_2_1, 0.25, na.rm = TRUE), 3), 
      "–", 
      round(quantile(lambda_2_1, 0.75, na.rm = TRUE), 3), 
      ")"
    )
  ) %>%
  rename(
    `Acquisition Rate` = lambda_1_2_summary,
    `Decolonisation Rate` = lambda_2_1_summary
  )%>%
  mutate(Time_Period = factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))


grid.table(lambda_summary_stats[c(1,4,2,3),])


# Plot λ₁₂ (Acquisition Rate)
p3 = ggplot(lambda_summary, aes(x = Week, y = Mean_lambda_1_2, color = Time_Period)) +
  geom_point(size= 3) +
  labs(title = expression(paste("Weekly Mean ", lambda[1*2], " (Acquisition Rate)")),
       x = "Week",
       y = expression(lambda[1*2])) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_color_brewer(palette = "Set1")

p3.2 = ggplot(lambda_summary, aes(x = Time_Period, y = Mean_lambda_1_2, fill = Time_Period)) +
  geom_boxplot() +
  labs(title = expression(paste("Weekly Mean ", lambda[1*2], " (Acquisition Rate)")),
       x = "Week",
       y = expression(lambda[2*1])) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  #scale_x_date(date_breaks = "4 week", date_labels = "%b %d") +
  scale_fill_brewer(palette = "Set1")
p3.2

# Plot λ₂₁ (Decolonization Rate)
p4 = ggplot(lambda_summary, aes(x = Week, y = Mean_lambda_2_1, color = Time_Period)) +
  geom_point(size= 3) +
  labs(title = expression(paste("Weekly Mean ", lambda[2*1], " (Decolonization Rate)")),
       x = "Week",
       y = expression(lambda[2*1])) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_color_brewer(palette = "Set1")

p4.2 = ggplot(lambda_summary, aes(x = Time_Period, y = Mean_lambda_2_1, fill = Time_Period)) +
  geom_boxplot() +
  labs(title = expression(paste("Weekly Mean ", lambda[2*1], " (Decolonization Rate)")),
       x = "Week",
       y = expression(lambda[2*1])) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  #scale_x_date(date_breaks = "4 week", date_labels = "%b %d") +
  scale_fill_brewer(palette = "Set1")
p4.2

# Acquisition rate - simulated vs expected
p5 = ggplot(sim_exp_dataset, aes(x = lambda_1_2, y = expected_lambda_1_2, col = Time_Period)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~Round) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Simulated λ12 (used in Stan simulation)",
    y = "Expected λ12 (recomputed from true_params)",
    title = "Check: Acquisition Rate λ12 – Expected vs Simulated"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
print(p5)

# Decolonisation rate - simulated vs expected
# Only off for round 2/month february. This likely has to do with how in the expected function, I am not splitting 
# between pre-post intervention when the subinterval covers both
p6 = ggplot(sim_exp_dataset, aes(x = lambda_2_1, y = expected_lambda_2_1, col=Time_Period)) + 
  geom_point(alpha = 0.2) +
  facet_wrap(~Round) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Simulated λ21 (used in Stan simulation)",
    y = "Expected λ21 (recomputed from true_params)",
    title = "Check: Decolonisation Rate λ21 – Expected vs Simulated"
  ) +
  theme_minimal()+
  scale_color_brewer(palette = "Set1")
p6
# + 
#guides(color=F) # Something is off for a particular date


# check_d = (sim_exp_dataset %>%
#              filter(Round == 2, Intervention == 1, round(diff_lambda_2_1,7)!=0) %>%
#              mutate(
#                t_obs = as.numeric(Date - as.Date(global_interval_start)),
#                t_star1 = as.numeric(Intervention_Date1 - as.Date(global_interval_start)),
#                diff_days = t_obs-t_star1,
#                int1_21_applied = t_obs >= t_star1,
#                first_sub_diff_applied = t_star1 < First_Subinterval_Sim, #   if (t0 < t_star1 && t_star1 < t1 && intervention[n] == 1) then do not apply intervention
#                household_u = map_dbl(HouseID, ~ true_params$u[.x]),
#                seasonal_t = as.numeric(Date - as.Date(global_interval_start)) / 365
#              ) %>%
#              select(Round, Date, seasonal_t,VillageID, HouseID,MenagememberID, t_obs,t_star1, int1_21_applied, first_sub_diff_applied,lambda_2_1, expected_lambda_2_1, diff_lambda_2_1, HouseID, household_u, diff_days))
# 
# View(check_d)
# length(unique(check_d$VillageID))
# unique(check_d$seasonal_t); min(check_d$seasonal_t); max(check_d$seasonal_t)
# unique(X)
# 
# Plot
# The number of acquisitions are lower in the intervention group, but do go down in the second two rounds. Could this be as people are still colonised?

acq_summary_weekly <- simulated_dataset %>%
  filter(Round>1) %>%
  mutate(Week = floor_date(Date, unit = "week")) %>%
  group_by(Round,Intervention, Time_Period, Week) %>%
  summarise(
    Acquisitions = sum(Acquisition),
    Decolonisations = sum(Decolonisations),
    AtRisk_Acquisitions = sum(AtRisk_Acquisitions),
    AcquisitionRate = Acquisitions / AtRisk_Acquisitions,
    AtRisk_Decolonisations = sum(AtRisk_Decolonisations),
    DecolonisationRate = Decolonisations / AtRisk_Decolonisations,
    Positive = sum(Observed_State_Sim-1),
    Sampled = n(),  # total individuals sampled in that week
    .groups = "drop"
  )%>%
  mutate(Time_Period = factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))


acq_summary_monthly <- simulated_dataset %>%
  filter(Round>1) %>%
  mutate(Month = floor_date(Date, unit = "month")) %>%
  group_by(Round,Intervention, Time_Period, Month) %>%
  summarise(
    Acquisitions = sum(Acquisition),
    Decolonisations = sum(Decolonisations),
    AtRisk_Acquisitions = sum(AtRisk_Acquisitions),
    AcquisitionRate = Acquisitions / AtRisk_Acquisitions,
    AtRisk_Decolonisations = sum(AtRisk_Decolonisations),
    DecolonisationRate = Decolonisations / AtRisk_Decolonisations,
    Positive = sum(Observed_State_Sim-1),
    Sampled = n(),  # total individuals sampled in that week
    .groups = "drop"
  )%>%
  mutate(Time_Period = factor(Time_Period, levels = c("Control", "Pre-Intervention", "Post-Intervention 1", "Post-Intervention 2")))


p7 = ggplot(acq_summary_weekly, aes(x = Week, y = Sampled, fill = Time_Period)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Intervention, labeller = as_labeller(c(`0` = "Control", `1` = "Intervention"))) +
  labs(
    title = "Weekly Sampled by Intervention Period",
    x = "Week",
    y = "Number of Sampled individuals",
    fill = "Time Period"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")

p8 =ggplot(acq_summary_weekly, aes(x = Week, y = AcquisitionRate, fill = Time_Period)) +
  geom_col() +
  facet_wrap(~ Intervention, scales = "free_x") +
  geom_text(aes(label = AtRisk_Acquisitions), vjust = -0.5, size = 2.5) +
  labs(title = "Weekly Acquisition rate",
       y = "Weekly acquisition rate", x = "Week")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p8

p8.2 =ggplot(acq_summary_weekly, aes(x = factor(Round), y = AcquisitionRate, fill = factor(Round))) +
  geom_boxplot() +
  facet_wrap(~ Intervention, labeller = as_labeller(c(`0` = "Control", `1` = "Intervention"))) +
  #geom_text(aes(label = AtRisk_Acquisitions), vjust = -0.5, size = 2.5) +
  labs(title = "Weekly Acquisition rate",
       y = "Weekly acquisition rate", x = "Round")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  #scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p8.2

p9 = ggplot(acq_summary_weekly, aes(x = Week, y = Positive/Sampled, fill = Time_Period)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Intervention, labeller = as_labeller(c(`0` = "Control", `1` = "Intervention"))) +
  labs(
    title = "Weekly Prevalence by Intervention Period",
    x = "Week",
    y = "Prevalence",
    fill = "Time Period"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p9


p9.2 = ggplot(acq_summary_weekly, aes(x = factor(Round), y = Positive/Sampled, fill = factor(Round))) +
  geom_boxplot() +
  facet_wrap(~ Intervention, labeller = as_labeller(c(`0` = "Control", `1` = "Intervention"))) +
  labs(
    title = "Weekly Prevalence by Intervention Period",
    x = "Round",
    y = "Prevalence",
    fill = "Time Period"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  # scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p9.2

p11 =ggplot(acq_summary_weekly, aes(x = Week, y = DecolonisationRate, fill = Time_Period)) +
  geom_col() +
  facet_wrap(~ Intervention, scales = "free_x") +
  geom_text(aes(label = AtRisk_Decolonisations), vjust = -0.5, size = 2.5) +
  labs(title = "Weekly decolonisation rate",
       y = "Weekly decolonisation rate", x = "Week")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p11


p11.2 =ggplot(acq_summary_weekly, aes(x = factor(Round), y = DecolonisationRate, fill = factor(Round))) +
  geom_boxplot() +
  facet_wrap(~ Intervention, labeller = as_labeller(c(`0` = "Control", `1` = "Intervention"))) +
  #geom_text(aes(label = AtRisk_Acquisitions), vjust = -0.5, size = 2.5) +
  labs(title = "Weekly Decolonisation rate",
       y = "Weekly decolonisation rate", x = "Round")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  #scale_x_date(date_breaks = "4 week", date_labels = "%b %Y") +
  scale_fill_brewer(palette = "Set1")
p11.2


# Loop over a list of plots
grid.table(lambda_summary_stats[c(1,4,2,3),])
plots <- list(p, p2, p3,p3.2, p4,p4.2, p5, p6, p7, p8,p8.2, p9, p9.2, p11, p11.2)
for (plt in plots) {
  print(plt)
}

dev.off()




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
#####################

# Ensure date columns are in Date format
if (!inherits(simulated_dataset$Date, "Date")) {
  simulated_dataset$Date <- as.Date(simulated_dataset$Date)
}

if (!inherits(simulated_dataset$Intervention_Date1, "Date")) {
  simulated_dataset$Intervention_Date1 <- as.Date(simulated_dataset$Intervention_Date1, origin = "1970-01-01")
}

if (!inherits(simulated_dataset$Intervention_Date2, "Date")) {
  simulated_dataset$Intervention_Date2 <- as.Date(simulated_dataset$Intervention_Date2, origin = "1970-01-01")
}

# Create intervention activity flags
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

# Check that Intervention 2 never activates before Intervention 1
if (any(simulated_dataset$Intervention2_active == 1 & simulated_dataset$Intervention1_active == 0)) {
  warning("Some rows have Intervention 2 active before Intervention 1!")
}
