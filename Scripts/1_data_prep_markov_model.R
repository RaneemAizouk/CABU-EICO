####################################################################################
# THIS CODE CONVERTS THE OBSERVED DATA IN THE CABU EICO STUDY TO A STAN DATA FORMAT
####################################################################################

# Last updated: 15 December 2025

# Clear environment
rm(list = ls())

# Libraries 
pacman::p_load(ggplot2, rstan, tidyr, bayesplot, gridExtra, dplyr, GGally, lubridate, readr)

#----------------------------------------------------------------------
# load dataset
#----------------------------------------------------------------------

#load("./Data/BF/clean/use_in_analyses/bf_esbl0123_long_all.rda")
#timing_interventions <- read.csv("./Data/BF/timing_interventions.csv")

#----------------------------------------------------------------------
# Data pre-processing 
#----------------------------------------------------------------------

# Lines 37-230 have been used to anonimise the data for publication therefore are canceled out
#------------------------------------------------------------------------------------------------


# Create an explicit "observed_state" variable from the original carriage indicator
# This makes it easier to recode or rename later without touching the raw "esble" column.
# data$observed_state <- data$esble
# data <- data %>%
#   relocate(observed_state, .after = esble)
# 
# # Assign each household a unique integer ID
# # factor() converts the string menage_id into a factor with levels in the order encountered,
# # and as.numeric() then maps those levels 1,2,3,… to the HouseID column.
# data$HouseID <- as.numeric(factor(data$menage_id))
# data <- data %>%
#   relocate(HouseID, .after = menage_id)
# 
# # Extract each person’s within-household number from menage_id_member
# # The pattern ".*-" strips off everything up to the dash, leaving the numeric suffix.
# data$individual_number <- as.numeric(sub(".*-", "", data$menage_id_member))
# 
# # Combine HouseID and individual_number into a single integer identifier:
# # multiply HouseID by 100 (reserving two digits) then add the individual number.
# # For example, HouseID=1 & individual_number=5 → 100*1 + 5 = 105
# multiplier <- 100
# data$menage_id_member_orig <- data$menage_id_member
# data$menage_id_member <- data$HouseID * multiplier + data$individual_number
# 
# # Parse the survey round from the RedCap event name
# # Extract the number after "round_" and before "_arm_1", then add 1 so that "round_1" → 2, etc.
# data$round <- as.integer(
#   gsub("round_(\\d+)_arm_1", "\\1", data$redcap_event_name)
# ) + 1
# 
# data <- data %>%
#   relocate(round, .after = redcap_event_name)
# 
# # store dataset for publication
# dpub <- data[, c(1:38)]
# 
# names(dpub)
# 
# dpub2 <- dpub %>% select(!c(redcap_event_name,time, menage_id, date.consent,date.stool.collection, date, date_conserv, germe_c, data.row,date.enquete,
#                             cleaning.water.storage.binary))
# 
# 
# #----------------------------------------------------------------------
# # Add intervention date
# #----------------------------------------------------------------------
# 
# # Clean and prepare timing_interventions
# timing_interventions <- timing_interventions %>%
#   mutate(
#     village_name = tolower(trimws(Village_name)),
#     Intervention_round = as.integer(Round),
#     Intervention_start_date = dmy(Intervention_start_date)  # handles dd/mm/yyyy safely
#   ) %>%
#   #mutate(data_round = Round + 1) %>%  # shift rounds by +1 to match data_complete
#   select(village_name, Intervention_round, Intervention_start_date)
# 
# # make wide
# timing_interventions_wide <- timing_interventions %>%
#   pivot_wider(
#     id_cols = village_name,
#     names_from = Intervention_round,
#     values_from = Intervention_start_date,
#     names_prefix = "intervention_round"
#   )
# 
# # Clean and prepare data_complete
# dpub2 <- dpub2 %>%
#   mutate(
#     village_name = tolower(trimws(village_name)),
#     round = as.integer(round)
#   )
# unique(dpub2$village_name)
# unique(timing_interventions$village_name)
# 
# 
# # Join intervention dates to data_complete by village and adjusted round
# dpub2 <- dpub2 %>%
#   left_join(
#     timing_interventions_wide,
#     by = c("village_name" = "village_name")
#   )
# 
# table(dpub2$intervention_round1[dpub2$intervention.text==1], useNA="always") # no NAs
# 
# # Fill Intervention_start_date with first intervention round start
# dpub2 <- dpub2 %>%
#   mutate(
#     Intervention_start_date = intervention_round1,
#     intervention_date = intervention_round1, # in the baseline scenario, we assume intervention rounds 1+2 have the same effect and intervention round 3 has an additional effect
#     intervention_date2 = intervention_round3 #
#   )
# 
# 
# # Relocate the column after 'round'
# dpub2 <- dpub2 %>%
#   relocate(Intervention_start_date, intervention_date, intervention_date2, .after = round)
# 
# # Check the result
# dpub2 %>%
#   select(menage_id_member, round, Intervention_start_date, intervention_date, intervention_date2) %>%
#   distinct() %>%
#   head(10)
# 
# sum(is.na(dpub2$intervention_date2))
# sum(is.na(dpub2$Intervention_start_date))
# class(dpub2$intervention_date)
# class(dpub2$intervention_date2)
# 
# unique(dpub2$intervention_date)
# unique(dpub2$intervention_date2)
# 
# # fill Na values in the intervention date  for control group with date after the study
# #------------------------------------------------------------------------------
# # 1) Compute the last sampling date across everyone:
# # 1) Compute the last sampling date across all observations:
# max_obs_date <- max(dpub2$date.use, na.rm = TRUE)
# 
# interval_length <- 28
# 
# # 2) Fill NA in intervention_date and intervention_date2 (both are Date already):
# dpub2 <- dpub2 %>%
#   mutate(
#     intervention_date  = if_else(
#       is.na(intervention_date),
#       max_obs_date + interval_length,
#       intervention_date
#     ),
#     intervention_date2 = if_else(
#       is.na(intervention_date2),
#       max_obs_date + interval_length,
#       intervention_date2
#     ),
#     Intervention_start_date = intervention_date
#   ) %>%
#   # 3) Sanity‐check: no NAs should remain
#   { stopifnot(!any(is.na(.$intervention_date))); . } %>%
#   { stopifnot(!any(is.na(.$intervention_date2))); . } #%>%
#   # 4) Convert back to numeric days‐since‐1970 for Stan
#   #mutate(
#   #  intervention_date  = as.numeric(intervention_date),
#   #  intervention_date2 = as.numeric(intervention_date2)
#   #)
# 
# sum(is.na(dpub2$intervention_date))    # should print 0
# sum(is.na(dpub2$intervention_date2))   # should print 0
# sum(is.na(dpub2$Intervention_start_date))  # should print 0
# 
# 
# # Get unique villages
# villages <- unique(dpub2$village_name)
# 
# # Sort alphabetically (ignore case)
# villages_sorted <- sort(villages)
# 
# # Create numeric codes 1:22
# village_codes <- setNames(seq_along(villages_sorted), villages_sorted)
# 
# # Recode into a new variable
# dpub2$village_code <- village_codes[dpub2$village_name]
# 
# 
# dpub3 <- dpub2 %>% select(-c(village, village_name))
# #data <- data[,c(1:31,113:116)]
# #names(data)
# 
# data = dpub3
# 
# #-----------------------------------------------------------------------
# # Convert variables to binary
# #-----------------------------------------------------------------------
# 
# # Create a binary variable for intervention villages with date condition
# data <- data %>%
#   mutate(
#     intervention_village = case_when(
#       intervention.text == "intervention"  ~ 1,
#       TRUE ~ 0
#     )
#   )
# 
# # Convert sexe and esble =observed_state to binary 0,1 for male and female respectively.
# # First, ensure 'sexe' is a character type
# data$sexe <- as.character(data$sexe)
# 
# # Now, correctly convert 'sexe' to binary (1 = Female, 0 = Male)
# data$sexe <- ifelse(data$sexe == "Female", 1, 0)
# 
# # Confirm conversion worked
# unique(data$sexe)
# 
# dfls0 = data
# 
# # Store data for paper
# saveRDS(dfls0, file="./Data/Manuscript/bf_esbl0123_long_all.rds")
# write.csv(dfls0, "./Data/Manuscript/bf_esbl0123_long_all.csv")


# LOAD IN DATA
#----------------------------------------------------------------------

data = readRDS("./Public_data/Observed/bf_esbl0123_long_all.rds")

#----------------------------------------------------------------------
# Handle missing data
#----------------------------------------------------------------------

# Removes individuals with all missing age or sexe values, then fills remaining missing values within individuals using fill(.direction = "downup").
# Identify individuals with all‐missing age or sex
missing_data_individuals <- data %>%
  group_by(menage_id_member) %>%  # for each person…
  summarise(
    missing_age = all(is.na(age)),    #   is age always missing?
    missing_sexe = all(is.na(sexe))  #   is sex always missing?
  ) %>%
  filter(missing_age | missing_sexe)  # keep only those with a fully missing field

# Drop those completely uninformative individuals(we dont have many of them around 4)
data_filtered <- data %>%
  filter(!menage_id_member %in% missing_data_individuals$menage_id_member)

# Carry forward/backward any remaining sporadic NAs
# (we assume age/sex are constant within a person)
data_filled <- data_filtered %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()

# Check that no NAs remain for age or sex
missing_after_fill <- data_filled %>%
  filter(is.na(age) | is.na(sexe))

if (nrow(missing_after_fill) > 0) {
  print("Rows with remaining missing values:")
  print(missing_after_fill)
} else {
  print("All missing values for age and sexe have been filled.")
}

# Report how many individuals survived this filtering
num_individuals_after_filling <- n_distinct(data_filled$menage_id_member)
print(paste("Number of individuals after filling missing age and sexe:", num_individuals_after_filling)) #1220

unique(data_filled$sexe)  # Should print only 0 and 1


#----------------------------------------------------------------------
# Filter on individuals with complete data or not
#----------------------------------------------------------------------

# If we want to use the full data set just filter out this section.
# Filter individuals who have participated in all rounds

# Count how many rounds each person actually appears in
round_counts <- data_filled %>%
  group_by(menage_id_member) %>%
  summarise(rounds_participated = n_distinct(round))

# Find the true maximum of that count
max_participation <- max(round_counts$rounds_participated)

# Keep only individuals who hit that maximum
individuals_in_all_rounds <- round_counts %>%
  filter(rounds_participated == max_participation) %>%
  pull(menage_id_member)

# Keep individuals with >= 2 observations
data_complete <- data_filled %>%
  group_by(menage_id_member) %>%
  filter(n() >= 2) %>%
  ungroup()

# Remove NA in observed_state
data_complete <- data_complete %>%
  filter(!is.na(esble) & esble %in% c(0, 1)) %>%
  mutate(observed_state = ifelse(esble == 0, 1, 2))


# Check that each remaining person appears exactly max_participation times
table(data_complete$menage_id_member) %>% 
  table()   # counts how many individuals have 1 obs, 2 obs, etc; 102 = 2; 306 = 3; 743 = 4

# Impute the first rounds for this that have no first round (17 in total)
#--------------------------------------------------------------------------------

# Impute missing Round 1 using population prevalence

# Who needs Round 1?
ids_missing_r1 <- data_complete %>%
  group_by(menage_id_member) %>%
  summarise(has_r1 = any(round == 1L), .groups = "drop") %>%
  filter(!has_r1) %>%
  pull(menage_id_member)

length(ids_missing_r1)  # should be 17

# Prevalence to use for imputation
#    Prefer baseline (Round 1) prevalence among those who actually have Round 1.
p_baseline <- data_complete %>%
  filter(round == 1L) %>%
  summarise(p = mean(observed_state == 2, na.rm = TRUE)) %>%
  pull(p)

# fallback to overall prevalence if no Round 1 available
# if (is.na(p_baseline) || length(p_baseline) == 0) {
#   p_baseline <- mean(data_complete$observed_state == 2, na.rm = TRUE)
# }

# Stratify prevalence by intervention_village at Round 1
p_by_stratum <- data_complete %>%
  filter(round == 1L) %>%
  group_by(intervention_village) %>%
  summarise(p = mean(observed_state == 2, na.rm = TRUE), .groups = "drop")

# Build rows to add (clone each person’s earliest row; set round=1; draw observed_state ~ Bernoulli(p))
set.seed(123)  # reproducible imputation

rows_to_add <- data_complete %>%
  filter(menage_id_member %in% ids_missing_r1) %>%
  arrange(menage_id_member, round) %>%
  group_by(menage_id_member) %>%
  slice(1) %>%
  ungroup() %>%
  # attach stratum prevalence and fallback to baseline if missing
  left_join(p_by_stratum, by = c("intervention_village")) %>%
  mutate(p = ifelse(is.na(p), p_baseline, p)) %>%
  mutate(
    round = 1L,
    # Draw colonisation status using prevalence p:
    # observed_state: 1 = negative, 2 = positive
    observed_state = ifelse(rbinom(n(), size = 1, prob = p) == 1, 2L, 1L),
    esble = ifelse(observed_state == 2L, 1L, 0L),
    imputed_round1 = TRUE
    # If you have a round->date map, also set date_use here, e.g.:
    # date_use = round_date_map[1]
  )

# Keep column order tidy (add flag to existing data, then bind)
if (!"imputed_round1" %in% names(data_complete)) {
  data_complete$imputed_round1 <- FALSE
}

# Append imputed Round 1 rows and de-duplicate
data_complete <- bind_rows(data_complete, rows_to_add) %>%
  arrange(menage_id_member, round) %>%
  distinct(menage_id_member, round, .keep_all = TRUE)

# Quick checks
table(data_complete$round)                      # now includes Round 1 everywhere
sum(data_complete$imputed_round1, na.rm = TRUE) # number of imputed R1 rows


#-----------------------------
# Create transition data, this constructs a transition dataset, a row for each pair of consecutive observations within an individual.
data_complete  <- data_complete %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    state_prev = lag(observed_state),
    state_next = observed_state
  ) 
  


##########################################################
# Global Interval Setup
# Defines a global time interval (Oct 2022 to Feb 2024) divided into 28-day intervals.
# Computes midpoints for spline and step function evaluation, adjusting for short final intervals.
#  Scales time (X) to years, centered at the median.
#  Sets up 5 knots for cubic B-splines (spline_degree = 3).

 
# Define global interval structure
global_interval_start <- min(data_complete$date.use) # start date of the study 03-10-2022
global_interval_end <- max(data_complete$date.use) # end date of the study 19-02-2024

global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric <- as.numeric(global_interval_end)

interval_length <- 28

# Compute number of full intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length) # 18
# Pass only the maximum possible number of full middle‐intervals to Stan
max_middle <- num_intervals - 1 # 17

# Compute global interval starts,seq(from, to, by)
interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)

# Compute global interval end, we dont need to use seq() as it will create a seqence of values for each start point while this way is vectorised and will create a vector for all the points the start intervention points sequence.
interval_ends <- interval_starts + (interval_length - 1)

# Compute seasonal time points (X) — use midpoint of each full global interval,
# but if the last interval is short, use the start of that interval instead
X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2

# For short last interval (e.g., <28 days), use start instead of midpoint
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
spline_degree <- 3  # Cubic splines

# Number of basis functions
num_basis <- num_knots + spline_degree - 1  # Correctly computed basis functions

# Print to verify
print(knots)
print(num_basis) # 7

#-------------------------------------------------------------------
# Conversion to stan units 
#-------------------------------------------------------------------

# Ensure HouseID is correctly converted to an integer sequence
data_complete$HouseID <- as.integer(factor(data_complete$HouseID))
# # Scale age before creating the stan_data list
# data_complete$age_scaled <- as.integer(data_complete$age - mean(data_complete$age)) / sd(data_complete$age)

data_complete$month <- as.integer(data_complete$month)

# Ensure all necessary conversions are done

data_complete$date.use <- as.Date(data_complete$date.use)  # Ensure it's in Date format
date_use <- as.numeric(data_complete$date.use)  # Convert to numeric

# Convert categorical variables
data_complete$round <- as.integer(data_complete$round)
data_complete$intervention_village <- as.integer(data_complete$intervention_village)

summary(data_complete$date.use)


#--------------------------------------------------------------------------
# Interval setup
#--------------------------------------------------------------------------

split_intervals <- function(start_day, end_day, interval_length) {
  total_days <- end_day - start_day
  if (total_days <= 0) return(c(0, 0, 0))
  
  first_end <- ceiling(start_day / interval_length) * interval_length
  first <- min(first_end - start_day, total_days)
  
  remaining <- total_days - first
  middle <- floor(remaining / interval_length)
  last <- remaining %% interval_length
  
  c(first, middle, last)
}


N <- nrow(data_complete)
id <- data_complete$menage_id_member

if (exists("interval_length_by_id")) {
  ilen_row <- unname(interval_length_by_id[as.character(id)])
} else if ("interval_length" %in% names(data_complete)) {
  ilen_row <- data_complete$interval_length
} else {
  ilen_row <- rep(interval_length, N)
}
stopifnot(length(ilen_row) == N)

# Convert dates to numeric (as you had)
date_use_num <- as.numeric(data_complete$date.use)
global_interval_start_num <- as.numeric(global_interval_start)
global_interval_end_num   <- as.numeric(global_interval_end)
total_days <- global_interval_end_num - global_interval_start_num + 1

# per-row number of intervals and max_middle for matrix width
num_intervals_row <- ceiling(total_days / ilen_row)
max_middle <- max(1L, max(pmax(0L, num_intervals_row - 2L)))

# initialise
first_subinterval_sim        <- rep(NA_integer_, N)
num_middle_subintervals_sim  <- rep(0L, N)
last_subinterval_sim         <- rep(NA_integer_, N)

idx_first_sim  <- rep(1L, N)
idx_last_sim   <- rep(1L, N)
idx_middle_sim <- matrix(1L, nrow = N, ncol = max_middle)
global_interval_index_start <- rep(1L, N)

## first/middle/last per transition (by id; any # of obs) 
idx_by_id <- split(seq_along(id), id)

for (idx in idx_by_id) {
  # ensure chronological within person
  idx <- idx[order(date_use_num[idx])]
  if (length(idx) < 2) next
  
  for (k in 2:length(idx)) {
    i_prev <- idx[k - 1]
    i_curr <- idx[k]
    
    d_prev <- date_use_num[i_prev]
    d_curr <- date_use_num[i_curr]
    ilen   <- ilen_row[i_curr]  # per-person (row) interval size
    
    out <- split_intervals(d_prev, d_curr, ilen)
    
    first_subinterval_sim[i_curr]       <- out[1]
    num_middle_subintervals_sim[i_curr] <- out[2]
    last_subinterval_sim[i_curr]        <- out[3]
  }
}

# indices of middle subintervals (same logic, per-row ilen)
for (n in 2:N) {
  if (id[n] == id[n - 1] && num_middle_subintervals_sim[n] > 0) {
    ilen_n <- ilen_row[n]
    num_int_n <- num_intervals_row[n]
    
    t_prev <- date_use_num[n - 1]
    raw_start_idx <- floor((t_prev - global_interval_start_num) / ilen_n)
    
    M <- min(num_middle_subintervals_sim[n], max_middle)
    for (m in 1:M) {
      interval_start_m <- global_interval_start_num + (raw_start_idx + m) * ilen_n
      interval_end_m   <- min(interval_start_m + ilen_n - 1, global_interval_end_num)
      interval_mid_m   <- (interval_start_m + interval_end_m) / 2
      raw_mid_idx <- floor((interval_mid_m - global_interval_start_num) / ilen_n) + 1L
      idx_middle_sim[n, m] <- max(1L, min(num_int_n, raw_mid_idx))
    }
    # (optional) leave remaining columns as-in; or set to 1L/NA if you prefer
  }
}

# idx_first_sim and idx_last_sim (same structure, per-row ilen) 
for (n in 1:N) {
  ilen_n <- ilen_row[n]
  num_int_n <- num_intervals_row[n]
  
  if (n == 1 || id[n] != id[n - 1]) {
    idx_first_sim[n] <- floor((date_use_num[n] - global_interval_start_num) / ilen_n) + 1L
    idx_last_sim[n]  <- idx_first_sim[n]
  } else {
    idx_first_sim[n] <- floor((date_use_num[n - 1] - global_interval_start_num) / ilen_n) + 1L
    idx_last_sim[n]  <- floor((date_use_num[n]   - global_interval_start_num) / ilen_n) + 1L
  }
  idx_first_sim[n] <- max(1L, min(num_int_n, idx_first_sim[n]))
  idx_last_sim[n]  <- max(1L, min(num_int_n, idx_last_sim[n]))
}

# global_interval_index_start (same structure, per-row ilen) 
ilen_1 <- ilen_row[1]; num_int_1 <- num_intervals_row[1]
global_interval_index_start[1] <- floor((date_use_num[1] - global_interval_start_num) / ilen_1) + 1L
global_interval_index_start[1] <- max(1L, min(num_int_1, global_interval_index_start[1]))

for (n in 2:N) {
  ilen_n <- ilen_row[n]; num_int_n <- num_intervals_row[n]
  if (id[n] != id[n - 1]) {
    global_interval_index_start[n] <- floor((date_use_num[n] - global_interval_start_num) / ilen_n) + 1L
  } else {
    t_prev <- date_use_num[n - 1]
    global_interval_index_start[n] <- floor((t_prev - global_interval_start_num) / ilen_n) + 1L
  }
  global_interval_index_start[n] <- max(1L, min(num_int_n, global_interval_index_start[n]))
}

#--------------------------------------------------------------------
# Check
#--------------------------------------------------------------------
# build id mapping (preserve order of first appearance)
id <- data_complete$menage_id_member
idx_by_id <- split(seq_len(N), factor(id, levels = unique(id)))

# rounds = actual rows per person, date-sorted
get_round_indices <- function(i) {
  idx <- idx_by_id[[i]]
  idx[order(date_use_num[idx])]
}

# loop over people (not fixed 4)
n_check <- length(idx_by_id)

diagnostic_table <- data.frame(stringsAsFactors = FALSE)

for (i in seq_len(n_check)) {
  idx <- get_round_indices(i)
  L <- length(idx)
  if (L < 2) next
  
  for (r in 2:L) {
    n      <- idx[r]
    n_prev <- idx[r - 1]
    
    date_n    <- as.Date(date_use_num[n], origin = "1970-01-01")
    date_prev <- as.Date(date_use_num[n_prev], origin = "1970-01-01")
    days_between <- as.numeric(date_n - date_prev)
    
    ## interval breakdown (use per-row interval length)
    ilen   <- ilen_row[n]
    first  <- first_subinterval_sim[n]
    middle <- num_middle_subintervals_sim[n] * ilen
    last   <- last_subinterval_sim[n]
    total_from_subs <- first + middle + last
    
    ## middle interval indices (guard width)
    M <- min(num_middle_subintervals_sim[n], ncol(idx_middle_sim))
    mids <- if (M > 0) {
      paste(idx_middle_sim[n, seq_len(M)], collapse = ",")
    } else {
      NA_character_
    }
    
    row <- data.frame(
      Individual      = i,
      Round           = r,
      Date_prev       = date_prev,
      Date_n          = date_n,
      Days_Between    = days_between,
      First_Sub       = first,
      Middle_Sub      = middle,
      Last_Sub        = last,
      Total_From_Subs = total_from_subs,
      Match           = (days_between == total_from_subs),
      Idx_First       = idx_first_sim[n],
      Idx_Last        = idx_last_sim[n],
      Num_Middle      = num_middle_subintervals_sim[n],
      Middle_Idx      = mids,
      stringsAsFactors = FALSE
    )
    
    diagnostic_table <- rbind(diagnostic_table, row)
  }
}

#write.csv(diagnostic_table, "./Data/BF/clean/use_in_analyses/bf_check_stan_data.csv")

print(diagnostic_table[c(1:20),])


#------------------------------------------------------
# Create stan_data frame
#------------------------------------------------------

stan_data <- list(
  N = nrow(data_complete),
  menage_id_member = data_complete$menage_id_member,
  HouseID = data_complete$HouseID,
  #VillageID = data_complete$village_name,
  VillageID = data_complete$village_code,
  H = length(unique(data_complete$HouseID)),
  age = as.integer(data_complete$age),
  round = as.integer(data_complete$round),
  sexe = as.integer(data_complete$sexe),
  observed_state = data_complete$state_next,
  date_use = as.numeric(as.Date(data_complete$date.use)),
  
  intervention = as.integer(data_complete$intervention_village),
  Intervention_start_date = as.numeric(data_complete$Intervention_start_date),
  intervention_date  = as.numeric(data_complete$intervention_date),
  intervention_date2 = as.numeric(data_complete$intervention_date2),
  
  # Interval length and date range
  global_interval_start = global_interval_start_numeric,
  global_interval_end = global_interval_end_numeric,
  interval_length = interval_length,
  num_data = length(X), 
  X = X,
  num_intervals = num_intervals, 
  max_middle = max_middle,
  idx_first_sim = idx_first_sim,
  idx_last_sim = idx_last_sim,
  idx_middle_sim = idx_middle_sim,
  first_subinterval_sim = first_subinterval_sim,
  last_subinterval_sim = first_subinterval_sim,
  num_middle_subintervals_sim = num_middle_subintervals_sim,
  
  # Spline settings (now based on global intervals)
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  num_basis = num_basis
)

# Convert Intervention_start_date to numeric (days since 1970-01-01)

# Replace NA in intervention dates with -999 for control group
stan_data$intervention_date[data_complete$intervention_village == 0] <- -999
stan_data$intervention_date2[data_complete$intervention_village == 0] <- -999
stan_data$intervention <- data_complete$intervention_village

str(stan_data)

# Save data
#-------------------------------------------------------------

# saveRDS(stan_data, file="./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")
# write.csv(data_complete, file = "./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")

saveRDS(stan_data, file="./Public_data/Observed/bf_stan_data_all.rds")
#write.csv(data_complete, file = "./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")
