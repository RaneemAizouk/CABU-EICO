# Clear environment
rm(list = ls())
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

# # old 4) Subset your main data, this part will filter for 4 observations per individual 
# data_complete <- data_filled %>%
#   filter(menage_id_member %in% individuals_in_all_rounds)
# 4) Keep individuals with >= 2 observations
data_complete <- data_filled %>%
  group_by(menage_id_member) %>%
  filter(n() >= 2) %>%
  ungroup()
#Remove NA in observed_state
data_complete <- data_complete %>%
  filter(!is.na(esble) & esble %in% c(0, 1)) %>%
  mutate(observed_state = ifelse(esble == 0, 1, 2))



# 5) Check that each remaining person appears exactly max_participation times
table(data_complete$menage_id_member) %>% 
  table()   # counts how many individuals have 1 obs, 2 obs, etc



##### 6)Create transition data,This constructs a transition dataset, a row for each pair of consecutive observations within an individual.
data_complete  <- data_complete %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    state_prev = lag(observed_state),
    state_next = observed_state
  ) %>%
  filter(!is.na(state_prev) & !is.na(state_next)) %>%
  ungroup()



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

# # Recode esble to 1 or 2 only (assuming it's 0/1)
# data_complete$observed_state <- ifelse(data_complete$esble == 0, 1,
#                                        ifelse(data_complete$esble == 1, 2, NA))


##########################################################


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
# # Scale age before creating the stan_data list
# data_complete$age_scaled <- as.integer(data_complete$age - mean(data_complete$age)) / sd(data_complete$age)

data_complete$month <- as.integer(data_complete$month)

# Ensure all necessary conversions are done

data_complete$date.use <- as.Date(data_complete$date.use)  # Ensure it's in Date format
date_use <- as.numeric(data_complete$date.use)  # Convert to numeric

# Convert categorical variables
data_complete$round <- as.integer(data_complete$round)
data_complete$intervention_village <- as.integer(data_complete$intervention_village)

summary(data_complete$ date.use)






#########################################
stan_data <- list(
  N = nrow(data_complete),
  observed_state = data_complete $state_next,
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
  X = X,# Time variable based on global intervals
  num_data = length(X), # Number of global intervals
  max_middle = max_middle,
  Intervention_start_date = as.numeric(data_complete$Intervention_start_date),
  intervention_date  = as.numeric(data_complete$intervention_date),
  intervention_date2 = data_complete$intervention_date2
)

# Convert Intervention_start_date to numeric (days since 1970-01-01)

# Replace NA in intervention dates with -999 for control group
stan_data$intervention_date[data_complete$intervention_village == 0] <- -999
stan_data$intervention_date2[data_complete$intervention_village == 0] <- -999
stan_data$intervention <- data_complete$intervention_village

str(stan_data)