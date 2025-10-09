####################################################################################################################
# THIS CODE IMPLEMENTS A COX-PROPORTIONAL HAZARD MODEL TO ASSESS THE ASSOCATION BETWEEN WASH AND ESBL-E ACQUISITION
####################################################################################################################
# Author: Raneem Aizouk
# Date: April 2025
# Last updated: 6 October 2025

rm(list = ls())

pacman::p_load(survival, dplyr, tidyr, coxme, splines, mgcv, ggplot2, scales, cowplot, forcats, lubridate, purrr, splines)

# Load dataset
#setwd("/Users/raizouk/Desktop")
#getwd()

#load("/Users/raizouk/Desktop/bf_esbl0123_long_all.rda")

# SET DIRECTORY
DirectoryData <- "./Data/BF/clean"

load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_all.rda"))
load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_completecases.rda"))

data <- dfls0

##############
# Convert round from redcap_event_name
data$round <- as.integer(gsub("round_(\\d+)_arm_1", "\\1", data$redcap_event_name)) + 1
data <- data %>% relocate(round, .after = redcap_event_name)

#-----------------------------------------------------------
# Check number of individuals
#-----------------------------------------------------------
num_unique_members <- data %>%
  filter(intervention.text == "control") %>%
  summarise(unique_members = n_distinct(menage_id_member)) 

print(num_unique_members) # 617

#-----------------------------------------------------------
#Find control individuals with 4 rounds of data
#-----------------------------------------------------------

four_sample_members <- data %>%
  filter(intervention.text == "control") %>%
  group_by(menage_id_member) %>%
  summarise(n_rounds = n_distinct(round)) %>%
  filter(n_rounds == 4)

# Number of individuals
n_four_sample_members <- nrow(four_sample_members)
print(n_four_sample_members) # 390

#-----------------------------------------------------------
# Prevalence at baseline
#-----------------------------------------------------------

esbl_baseline <- data %>%
  filter(intervention.text == "control") %>%
  filter(menage_id_member %in% four_sample_members$menage_id_member) %>%
  summarise(pos_rate = mean(esble == 1) * 100)
print(esbl_baseline)

# Handle "Missing" data
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

# Filter individuals who have participated in all rounds
round_counts <- data_filled %>%
  group_by(menage_id_member) %>%
  summarise(rounds_participated = n_distinct(round))
print(round_counts)

#odd results 3 rounds not 4 
max_rounds <- max(data_filled$round)
individuals_in_all_rounds <- round_counts %>%
  filter(rounds_participated == max_rounds) %>%
  pull(menage_id_member)
print (max_rounds)

data_complete <- data_filled 


# Create HouseID first
data_complete$HouseID <- as.numeric(factor(data_complete$menage_id))

# Define sex column where Female = 1 and Male = 0
data_complete$sexe <- as.integer(ifelse(data_complete$sexe == "Female", 1, 0))

data_complete$month <- as.integer(data_complete$month)

# 1) Set calendar time reference point
# -------------------------------------------------------------
# We define a fixed calendar date `t` as the reference start time for the study: 
t <- as.Date("2022-10-03")

# This date represents "time zero" (t0 = the study baseline) for our survival analysis.
# However, we do NOT actually observe ESBL carriage status at this exact date.

# The earliest observed data on colonization status is from round 1.
# So we *assume that each individual's ESBL status at round 1 reflects their status at time t*.
# This assumption allows us to compute acquisition (ACQ) events in subsequent rounds (e.g., from t1 to t2).

# We do not use time differences between rounds (e.g., round 2 - round 1) because:
#   The follow-up intervals are not consistent across participants.
#   The hazard of acquisition is likely **not constant** over calendar time.
# Therefore, we rely on absolute calendar time (using actual visit dates) instead of relative time between rounds.
# Convert all date fields to proper Date objects.
data_complete$date.use <- as.Date(data_complete$date.use)

# data_complete <- data_complete %>%
#   mutate(date.use = case_when(
#     round == 1 ~ t,  # I don't understand this. Why not just set the first observation at 0. As now the interval between t1 and t1 is not correctly computed
#     TRUE ~ as.Date(date.use)
#   ))

#By setting date.use = t for round 1, the time difference (dt_n) for round 1 becomes 0, and events (e.g. acquisition) are measured only from round 1 onward — avoiding incorrect inference from unobserved baseline.

# Calculate dt_n (time to event variable)
# dt_n is only relevant from round 2 onward (time since baseline `t`)

global_origin <- t

data_complete <- data_complete %>%
  arrange(menage_id_member, date.use) %>%                      # make sure rows are ordered
  group_by(menage_id_member) %>%
  mutate(
    first_date   = first(date.use),              # per-person origin
    # Per-individual counting-process time scale
    tstop        = as.numeric(difftime(date.use, first_date, units = "days")),
    tstart       = lag(tstop, default = 0),
    # Calendar-time scale (absolute days since a fixed origin)
    cal_stop     = as.numeric(difftime(date.use, global_origin, units = "days")),
    cal_start    = lag(cal_stop, default = cal_stop[1]),
    # Interval length check
    dt_days      = tstop - tstart
  ) %>%
  ungroup()

data_complete <- data_complete %>%
  mutate(
    dt_n = case_when(
      round %in% c(2, 3, 4) ~ as.numeric(difftime(date.use, t, units = "days")),
      TRUE ~ NA_real_
    )
  )

data_complete<- data_complete  %>% 
  relocate(dt_n , .after= esble)

data_complete<- data_complete  %>% 
  relocate(dt_days, .after= esble)

# Filter to control group
data_complete <- data_complete %>% filter(intervention.text == "control") 

# Instead of using a start date that is the same for all, use days (add on Esther 6 October 2025)
#data_complete$dt_n = data_complete$days

##################

# Generate IDs
data_complete$HouseID <- as.numeric(factor(data_complete$menage_id))
data_complete$individual_number <- as.numeric(sub(".*-", "", data_complete$menage_id_member))
data_complete$menage_id_member <- data_complete$HouseID * 100 + data_complete$individual_number

# 3) Define outcome (ACQ)

# We track ESBL acquisition (ACQ) and clearance events across calendar time. --> This is then no longer calendar time (evk 6/10/2025)
# Each individual can contribute up to 3 observed transitions: 
#      From round 1 -> round 2
#      From round 2 -> round 3
#      From round 3 -> round 4
#
# These transitions are aligned to real calendar time (dt_n) starting from 
# reference date t = "2022-10-03"` (assigned to round 1 for all).
#
# The ACQ and clearance outcomes are defined between two visits:
#      ACQ = 1 if a person changed from uncolonized (0) -> colonized (1)
#     clearance = 1 if a person changed from colonized (1) -> uncolonized (0)
#
# Since we have 4 rounds, this setup gives up to 3 transitions (and thus 3 time intervals)
# and 3 possible acquisition/clearance events per individual.
#
# Round 1 is treated as baseline (t0), so no ACQ or clearance is defined  at t0.

data_complete <- data_complete %>% arrange(menage_id_member, round)
data_complete$ACQ <- NA
for (i in 2:nrow(data_complete)) {
  if (data_complete$menage_id_member[i] == data_complete$menage_id_member[i - 1]) {
    if (data_complete$esble[i - 1] == 0 && data_complete$esble[i] == 1) {
      data_complete$ACQ[i] <- 1
    } else {
      data_complete$ACQ[i] <- 0
    }
  }
}
data_complete<- data_complete  %>% 
  relocate(  ACQ , .after= dt_n)

##############
# Define clearance (loss of colonization)
data_complete$clearance <- NA
n <- nrow(data_complete)

for (i in 2:n) {
  if (data_complete$menage_id_member[i] == data_complete$menage_id_member[i - 1]) {
    if (data_complete$esble[i - 1] == 1) {
      if (data_complete$esble[i] == 0) {
        data_complete$clearance[i] <- 1  # clearance happened
      } else {
        data_complete$clearance[i] <- 0  # still colonized
      }
    } else {
      data_complete$clearance[i] <- 0  # no clearance possible if previously not colonised
    }
  }
}


data_complete<- data_complete  %>% 
  relocate( clearance , .after= ACQ)
##########
# ------------------------------------------------------------------------------
# Step 4: WASH Variable Harmonization Across Rounds
#
# Rationale:
# Many WASH indicators (and SES) are only collected or meaningfully answered in rounds 1 and 4.
# To ensure complete covariate information for modeling time-varying effects:
#   -> Copy values from round 1 to round 2 (assume no change in short term)
#   -> Copy values from round 4 to round 3
#
# This assumes WASH characteristics are stable between those rounds for individuals,
# which is reasonable for structural household-level indicators.
#
# Variables included:
# - ses.quintile
# - main.drinking.water.dry.binary
# - main.drinking.water.rainy.binary
# - cleaning.water.storage.binary
# - correct.handwashing.binary
# - improved.sanitation.binary
# - livestock.access.house.binary
# - animal.excrement.floor.binary
# - filling round 3 from round 1 when round 4 is missing
# - fill round 2 from round 4 if round 1 is missing.
# ------------------------------------------------------------------------------

wash_vars <- c(
  #"ses.quintile",
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)

# Loop through each WASH variable and fill round 2 from 1, and round 3 from 4
for (var in wash_vars) {
  data_complete <- data_complete %>%
    group_by(menage_id_member) %>%
    group_modify(~ {
      df <- .x
      round1_value <- first(df[[var]][df$round == 1], default = NA)
      round4_value <- first(df[[var]][df$round == 4], default = NA)
      
      df[[var]] <- dplyr::case_when(
        # Fill round 2: use round 1 first, then round 4 if round 1 is missing
        df$round == 2 & is.na(df[[var]]) & !is.na(round1_value) ~ round1_value,
        df$round == 2 & is.na(df[[var]]) & is.na(round1_value) & !is.na(round4_value) ~ round4_value,
        
        # Fill round 3: use round 4 first, then round 1 if round 4 is missing
        df$round == 3 & is.na(df[[var]]) & !is.na(round4_value) ~ round4_value,
        df$round == 3 & is.na(df[[var]]) & is.na(round4_value) & !is.na(round1_value) ~ round1_value,
        
        TRUE ~ df[[var]]
      )
      df
    }) %>%
    ungroup()
}



# Check one example variable before and after
data_complete %>%
  filter(round %in% 1:4) %>%
  select(menage_id_member, round, main.drinking.water.dry.binary) %>%
  arrange(menage_id_member, round) %>%
  head(12)

#Move all WASH variables right after the "round" column
data_complete <- data_complete %>%
  relocate(all_of(wash_vars), .after = round)
data_complete %>%
  filter(menage_id_member == 202) %>%   # example from your screenshot
  select(menage_id_member, round, main.drinking.water.dry.binary)

table(data_complete$correct.handwashing.binary)

###############

## 5) Add rainy_effect/seasonal
# ------------------------------------------------------------------------------
# Create a "rainy_effect" time-varying covariate to capture rainy season exposure
#
# Goal: For each interval (between two consecutive observations of the same individual),
#       determine whether the time was predominantly during rainy months.
#
# Step-by-step:
# 1. Extract all dates where the "rainy" field is marked "yes".
# 2. Convert those dates into a unique list of rainy months ("YYYY-MM").
# 3. For each interval between two rounds per individual:
#     - Extract the calendar months that span the interval.
#     - Count how many of those are rainy vs. non-rainy months.
#     -If rainy months are equal or more frequent, set rainy_effect = 1 (rainy-dominant period),
#       otherwise set rainy_effect = 0.
#
# This variable is used later in the Cox model to assess whether acquisition risk
# is associated with being observed during a predominantly rainy period.
# ------------------------------------------------------------------------------

rainy_dates <- unique(data_complete[data_complete[["rainy"]] == "yes", "date.use", drop = TRUE])
print("Unique Rainy Dates (Before Conversion):")
print(rainy_dates)

rainy_dates <- as.Date(rainy_dates, format = "%Y-%m-%d")  # Adjusted format
print("Rainy Dates (After Conversion to Date):")
print(rainy_dates)

# Format rainy dates to "YYYY-MM" for consistency
rainy_months_formatted <- unique(format(rainy_dates, "%Y-%m"))
print(rainy_months_formatted)

data_complete$rainy_effect <- 0

for (i in 2:nrow(data_complete)) {
  if (data_complete$menage_id_member[i] == data_complete$menage_id_member[i - 1]) {
    start_date <- as.Date(data_complete$date.use[i - 1])
    end_date <- as.Date(data_complete$date.use[i])
    
    # Skip NA or reversed dates
    if (!is.na(start_date) && !is.na(end_date) && start_date < end_date) {
      # Generate monthly sequence
      months_seq <- seq(from = as.Date(format(start_date, "%Y-%m-01")),
                        to = as.Date(format(end_date, "%Y-%m-01")),
                        by = "month")
      month_labels <- format(months_seq, "%Y-%m")
      
      rainy_count <- sum(month_labels %in% rainy_months_formatted)
      non_rainy_count <- length(month_labels) - rainy_count
      
      data_complete$rainy_effect[i] <- ifelse(rainy_count >= non_rainy_count, 1, 0)
    }
  }
}

data_complete <- data_complete %>%
  relocate(rainy_effect, .after = rainy)
table(data_complete$rainy_effect)
#table(data_complete$rainy)

# ------------------------------------------------------------------------------
# Handle Missingness and Factor Levels for WASH and SES Variables
#
# Purpose:
# Prepare categorical variables (WASH indicators + SES) for modeling by:
#   - Converting each variable to a factor` (categorical type),
#   - Explicitly creating a "Missing" category for NA values,
#   - Defining meaningful level orderings (e.g., Improved > Unimproved > Missing).
#
# Why this matters:
#   - The coxme() model (or any survival model) cannot handle NA values directly.
#   - Simply removing rows with missing WASH data may introduce selection bias.
#   - Instead, we retain these observations and treat "Missing" as an informative category.
#     (this also helps detect if missingness is associated with the outcome).
#
# Included variables:
# - ses.quintile  socioeconomic status (ordinal from "lowest" to "highest")
# - WASH indicators (e.g., water source, sanitation, hygiene, animal exposure)
#
# After this block:
# - All modeling variables will be factors with 3 levels:
#     e.g., "Improved", "Unimproved", "Missing"
# ------------------------------------------------------------------------------

# Factor WASH variables and include missing as a level
data_complete <- data_complete %>%
  mutate(
    # ses.quintile = factor(
    #   ifelse(is.na(ses.quintile), "Missing", ses.quintile),
    #   levels = c("lowest", "second", "third", "fourth", "highest", "Missing")
    # ),
    main.drinking.water.dry.binary = factor(
      ifelse(is.na(main.drinking.water.dry.binary), "Missing", main.drinking.water.dry.binary),
      levels = c("Improved", "Unimproved", "Missing")
    ),
    main.drinking.water.rainy.binary = factor(
      ifelse(is.na(main.drinking.water.rainy.binary), "Missing", main.drinking.water.rainy.binary),
      levels = c("Improved", "Unimproved", "Missing")
    ),
    cleaning.water.storage.binary = factor(
      ifelse(is.na(cleaning.water.storage.binary), "Missing",
             ifelse(cleaning.water.storage.binary == "Yes", "Treated", "Not treated")),
      levels = c("Treated", "Not treated", "Missing")
    ),
    correct.handwashing.binary = factor(
      ifelse(is.na(correct.handwashing.binary), "Missing",
             ifelse(correct.handwashing.binary == "Yes", "Correct", "Not Correct")),
      levels = c("Correct", "Not Correct", "Missing")
    ),
    improved.sanitation.binary = factor(
      ifelse(is.na(improved.sanitation.binary), "Missing",
             ifelse(improved.sanitation.binary == "Yes", "Improved", "Unimproved")),
      levels = c("Improved", "Unimproved", "Missing")
    ),
    livestock.access.house.binary = factor(
      ifelse(is.na(livestock.access.house.binary), "Missing",
             ifelse(livestock.access.house.binary == "No", "No Access", "Access")),
      levels = c("No Access", "Access", "Missing")
    ),
    animal.excrement.floor.binary = factor(
      ifelse(is.na(animal.excrement.floor.binary), "Missing",
             ifelse(animal.excrement.floor.binary == "No", "Not Exposed", "Exposed")),
      levels = c("Not Exposed", "Exposed", "Missing")
    )
  )


####################

lapply(data_complete[c("main.drinking.water.dry.binary", "cleaning.water.storage.binary", "main.drinking.water.rainy.binary","improved.sanitation.binary",
              "correct.handwashing.binary", "livestock.access.house.binary", "animal.excrement.floor.binary")], table, useNA = "always")

summary(data_complete[, c(
  "age",
  "sexe",
  "main.drinking.water.dry.binary",
  "cleaning.water.storage.binary",
  "main.drinking.water.rainy.binary",
  "correct.handwashing.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary",
  "improved.sanitation.binary"
)])

################
# Drop unused factor levels to avoid model convergence issues
# Filter out "Missing" ses.quintile
# data_complete_model <- data_complete_model %>%
#   filter(ses.quintile != "Missing")

# Drop unused levels across all factor variables again
data_complete <- data_complete %>%
  mutate(
    rainy_effect = as.factor(rainy_effect),  
    #ses.quintile = droplevels(ses.quintile),  # include this if used
    main.drinking.water.dry.binary = droplevels(main.drinking.water.dry.binary),
    main.drinking.water.rainy.binary = droplevels(main.drinking.water.rainy.binary),
    cleaning.water.storage.binary = droplevels(cleaning.water.storage.binary),
    correct.handwashing.binary = droplevels(correct.handwashing.binary),
    improved.sanitation.binary = droplevels(improved.sanitation.binary),
    livestock.access.house.binary = droplevels(livestock.access.house.binary),
    animal.excrement.floor.binary = droplevels(animal.excrement.floor.binary)
  )


t1 <- table(data_complete$livestock.access.house.binary)
print(t1)

#################
###Debugging
lapply(data_complete[, c(
  "rainy_effect",
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "improved.sanitation.binary",
  "correct.handwashing.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)], function(x) table(x))

table(data_complete$main.drinking.water.dry.binary, useNA = "always")
# Vector of all WASH binary indicators
wash_vars <- c(
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)

# Remove rows where any WASH variable is "Missing"
data_model <- data_complete %>%
  filter(if_all(all_of(wash_vars), ~ . != "Missing")) %>%
  droplevels()

##last)# Remove round 1 rows — they don't represent any observed transition period, and having NA values in dt_n, ACQ, and clearance for round 1 may cause confusion and could lead to errors or warnings.
data_model_tv <- data_model # For timevarying dataset keep round one

data_model <- data_model %>% filter(round != 1)


#########################################################
# FIT MODEL USING RAINY VARIABLE
#########################################################
table(data_model$correct.handwashing.binary, data_model$ACQ)
table(data_model$main.drinking.water.rainy.binary, data_model$ACQ)
table(data_model$main.drinking.water.dry.binary, data_model$ACQ)
table(data_model$improved.sanitation.binary, data_model$ACQ)
table(data_model$rainy, data_model$ACQ)
table(data_model$rainy_effect, data_model$ACQ)

# Fit model with age and sex - time-to-event
# coxme_mod <- coxme(Surv(days, ACQ) ~ age + sexe + rainy_effect +
#                      main.drinking.water.dry.binary +
#                      main.drinking.water.rainy.binary +
#                      cleaning.water.storage.binary +
#                      correct.handwashing.binary +
#                      improved.sanitation.binary +
#                      livestock.access.house.binary +
#                      animal.excrement.floor.binary +
#                      (1 | HouseID),
#                    data = data_model)
# 
# # Fit model without age and sex
# coxme_mod_nd <- coxme(Surv(days, ACQ) ~ rainy_effect +
#                      main.drinking.water.dry.binary +
#                      main.drinking.water.rainy.binary +
#                      cleaning.water.storage.binary +
#                      correct.handwashing.binary +
#                      improved.sanitation.binary +
#                      livestock.access.house.binary +
#                      animal.excrement.floor.binary +
#                      (1 | HouseID),
#                    data = data_model)
# 
# AIC(coxme_mod, coxme_mod_nd)
# 
# # Summarize
# coxm_summary <- summary(coxme_mod_nd)
# coxm_summary
# exp(coxm_summary$coefficients)
# 
# # Clean labels (removed age & sex)
# pretty_labels <- c(
#   "rainy_effect1" = "Seasonality: Rainy season exposure",
#   "main.drinking.water.dry.binaryUnimproved" = "W1a: Drinking water (dry): Unimproved",
#   "main.drinking.water.rainy.binaryUnimproved" = "W1b: Drinking water (rainy): Unimproved",
#   "cleaning.water.storage.binaryNot treated" = "W2: Stored water: Not treated",
#   "correct.handwashing.binaryNot Correct" = "W3: Handwashing: Incorrect",
#   "improved.sanitation.binaryUnimproved" = "W4: Sanitation: Unimproved",
#   "livestock.access.house.binaryAccess" = "W5: Livestock access inside",
#   "animal.excrement.floor.binaryExposed" = "W6: Animal excrement on floor"
# )
# 
# ordered_terms <- names(pretty_labels)
# 
# # Build dataframe for forest plot
# df_forest <- data.frame(
#   term = rownames(coxm_summary$coefficients),
#   estimate = coxm_summary$coefficients[, "coef"],
#   std.error = coxm_summary$coefficients[, "se(coef)"]
# ) %>%
#   filter(term %in% ordered_terms) %>%
#   mutate(
#     HR = exp(estimate),
#     HR_50_low = exp(estimate - 0.674 * std.error),
#     HR_50_high = exp(estimate + 0.674 * std.error),
#     HR_80_low = exp(estimate - 1.28 * std.error),
#     HR_80_high = exp(estimate + 1.28 * std.error),
#     HR_95_low = exp(estimate - 1.96 * std.error),
#     HR_95_high = exp(estimate + 1.96 * std.error),
#     term = factor(term, levels = ordered_terms),
#     label = pretty_labels[as.character(term)]
#   )
# df_forest
# # Forest plot
# base_size = 16
# p_rainy <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
#   geom_rect(aes(xmin = HR_95_low, xmax = HR_95_high,
#                 ymin = as.numeric(fct_rev(label)) - 0.3,
#                 ymax = as.numeric(fct_rev(label)) + 0.3,
#                 fill = "95% CI"), alpha = 0.4) +
#   geom_rect(aes(xmin = HR_80_low, xmax = HR_80_high,
#                 ymin = as.numeric(fct_rev(label)) - 0.2,
#                 ymax = as.numeric(fct_rev(label)) + 0.2,
#                 fill = "80% CI"), alpha = 0.6) +
#   geom_rect(aes(xmin = HR_50_low, xmax = HR_50_high,
#                 ymin = as.numeric(fct_rev(label)) - 0.1,
#                 ymax = as.numeric(fct_rev(label)) + 0.1,
#                 fill = "50% CI"), alpha = 0.8) +
#   geom_point(size = 3, color = "black") +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
#   scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.25)) +
#   scale_fill_manual(
#     name = "Shading Level",
#     values = c("50% CI" = "#7ca982", "80% CI" = "#cde5b2", "95% CI" = "#eeeec8")
#   ) +
#   labs(
#     title = "",
#     x = "Hazard Ratio ",
#     y = NULL
#   ) +
#   theme_cowplot() +
#   theme(
#     legend.position   = "bottom",
#     legend.box        = "horizontal",
#     legend.box.just   = "center",
#     legend.justification = "center",
#     strip.text.y = element_text(),
#     plot.title    = element_text(size = base_size + 4, face = "bold"),
#     plot.subtitle = element_text(size = base_size + 2),
#     axis.title    = element_text(size = base_size + 2),
#     axis.text     = element_text(size = base_size + 2),
#     plot.caption  = element_text(size = base_size - 3),
#     legend.text   = element_text(size = 16),
#     legend.title  = element_blank(),
#     panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
#     panel.grid.minor = element_line(linetype = "dotted", color = "grey90")
#   )
# 
# print(p_rainy)

# Save forest plot as high-resolution TIFF
# ggsave(
#   filename = "FigureX_forestplot_acquisition_no_age_sex.tiff", 
#   plot = p_rainy, 
#   device = "tiff",
#   dpi = 300,           # ensures print quality
#   width = 8, height = 6, units = "in", 
#   compression = "lzw"  # lossless compression
# )
# getwd()
# list.files()


#-----------------------------------------------------------------------------
# Calendar-time - assuming all individuals entered the study at the same time
#-----------------------------------------------------------------------------

# Fit model with age and sex - time-to-event
coxme_mod_ct <- coxme(Surv(dt_n, ACQ) ~ age + sexe + rainy_effect +
                     main.drinking.water.dry.binary +
                     main.drinking.water.rainy.binary +
                     cleaning.water.storage.binary +
                     correct.handwashing.binary +
                     improved.sanitation.binary +
                     livestock.access.house.binary +
                     animal.excrement.floor.binary +
                     (1 | HouseID) + (1 | menage_id_member),,
                   data = data_model)

# Fit model without age and sex
coxme_mod_nd_ct <- coxme(Surv(dt_n, ACQ) ~ rainy_effect +
                        main.drinking.water.dry.binary +
                        main.drinking.water.rainy.binary +
                        cleaning.water.storage.binary +
                        correct.handwashing.binary +
                        improved.sanitation.binary +
                        livestock.access.house.binary +
                        animal.excrement.floor.binary +
                        (1 | HouseID) + (1 | menage_id_member), # Earlier version did not incorporate a clustering at household member id
                      data = data_model)

anova(coxme_mod_ct, coxme_mod_nd_ct) # Similar model fit
AIC(coxme_mod_ct, coxme_mod_nd_ct)

# Summarize
coxm_summary <- summary(coxme_mod_ct)
coxm_summary
exp(coxm_summary$coefficients)

# Clean labels (removed age & sex)
pretty_labels <- c(
  "rainy_effect1" = "Seasonality: Rainy season exposure",
  "main.drinking.water.dry.binaryUnimproved" = "W1a: Drinking water (dry): Unimproved",
  "main.drinking.water.rainy.binaryUnimproved" = "W1b: Drinking water (rainy): Unimproved",
  "cleaning.water.storage.binaryNot treated" = "W2: Stored water: Not treated",
  "correct.handwashing.binaryNot Correct" = "W3: Handwashing: Incorrect",
  "improved.sanitation.binaryUnimproved" = "W4: Sanitation: Unimproved",
  "livestock.access.house.binaryAccess" = "W5: Livestock access inside",
  "animal.excrement.floor.binaryExposed" = "W6: Animal excrement on floor"
)

ordered_terms <- names(pretty_labels)

# Build dataframe for forest plot
df_forest <- data.frame(
  term = rownames(coxm_summary$coefficients),
  estimate = coxm_summary$coefficients[, "coef"],
  std.error = coxm_summary$coefficients[, "se(coef)"]
) %>%
  filter(term %in% ordered_terms) %>%
  mutate(
    HR = exp(estimate),
    HR_50_low = exp(estimate - 0.674 * std.error),
    HR_50_high = exp(estimate + 0.674 * std.error),
    HR_80_low = exp(estimate - 1.28 * std.error),
    HR_80_high = exp(estimate + 1.28 * std.error),
    HR_95_low = exp(estimate - 1.96 * std.error),
    HR_95_high = exp(estimate + 1.96 * std.error),
    term = factor(term, levels = ordered_terms),
    label = pretty_labels[as.character(term)]
  )
df_forest
# Forest plot
base_size = 16
p_rainy <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
  geom_rect(aes(xmin = HR_95_low, xmax = HR_95_high,
                ymin = as.numeric(fct_rev(label)) - 0.3,
                ymax = as.numeric(fct_rev(label)) + 0.3,
                fill = "95% CI"), alpha = 0.4) +
  geom_rect(aes(xmin = HR_80_low, xmax = HR_80_high,
                ymin = as.numeric(fct_rev(label)) - 0.2,
                ymax = as.numeric(fct_rev(label)) + 0.2,
                fill = "80% CI"), alpha = 0.6) +
  geom_rect(aes(xmin = HR_50_low, xmax = HR_50_high,
                ymin = as.numeric(fct_rev(label)) - 0.1,
                ymax = as.numeric(fct_rev(label)) + 0.1,
                fill = "50% CI"), alpha = 0.8) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.25)) +
  scale_fill_manual(
    name = "Shading Level",
    values = c("50% CI" = "#7ca982", "80% CI" = "#cde5b2", "95% CI" = "#eeeec8")
  ) +
  labs(
    title = "",
    x = "Hazard Ratio ",
    y = NULL
  ) +
  theme_cowplot() +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.justification = "center",
    strip.text.y = element_text(),
    plot.title    = element_text(size = base_size + 4, face = "bold"),
    plot.subtitle = element_text(size = base_size + 2),
    axis.title    = element_text(size = base_size + 2),
    axis.text     = element_text(size = base_size + 2),
    plot.caption  = element_text(size = base_size - 3),
    legend.text   = element_text(size = 16),
    legend.title  = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor = element_line(linetype = "dotted", color = "grey90")
  )

print(p_rainy)

# Save forest plot as high-resolution TIFF
# ggsave(
#   filename = "FigureX_forestplot_acquisition_no_age_sex.tiff", 
#   plot = p_rainy, 
#   device = "tiff",
#   dpi = 300,           # ensures print quality
#   width = 8, height = 6, units = "in", 
#   compression = "lzw"  # lossless compression
# )
# getwd()
# list.files()

ggsave(
  filename = "./Output/Figures_and_tables/Paper/Figure4_CoxModel_WASH_effect_timetoevent.png",
  plot = p_rainy, 
  device = "png",
  dpi = 300,           # ensures print quality
  width = 14, height = 8, units = "in", 
  # compression = "lzw",  # lossless compression 
  bg = "white"
)

write.csv(df_forest, file="./Output/Model_results/Model_summaries/CoxPH_calendar_time.csv")

#---------------------------------------------------------------------------------------------------
# Calenar-time alternative - left-truncuated staggered entry (accurate representation of the study)
#---------------------------------------------------------------------------------------------------

period <- 365  # yearly seasonality; set 366 if you insist on leap-year scale

d_long <- data_model_tv %>%
  mutate(date.use = as.Date(date.use)) %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    start_date = cal_start,
    end_date   = cal_stop,
    dt_start   = lag(dt_days),
    dt_end     = dt_days,
    event      = ACQ
  ) %>%
  ungroup() %>%
  # keep proper intervals only
  filter(!is.na(start_date), !is.na(dt_start), end_date > start_date, dt_end > dt_start) %>%
  select(-c(rainy_effect))

d_long <- d_long[,c(1,3:19,31,116:128)]

# Interval id per person
d_long_id <- d_long %>%
  group_by(menage_id_member) %>%
  mutate(interval_id = row_number()) %>%
  ungroup()

d_long_id <- d_long_id %>%
  mutate(
    start_date = global_origin + start_date,  # start_date is numeric days
    end_date   = global_origin + end_date
  )

# Rainy months vector like "YYYY-MM" from your rainy_dates
rainy_months_formatted <- unique(format(as.Date(rainy_dates), "%Y-%m"))

# Generate one row per interval–month (end exclusive: [start, end))
monthly_grid <- d_long_id %>%
  mutate(
    start_date = as.Date(start_date),
    end_date   = as.Date(end_date)
  ) %>%
  transmute(
    menage_id_member, interval_id, start_date, end_date,
    month_start = map2(start_date, end_date, ~{
      if (is.na(.x) || is.na(.y)) return(as.Date(character()))
      from <- floor_date(.x, "month")
      # subtract 1 day so the end month is included when end_date is the first of a month
      to   <- floor_date(.y - days(1), "month")
      if (from > to) return(as.Date(character()))
      seq(from, to, by = "month")
    })
  ) %>%
  unnest(cols = month_start) %>%
  mutate(month_label = format(month_start, "%Y-%m"))


# Majority-by-months summary per interval
rainy_by_interval <- monthly_grid %>%
  group_by(menage_id_member, interval_id) %>%
  summarise(
    total_months  = n_distinct(month_label),
    rainy_months  = sum(month_label %in% rainy_months_formatted),
    rainy_effect  = as.integer(rainy_months > total_months / 2), # use >= for ties-as-rainy
    .groups = "drop"
  )

# Join back
d_long_with_rain <- d_long_id %>%
  left_join(rainy_by_interval, by = c("menage_id_member", "interval_id"))


vars_fixed <- c(
  "age","sexe","rainy_effect",
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)

# Count unique non-NA values
nuniq <- function(x) length(unique(x[!is.na(x)]))

# check variation (handles numeric, logical, factors)
check_tbl <- data.frame(
  var = vars_fixed,
  n_unique = sapply(d_long_with_rain[vars_fixed], nuniq),
  stringsAsFactors = FALSE
) |>
  transform(is_constant = n_unique <= 1)

check_tbl

fit_season_ct <- coxme(
  Surv(cal_start, cal_stop, event) ~
    age + sexe +
    rainy_effect +
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) +
    (1 | HouseID),
  data = d_long_with_rain
)

# Not accounting for dependent observations at individual-level
fit_season_ct_ni <- coxme(
  Surv(cal_start, cal_stop, event) ~
    age + sexe +
    rainy_effect +
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | HouseID),
  data = d_long_with_rain
)

fit_season_nd_ct <- coxme(
  Surv(cal_start, cal_stop, event) ~
    rainy_effect +
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) +
    (1 | HouseID),
  data = d_long_with_rain
)

summary(fit_season_ct)
summary(fit_season_ct_ni)
summary(fit_season_nd_ct)

AIC(fit_season_ct, fit_season_nd_ct, fit_season_ct_ni) # Hardly difference in fit, keep with age and sex

# Summarise model
coxm_summary <- summary(fit_season_ct)

# Pretty labels (edit as you like)
pretty_labels <- c(
  "rainy_effect" = "Seasonality: Rainy season exposure",
  "main.drinking.water.dry.binaryUnimproved"   = "W1a: Drinking water (dry): Unimproved",
  "main.drinking.water.rainy.binaryUnimproved" = "W1b: Drinking water (rainy): Unimproved",
  "cleaning.water.storage.binaryNot treated"   = "W2: Stored water: Not treated",
  "correct.handwashing.binaryNot Correct"      = "W3: Handwashing: Incorrect",
  "improved.sanitation.binaryUnimproved"       = "W4: Sanitation: Unimproved",
  "livestock.access.house.binaryAccess"        = "W5: Livestock access inside",
  "animal.excrement.floor.binaryExposed"       = "W6: Animal excrement on floor"
)

# Order to display (top→bottom)
ordered_terms <- names(pretty_labels)

# Build forest-plot dataframe
tb <- as.data.frame(coxm_summary$coefficients)
tb$term <- rownames(coxm_summary$coefficients)

df_forest <- tb %>%
  filter(term %in% ordered_terms) %>%
  rename(estimate = coef, std.error = `se(coef)`) %>%
  mutate(
    HR         = exp(estimate),
    HR_50_low  = exp(estimate - 0.674 * std.error),
    HR_50_high = exp(estimate + 0.674 * std.error),
    HR_80_low  = exp(estimate - 1.2816 * std.error),
    HR_80_high = exp(estimate + 1.2816 * std.error),
    HR_95_low  = exp(estimate - 1.96 * std.error),
    HR_95_high = exp(estimate + 1.96 * std.error),
    term       = factor(term, levels = ordered_terms),
    label      = pretty_labels[as.character(term)]
  )

# Plot (shaded CI rectangles + point estimate)
base_size = 16
p <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
  geom_rect(aes(xmin = HR_95_low, xmax = HR_95_high,
                ymin = as.numeric(fct_rev(label)) - 0.3,
                ymax = as.numeric(fct_rev(label)) + 0.3,
                fill = "95% CI"), alpha = 0.4) +
  geom_rect(aes(xmin = HR_80_low, xmax = HR_80_high,
                ymin = as.numeric(fct_rev(label)) - 0.2,
                ymax = as.numeric(fct_rev(label)) + 0.2,
                fill = "80% CI"), alpha = 0.6) +
  geom_rect(aes(xmin = HR_50_low, xmax = HR_50_high,
                ymin = as.numeric(fct_rev(label)) - 0.1,
                ymax = as.numeric(fct_rev(label)) + 0.1,
                fill = "50% CI"), alpha = 0.8) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.25)) +
  scale_fill_manual(
    name = "Shading Level",
    values = c("50% CI" = "#7ca982", "80% CI" = "#cde5b2", "95% CI" = "#eeeec8")
  ) +
  labs(
    title = "",
    x = "Hazard Ratio ",
    y = NULL
  ) +
  theme_cowplot() +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.justification = "center",
    strip.text.y = element_text(),
    plot.title    = element_text(size = base_size + 4, face = "bold"),
    plot.subtitle = element_text(size = base_size + 2),
    axis.title    = element_text(size = base_size + 2),
    axis.text     = element_text(size = base_size + 2),
    plot.caption  = element_text(size = base_size - 3),
    legend.text   = element_text(size = 16),
    legend.title  = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor = element_line(linetype = "dotted", color = "grey90")
  )

print(p)

ggsave(
  filename = "./Output/Figures_and_tables/Paper/Figure4_CoxModel_WASH_effect_calendar_time_alt.png",
  plot = p, 
  device = "png",
  dpi = 300,           # ensures print quality
  width = 14, height = 8, units = "in", 
  # compression = "lzw",  # lossless compression 
  bg = "white"
)

write.csv(df_forest, file="./Output/Model_results/Model_summaries/CoxPH_calendar_time.alt.csv")



#--------------------------------------------------
# Time-to-event 
#--------------------------------------------------
# (see https://doi.org/10.1186/s12874-016-0199-y)
# Fitting two time scales may be tricky. We are doing that somewhat by including seasonality as a binary variable based on 
# The calendar dates between interval. Just fitting to date.use would not make sense as 

# Build start–stop intervals per individual 
# Assumes columns:
#   menage_id        = individual or household ID (use your preferred random-effect ID)
#   date.consent     = visit date (Date or character ISO)
#   dt_n             = time since first sample (days) at the visit
#   ACQ              = event indicator at END of the interval (1 if acquired by this visit)
#   age, sexe, ...   = covariates you want in the model
# If your ID is ‘menage_id_member’ or random effect is ‘HouseID’, just swap names below
# 
period <- 365  # yearly seasonality; set 366 if you insist on leap-year scale

d_long <- data_model_tv %>%
  mutate(date.use = as.Date(date.use)) %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    start_date = lag(date.use),
    end_date   = date.use,
    dt_start   = lag(days),
    dt_end     = days,
    event      = ACQ
  ) %>%
  ungroup() %>%
  # keep proper intervals only
  filter(!is.na(start_date), !is.na(dt_start), end_date > start_date, dt_end > dt_start) %>%
  select(-c(rainy_effect))

d_long <- d_long[,c(1,3:19,31,116:128)]


# Interval id per person
d_long_id <- d_long %>%
  group_by(menage_id_member) %>%
  mutate(interval_id = row_number()) %>%
  ungroup()


# Rainy months vector like "YYYY-MM" from your rainy_dates
rainy_months_formatted <- unique(format(as.Date(rainy_dates), "%Y-%m"))

# Generate one row per interval–month (end exclusive: [start, end))
monthly_grid <- d_long_id %>%
  transmute(
    menage_id_member, interval_id, start_date, end_date,
    month_start = purrr::map2(
      start_date, end_date,
      ~ seq(lubridate::floor_date(.x, "month"), lubridate::floor_date(.y - days(1), "month"), by = "month")
    )
  ) %>%
  unnest(month_start) %>%
  mutate(month_label = format(month_start, "%Y-%m"))

# Majority-by-months summary per interval
rainy_by_interval <- monthly_grid %>%
  group_by(menage_id_member, interval_id) %>%
  summarise(
    total_months  = n_distinct(month_label),
    rainy_months  = sum(month_label %in% rainy_months_formatted),
    rainy_effect  = as.integer(rainy_months > total_months / 2), # use >= for ties-as-rainy
    .groups = "drop"
  )

# Join back
d_long_with_rain <- d_long_id %>%
  left_join(rainy_by_interval, by = c("menage_id_member", "interval_id"))


vars_fixed <- c(
  "age","sexe","rainy_effect",
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)

# Count unique non-NA values
nuniq <- function(x) length(unique(x[!is.na(x)]))

# check variation (handles numeric, logical, factors)
check_tbl <- data.frame(
  var = vars_fixed,
  n_unique = sapply(d_long_with_rain[vars_fixed], nuniq),
  stringsAsFactors = FALSE
) |>
  transform(is_constant = n_unique <= 1)

check_tbl

fit_season <- coxme(
  Surv(dt_start, dt_end, event) ~
    age + sexe + 
    rainy_effect +
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) +
    (1 | HouseID),
  data = d_long_with_rain
)

summary(fit_season)

# Summarise model
coxm_summary <- summary(fit_season)

# Pretty labels (edit as you like)
pretty_labels <- c(
  "rainy_effect" = "Seasonality: Rainy season exposure",
  "main.drinking.water.dry.binaryUnimproved"   = "W1a: Drinking water (dry): Unimproved",
  "main.drinking.water.rainy.binaryUnimproved" = "W1b: Drinking water (rainy): Unimproved",
  "cleaning.water.storage.binaryNot treated"   = "W2: Stored water: Not treated",
  "correct.handwashing.binaryNot Correct"      = "W3: Handwashing: Incorrect",
  "improved.sanitation.binaryUnimproved"       = "W4: Sanitation: Unimproved",
  "livestock.access.house.binaryAccess"        = "W5: Livestock access inside",
  "animal.excrement.floor.binaryExposed"       = "W6: Animal excrement on floor"
)

# Order to display (top→bottom)
ordered_terms <- names(pretty_labels)

# Build forest-plot dataframe
tb <- as.data.frame(coxm_summary$coefficients)
tb$term <- rownames(coxm_summary$coefficients)

df_forest <- tb %>%
  filter(term %in% ordered_terms) %>%
  rename(estimate = coef, std.error = `se(coef)`) %>%
  mutate(
    HR         = exp(estimate),
    HR_50_low  = exp(estimate - 0.674 * std.error),
    HR_50_high = exp(estimate + 0.674 * std.error),
    HR_80_low  = exp(estimate - 1.2816 * std.error),
    HR_80_high = exp(estimate + 1.2816 * std.error),
    HR_95_low  = exp(estimate - 1.96 * std.error),
    HR_95_high = exp(estimate + 1.96 * std.error),
    term       = factor(term, levels = ordered_terms),
    label      = pretty_labels[as.character(term)]
  )

# Plot (shaded CI rectangles + point estimate)
base_size = 16
p <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
  geom_rect(aes(xmin = HR_95_low, xmax = HR_95_high,
                ymin = as.numeric(fct_rev(label)) - 0.3,
                ymax = as.numeric(fct_rev(label)) + 0.3,
                fill = "95% CI"), alpha = 0.4) +
  geom_rect(aes(xmin = HR_80_low, xmax = HR_80_high,
                ymin = as.numeric(fct_rev(label)) - 0.2,
                ymax = as.numeric(fct_rev(label)) + 0.2,
                fill = "80% CI"), alpha = 0.6) +
  geom_rect(aes(xmin = HR_50_low, xmax = HR_50_high,
                ymin = as.numeric(fct_rev(label)) - 0.1,
                ymax = as.numeric(fct_rev(label)) + 0.1,
                fill = "50% CI"), alpha = 0.8) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.25)) +
  scale_fill_manual(
    name = "Shading Level",
    values = c("50% CI" = "#7ca982", "80% CI" = "#cde5b2", "95% CI" = "#eeeec8")
  ) +
  labs(
    title = "",
    x = "Hazard Ratio ",
    y = NULL
  ) +
  theme_cowplot() +
  theme(
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.justification = "center",
    strip.text.y = element_text(),
    plot.title    = element_text(size = base_size + 4, face = "bold"),
    plot.subtitle = element_text(size = base_size + 2),
    axis.title    = element_text(size = base_size + 2),
    axis.text     = element_text(size = base_size + 2),
    plot.caption  = element_text(size = base_size - 3),
    legend.text   = element_text(size = 16),
    legend.title  = element_blank(),
    panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor = element_line(linetype = "dotted", color = "grey90")
  )

print(p)

ggsave(
  filename = "./Output/Figures_and_tables/Paper/Figure4_CoxModel_WASH_effect_timetoevent.png",
  plot = p, 
  device = "png",
  dpi = 300,           # ensures print quality
  width = 14, height = 8, units = "in", 
  # compression = "lzw",  # lossless compression 
  bg = "white"
)

write.csv(df_forest, file="./Output/Model_results/Model_summaries/CoxPH_time_to_event.csv")

