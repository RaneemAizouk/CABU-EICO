###############################################################################
# Analysis_wash_coxme_midpoint.R
###############################################################################
# Title:WASH and ESBL-E acquisition – Cox mixed-effects model (midpoint)
#
# Purpose:
#   Fit Cox proportional hazards mixed-effects models (coxme) to estimate the
#   association between household WASH indicators and ESBL-E acquisition.
#
# Observed data:
# Outcome: ESBL-E acquisition (transition 0 → 1), defined as a change from
#   ESBL-negative to ESBL-positive status between two consecutive stool sampling rounds.
#   ESBL-E clearance (1 → 0) is not modelled here; clearance is handled separately
#   in multi-state / MSM models.
#
# Model:
#  Cox proportional hazards mixed-effects models (coxme).
# Random intercepts (frailty-like) for:
#     (1) individual (menage_id_member)
#     (2) household (HouseID)
#   to account for within-person and within-household correlation.
#  WASH indicators are treated as time-varying covariates, updated at each
#   survey round and carried forward/backward within individuals.
#
# Time scale:
#  Calendar time (days since global origin: 2022-10-03).
#  Each observation interval corresponds to the time between two consecutive
#   stool sampling dates for the same individual.
#
# Event timing (interval-censored acquisition):
# Exact acquisition times are unknown and occur between two sampling dates.
#  Midpoint imputation is used for the analysis:
#     event_time = (start_date + end_date) / 2
#  Non-acquisition intervals are right-censored at end_date.
#  Left truncation is considered.
#
# Assumptions:
#  Piecewise-constant hazard within observation intervals.
#  Proportional hazards for covariate effects.
#  Primary analysis excludes rows with missing WASH indicators.
#
# Outputs:
# Fitted coxme model objects (.rds)
# Text summaries of model results
# Tables of hazard ratios (CSV)
#
# Author:   Raneem Aizouk
# Created:   April 2025
# Updated:    Dec  2025
###############################################################################

rm(list = ls())

#------------------------------------------------------------------------------
# Command-line arguments / environment variables
#------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
scenario    <- if (length(args) >= 1) args[[1]] else Sys.getenv("scenario", "Primary_midpoint")
data_source <- if (length(args) >= 2) args[[2]] else Sys.getenv("data_source", "observed")

cat("Scenario:", scenario, "\n")
cat("Data source:", data_source, "\n")

#------------------------------------------------------------------------------
# Directories
#------------------------------------------------------------------------------

if (data_source == "observed") {
  data_dir   <- "./Data/BF/clean/use_in_analyses"
  output_dir <- "./Output/Modeling/WASH_cox/Observed_data"
} else {
  stop("This script is intended for observed data.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------------------------
# Libraries
#------------------------------------------------------------------------------

pacman::p_load(
  survival,
  coxme,
  dplyr,
  tidyr,
  lubridate,
  zoo
)

#------------------------------------------------------------------------------
# Load cleaned dataset
#------------------------------------------------------------------------------

load(file.path(data_dir, "bf_esbl0123_long_all.rda"))  # loads dfls0
data <- dfls0
rm(dfls0)

cat("Rows:", nrow(data),
    "| Individuals:", n_distinct(data$menage_id_member), "\n")

#------------------------------------------------------------------------------
# Basic preparation
#------------------------------------------------------------------------------

# Round index
data <- data %>%
  mutate(round = as.integer(gsub("round_(\\d+)_arm_1", "\\1", redcap_event_name)) + 1)

# Remove individuals with completely missing age or sex
bad_ids <- data %>%
  group_by(menage_id_member) %>%
  summarise(miss_age = all(is.na(age)),
            miss_sex = all(is.na(sexe)),
            .groups = "drop") %>%
  filter(miss_age | miss_sex) %>%
  pull(menage_id_member)

data <- data %>% filter(!menage_id_member %in% bad_ids)

# Fill remaining missing age/sex within individuals
data <- data %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()

# IDs and dates
data$HouseID <- as.numeric(factor(data$menage_id))
data$sexe    <- ifelse(data$sexe %in% c("Female", "female"), 1, 0)
data$date.use <- as.Date(data$date.use)

#------------------------------------------------------------------------------
# Time variables (calendar time)
#------------------------------------------------------------------------------

global_origin <- as.Date("2022-10-03")

data <- data %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    cal_stop  = as.numeric(date.use - global_origin),
    cal_start = lag(cal_stop, default = cal_stop[1]),
    dt_days   = cal_stop - cal_start
  ) %>%
  ungroup()

#------------------------------------------------------------------------------
# Define acquisition outcome
#------------------------------------------------------------------------------

data <- data %>% arrange(menage_id_member, round)
data$ACQ <- NA_integer_

for (i in 2:nrow(data)) {
  if (data$menage_id_member[i] == data$menage_id_member[i - 1]) {
    data$ACQ[i] <- ifelse(data$esble[i - 1] == 0 & data$esble[i] == 1, 1, 0)
  }
}

#------------------------------------------------------------------------------
# WASH harmonisation
#------------------------------------------------------------------------------

wash_vars <- c(
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)

data <- data %>%
  group_by(menage_id_member) %>%
  arrange(round) %>%
  mutate(across(all_of(wash_vars), ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
  mutate(across(all_of(wash_vars), ~ zoo::na.locf(.x, fromLast = TRUE, na.rm = FALSE))) %>%
  ungroup()

# Drop missing WASH for primary analysis
data_model <- data %>%
  filter(if_all(all_of(wash_vars), ~ !is.na(.)))

#------------------------------------------------------------------------------
# Long format for Cox (midpoint)
#------------------------------------------------------------------------------

d_long <- data_model %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    start_date = cal_start,
    end_date   = cal_stop,
    event      = ACQ
  ) %>%
  ungroup() %>%
  filter(!is.na(event) & end_date > start_date)

d_long <- d_long %>%
  mutate(
    event_time = ifelse(event == 1,
                        (start_date + end_date) / 2,
                        end_date)
  )

#------------------------------------------------------------------------------
# Cox mixed-effects model (primary)
#------------------------------------------------------------------------------

cox_formula <- Surv(start_date, event_time, event) ~
  age + sexe +
  main.drinking.water.dry.binary +
  main.drinking.water.rainy.binary +
  cleaning.water.storage.binary +
  correct.handwashing.binary +
  improved.sanitation.binary +
  livestock.access.house.binary +
  animal.excrement.floor.binary +
  (1 | menage_id_member) + (1 | HouseID)

fit_coxme_midpoint <- coxme(cox_formula, data = d_long)

#------------------------------------------------------------------------------
# Save outputs
#------------------------------------------------------------------------------

saveRDS(
  fit_coxme_midpoint,
  file = file.path(output_dir, paste0("coxme_midpoint_", scenario, ".rds"))
)

sink(file = file.path(output_dir, paste0("coxme_midpoint_summary_", scenario, ".txt")))
print(summary(fit_coxme_midpoint))
sink()

# Hazard ratios table
coef_df <- as.data.frame(summary(fit_coxme_midpoint)$coefficients)
coef_df$HR <- exp(coef_df$coef)
coef_df$HR_low <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
coef_df$HR_high <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)

write.csv(
  coef_df,
  file = file.path(output_dir, paste0("coxme_midpoint_HR_", scenario, ".csv"))
)

cat("Cox mixed-effects midpoint analysis completed.\n")
