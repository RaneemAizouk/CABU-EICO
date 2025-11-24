####################################################################################################################
# THIS CODE IMPLEMENTS A COX-PROPORTIONAL HAZARD MODEL TO ASSESS THE ASSOCIATION BETWEEN WASH AND ESBL-E ACQUISITION/ sensitivity analyses
####################################################################################################################
# Author: Raneem Aizouk
# Date: April 2025
# Last updated: 12 November 2025

rm(list = ls())

pacman::p_load(survival, dplyr, tidyr, coxme, splines, mgcv, ggplot2,
               scales, cowplot, forcats, lubridate, purrr, splines)

#-----------------------------------------
# Load dataset
#-----------------------------------------
load("/Users/raizouk/Desktop/bf_esbl0123_long_all.rda")
data <- dfls0

##############
# Convert round from redcap_event_name to numeric round index
data$round <- as.integer(gsub("round_(\\d+)_arm_1", "\\1", data$redcap_event_name)) + 1
data <- data %>% relocate(round, .after = redcap_event_name)

#-----------------------------------------------------------
# Check number of individuals
#-----------------------------------------------------------
num_unique_members <- data %>%
  filter(intervention.text == "control") %>%
  summarise(unique_members = n_distinct(menage_id_member))
print(num_unique_members)

#-----------------------------------------------------------
# Find control individuals with 4 rounds of data
#-----------------------------------------------------------
four_sample_members <- data %>%
  filter(intervention.text == "control") %>%
  group_by(menage_id_member) %>%
  summarise(n_rounds = n_distinct(round)) %>%
  filter(n_rounds == 4)
n_four_sample_members <- nrow(four_sample_members)
print(n_four_sample_members)

#-----------------------------------------------------------
# Prevalence of ESBL at baseline (round 1, control group)
#-----------------------------------------------------------
esbl_baseline <- data %>%
  filter(intervention.text == "control",
         menage_id_member %in% four_sample_members$menage_id_member) %>%
  summarise(pos_rate = mean(esble == 1) * 100)
print(esbl_baseline)

#-----------------------------------------------------------
# Handle individuals with completely missing age or sex
#-----------------------------------------------------------
missing_data_individuals <- data %>%
  group_by(menage_id_member) %>%
  summarise(
    missing_age = all(is.na(age)),
    missing_sexe = all(is.na(sexe))
  ) %>%
  filter(missing_age | missing_sexe)

# Remove those individuals entirely
data_filtered <- data %>%
  filter(!menage_id_member %in% missing_data_individuals$menage_id_member)

# Fill in missing age and sex values within each individual
data_filled <- data_filtered %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()

# Verify all missing age/sexe are filled
missing_after_fill <- data_filled %>%
  filter(is.na(age) | is.na(sexe))
if (nrow(missing_after_fill) > 0) {
  print("Rows with remaining missing values:")
  print(missing_after_fill)
} else {
  print("All missing values for age and sexe have been filled.")
}
print(paste("Number of individuals after filling missing age and sexe:",
            n_distinct(data_filled$menage_id_member)))

#-----------------------------------------------------------
# Identify individuals in all rounds
#-----------------------------------------------------------
round_counts <- data_filled %>%
  group_by(menage_id_member) %>%
  summarise(rounds_participated = n_distinct(round))
print(round_counts)
max_rounds <- max(data_filled$round)
individuals_in_all_rounds <- round_counts %>%
  filter(rounds_participated == max_rounds) %>%
  pull(menage_id_member)
print(max_rounds)

data_complete <- data_filled

#-----------------------------------------------------------
# Create numeric household ID and encode sex
#-----------------------------------------------------------
data_complete$HouseID <- as.numeric(factor(data_complete$menage_id))
data_complete$sexe <- as.integer(ifelse(data_complete$sexe == "Female", 1, 0))
data_complete$month <- as.integer(data_complete$month)

#-----------------------------------------------------------
# Define calendar baseline date
#-----------------------------------------------------------
t <- as.Date("2022-10-03")
global_origin <- t

# Ensure date is Date class
data_complete$date.use <- as.Date(data_complete$date.use)

#-----------------------------------------------------------
# Compute time variables (per individual and calendar time)
#-----------------------------------------------------------
data_complete <- data_complete %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    first_date = first(date.use),
    tstop    = as.numeric(difftime(date.use, first_date, units = "days")),
    tstart   = lag(tstop, default = 0),
    cal_stop  = as.numeric(difftime(date.use, global_origin, units = "days")),
    cal_start = lag(cal_stop, default = cal_stop[1]),
    dt_days  = tstop - tstart
  ) %>%
  ungroup()

# Filter to control group data
data_complete <- data_complete %>% filter(intervention.text == "control")

# Generate numeric IDs for household and member
data_complete$HouseID <- as.numeric(factor(data_complete$menage_id))
data_complete$individual_number <- as.numeric(sub(".*-", "", data_complete$menage_id_member))
data_complete$menage_id_member <- data_complete$HouseID * 100 + data_complete$individual_number

#-----------------------------------------------------------
# Define acquisition outcome (ACQ) for transitions between rounds
#-----------------------------------------------------------
data_complete <- data_complete %>% arrange(menage_id_member, round)
data_complete$ACQ <- NA
for (i in 2:nrow(data_complete)) {
  if (data_complete$menage_id_member[i] == data_complete$menage_id_member[i - 1]) {
    data_complete$ACQ[i] <- ifelse(data_complete$esble[i - 1] == 0 && data_complete$esble[i] == 1, 1, 0)
  }
}
data_complete <- data_complete %>% relocate(ACQ, .after = dt_days)

#-----------------------------------------------------------
# Define clearance outcome for transitions
#-----------------------------------------------------------
data_complete$clearance <- NA
for (i in 2:nrow(data_complete)) {
  if (data_complete$menage_id_member[i] == data_complete$menage_id_member[i - 1]) {
    if (data_complete$esble[i - 1] == 1) {
      data_complete$clearance[i] <- ifelse(data_complete$esble[i] == 0, 1, 0)
    } else {
      data_complete$clearance[i] <- 0
    }
  }
}
data_complete <- data_complete %>% relocate(clearance, .after = ACQ)

#-----------------------------------------------------------
# WASH variable harmonization (carry values between rounds)
#-----------------------------------------------------------
wash_vars <- c("main.drinking.water.dry.binary", "main.drinking.water.rainy.binary",
               "cleaning.water.storage.binary", "correct.handwashing.binary",
               "improved.sanitation.binary", "livestock.access.house.binary",
               "animal.excrement.floor.binary")

for (var in wash_vars) {
  data_complete <- data_complete %>%
    group_by(menage_id_member) %>%
    group_modify(~ {
      df <- .x
      round1_value <- first(df[[var]][df$round == 1], default = NA)
      round4_value <- first(df[[var]][df$round == 4], default = NA)
      df[[var]] <- case_when(
        df$round == 2 & is.na(df[[var]]) & !is.na(round1_value) ~ round1_value,
        df$round == 2 & is.na(df[[var]]) & is.na(round1_value) & !is.na(round4_value) ~ round4_value,
        df$round == 3 & is.na(df[[var]]) & !is.na(round4_value) ~ round4_value,
        df$round == 3 & is.na(df[[var]]) & is.na(round4_value) & !is.na(round1_value) ~ round1_value,
        TRUE ~ df[[var]]
      )
      df
    }) %>%
    ungroup()
}

data_complete <- data_complete %>% relocate(all_of(wash_vars), .after = round)

#-----------------------------------------------------------
# Convert WASH binary variables to factors
#-----------------------------------------------------------
data_complete <- data_complete %>%
  mutate(
    main.drinking.water.dry.binary = factor(
      ifelse(is.na(main.drinking.water.dry.binary), "Missing", main.drinking.water.dry.binary),
      levels = c("Unimproved", "Improved", "Missing")
    ),
    main.drinking.water.rainy.binary = factor(
      ifelse(is.na(main.drinking.water.rainy.binary), "Missing", main.drinking.water.rainy.binary),
      levels = c("Unimproved", "Improved", "Missing")
    ),
    cleaning.water.storage.binary = factor(
      ifelse(is.na(cleaning.water.storage.binary), "Missing",
             ifelse(cleaning.water.storage.binary == "Yes", "Treated", "Not treated")),
      levels = c("Not treated", "Treated", "Missing")
    ),
    correct.handwashing.binary = factor(
      ifelse(is.na(correct.handwashing.binary), "Missing",
             ifelse(correct.handwashing.binary == "Yes", "Correct", "Not Correct")),
      levels = c("Not Correct", "Correct", "Missing")
    ),
    improved.sanitation.binary = factor(
      ifelse(is.na(improved.sanitation.binary), "Missing",
             ifelse(improved.sanitation.binary == "Yes", "Improved", "Unimproved")),
      levels = c("Unimproved", "Improved", "Missing")
    ),
    livestock.access.house.binary = factor(
      ifelse(is.na(livestock.access.house.binary), "Missing",
             ifelse(livestock.access.house.binary == "No", "No Access", "Access")),
      levels = c("Access", "No Access", "Missing")
    ),
    animal.excrement.floor.binary = factor(
      ifelse(is.na(animal.excrement.floor.binary), "Missing",
             ifelse(animal.excrement.floor.binary == "No", "Not Exposed", "Exposed")),
      levels = c("Exposed", "Not Exposed", "Missing")
    )
  )

#-----------------------------------------------------------
# Drop unused factor levels
#-----------------------------------------------------------
data_complete <- data_complete %>%
  mutate(across(all_of(wash_vars), droplevels))

#-----------------------------------------------------------
# Prepare modeling datasets
#-----------------------------------------------------------
data_model <- data_complete %>%
  filter(if_all(all_of(wash_vars), ~ . != "Missing")) %>%
  droplevels()
data_model_tv <- data_model
data_model <- data_model %>% filter(round != 1)


#  Calendar time Cox model (midpoint event timing sensitivity)

d_long <- data_model_tv %>%
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
  filter(!is.na(dt_start) & end_date > start_date & dt_end > dt_start)

# midpoint adjustment for event timing
d_long <- d_long %>%
  mutate(
    event_time_raw = (start_date + end_date) / 2,
    event_time     = ifelse(event == 1, event_time_raw, end_date),
    event_time     = pmax(event_time, start_date + 1e-6)
  )

#-----------------------------------------------------------
# Fit Cox model (with midpoint-adjusted event times)
#-----------------------------------------------------------
fit_ct_midpoint <- coxme(
  Surv(start_date, event_time, event) ~
    age + sexe +
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) + (1 | HouseID),
  data = d_long
)

#-----------------------------------------------------------
# Summarize model results
#-----------------------------------------------------------
cox_summary <- summary(fit_ct_midpoint)
print(cox_summary)
print(exp(cox_summary$coefficients))

#-----------------------------------------------------------
# (Optional) Forest plot of hazard ratios
#-----------------------------------------------------------
pretty_labels <- c(
  "main.drinking.water.dry.binaryImproved"   = "W1a: Drinking water (dry): Improved",
  "main.drinking.water.rainy.binaryImproved" = "W1b: Drinking water (rainy): Improved",
  "cleaning.water.storage.binaryTreated"     = "W2: Stored water: Treated",
  "correct.handwashing.binaryCorrect"        = "W3: Handwashing: Correct",
  "improved.sanitation.binaryImproved"       = "W4: Sanitation: Improved",
  "livestock.access.house.binaryNo Access"   = "W5: Livestock access inside: No",
  "animal.excrement.floor.binaryNot Exposed" = "W6: Animal excrement on floor: No"
)

coef_df <- as.data.frame(cox_summary$coefficients)
coef_df$term <- rownames(cox_summary$coefficients)
df_forest <- coef_df %>%
  filter(term %in% names(pretty_labels)) %>%
  rename(estimate = coef, std.error = `se(coef)`) %>%
  mutate(
    HR = exp(estimate),
    HR_95_low = exp(estimate - 1.96 * std.error),
    HR_95_high = exp(estimate + 1.96 * std.error),
    term = factor(term, levels = names(pretty_labels)),
    label = pretty_labels[term]
  )

# Smaller font for better balance
base_size <- 11
p <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
  geom_errorbarh(aes(xmin = HR_95_low, xmax = HR_95_high), height = 0.2) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  labs(x = "Hazard Ratio", y = NULL) +
  theme_cowplot() +
  theme(
    axis.text = element_text(size = base_size - 2),
    axis.title = element_text(size = base_size + 2),
    legend.position = "none"
  )
print(p)


# Fit Cox proportional hazards model (no age, no sex, no seasonality; midpoint timing for events)
fit_ct_midpoint <- coxme(
  Surv(start_date, end_date, event) ~ 
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) + (1 | HouseID),
  data = d_long
)

# Summarize model results
cox_summary <- summary(fit_ct_midpoint)
print(cox_summary)
print(exp(cox_summary$coefficients))  # Exponentiated coefficients (hazard ratios)

# Prepare data for forest plot
pretty_labels <- c(
  "main.drinking.water.dry.binaryImproved"   = " Drinking water (dry): Improved",
  "main.drinking.water.rainy.binaryImproved" = " Drinking water (rainy): Improved",
  "cleaning.water.storage.binaryTreated"     = " Stored water: Treated",
  "correct.handwashing.binaryCorrect"        = " Handwashing: Correct",
  "improved.sanitation.binaryImproved"       = " Sanitation: Improved",
  "livestock.access.house.binaryNo Access"   = " Livestock access inside: No",
  "animal.excrement.floor.binaryNot Exposed" = " Animal excrement on floor: No"
)
ordered_terms <- names(pretty_labels)

coef_df <- as.data.frame(cox_summary$coefficients)
coef_df$term <- rownames(cox_summary$coefficients)
df_forest <- coef_df %>%
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
    label      = pretty_labels[term]
  )

base_size <- 10
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
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  scale_fill_manual(values = c("50% CI" = "#7ca982", "80% CI" = "#cde5b2", "95% CI" = "#eeeec8")) +
  labs(x = "Hazard Ratio", y = NULL) +
  theme_cowplot() +
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.title         = element_blank(),
    legend.text          = element_text(size = base_size - 2),
    axis.text            = element_text(size = base_size),
    axis.title           = element_text(size = base_size + 2),
    panel.grid.major     = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor     = element_line(linetype = "dotted", color = "grey90")
  )
print(p)
##########################

# Fit Cox proportional hazards model (no age, no sex, no seasonality; midpoint timing for events)
fit_ct_midpoint <- coxme(
  Surv(start_date, end_date, event) ~ 
    main.drinking.water.dry.binary +
    main.drinking.water.rainy.binary +
    cleaning.water.storage.binary +
    correct.handwashing.binary +
    improved.sanitation.binary +
    livestock.access.house.binary +
    animal.excrement.floor.binary +
    (1 | menage_id_member) + (1 | HouseID),
  data = d_long
)

# Summarize and build coefficient dataframe
cox_summary <- summary(fit_ct_midpoint)
coef_df <- as.data.frame(cox_summary$coefficients)
coef_df$term <- rownames(cox_summary$coefficients)

# Pretty labels mapping
pretty_labels <- c(
  "main.drinking.water.dry.binaryImproved"   = " Drinking water (dry): Improved",
  "main.drinking.water.rainy.binaryImproved" = " Drinking water (rainy): Improved",
  "cleaning.water.storage.binaryTreated"     = " Stored water: Treated",
  "correct.handwashing.binaryCorrect"        = " Handwashing: Correct",
  "improved.sanitation.binaryImproved"       = " Sanitation: Improved",
  "livestock.access.house.binaryNo Access"   = " Livestock access inside: No",
  "animal.excrement.floor.binaryNot Exposed" = " Animal excrement on floor: No"
)

ordered_terms <- names(pretty_labels)

# Prepare data for forest plot
df_forest <- coef_df %>%
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
    label      = pretty_labels[term],
    
    # Create formatted hazard ratio text
    HR_label = sprintf("HR = %.2f (%.2f–%.2f)", HR, HR_95_low, HR_95_high)
    
  )

# --------- PLOT ---------

base_size <- 10

p <- ggplot(df_forest, aes(x = HR, y = fct_rev(label))) +
  
  # 95% CI bar
  geom_rect(aes(xmin = HR_95_low, xmax = HR_95_high,
                ymin = as.numeric(fct_rev(label)) - 0.3,
                ymax = as.numeric(fct_rev(label)) + 0.3,
                fill = "95% CI"), alpha = 0.4) +
  
  # 80% CI bar
  geom_rect(aes(xmin = HR_80_low, xmax = HR_80_high,
                ymin = as.numeric(fct_rev(label)) - 0.2,
                ymax = as.numeric(fct_rev(label)) + 0.2,
                fill = "80% CI"), alpha = 0.6) +
  
  # 50% CI bar
  geom_rect(aes(xmin = HR_50_low, xmax = HR_50_high,
                ymin = as.numeric(fct_rev(label)) - 0.1,
                ymax = as.numeric(fct_rev(label)) + 0.1,
                fill = "50% CI"), alpha = 0.8) +
  
  # Point for HR
  geom_point(size = 3, color = "black") +
  
  # Vertical reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  
  # ---- ADD HAZARD LABELS HERE ----
geom_text(
  aes(x = HR_95_high, y = fct_rev(label), label = HR_label),
  hjust = -0.2, size = 3.5, color = "black"
) +
  
  # Axis scale
  scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  
  # Colors for CI shading
  scale_fill_manual(values = c("50% CI" = "#7ca982",
                               "80% CI" = "#cde5b2",
                               "95% CI" = "#eeeec8")) +
  
  labs(x = "Hazard Ratio", y = NULL) +
  theme_cowplot() +
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.title         = element_blank(),
    legend.text          = element_text(size = base_size - 2),
    axis.text            = element_text(size = base_size),
    axis.title           = element_text(size = base_size + 2),
    panel.grid.major     = element_line(linetype = "dotted", color = "grey70"),
    panel.grid.minor     = element_line(linetype = "dotted", color = "grey90")
  )

print(p)

