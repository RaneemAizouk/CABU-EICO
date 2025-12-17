########################################################################################
# THIS CODE ESTIMATES CHANGE IN WASH 
########################################################################################
# Author: Esther van Kleef
# Date: 4 October 2025
# Last updated: 14 October 2025

rm(list = ls())

pacman::p_load(survey, dplyr, tidyr, AER, stringr, broom, readr, purrr, tibble,forcats, cowplot, metafor, lubridate)

# -------------------------------------------------------------------------------
# Directory and load in data
# -------------------------------------------------------------------------------
DirectoryData <- "./Data/BF/clean"

# HH WaSH (all HH) and stool-HH subset
wash        <- read.csv(file.path(DirectoryData, "FINAL_FOR_SHARING/Household_WASH_BF.csv"))
wash_stool  <- read.csv(file.path(DirectoryData, "FINAL_FOR_SHARING/Household_stool_WASH_BF.csv"))

# Provider × village antibiotic use
watch <- read.csv(file.path(DirectoryData, "use_in_analyses/watch_acute_offset.csv"))

# MSM village ESBL results (calendar-time adjusted) — long with village, Group, phase_label, median, var, (lo/hi, n_draws)
msm_draws <- read.csv("./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario1/inc_s1.csv")

# Villages (that are the clusters) of CABU-EICO
villages = readxl::read_xlsx(paste0("./Data/BF/Raw/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","ajouter")

village_size = readxl::read_xlsx(paste0("./Data/BF/Raw/Sélection et randomisation grappes CABU_B Nanoro.xlsx"), sheet=1)
names(village_size)
village_size = village_size %>% select(c(`Nom des villages`, ménages, `population active (Mar 2023 update)`))%>%
  rename(
    village_name = `Nom des villages`,
    n.households =  ménages,
    n.population = `population active (Mar 2023 update)`
  ) %>%
  filter(!is.na(village_name))
head(village_size)

villages = left_join(villages, village_size, by="village_name")
villages$village_name = tolower(villages$village_name)

# Inform the number of households based on average household size
m_size = median(villages$n.population[villages$village_name!="pella"]/villages$n.households[villages$village_name!="pella"])

wash$n.households[wash$village_name=="pella"] = wash$n.population[wash$village_name=="pella"]/m_size

#wash = left_join(wash, villages%>%select(c(village_name, n.population, n.households)))

# -------------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------------
pick_col <- function(df, candidates, default = NULL) {
  nm <- candidates[candidates %in% names(df)][1]
  if (is.na(nm) || is.null(nm)) default else nm
}

to_arm <- function(x) {
  case_when(
    str_detect(tolower(x), "inter") ~ "intervention",
    str_detect(tolower(x), "contr") ~ "control",
    TRUE ~ as.character(x)
  )
}

# -------------------------------------------------------------------------------
# HH-level WaSH indicators (good = 1), harmonize labels
# -------------------------------------------------------------------------------
prep_wash <- function(df) {
  round_nm <- pick_col(df, c("redcap_event_name","round"))
  arm_nm   <- pick_col(df, c("intervention.text","arm","intervention","intervention_text"))
  vill_nm  <- pick_col(df, c("village_name","village","village.cluster","VillageID"))
  date <- pick_col(df, c("date.enquete"))
  
  out <- df %>%
    transmute(
      menage_id = .data[["menage_id"]],
      n.householdmember = .data[["n.householdmember"]],
      n.population = .data[["n.population"]],
      n.households = .data[["n.households"]],
      village   = .data[[vill_nm]],
      arm_raw   = .data[[arm_nm]],
      round_raw = .data[[round_nm]],
      date = .data[[date]],
      date = as.Date(date, format = "%Y-%m-%d"),
      month = month(date),
      rainy = ifelse(month %in% c(6:11), 1, 0),
      main.drinking.water.dry.binary    = case_when(main.drinking.water.dry.binary == "Improved" ~ 0, 
                              main.drinking.water.dry.binary == "Unimproved" ~ 1,
                              TRUE ~ NA_real_),
      main.drinking.water.rainy.binary  = case_when(main.drinking.water.rainy.binary == "Improved" ~ 0, 
                              main.drinking.water.rainy.binary == "Unimproved" ~ 1,
                              TRUE ~ NA_real_),
      no.cleaning.water.storage.binary = case_when(cleaning.water.storage.binary== "Yes" ~ 0, 
                                cleaning.water.storage.binary == "No" ~ 1,
                                TRUE ~ NA_real_),
      no.correct.handwashing.binary  = case_when(correct.handwashing.binary == "Yes" ~ 0,
                                correct.handwashing.binary == "No" ~ 1,
                                TRUE ~ NA_real_),
      no.improved.sanitation.binary = case_when(improved.sanitation.binary == "Yes" ~  0,
                                improved.sanitation.binary== "No" ~ 1,
                                TRUE ~ NA_real_),
      livestock.access.house.binary = case_when(livestock.access.house.binary == "Yes" ~ 1,
                                        livestock.access.house.binary == "No" ~ 0,
                                        TRUE ~ NA_real_),
      animal.excrement.floor.binary  = case_when(animal.excrement.floor.binary == "Yes" ~ 1,
                                               animal.excrement.floor.binary == "No" ~ 0,
                                               TRUE ~ NA_real_)
    ) %>%
    mutate(
      arm = factor(to_arm(arm_raw), levels = c("control","intervention")),
      round = factor(case_when(
        round_raw == "round_0_arm_1" ~ "baseline",
        round_raw == "round_3_arm_1" ~ "post",
        TRUE ~ as.character(round_raw)
      ), levels = c("baseline","post")),
      village = tolower(village)
    )
  out
}

# While waiting for sampling date of second round, impute the season with data from the first round (as is a year later)
wash_hh_all   <- prep_wash(wash)
wash_hh_all$rainy = ifelse(is.na(wash_hh_all$date), NA, wash_hh_all$rainy) # This as date is missing for all households in the 2nd round, checked with Franck, waiting for reply

wash_hh_all <- wash_hh_all %>%
  group_by(menage_id) %>%
  mutate(rainy = first(na.omit(rainy))) %>%
  ungroup()


# Select only those where stool samples were taken
#wash_hh_stool <- prep_wash(wash_stool)
wash_hh_stool <- wash_hh_all %>% filter(menage_id %in% unique(wash_stool$menage_id))
length(unique(wash_hh_stool$menage_id[wash_hh_stool$round=="baseline"])) #266
length(unique(wash_hh_stool$menage_id[wash_hh_stool$round=="post"])) # 258

# Ensure only one observation per hh (currently at individual level for wash + stool)
indicators <- c(
  "main.drinking.water.dry.binary","main.drinking.water.rainy.binary",#"no.cleaning.water.storage.binary",
  "no.correct.handwashing.binary","no.improved.sanitation.binary",
  "livestock.access.house.binary","animal.excrement.floor.binary"
)

length(unique(wash_hh_all$menage_id))
length(unique(wash_hh_all$menage_id[wash_hh_all$round=="post"]))

# -------------------------------------------------------------------------------
# HH weights for SURVEY (population / sampled HHs in each village-round)
# -------------------------------------------------------------------------------

# Compute population-based survey weights:
# Each household gets a weight = village population / number of sampled households
# so that results are more representative of the true population.
# Larger villages thus contribute proportionally more, even if all have equal sample sizes.

# Extract the total population per village from the WASH dataset
pop_nm_wash  <- pick_col(wash,  c("n.population"))

# Create a reference table: one row per village with its total population
pop_ref <- wash_hh_all %>% select(c(village, n.population)) %>% unique()

# Number of sampled households per village per round
n_hh_vr <- wash_hh_all %>% count(village, round, name = "n_hh_round")

# Combine population and sample counts to compute survey weights
weights_hh <- n_hh_vr %>%
  left_join(pop_ref, by = "village") %>%
  mutate(weight = ifelse(is.na(n.population), 1, n.population / pmax(n_hh_round, 1)))

# Merge weights back into the main household- and stool-level datasets
wash_hh_all   <- wash_hh_all   %>% left_join(weights_hh)
wash_hh_stool <- wash_hh_stool %>% left_join(weights_hh %>% select(village, round, weight))

# Normalise the weights
norm_w <- function(w) w / mean(w, na.rm = TRUE)
wash_hh_all$weight   <- ifelse(is.na(wash_hh_all$weight),   1, norm_w(wash_hh_all$weight))
wash_hh_stool$weight <- ifelse(is.na(wash_hh_stool$weight), 1, norm_w(wash_hh_stool$weight))

# CHECK DIFFERENCE IN RAINY VS DRY 
#--------------------------------------------------------------------------------

table(wash_hh_all$rainy, wash_hh_all$round, useNA="always") # The second round is not in the rainy season at all
table(wash_hh_all$rainy, wash_hh_all$arm, useNA="always") # 249/(249+160) = 61% of control group in rainy season; 189/(189+210) this was 47% for the intervention group

table(wash_hh_all$rainy, wash_hh_all$round, wash_hh_all$animal.excrement.floor.binary)

# -------------------------------------------------------------------------------
# SURVEY-WEIGHTED WASH change (ALL HH): arm × round prevalences + DiD PR
# -------------------------------------------------------------------------------

# In the below, I have used a quasipoisson as the binomial model did not converge for all variables
# I have followed what is described here https://doi.org/10.1093/aje/kwh090


# Prepare long dataset
# -------------------------------------------------------------------

# !!! USE ALL DATA (wash_hh_all) OR ONLY OF THOSE HOUSEHOLDS WHERE STOOL IS COLLECTED (wash_hh_stool) !!!

# Below is with all hh data, have for now run it twice, changing the dataset to wash_hh_stool and store seperately to get the two tables stored
wash_long_all <- wash_hh_stool%>%
  pivot_longer(all_of(indicators), names_to = "indicator", values_to = "y") %>%
  mutate(
    y = as.numeric(y),
    indicator = factor(indicator),
    arm = factor(arm, levels = c("control","intervention")),
    round = factor(round, levels = c("baseline","post"))
  )
 
inds <- sort(unique(as.character(wash_long_all$indicator)))

# Function to fit model and extract summary table for each indicator
# -------------------------------------------------------------------


fit_one_indicator <- function(ind) {
  print(ind)
  d <- wash_long_all %>% filter(indicator == ind) 
  # Fit survey-weighted log-binomial model
  des <- svydesign(
    id = ~village,
    weights = ~weight,
    data = d,
    nest = TRUE
  )
  
  # Regression model accounting for survey design and seasonality
  fit <- svyglm(y ~ arm * round + rainy, # + rainy, can not account for rainy season as is date is mssing in the last round
                design = des,
                family = quasipoisson(link = "log"))
  
  # Predicted adjusted prevalences (for rainy season)
  # adjusted prevalences (set rainy = 0, this as most surveys were taken in the rainy season so most representative)
  newdat <- expand.grid(
    arm   = levels(d$arm),
    round = levels(d$round),
    rainy = 0
  )
  pred <- as.data.frame(predict(fit, newdata = newdat, type = "response", se.fit = TRUE))
  
  newdat <- newdat %>%
    mutate(
      estimate  = pred$response,
      conf.low  = pmax(0, pred$response - 1.96 * pred$SE),
      conf.high = pred$response + 1.96 * pred$SE,
      prev_lab  = sprintf("%.1f (%.1f–%.1f)", estimate*100, conf.low*100, conf.high*100)
    ) %>%
    select(arm, round, prev_lab) %>%
    pivot_wider(names_from = c(round, arm), values_from = prev_lab)
  
  # DiD PR (interaction term)
  did_effect <- tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "armintervention:roundpost") %>%
    mutate(DiD_PR = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high)) %>%
    pull(DiD_PR)
  did_effect
  
  # Baseline-adjusted PR (post only; intervention vs control), delta method on link scale
  newdat_post <- expand.grid(
    arm   = c("control","intervention"),
    round = "post",
    rainy = 0
  )
  # Factors with the same levels as used in the model
  newdat_post$arm   <- factor(newdat_post$arm,   levels = levels(d$arm))
  newdat_post$round <- factor(newdat_post$round, levels = levels(d$round))
  
  # linear predictor and SE
  pred_post <- as.data.frame(predict(fit, newdata = newdat_post, type = "link", se.fit = TRUE))
  
  X <- model.matrix(~ arm * round + rainy, newdat_post)
  #X <- model.matrix(~ arm * round, newdat_post)
  
  # Concluded that the post only PR is not meaningfull in our context, as the baseline differences between the groups are quite high, so report the DiD only
  # log(PR) = eta_int - eta_ctrl ; Var = (x2-x1)' V (x2-x1)
  var_est <- vcov(fit)
  diff_x <- X[2,] - X[1,]
  log_pr <- sum(diff_x * coef(fit))
  se_log_pr <- as.numeric(sqrt(t(diff_x) %*% var_est %*% diff_x))
  pr_post_only <- exp(log_pr)
  ci_low  <- exp(log_pr - 1.96 * se_log_pr)
  ci_high <- exp(log_pr + 1.96 * se_log_pr)
  post_only_pr <- sprintf("%.2f (%.2f–%.2f)", pr_post_only, ci_low, ci_high)
  
  newdat %>%
    mutate(
      indicator = ind,
      `Post-only PR` = post_only_pr,
      `Difference-in-differences PR` = did_effect
    ) %>%
    select(
      indicator,
      Baseline_Control       = baseline_control,
      Baseline_Intervention  = baseline_intervention,
      Post_Control           = post_control,
      Post_Intervention      = post_intervention,
      `Post-only PR`,
      `Difference-in-differences PR`
    )
}

# -------------------------------------------------------------------
# Run for all indicators
# -------------------------------------------------------------------

results_table <- map_dfr(inds, fit_one_indicator)

# -------------------------------------------------------------------
# Formatted output
# -------------------------------------------------------------------

results_table <- results_table %>%
  rename(`Practice/condition` = indicator,
         `Prevalence ratio (Post-only)` =  `Post-only PR`,
         `Prevalence ratio (did)` =  `Difference-in-differences PR`) %>%
  mutate(`Practice/condition` = forcats::fct_recode(`Practice/condition`,
                                                    "Use of unimproved drinking water source (dry)" = "main.drinking.water.dry.binary",
                                                    "Use of unimproved drinking water source (rainy)" = "main.drinking.water.rainy.binary",
                                                    #"No cleaning of drinking water storage containers" = "no.cleaning.water.storage.binary",
                                                    "Incorrect handwashing" = "no.correct.handwashing.binary",
                                                    "Use of unimproved sanitary facility" = "no.improved.sanitation.binary",
                                                    "Livestock animals access the house" = "livestock.access.house.binary",
                                                    "Animal excrement on the house floor" = "animal.excrement.floor.binary"
  )) 

results_table = results_table[c(3,4,5,6,2,1),]

# View final table
results_table

# STORE OUTPUT FOR SI
write_xlsx(results_table, path = "./Output/Figures_and_tables/Paper/Final/Table2_change_wash_unformatted.xlsx")

# STORE OUTPUT FOR SI - STOOL ONLY
write_xlsx(results_table, path = "./Output/Figures_and_tables/Paper/Final/Table2_change_wash_STOOL_hh_unformatted.xlsx")

# Descriptive table WASH
#table(wash_hh_all$no_animal_excrement_on_floor[wash_hh_all$round=="baseline"],wash_hh_all$arm[wash_hh_all$round=="baseline"])
