########################################################################################
# THIS CODE ASSOCIATES CHANGE IN VILLAGE LEVEL ABX AND WASH WITH ESBL-E ACQUISITION
########################################################################################
# Author: Esther van Kleef
# Date: 4 October 2025

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
      village   = .data[[vill_nm]],
      arm_raw   = .data[[arm_nm]],
      round_raw = .data[[round_nm]],
      date = .data[[date]],
      date = as.Date(date, format = "%Y-%m-%d"),
      month = month(date),
      rainy = ifelse(month %in% c(6:11), 1, 0),
      safe_dry    = case_when(main.drinking.water.dry.binary == "Improved" ~ 1, 
                              main.drinking.water.dry.binary == "Unimproved" ~ 0,
                              TRUE ~ NA_real_),
      safe_rainy  = case_when(main.drinking.water.rainy.binary == "Improved" ~ 1, 
                              main.drinking.water.rainy.binary == "Unimproved" ~ 0,
                              TRUE ~ NA_real_),
      clean_storage = case_when(cleaning.water.storage.binary== "Yes" ~ 1, 
                                cleaning.water.storage.binary == "No" ~ 0,
                                TRUE ~ NA_real_),
      correct_hw    = case_when(correct.handwashing.binary == "Yes" ~ 1,
                                correct.handwashing.binary == "No" ~ 0,
                                TRUE ~ NA_real_),
      improved_san  = case_when(improved.sanitation.binary == "Yes" ~  1,
                                improved.sanitation.binary== "No" ~ 0,
                                TRUE ~ NA_real_),
      no_livestock_in_house = case_when(livestock.access.house.binary == "Yes" ~ 1,
                                        livestock.access.house.binary == "No" ~ 0,
                                        TRUE ~ NA_real_),
      no_animal_excrement_on_floor = case_when(animal.excrement.floor.binary == "Yes" ~ 1,
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

wash_hh_all   <- prep_wash(wash)
wash_hh_stool <- prep_wash(wash_stool)

indicators <- c(
  "safe_dry","safe_rainy","clean_storage","correct_hw","improved_san",
  "no_livestock_in_house","no_animal_excrement_on_floor"
)

# -------------------------------------------------------------------------------
# HH weights for SURVEY (population / sampled HHs in each village-round)
# -------------------------------------------------------------------------------
pop_nm_wash  <- pick_col(wash,  c("n.population"))
pop_ref <- wash_hh_all %>% select(c(village, n.population)) %>% unique()

# Number of households per village per round
n_hh_vr <- wash_hh_all %>% count(village, round, name = "n_hh_round")

weights_hh <- n_hh_vr %>%
  left_join(pop_ref, by = "village") %>%
  mutate(weight = ifelse(is.na(n.population), 1, n.population / pmax(n_hh_round, 1)))

wash_hh_all   <- wash_hh_all   %>% left_join(weights_hh)
wash_hh_stool <- wash_hh_stool %>% left_join(weights_hh %>% select(village, round, weight))

norm_w <- function(w) w / mean(w, na.rm = TRUE)
wash_hh_all$weight   <- ifelse(is.na(wash_hh_all$weight),   1, norm_w(wash_hh_all$weight))
wash_hh_stool$weight <- ifelse(is.na(wash_hh_stool$weight), 1, norm_w(wash_hh_stool$weight))

# -------------------------------------------------------------------------------
# SURVEY-WEIGHTED WaSH change (ALL HH): arm×round prevalences + DiD PR
# -------------------------------------------------------------------------------

# In the below, I have used a quasipoisson as the binomial model did not converge for all variables
# I have followed what is described here https://doi.org/10.1093/aje/kwh090
# 

# Prepare long dataset
# -------------------------------------------------------------------

# USE ALL DATA OR ONLY OF THOSE HOUSEHOLDS WHERE STOOL IS COLLECTED
wash_long_all <- wash_hh_stool %>%
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
    id = ~village + menage_id,
    weights = ~weight,
    data = d,
    nest = TRUE
  )
  
  # Regression model accounting for survey design and seasonality
  fit <- svyglm(y ~ arm * round + rainy,
                design = des,
                family = quasipoisson(link = "log"))
  
  # Predicted adjusted prevalences (for rainy season)
  # adjusted prevalences (set rainy = 1; keep your predict() syntax)
  newdat <- expand.grid(
    arm   = levels(d$arm),
    round = levels(d$round),
    rainy = 1
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
    rainy = 1
  )
  # make them factors with the SAME levels as used in the model
  newdat_post$arm   <- factor(newdat_post$arm,   levels = levels(d$arm))
  newdat_post$round <- factor(newdat_post$round, levels = levels(d$round))
  
  # linear predictor & SE
  pred_post <- as.data.frame(predict(fit, newdata = newdat_post, type = "link", se.fit = TRUE))
  
  X <- model.matrix(~ arm * round + rainy, newdat_post)
  
  # log(PR) = eta_int - eta_ctrl ; Var = (x2-x1)' V (x2-x1)
  var_est <- vcov(fit)
  diff_x <- X[2,] - X[1,]
  log_pr <- sum(diff_x * coef(fit))
  se_log_pr <- as.numeric(sqrt(t(diff_x) %*% var_est %*% diff_x))
  pr_baseline_adj <- exp(log_pr)
  ci_low  <- exp(log_pr - 1.96 * se_log_pr)
  ci_high <- exp(log_pr + 1.96 * se_log_pr)
  baseline_adj_pr <- sprintf("%.2f (%.2f–%.2f)", pr_baseline_adj, ci_low, ci_high)
  
  newdat %>%
    mutate(
      indicator = ind,
      `Baseline-adjusted PR` = baseline_adj_pr,
      `Difference-in-differences PR` = did_effect
    ) %>%
    select(
      indicator,
      Baseline_Control       = baseline_control,
      Baseline_Intervention  = baseline_intervention,
      Post_Control           = post_control,
      Post_Intervention      = post_intervention,
      `Baseline-adjusted PR`,
      `Difference-in-differences PR`
    )
}

# -------------------------------------------------------------------
# Run for all indicators
# -------------------------------------------------------------------

results_table <- purrr::map_dfr(inds, fit_one_indicator)

# -------------------------------------------------------------------
# Formatted output
# -------------------------------------------------------------------

results_table <- results_table %>%
  rename(`Practice/condition` = indicator,
         `Prevalence ratio (baseline-adjusted)` =  `Baseline-adjusted PR`,
         `Prevalence ratio (did)` =  `Difference-in-differences PR`) %>%
  mutate(`Practice/condition` = forcats::fct_recode(`Practice/condition`,
                                                    "Access to safe drinking water (dry season)" = "safe_dry",
                                                    "Access to safe drinking water (rainy season)" = "safe_rainy",
                                                    "Cleaning drinking water storage containers before reuse" = "clean_storage",
                                                    "Correct handwashing" = "correct_hw",
                                                    "Using improved sanitary facility" = "improved_san",
                                                    "Livestock animals access the house" = "no_livestock_in_house",
                                                    "Animal excrement on the house floor" = "no_animal_excrement_on_floor"
  )) 

results_table = results_table[c(6,7,1,2,3,4,5),]

# View final table
results_table

# STORE OUTPUT FOR SI
writexl::write_xlsx(results_table, path = "./Output/Figures_and_tables/Paper/TableS3_change_wash.xlsx")

# STORE OUTPUT FOR SI - STOOL ONLY
writexl::write_xlsx(results_table, path = "./Output/Figures_and_tables/Paper/TableS3_change_wash_STOOL_hh.xlsx")
