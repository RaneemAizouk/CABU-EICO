
###############################################################################
# Data_prep_coxme.R
# Data preparation for Cox mixed-effects modeling (coxme).
# This script prepares longitudinal cohort data for Cox mixed-effects
# analysis of WASH indicators and ESBL-E acquisition, including:
#    cleaning and harmonisation of covariates
#    definition of ESBL-E acquisition events (0 → 1 transitions)
#    construction of start–stop time intervals
#    midpoint imputation of event times
#    anonymisation for data sharing
#
# Outputs:
#    An anonymised CSV for sharing (no original IDs)
#    An internal modeling RDS (contains original IDs)
#   Private mapping files linking original IDs to anonymised IDs
#
# Related analysis script:
#    5_CoxModel_Wash.R (fits the Cox mixed-effects models)
# Author: adapted from Raneem Aizouk
# Created: April 2025
# Updated: Dec 2025

###############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(lubridate)
  library(zoo)
  library(tidyr)
  library(stringr)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "./Data/BF/clean/use_in_analyses/bf_esbl0123_long_all.rda",
              help = "Path to .rda containing object dfls0 (default ./Data/BF/clean/use_in_analyses/bf_esbl0123_long_all.rda)"),
  make_option(c("-o", "--out"), type = "character",
              default = "./Data/BF/clean/FINAL_FOR_SHARING/",
              help = "Output directory for shareable outputs (default ./Data/BF/clean/FINAL_FOR_SHARING/)"),
  make_option(c("-p", "--private_out"), type = "character",
              default = "../private_mappings/",
              help = "Private output dir for sensitive mapping files (default ../private_mappings/). MUST be outside repo."),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Verbose output")
)
opt <- parse_args(OptionParser(option_list = option_list))

INPUT_RDA         <- opt$input
OUTPUT_DIR        <- opt$out
PRIVATE_OUTPUT_DIR<- opt$private_out
VERBOSE           <- opt$verbose

msg <- function(...) if (VERBOSE) cat(..., "\n")

# Create folders
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PRIVATE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Safety warning: if private dir is a subdir of OUTPUT_DIR, warn the user
normalize_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
if (startsWith(normalize_path(PRIVATE_OUTPUT_DIR), normalize_path(OUTPUT_DIR))) {
  warning("PRIVATE_OUTPUT_DIR is inside OUTPUT_DIR. This increases risk of accidentally committing mapping files. Consider setting --private_out to a path outside the repo.")
  msg("PRIVATE_OUTPUT_DIR currently:", PRIVATE_OUTPUT_DIR)
}

# ---- Load input ------------------------------------------------------------
if (!file.exists(INPUT_RDA)) stop("Input file not found: ", INPUT_RDA)
load(INPUT_RDA)    # expects object dfls0
if (!exists("dfls0")) stop("Object `dfls0` not found inside the input .rda")
data <- dfls0
rm(dfls0)
msg("Loaded input. Rows:", nrow(data), "Unique persons:", n_distinct(data$menage_id_member))

# ---- Basic sanity checks ---------------------------------------------------
required_vars <- c("menage_id_member", "menage_id", "redcap_event_name", "date.use")
miss_req <- setdiff(required_vars, names(data))
if (length(miss_req) > 0) stop("Missing required variables in dfls0: ", paste(miss_req, collapse = ", "))

# Ensure date.use is Date
data <- data %>% mutate(date.use = as.Date(date.use))

# ---- Normalise outcome and sex variables (canonical names) -----------------
# Create a canonical `esble` variable
if ("esble" %in% names(data)) {
  data$esble <- data$esble
} else if ("esbl_pos" %in% names(data)) {
  data$esble <- data$esbl_pos
} else {
  stop("No ESBL status variable found (esble or esbl_pos).")
}
# Normalise sexe to integer: Female -> 1, else 0
data$sexe <- as.character(data$sexe)
data$sexe <- ifelse(tolower(data$sexe) %in% c("female", "f"), 1L, 0L)

# Ensure ID columns have sensible types
data$menage_id_member <- as.character(data$menage_id_member)
data$menage_id <- as.character(data$menage_id)

# ---- Derive round index (same convention as analysis) ----------------------
data <- data %>%
  mutate(round = as.integer(gsub("round_(\\d+)_arm_1", "\\1", redcap_event_name)) + 1)

# ---- Remove individuals with completely missing age or sex -----------------
bad_ids <- data %>%
  group_by(menage_id_member) %>%
  summarise(miss_age = all(is.na(age)), miss_sex = all(is.na(sexe)), .groups = "drop") %>%
  filter(miss_age | miss_sex) %>%
  pull(menage_id_member)

msg("Dropping individuals with all-missing age or sex:", length(bad_ids))
data <- data %>% filter(!menage_id_member %in% bad_ids)

# Fill remaining missing age/sex within individuals (assume constant over time)
data <- data %>%
  group_by(menage_id_member) %>%
  fill(age, sexe, .direction = "downup") %>%
  ungroup()

# ---- IDs and coding used in analysis -------------------------------------
# Numeric household ID for random effect (but keep string menage_id available)
data$HouseID <- as.numeric(factor(data$menage_id))

# ---- Calendar time variables (days since global origin) --------------------
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

# ---- Define acquisition outcome ACQ (0 -> 1 transitions), NA-safe -----------
# prefer variable `esble` (created above)
data <- data %>% arrange(menage_id_member, round)
data$ACQ <- NA_integer_
for (i in 2:nrow(data)) {
  if (!is.na(data$menage_id_member[i]) && data$menage_id_member[i] == data$menage_id_member[i - 1]) {
    prev <- data$esble[i - 1]
    cur  <- data$esble[i]
    if (!is.na(prev) && !is.na(cur)) {
      data$ACQ[i] <- ifelse(prev == 0 & cur == 1, 1L, 0L)
    } else {
      data$ACQ[i] <- NA_integer_
    }
  }
}
msg("ACQ transitions computed (NA-safe).")

# ---- WASH harmonisation: LOCF then fromLast LOCF ----------------------------
wash_vars <- c(
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "correct.handwashing.binary",
  "improved.sanitation.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)
wash_vars <- intersect(wash_vars, names(data))
if (length(wash_vars) == 0) warning("No WASH variables found in data; check names.")

data <- data %>%
  group_by(menage_id_member) %>%
  arrange(round) %>%
  mutate(across(all_of(wash_vars), ~ zoo::na.locf(.x, na.rm = FALSE))) %>%
  mutate(across(all_of(wash_vars), ~ zoo::na.locf(.x, fromLast = TRUE, na.rm = FALSE))) %>%
  ungroup()

# ---- Filter to rows with complete WASH (primary analysis) ------------------
data_model <- data %>%
  filter(if_all(all_of(wash_vars), ~ !is.na(.)))
msg("Rows after dropping missing WASH:", nrow(data_model), "Unique persons:", n_distinct(data_model$menage_id_member))

# ---- Build long dataset for Cox (midpoint timing) -------------------------
d_long <- data_model %>%
  arrange(menage_id_member, date.use) %>%
  group_by(menage_id_member) %>%
  mutate(
    start_date = cal_start,
    end_date   = cal_stop,
    event      = ACQ
  ) %>%
  ungroup() %>%
  filter(!is.na(event) & end_date > start_date) %>%
  mutate(
    event_time = ifelse(event == 1, (start_date + end_date) / 2, end_date),
    event = as.integer(event)
  )

msg("d_long prepared. Rows:", nrow(d_long), "Unique persons:", n_distinct(d_long$menage_id_member))

# ---- Anonymisation: create HouseAnon / PersonAnon (save mapping off-repo) ---
anon_house <- d_long %>% distinct(menage_id) %>% arrange(menage_id) %>%
  mutate(HouseAnon = sprintf("H%03d", row_number()))
anon_person <- d_long %>% distinct(menage_id_member) %>% arrange(menage_id_member) %>%
  mutate(PersonAnon = sprintf("P%04d", row_number()))

# Merge anonymised IDs into modeling df
d_long <- d_long %>%
  left_join(anon_house, by = "menage_id") %>%
  left_join(anon_person, by = "menage_id_member")

# ---- Save mapping files to PRIVATE_OUTPUT_DIR (DO NOT PUBLISH) -------------
mapping_person_path <- file.path(PRIVATE_OUTPUT_DIR, "mapping_personid_to_anon__DO_NOT_PUBLISH.csv")
mapping_house_path  <- file.path(PRIVATE_OUTPUT_DIR, "mapping_houseid_to_anon__DO_NOT_PUBLISH.csv")
write.csv(anon_person, mapping_person_path, row.names = FALSE)
write.csv(anon_house, mapping_house_path, row.names = FALSE)
msg("Wrote mapping files to PRIVATE_OUTPUT_DIR (do NOT commit/push these):", mapping_person_path, mapping_house_path)

# ---- Write anonymised minimal CSV for sharing (no original IDs) ------------
cox_anonymised <- d_long %>%
  select(PersonAnon, HouseAnon, round, date.use, start_date, end_date, event, event_time, age, sexe, all_of(wash_vars)) %>%
  # Safety: drop any accidentally included original ID columns (just in case)
  select(-any_of(c("menage_id", "menage_id_member", "cs.id.individu", "record_id")))

anon_out_path <- file.path(OUTPUT_DIR, "cox_wash_anonymised_minimal.csv")
write.csv(cox_anonymised, anon_out_path, row.names = FALSE)
msg("Wrote anonymised minimal CSV (shareable):", anon_out_path)

# ---- Save full modeling dataset (contains original IDs; keep local/private) -
model_rds_path <- file.path(OUTPUT_DIR, "cox_d_long_for_modeling.rds")
saveRDS(d_long, model_rds_path)
msg("Saved full d_long for modeling (RDS):", model_rds_path)

# ---- Save session info for provenance -------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo_data_prep_coxme.txt"))

# ---- Final summary --------------------------------------------------------
msg("Data preparation complete.")
msg("Unique individuals in d_long:", n_distinct(d_long$menage_id_member))
msg("Rows in d_long:", nrow(d_long))
msg("Anonymised CSV rows:", nrow(cox_anonymised))
msg("Shareable outputs written to:", normalize_path(OUTPUT_DIR))
msg("Private mapping outputs written to:", normalize_path(PRIVATE_OUTPUT_DIR))

# ---- .gitignore snippet (add this to your repository's .gitignore) ---------
cat("\nIMPORTANT: Add the following to your repo .gitignore to avoid leaking\nprivate mapping files and raw data.\n\n")
cat("### .gitignore suggested entries (append to .gitignore) ###\n")
cat("/Data/\n")
cat("/Output/\n")
cat("/Data/BF/clean/FINAL_FOR_SHARING/\n")
cat(sprintf("%s\n", paste0(normalize_path(PRIVATE_OUTPUT_DIR), collapse = "")))
cat("*__DO_NOT_PUBLISH.csv\n")
cat("mapping_personid_to_anon__DO_NOT_PUBLISH.csv\n")
cat("mapping_houseid_to_anon__DO_NOT_PUBLISH.csv\n")
cat("### End .gitignore snippet ###\n\n")
