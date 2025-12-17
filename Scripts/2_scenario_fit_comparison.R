###############################################################################
# THIS CODE ALLOWS TO COMPARE SCENARIOS AND PLOT THEIR RESULTS
###############################################################################

# Date created: 1 October 2025
# Date last updated: 1 October 2025
# Author: Esther van Kleef

rm(list = ls())

# ----------------------------
# Packages
# ----------------------------
pacman::p_load(
  posterior, bayesplot, dplyr, ggplot2, tibble, tidyr, purrr, stringr,
  readr, forcats, cowplot, rstan, loo, flextable, gridExtra, grid, furrr, matrixStats
)

# ----------------------------
# USER SETTINGS
# ----------------------------
# Where the fitted stan objects live (as .rds)
dir_sim  <- "Simulated_data"
dir_obs  <- "Observed_data"      
dir_fits <- "./Output/Model_results"

# ----------------------------
# Scenario registry
# ----------------------------
# Define all scenarios to compare
# Columns match how the file names are composed: paste0(intervention, seasonality, intervention_effect, scen, data_source, ".rds")
scenarios <- tribble(
  ~data_source, ~intervention, ~seasonality, ~intervention_effect, ~scen,
  # observed
  "observed",   "Two_step_",   "sine_",      "",                   "seasonality",
  "observed",   "Two_step_",   "spline_",    "",                   "seasonality",
  "observed",   "Two_step_",   "sine_",      "NonAdd_NonCol_",      "seasonality",
  "observed",   "Two_step_",   "spline_",    "NonAdd_NonCol_",      "seasonality",
  "observed",   "One_step_",   "sine_",      "",                   "seasonality",
  "observed",   "One_step_",   "spline_",    "",                   "seasonality",
  "observed",   "One_step_",   "sine_",      "NonAdd_NonCol_",      "seasonality",
  "observed",   "One_step_",   "spline_",    "NonAdd_NonCol_",      "seasonality",
) %>%
  mutate(
    scenario_label = paste0(intervention, seasonality, intervention_effect, scen, " (", data_source, ")"),
    # build the filename stem used by your save calls
    file_stem = paste0(intervention, seasonality, intervention_effect, scen, data_source, ".rds"),
    fit_path  = file.path(dir_fits, ifelse(data_source == "simulated", dir_sim, dir_obs), file_stem)
  )


# -------------------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------------------
safe_read_rds <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    warning(sprintf("Could not read %s: %s", path, e$message))
    NULL
  })
}

# LOO/diag function
loo_from_stanfit <- function(fit, log_lik_param = "log_lik") {
  ll <- tryCatch(
    loo::extract_log_lik(fit, parameter_name = log_lik_param, merge_chains = FALSE),
    error = function(e) NULL
  )
  if (is.null(ll)) return(NULL)
  
  # drop obs that are NA for all draws/chains
  drop_idx <- which(vapply(seq_len(dim(ll)[3]),
                           function(n) all(is.na(ll[,,n])),
                           logical(1)))
  ll2 <- if (length(drop_idx)) ll[, , -drop_idx, drop = FALSE] else ll
  if (length(dim(ll2)) != 3 || dim(ll2)[3] == 0) return(NULL)
  
  tryCatch(
    withCallingHandlers(
      loo.array(ll2),
      warning = function(w) tryInvokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
}

draws_df <- function(fit) as_draws_df(fit)

# Summarise named parameters into median & 95% CrI on LOG and HR scales
sum_param <- function(dd, par_name) {
  ix <- grep(paste0("^", par_name, "(\\[|$)"), names(dd))
  if (length(ix) == 0) return(NULL)
  M <- as.matrix(dd[, ix, drop = FALSE])
  qs <- apply(M, 2, function(x) stats::quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE))
  tibble(
    variable = colnames(M),
    log_lo = qs[1,], log_med = qs[2,], log_hi = qs[3,]
  ) %>%
    mutate(
      HR_med = exp(log_med),
      HR_lo  = exp(log_lo),
      HR_hi  = exp(log_hi),
      base   = sub("\\[.*$", "", variable),
      index  = suppressWarnings(as.integer(sub("^.*\\[(\\d+)\\].*$", "\\1", variable)))
    )
}

intervention_pars <- c("beta_int1_1", "beta_int1_2", "beta_int2_1", "beta_int2_2")

# count divergences (works for rstan brmsfit and cmdstanr)
# count_divergences <- function(fit) {
#   # try rstan path (brmsfit$fit is stanfit)
#   divs <- tryCatch({
#     if (!is.null(fit$fit)) {
#       sp <- get_sampler_params(fit$fit, inc_warmup = FALSE)
#       return(sum(vapply(sp, function(ch) sum(ch[, "divergent__"]), numeric(1))))
#     }
#     NA_real_
#   }, error = function(e) NA_real_)
#   if (!is.na(divs)) return(as.integer(divs))
#   
#   # fallback via posterior draws df (cmdstanr or when available)
#   divs2 <- tryCatch({
#     ddf <- as_draws_df(fit)
#     if (".divergent" %in% names(ddf)) sum(as.integer(ddf$.divergent), na.rm = TRUE)
#     else if ("divergent__" %in% names(ddf)) sum(as.integer(ddf$divergent__), na.rm = TRUE)
#     else NA_real_
#   }, error = function(e) NA_real_)
#   if (is.na(divs2)) return(NA_integer_) else as.integer(divs2)
# }

count_divergences_stanfit <- function(sf) {
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  as.integer(sum(vapply(sp, function(ch) sum(ch[, "divergent__"], na.rm = TRUE), 0)))
}

# Rhat min/max using posterior
rhat_min_max <- function(fit) {
  tryCatch({
    dr <- as_draws(fit)
    sdf <- posterior::summarise_draws(dr, "rhat")
    c(min = suppressWarnings(min(sdf$rhat, na.rm = TRUE)),
      max = suppressWarnings(max(sdf$rhat, na.rm = TRUE)))
  }, error = function(e) c(min = NA_real_, max = NA_real_))
}

# -------------------------------------------------------------------------------------
# Load all fits + compute diagnostics/LOO + extract intervention params
# -------------------------------------------------------------------------------------
results <- scenarios %>%
  mutate(
    scenario_id = row_number(),
    fit = imap(fit_path, ~ {
      message(">>> Loading scenario ", .y, " of ", nrow(scenarios), " : ", .x)
      safe_read_rds(.x)
    }),
    ok = map_lgl(fit, ~ !is.null(.x))
  ) %>%
  filter(ok) %>%
  mutate(
    # diagnostics columns in the output tibble
    # divergences = map_int(fit, ~ {
    #   d <- count_divergences_stanfit(.x)
    #   d
    # }),
    rhat_stats = map(fit, ~ {
      rr <- rhat_min_max(.x)
      tibble(rhat_min = rr["min"], rhat_max = rr["max"])
    }),
    rhat_min = map_dbl(rhat_stats, ~ .x$rhat_min),
    rhat_max = map_dbl(rhat_stats, ~ .x$rhat_max),
    
    # print a concise line per scenario
    print = pmap(list(scenario_id, divergences, rhat_min, rhat_max), function(id, divs, rmin, rmax) {
      message(
        sprintf("    Scenario %d: divergences = %s, Rhat min = %s, Rhat max = %s",
                id,
                ifelse(is.na(divs), "NA", as.character(divs)),
                ifelse(is.na(rmin), "NA", signif(rmin, 4)),
                ifelse(is.na(rmax), "NA", signif(rmax, 4)))
      )
      NULL
    }),
    
    loo_obj = imap(fit, ~ {
      message(">>> Computing LOO for scenario ", .y, " of ", n())
      loo_from_stanfit(.x, "log_lik")
    }),
    dd = imap(fit, ~ {
      message(">>> Extracting draws for scenario ", .y, " of ", n())
      draws_df(.x)
    }),
    ints = imap(dd, ~ {
      message(">>> Summarising intervention params for scenario ", .y, " of ", n())
      bind_rows(lapply(intervention_pars, \(p) sum_param(.x, p))) %>%
        arrange(base, index)
    })
  ) %>%
  select(-rhat_stats, -`print`)  # keep the clean columns needed

# add divergences column to your existing results tibble
results <- results %>%
  mutate(divergences = map_int(fit, count_divergences_stanfit))

# -------------------------------------------------------------------------------------
# Build LOO comparison table
# -------------------------------------------------------------------------------------

loo_table <- results %>%
  transmute(
    scenario_label,
    loo_obj
  )

# Run loo_compare on available LOO objects
# -------------------------------------------------------------------------------------

loo_list <- set_names(loo_table$loo_obj, loo_table$scenario_label)
lc <- tryCatch(loo_compare(loo_list), error = function(e) NULL)

# Tidy LOO comparison (ELPD, SE, ΔELPD)
loo_tidy <- if (!is.null(lc)) {
  as.data.frame(lc) %>%
    rownames_to_column("scenario_label") %>%
    as_tibble() %>%
    arrange(desc(ELPD)) %>%
    mutate(
      LOOIC    = round(-2 * ELPD,1),            # absolute LOOIC
      dELPD    = round(ELPD - max(ELPD, na.rm = TRUE),1),
      dLOOIC   = round(-2 * dELPD,1),           # delta LOOIC
      SE_LOOIC =  round(2 * SE,1)               # SE on LOOIC scale
    )
  }else{
  results %>%
    transmute(
      scenario_label,
      ELPD = map_dbl(loo_obj, ~ if (is.null(.x)) NA_real_ else .x$estimates["elpd_loo","Estimate"]),
      SE_ELPD   = map_dbl(loo_obj, ~ if (is.null(.x)) NA_real_ else .x$estimates["elpd_loo","SE"])
    )%>%
    arrange(desc(ELPD)) %>%
    mutate(
      LOOIC    = round(-2 * ELPD,1),            # absolute LOOIC
      dELPD    = round(ELPD - max(ELPD, na.rm = TRUE),1),
      dLOOIC   = round(-2 * dELPD,1),           # delta LOOIC
      SE_LOOIC =  round(2 * SE_ELPD,1),
      SE_ELPD = round(SE_ELPD, 1)
    )%>%
    select(scenario_label, ELPD, SE_ELPD, dELPD, LOOIC, dLOOIC, SE_LOOIC)
}

rhat_table <- results %>%
  transmute(
    scenario_label,
    rhat_min = round(rhat_min, 4),
    rhat_max = round(rhat_max, 4),
    divergences
  )

diagnostics_table <- left_join(loo_tidy, rhat_table)
diagnostics_table


# -------------------------------------------------------------------------------------
# Tidy intervention effects for plotting
# -------------------------------------------------------------------------------------

ints_long <- results %>%
  select(scenario_label, intervention, seasonality, intervention_effect, scen, data_source, ints) %>%
  unnest(ints) %>%
  mutate(
    pretty = case_when(
      base == "beta_int1_1" ~ "Intervention effect λ12 (Post 1)",
      base == "beta_int1_2" ~ "Intervention effect λ12 (Post 2)",
      base == "beta_int2_1" ~ "Intervention effect λ21 (Post 1)",
      base == "beta_int2_2" ~ "Intervention effect λ21 (Post 2)",
      TRUE ~ base
    ),
    scenario_label = fct_reorder(scenario_label, HR_med) 
  )

# -------------------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------------------

# LOO comparison (ΔLOOIC)
p_loo <- ggplot(loo_tidy, aes(x = fct_reorder(scenario_label, dLOOIC), y = dLOOIC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  coord_flip() +
  labs(title = "Model fit comparison across scenarios", subtitle = "ΔLOOIC (lower is better)",
       x = "", y = "ΔLOOIC vs best") +
  theme_minimal(base_size = 13)
p_loo

# Intervention effects (HR) side-by-side by scenario (one facet per parameter, CI = 95%)
p_int <- ggplot(
  ints_long %>% filter(!is.na(HR_med)),
  aes(x = scenario_label, y = HR_med)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.15) +
  geom_point(size = 2) +
  coord_flip() +
  facet_wrap(~ pretty, scales = "free_x") +
  scale_y_log10() +
  labs(
    title = "Intervention effects across scenarios",
    subtitle = "Posterior median HR with 95% CrI (log scale)",
    x = "", y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 12)
p_int

# break out by seasonality family for a cleaner look
p_int_byfam <- ggplot(
  ints_long %>% filter(!is.na(HR_med)) %>%
    mutate(family = ifelse(seasonality == "sine_", "Sine", "Spline")),
  aes(x = fct_reorder(scenario_label, HR_med), y = HR_med, color = family)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.15) +
  geom_point(size = 2) +
  coord_flip() +
  facet_wrap(~ pretty, scales = "free_x") +
  scale_y_log10() +
  labs(
    title = "Intervention effects by seasonal family",
    subtitle = "Posterior median HR with 95% CrI (log scale)",
    x = "", y = "Hazard Ratio", color = "Seasonality"
  ) +
  theme_minimal(base_size = 12)
p_int_byfam

# -------------------------------------------------------------------------------------
# Tables for the paper
# -------------------------------------------------------------------------------------

diag_ft <- flextable(
  diagnostics_table %>%
    select(c(scenario_label, LOOIC, SE_LOOIC, dLOOIC, rhat_min, rhat_max, divergences))
) |> autofit() |> set_caption("LOO model comparison across scenarios")
diag_ft

ints_ft <- flextable(
  ints_long %>%
    transmute(
      Scenario = scenario_label,
      Parameter = pretty,
      HR = sprintf("%.2f [%.2f–%.2f]", HR_med, HR_lo, HR_hi)
    )
) |> autofit() |> set_caption("Intervention effects (HR) across scenarios")

loo_ft; ints_ft

# -------------------------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------------------------

out_dir <- "./Output/Model_results/Model_summaries/Observed_data/Model_fit/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# CSVs
write_csv(diagnostics_table, file.path(out_dir, "loo_compare.csv"))
write_csv(
  ints_long %>% select(scenario_label, pretty, HR_med, HR_lo, HR_hi),
  file.path(out_dir, "intervention_effects_hr.csv")
)


# Combined PDF for quick review
pdf(file = file.path(out_dir, "scenario_compare.pdf"), width = 11.7, height = 8.3) # A4 landscape
grid.newpage(); grid.draw(ggplotGrob(p_loo))
grid.newpage(); grid.draw(ggplotGrob(p_int))
grid.newpage(); grid.draw(ggplotGrob(p_int_byfam))
grid.newpage(); grid.table(loo_tidy)
grid.newpage(); grid.table(
  ints_long %>% select(scenario_label, pretty, HR_med, HR_lo, HR_hi) %>%
    arrange(pretty, scenario_label)
)
dev.off()

# NOW STORE RESULTS
#---------------------------------------------------------------
saveRDS(results, file.path(out_dir, "summary_results.rds"))  

