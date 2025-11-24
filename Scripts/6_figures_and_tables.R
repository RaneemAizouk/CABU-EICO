###############################################################################
# THIS CODE VISUALISES THE MANUSCRIPT RESULTS
###############################################################################

# Date created: 1 October 2025
# Date last updated: 1 October 2025
# Author: Esther van Kleef

# Load in packages
pacman::p_load(dplyr, tibble, tidyr, stringr, ggplot2, posterior, rstan, readr, 
               tools, ggplot2, patchwork)

# SET DIRECTORY
DirectoryData <- "./Data/BF/clean"
load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_all.rda"))
load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_completecases.rda"))

#------------------------------------------------------------------------------
# Build data (draw-level incidences, summaries, HRs) for any fit
#------------------------------------------------------------------------------

# Quantile function
qsum <- function(x) c(
  med        = stats::median(x, na.rm = TRUE),
  "lo"  = stats::quantile(x, 0.025, na.rm = TRUE),
  "hi" = stats::quantile(x, 0.975, na.rm = TRUE)
)


# Build segments once from stan_data_fit 
build_segs <- function(stan_data_fit) {
  N  <- stan_data_fit$N
  L  <- stan_data_fit$interval_length
  g0 <- stan_data_fit$global_interval_start
  K  <- stan_data_fit$num_intervals
  
  obs <- tibble(
    n        = seq_len(N),
    house    = stan_data_fit$HouseID,
    id       = stan_data_fit$menage_id_member,
    age      = stan_data_fit$age,
    sexe     = stan_data_fit$sexe,
    t        = stan_data_fit$date_use,
    int_flag = as.integer(stan_data_fit$intervention > 0.5),
    t_star1  = stan_data_fit$intervention_date,
    t_star2  = stan_data_fit$intervention_date2
  ) %>% arrange(id, n) %>%
    mutate(t0 = lag(t), first_obs = is.na(t0) | (id != lag(id)))
  
  clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
  make_one_obs_subints <- function(i) {
    if (obs$first_obs[i]) return(NULL)
    t0 <- obs$t0[i]; t1 <- obs$t[i]
    raw_start <- floor((t0 - g0) / L); raw_end <- ceiling((t1 - g0) / L)
    gi_start <- clamp(raw_start + 1, 1, K); gi_end <- clamp(raw_end, 1, K)
    a1 <- g0 + (gi_start - 1) * L; b1 <- a1 + L - 1
    first_len <- min(c(b1 - t0 + 1, t1 - t0 + 1, L))
    first <- NULL
    if (first_len > 0) {
      mid1 <- (t0 + b1)/2; idx_first <- clamp(floor((mid1 - g0)/L) + 1, 1, K)
      first <- list(list(a = t0, b = t0 + first_len - 1, idx = idx_first))
    }
    n_mid <- gi_end - gi_start - 1; mids <- NULL
    if (n_mid > 0) {
      mids <- vector("list", n_mid)
      for (m in 1:n_mid) {
        a <- g0 + (gi_start + m - 1) * L; b <- a + L - 1
        midm <- (a + b)/2; idxm <- clamp(floor((midm - g0)/L) + 1, 1, K)
        mids[[m]] <- list(a=a,b=b,idx=idxm)
      }
    }
    last_len <- min(c(t1 - (g0 + (gi_end - 1) * L) + 1, L))
    last <- NULL
    if (last_len > 0) {
      al <- t1 - last_len + 1; bl <- t1
      midl <- ((g0 + (gi_end - 1) * L) + t1)/2; idx_last <- clamp(floor((midl - g0)/L) + 1, 1, K)
      last <- list(list(a=al,b=bl,idx=idx_last))
    }
    c(first, mids, last)
  }
  overlap_len <- function(a, b, Lb, Ub) {
    if (is.na(a) || is.na(b)) return(0)
    as.numeric(max(0, min(b, Ub) - max(a, Lb) + 1))
  }
  
  seg_list <- vector("list", N)
  for (i in seq_len(N)) {
    if (obs$first_obs[i]) next
    subints <- make_one_obs_subints(i)
    if (is.null(subints) || !length(subints)) next
    
    t1 <- obs$t_star1[i]; t2 <- obs$t_star2[i]
    grp <- ifelse(obs$int_flag[i] == 1L, "Intervention", "Control")
    
    for (si in subints) {
      a <- si$a; b <- si$b; idx <- si$idx
      if (grp == "Intervention") {
        d0 <- overlap_len(a,b, -Inf, t1-1)
        d1 <- overlap_len(a,b, t1,   t2-1)
        d2 <- overlap_len(a,b, t2,   Inf)
        rows <- NULL
        if (d0>0) rows <- rbind(rows, c(idx, obs$house[i], obs$age[i], obs$sexe[i], 0L, d0, grp))
        if (d1>0) rows <- rbind(rows, c(idx, obs$house[i], obs$age[i], obs$sexe[i], 1L, d1, grp))
        if (d2>0) rows <- rbind(rows, c(idx, obs$house[i], obs$age[i], obs$sexe[i], 2L, d2, grp))
        if (!is.null(rows)) seg_list[[i]] <- rbind(seg_list[[i]], rows)
      } else {
        seg_list[[i]] <- rbind(seg_list[[i]], c(idx, obs$house[i], obs$age[i], obs$sexe[i], NA, b-a+1, grp))
      }
    }
  }
  segs <- do.call(rbind, seg_list) %>% as_tibble()
  names(segs) <- c("idx","house","age","sexe","iperiod","dur","group")
  segs %>%
    mutate(
      idx=as.integer(idx), house=as.integer(house),
      iperiod=suppressWarnings(as.integer(iperiod)),
      dur=as.numeric(dur),
      group=factor(group, levels=c("Intervention","Control"))
    )
}

# Compute draw-level datasets for one saved fit
# compute_from_fit <- function(results_path, segs) {
#   fit   <- readRDS(results_path)
#   draws <- as_draws_df(fit)
#   
#   D <- nrow(draws); R <- nrow(segs)
#   
#   # seasonal term per bin
#   Y_cols <- grep("^Y_hat_1_2_out\\[\\d+\\]$", names(draws), value = TRUE)
#   Y_idx  <- as.integer(stringr::str_match(Y_cols, "\\[(\\d+)\\]")[,2])
#   Y_cols <- Y_cols[order(Y_idx)]
#   Y_mat  <- as.matrix(draws[, Y_cols, drop = FALSE])  # [D x K]
#   
#   # household RE
#   u_cols <- grep("^u\\[\\d+\\]$", names(draws), value = TRUE)
#   u_idx  <- as.integer(stringr::str_match(u_cols, "\\[(\\d+)\\]")[,2])
#   u_cols <- u_cols[order(u_idx)]
#   U_mat  <- as.matrix(draws[, u_cols, drop = FALSE])  # [D x H]
#   
#   # coefficients
#   q12_base <- draws[["q_1_2_base"]]
#   b_age    <- draws[["beta_1_2_age"]]
#   b_sex    <- draws[["beta_1_2_sexe"]]
#   d_int1   <- draws[["beta_int1_1"]]
#   d_int2   <- draws[["beta_int1_2"]]
#   
#   # linear predictor per segment row
#   Y_by_row   <- Y_mat[, segs$idx,   drop = FALSE]
#   U_by_row   <- U_mat[, segs$house, drop = FALSE]
#   AGE_by_row <- matrix(as.numeric(segs$age),  nrow = D, ncol = R, byrow = TRUE)
#   SEX_by_row <- matrix(as.numeric(segs$sexe), nrow = D, ncol = R, byrow = TRUE)
#   BASE       <- matrix(q12_base,  nrow = D, ncol = R)
#   LP         <- BASE + Y_by_row + U_by_row + (b_age * AGE_by_row) + (b_sex * SEX_by_row)
#   
#   # intervention effects only on I1/I2
#   int_eff <- matrix(0, nrow = D, ncol = R)
#   ix_I0 <- which(segs$group=="Intervention" & segs$iperiod==0L)
#   ix_I1 <- which(segs$group=="Intervention" & segs$iperiod==1L)
#   ix_I2 <- which(segs$group=="Intervention" & segs$iperiod==2L)
#   if (length(ix_I1)) int_eff[, ix_I1] <- matrix(d_int1, nrow = D, ncol = length(ix_I1))
#   if (length(ix_I2)) int_eff[, ix_I2] <- matrix(d_int2, nrow = D, ncol = length(ix_I2))
#   
#   lambda <- exp(LP + int_eff)
#   Dur    <- matrix(segs$dur, nrow = D, ncol = R, byrow = TRUE)
#   
#   sum_w <- function(mat, idx) if (length(idx)) rowSums(mat[, idx, drop=FALSE]) else rep(NA_real_, D)
#   rate_from_rows <- function(idx) {
#     if (!length(idx)) return(rep(NA_real_, D))
#     num <- sum_w(lambda * Dur, idx)
#     den <- sum(segs$dur[idx])
#     num / den
#   }
#   
#   # I0/I1/I2
#   R_I0 <- rate_from_rows(ix_I0)
#   R_I1 <- rate_from_rows(ix_I1)
#   R_I2 <- rate_from_rows(ix_I2)
#   
#   # calendar weights from intervention for time-matching controls
#   weights_from_I <- function(ix) {
#     tab <- tapply(segs$dur[ix], segs$idx[ix], sum)
#     idxs <- as.integer(names(tab)); w <- as.numeric(tab) / sum(tab)
#     list(idx = idxs, w = w)
#   }
#   w_I0 <- weights_from_I(ix_I0); w_I1 <- weights_from_I(ix_I1); w_I2 <- weights_from_I(ix_I2)
#   
#   # control per-bin hazards
#   ix_C_all  <- which(segs$group=="Control")
#   ctrl_bins <- sort(unique(segs$idx[ix_C_all]))
#   H_C <- sapply(ctrl_bins, function(j) {
#     jj  <- ix_C_all[segs$idx[ix_C_all] == j]
#     num <- sum_w(lambda * Dur, jj)
#     den <- sum(segs$dur[jj])
#     num / den
#   })
#   if (is.null(dim(H_C))) H_C <- matrix(H_C, nrow=D)
#   colnames(H_C) <- as.character(ctrl_bins)
#   
#   std_rate <- function(wlist) {
#     common <- intersect(colnames(H_C), as.character(wlist$idx))
#     if (!length(common)) return(rep(NA_real_, D))
#     H <- H_C[, match(common, colnames(H_C)), drop = FALSE]
#     w <- wlist$w[match(as.integer(common), wlist$idx)]
#     w <- w / sum(w)
#     as.vector(H %*% w)
#   }
#   Rt_C0 <- std_rate(w_I0)
#   Rt_C1 <- std_rate(w_I1)
#   Rt_C2 <- std_rate(w_I2)
#   
#   # per-draw incidence per 100 PD (for plotting)
#   make_df <- function(draw_vec, grp, phase_label) {
#     tibble(
#       draw         = seq_len(D),
#       Group        = grp,
#       phase_label  = phase_label,
#       ID12_per100  = 100 * as.numeric(draw_vec)
#     )
#   }
#   rates_draws <- bind_rows(
#     make_df(R_I0,  "Intervention", "Pre-intervention"),
#     make_df(R_I1,  "Intervention", "3-months\npost-intervention"),
#     make_df(R_I2,  "Intervention", "9-months\npost-intervention"),
#     make_df(Rt_C0, "Control",      "Pre-intervention"),
#     make_df(Rt_C1, "Control",      "3-months\npost-intervention"),
#     make_df(Rt_C2, "Control",      "9-months\npost-intervention")
#   ) %>%
#     filter(is.finite(ID12_per100)) %>%
#     mutate(
#       Group = factor(Group, levels = c("Intervention","Control")),
#       phase_label = factor(
#         phase_label,
#         levels = c("Pre-intervention","3-months\npost-intervention","9-months\npost-intervention")
#       )
#     )
#   
#   # tidy summary by Group x Phase
#   rates_cri <- rates_draws %>%
#     group_by(Group, phase_label) %>%
#     summarise(
#       Median = median(ID12_per100),
#       Lower  = quantile(ID12_per100, 0.025),
#       Upper  = quantile(ID12_per100, 0.975),
#       .groups = "drop"
#     )
#   
#   # HRs from betas with your indexing
#   hr_tbl <- tibble(
#     period = c("Phase 1","Phase 2"),
#     med    = c(qsum(exp(draws[["beta_int1_1"]]))["med"],
#                qsum(exp(draws[["beta_int1_2"]]))["med"]),
#     lo     = c(qsum(exp(draws[["beta_int1_1"]]))["lo.2.5%"],
#                qsum(exp(draws[["beta_int1_2"]]))["lo.2.5%"]),
#     hi     = c(qsum(exp(draws[["beta_int1_1"]]))["hi.97.5%"],
#                qsum(exp(draws[["beta_int1_2"]]))["hi.97.5%"])
#   ) %>%
#     mutate(
#       phase_label = recode(
#         period,
#         "Phase 1" = "Phase 1\n~3-months",
#         "Phase 2" = "Phase 2\n~9-months"
#       ),
#       hr_label = sprintf("HR %.2f (%.2f–%.2f)", med, lo, hi),
#       Group    = factor("Intervention", levels = c("Intervention","Control"))
#     )
#   
#   # y-positions for HR annotations (per facet & phase)
#   ypos_box <- rates_draws %>%
#     group_by(Group, phase_label) %>%
#     summarise(ypos = stats::quantile(ID12_per100, 0.975, na.rm = TRUE) * 1.08, .groups = "drop")
#   ann_box <- inner_join(hr_tbl, ypos_box, by = c("Group","phase_label")) %>%
#     select(Group, phase_label, hr_label, ypos)
#   
#   list(
#     rates_draws = rates_draws,
#     rates_cri   = rates_cri,
#     hr_tbl      = hr_tbl,
#     ann_box     = ann_box
#   )
# }

compute_from_fit <- function(results_path, segs) {
  fit   <- readRDS(results_path)
  draws <- as_draws_df(fit)
  
  D <- nrow(draws); R <- nrow(segs)
  
  # seasonal term per bin
  Y_cols <- grep("^Y_hat_1_2_out\\[\\d+\\]$", names(draws), value = TRUE)
  Y_idx  <- as.integer(stringr::str_match(Y_cols, "\\[(\\d+)\\]")[,2])
  Y_cols <- Y_cols[order(Y_idx)]
  Y_mat  <- as.matrix(draws[, Y_cols, drop = FALSE])  # [D x K]
  
  # household RE
  u_cols <- grep("^u\\[\\d+\\]$", names(draws), value = TRUE)
  u_idx  <- as.integer(stringr::str_match(u_cols, "\\[(\\d+)\\]")[,2])
  u_cols <- u_cols[order(u_idx)]
  U_mat  <- as.matrix(draws[, u_cols, drop = FALSE])  # [D x H]
  
  # coefficients
  q12_base <- draws[["q_1_2_base"]]
  b_age    <- draws[["beta_1_2_age"]]
  b_sex    <- draws[["beta_1_2_sexe"]]
  d_int1   <- draws[["beta_int1_1"]]  # early intv beta
  d_int2   <- draws[["beta_int1_2"]]  # late  intv beta
  
  # linear predictor per segment row
  Y_by_row   <- Y_mat[, segs$idx,   drop = FALSE]
  U_by_row   <- U_mat[, segs$house, drop = FALSE]
  AGE_by_row <- matrix(as.numeric(segs$age),  nrow = D, ncol = R, byrow = TRUE)
  SEX_by_row <- matrix(as.numeric(segs$sexe), nrow = D, ncol = R, byrow = TRUE)
  BASE       <- matrix(q12_base,  nrow = D, ncol = R)
  LP         <- BASE + Y_by_row + U_by_row + (b_age * AGE_by_row) + (b_sex * SEX_by_row)
  
  # intervention effects only on I1/I2
  int_eff <- matrix(0, nrow = D, ncol = R)
  ix_I0 <- which(segs$group=="Intervention" & segs$iperiod==0L)
  ix_I1 <- which(segs$group=="Intervention" & segs$iperiod==1L)
  ix_I2 <- which(segs$group=="Intervention" & segs$iperiod==2L)
  if (length(ix_I1)) int_eff[, ix_I1] <- matrix(d_int1, nrow = D, ncol = length(ix_I1))
  if (length(ix_I2)) int_eff[, ix_I2] <- matrix(d_int2, nrow = D, ncol = length(ix_I2))
  
  lambda <- exp(LP + int_eff)
  Dur    <- matrix(segs$dur, nrow = D, ncol = R, byrow = TRUE)
  
  sum_w <- function(mat, idx) {
    if (length(idx)) rowSums(mat[, idx, drop=FALSE]) else rep(NA_real_, D)
  }
  rate_from_rows <- function(idx) {
    if (!length(idx)) return(rep(NA_real_, D))
    num <- sum_w(lambda * Dur, idx)
    den <- sum(segs$dur[idx])
    num / den
  }
  
  ## Phase-specific rates (I0/I1/I2) for Intervention
  R_I0 <- rate_from_rows(ix_I0)  # Intervention baseline
  R_I1 <- rate_from_rows(ix_I1)  # Intervention early post
  R_I2 <- rate_from_rows(ix_I2)  # Intervention late post
  
  # calendar weights from intervention for time-matching controls
  weights_from_I <- function(ix) {
    if (!length(ix)) return(list(idx = integer(0), w = numeric(0)))
    tab <- tapply(segs$dur[ix], segs$idx[ix], sum)
    idxs <- as.integer(names(tab))
    w    <- as.numeric(tab) / sum(tab)
    list(idx = idxs, w = w)
  }
  w_I0 <- weights_from_I(ix_I0)
  w_I1 <- weights_from_I(ix_I1)
  w_I2 <- weights_from_I(ix_I2)
  
  # control per-bin hazards
  ix_C_all  <- which(segs$group=="Control")
  ctrl_bins <- sort(unique(segs$idx[ix_C_all]))
  H_C <- sapply(ctrl_bins, function(j) {
    jj  <- ix_C_all[segs$idx[ix_C_all] == j]
    num <- sum_w(lambda * Dur, jj)
    den <- sum(segs$dur[jj])
    num / den
  })
  if (is.null(dim(H_C))) H_C <- matrix(H_C, nrow=D)
  colnames(H_C) <- as.character(ctrl_bins)
  
  std_rate <- function(wlist) {
    common <- intersect(colnames(H_C), as.character(wlist$idx))
    if (!length(common)) return(rep(NA_real_, D))
    H <- H_C[, match(common, colnames(H_C)), drop = FALSE]
    w <- wlist$w[match(as.integer(common), wlist$idx)]
    w <- w / sum(w)
    as.vector(H %*% w)
  }
  Rt_C0 <- std_rate(w_I0)  # Control baseline (time-matched)
  Rt_C1 <- std_rate(w_I1)  # Control early post
  Rt_C2 <- std_rate(w_I2)  # Control late post
  
  ## Overall post-intervention rates as duration-weighted average of phase 1 & 2
  T1 <- sum(segs$dur[ix_I1])
  T2 <- sum(segs$dur[ix_I2])
  T_post <- T1 + T2
  
  R_I_overall  <- (R_I1 * T1 + R_I2 * T2) / T_post
  Rt_C_overall <- (Rt_C1 * T1 + Rt_C2 * T2) / T_post
  
  # Raw post-only IRR (Int vs Ctrl), per draw
  IRR_overall_raw_draws <- R_I_overall / Rt_C_overall
  
  # Difference-in-difference overall IRR: (postI/preI) / (postC/preC)
  DiD_overall_draws <- (R_I_overall / R_I0) / (Rt_C_overall / Rt_C0)
  
  overall_raw_irr <- tibble(
    effect = "Overall post-intervention IRR (Int vs Ctrl)",
    med    = qsum(IRR_overall_raw_draws)["med"],
    lo     = qsum(IRR_overall_raw_draws)["lo.2.5%"],
    hi     = qsum(IRR_overall_raw_draws)["hi.97.5%"]
  )
  
  overall_did_irr <- tibble(
    effect = "Overall post-intervention DiD-IRR",
    med    = qsum(DiD_overall_draws)["med"],
    lo     = qsum(DiD_overall_draws)["lo.2.5%"],
    hi     = qsum(DiD_overall_draws)["hi.97.5%"]
  )
  
  ## ----------------------------------------
  ## Per-draw incidence per 100 PD (for plotting)
  make_df <- function(draw_vec, grp, phase_label) {
    tibble(
      draw         = seq_len(D),
      Group        = grp,
      phase_label  = phase_label,
      ID12_per100  = 100 * as.numeric(draw_vec)
    )
  }
  rates_draws <- bind_rows(
    make_df(R_I0,        "Intervention", "Pre-intervention"),
    make_df(R_I1,        "Intervention", "3-months\npost-intervention"),
    make_df(R_I2,        "Intervention", "9-months\npost-intervention"),
    make_df(Rt_C0,       "Control",      "Pre-intervention"),
    make_df(Rt_C1,       "Control",      "3-months\npost-intervention"),
    make_df(Rt_C2,       "Control",      "9-months\npost-intervention"),
    make_df(R_I_overall, "Intervention", "Overall post-intervention"),
    make_df(Rt_C_overall,"Control",      "Overall post-intervention")
  ) %>%
    filter(is.finite(ID12_per100)) %>%
    mutate(
      Group = factor(Group, levels = c("Intervention","Control")),
      phase_label = factor(
        phase_label,
        levels = c(
          "Pre-intervention",
          "3-months\npost-intervention",
          "9-months\npost-intervention",
          "Overall post-intervention"
        )
      )
    )
  
  # tidy summary by Group x Phase
  rates_cri <- rates_draws %>%
    group_by(Group, phase_label) %>%
    summarise(
      Median = median(ID12_per100),
      Lower  = quantile(ID12_per100, 0.025),
      Upper  = quantile(ID12_per100, 0.975),
      .groups = "drop"
    )
  
  ## ----------------------------------------
  ## 1) HRs from betas (model-based)
  
  HR1_draws <- exp(d_int1)
  HR2_draws <- exp(d_int2)
  
  h1 <- qsum(HR1_draws)
  h2 <- qsum(HR2_draws)
  
  hr_tbl_HR <- tibble(
    period = c("Phase 1","Phase 2"),
    med    = c(h1["med"],      h2["med"]),
    lo     = c(h1["lo.2.5%"],  h2["lo.2.5%"]),
    hi     = c(h1["hi.97.5%"], h2["hi.97.5%"])
  )
  
  ## ----------------------------------------
  ## 2) DiD-IRRs from incidence (baseline-adjusted)
  
  IRR1_draws <- R_I1 / Rt_C1
  IRR2_draws <- R_I2 / Rt_C2
  
  DiD1_draws <- (R_I1 / R_I0) / (Rt_C1 / Rt_C0)  # early phase
  DiD2_draws <- (R_I2 / R_I0) / (Rt_C2 / Rt_C0)  # late phase
  
  s1 <- qsum(DiD1_draws)
  s2 <- qsum(DiD2_draws)
  
  irr_tbl_DiD <- tibble(
    period = c("Phase 1","Phase 2"),
    med    = c(s1["med"],      s2["med"]),
    lo     = c(s1["lo.2.5%"],  s2["lo.2.5%"]),
    hi     = c(s1["hi.97.5%"], s2["hi.97.5%"])
  )
  
  ## ----------------------------------------
  
  list(
    rates_draws              = rates_draws,
    rates_cri                = rates_cri,
    
    # HRs (model-based)
    HR1_draws                = HR1_draws,
    HR2_draws                = HR2_draws,
    hr_tbl_HR                = hr_tbl_HR,
    
    # IRRs & DiD-IRRs (incidence-based)
    IRR1_draws               = IRR1_draws,
    IRR2_draws               = IRR2_draws,
    DiD1_draws               = DiD1_draws,
    DiD2_draws               = DiD2_draws,
    irr_tbl_DiD              = irr_tbl_DiD,
    
    IRR_overall_raw_draws    = IRR_overall_raw_draws,
    DiD_overall_draws        = DiD_overall_draws,
    IRR_overall_raw_summary  = overall_raw_irr,
    IRR_overall_did_summary  = overall_did_irr
  )
}


#------------------------------------------------------------------------------
# PLOTTING
#------------------------------------------------------------------------------

base_size <- 14

# - hr_tbl needs: period ("Phase 1","Phase 2"), med, lo, hi
# - rates_draws needs: Group, phase_label, ID12_per100
# make_ann_box_from_hr_tbl <- function(rates_draws, hr_tbl, upper_cap = 50, collapse_post = FALSE) {
#   recode_phase <- function(x) case_when(
#     x %in% c("Pre-intervention","Baseline") ~ "Baseline",
#     x %in% c("3-months\npost-intervention","Post-Intervention 1","Phase 1") ~ "Post-Intervention 1",
#     x %in% c("9-months\npost-intervention","Post-Intervention 2","Phase 2") ~ "Post-Intervention 2",
#     TRUE ~ as.character(x)
#   )
#   
#   # HR labels (Intervention only)
#   hr_annot <- hr_tbl %>%
#     mutate(
#       phase_label = recode(
#         period,
#         "Phase 1" = "Post-Intervention 1",
#         "Phase 2" = "Post-Intervention 2"
#       ),
#       hr_label = sprintf("HR %.2f (%.2f–%.2f)", med, lo, hi),
#       Group    = factor("Intervention", levels = c("Intervention","Control"))
#     ) %>%
#     select(Group, phase_label, hr_label)
#   
#   # Per-phase y position from the box distribution (97.5% * 1.08), clamped
#   ypos_box <- rates_draws %>%
#     mutate(phase_label = recode_phase(as.character(phase_label))) %>%
#     group_by(Group, phase_label) %>%
#     summarise(
#       ypos = pmin(stats::quantile(ID12_per100, 0.975, na.rm = TRUE) * 1.08,
#                   upper_cap * 0.98),
#       .groups = "drop"
#     )
#   
#   ann <- inner_join(hr_annot, ypos_box, by = c("Group","phase_label"))
#   
#   if (collapse_post) {
#     # Collapse the two post phases to one label, placed at the higher of the two
#     ann <- ann %>%
#       mutate(phase_label = if_else(
#         phase_label %in% c("Post-Intervention 1","Post-Intervention 2"),
#         "Post-Intervention 2", # place label on the second row
#         phase_label
#       )) %>%
#       group_by(Group, phase_label) %>%
#       summarise(
#         hr_label = first(hr_label), # keep Phase 2 label (or replace with custom string)
#         ypos     = max(ypos, na.rm = TRUE),
#         .groups  = "drop"
#       )
#   }
#   
#   ann
# }

base_size <- 14

# - hr_tbl needs: period ("Phase 1","Phase 2"), med, lo, hi
# - rates_draws needs: Group, phase_label, ID12_per100
# hr_tbl needs: period (e.g. "Phase 1","Phase 2","Overall"),
# with med/lo/hi = DiD-IRR summaries.
# rates_draws needs: Group, phase_label, ID12_per100
make_ann_box_from_hr_tbl <- function(rates_draws, hr_tbl, IRR_overall_did_summary, upper_cap = 50, collapse_post = FALSE) {
  names(IRR_overall_did_summary) = c("period", "med", "lo", "hi")
  hr_tbl = rbind(hr_tbl, IRR_overall_did_summary)
  
  recode_phase <- function(x) dplyr::case_when(
    x %in% c("Pre-intervention", "Baseline") ~ "Baseline",
    x %in% c("3-months\npost-intervention", "Post-Intervention 1", "Phase 1") ~ "Post-Intervention 1",
    x %in% c("9-months\npost-intervention", "Post-Intervention 2", "Phase 2") ~ "Post-Intervention 2",
    x %in% c("Overall post-intervention", "Overall", "Overall post-intervention DiD-IRR") ~ "Post-intervention overall",
    TRUE ~ as.character(x)
  )
  
  # DiD-IRR labels (Intervention only) – assumes med/lo/hi in hr_tbl are DiD-IRRs
  hr_annot <- hr_tbl %>%
    dplyr::mutate(
      phase_label = dplyr::recode(
        period,
        "Phase 1" = "Post-Intervention 1",
        "Phase 2" = "Post-Intervention 2",
        "Overall" = "Post-intervention overall",
        "Overall post-intervention" = "Post-intervention overall",
        "Overall post-intervention DiD-IRR" = "Post-intervention overall"
      ),
      hr_label = sprintf("IRR %.2f (%.2f–%.2f)", med, lo, hi),
      Group    = factor("Intervention", levels = c("Intervention", "Control"))
    ) %>%
    dplyr::select(Group, phase_label, hr_label)
  
  # Per-phase y position from the box distribution (97.5% * 1.08), clamped
  ypos_box <- rates_draws %>%
    mutate(phase_label = recode_phase(as.character(phase_label))) %>%
    group_by(Group, phase_label) %>%
    summarise(
      ypos = pmin(
        stats::quantile(ID12_per100, 0.975, na.rm = TRUE) * 1.08,
        upper_cap * 0.98
      ),
      .groups = "drop"
    )
  
  ann <- dplyr::inner_join(hr_annot, ypos_box, by = c("Group", "phase_label"))
  
  if (collapse_post) {
    # Collapse the two post phases to one label, placed at the higher of the two
    ann <- ann %>%
      dplyr::mutate(
        phase_label = dplyr::if_else(
          phase_label %in% c("Post-Intervention 1", "Post-Intervention 2"),
          "Post-Intervention 2",  # place combined label on the second row
          phase_label
        )
      ) %>%
      dplyr::group_by(Group, phase_label) %>%
      dplyr::summarise(
        hr_label = dplyr::first(hr_label),
        ypos     = max(ypos, na.rm = TRUE),
        .groups  = "drop"
      )
  }
  
  ann
}



#---------------------------------------------------------------------------------------------
# PLOT INCIDENCE
#---------------------------------------------------------------------------------------------
plot_incidence_panel <- function(data_list,
                                 panel_title,
                                 panel_subtitle,
                                 upper_cap = 25,
                                 base_size = 16) {
  # from compute_from_fit(): rates + DiD-IRR summaries
  rates_draws <- data_list$rates_draws   # needs: Group, phase_label, ID12_per100
  #hr_tbl      <- data_list$hr_tbl        # here: DiD-IRR table with period = "Phase 1","Phase 2","Overall"
  hr_tbl      <- data_list$ann_box        # here: DiD-IRR table with period = "Phase 1","Phase 2","Overall"
  
  # Build annotation (DiD-IRR labels, incl. overall) from rates_draws + hr_tbl
  # ann_box <- make_ann_box_from_hr_tbl(
  #   rates_draws = rates_draws,
  #   hr_tbl      = irr_tbl_DiD,
  #   RR_overall_did_summary = IRR_overall_did_summary,
  #   upper_cap   = upper_cap,
  #   collapse_post = FALSE
  # )
  ann_box <- data_list$ann_box 

  
  # Display phases as Baseline / Early / Late / Overall
  recode_phase <- function(x) dplyr::case_when(
    x %in% c("Pre-intervention","Baseline") ~ "Baseline",
    x %in% c("3-months\npost-intervention","Post-Intervention 1","Phase 1") ~ 
      "Early intervention phase",
    x %in% c("9-months\npost-intervention","Post-Intervention 2","Phase 2") ~ 
      "Late intervention phase",
    x %in% c("Overall post-intervention","Overall", "Post-intervention overall") ~ 
      "Intervention overall",
    TRUE ~ as.character(x)
  )
  
  # Data to plot
  rates_draws <- rates_draws %>%
    dplyr::mutate(
      phase_disp = recode_phase(as.character(phase_label)),
      phase_disp = factor(
        phase_disp,
        levels = c(
          "Baseline",
          "Early intervention phase",
          "Late intervention phase",
          "Intervention overall"
        )
      ),
      Group = factor(Group, levels = c("Intervention","Control")),
      fill_key = interaction(Group, phase_disp, drop = TRUE, lex.order = TRUE)
    )
  
  # Annotation (DiD-IRR text)
  ann_box <- ann_box %>%
    dplyr::mutate(
      phase_disp = recode_phase(as.character(phase_label)),
      Group      = factor(Group, levels = c("Intervention","Control")),
      y_plot     = pmin(ypos, 0.98 * upper_cap)
    )
  
  
  # Coral/reds for Intervention; greys for Control
  # fill_vals <- c(
  #   "Intervention.Baseline" = "#FEE5D9",
  #   "Intervention.Early intervention phase\n(>1–6 months)" = "#CB181D",
  #   "Intervention.Late intervention phase\n(>6–9 months)" = "#CB181D",
  #   "Intervention.Intervention overall"    = "#CB181D",
  #   "Control.Baseline" = "#F0F0F0",
  #   "Control.Early intervention phase\n(>1–6 months)" = "#BDBDBD",
  #   "Control.Late intervention phase\n(>6–9 months)" = "#BDBDBD",
  #   "Control.Intervention overall"    = "#BDBDBD"
  # )
  
  ggplot(rates_draws, aes(x = phase_disp, y = ID12_per100, fill = fill_key)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_text(
      data = ann_box,
      aes(x = phase_disp, y = y_plot, label = hr_label),
      inherit.aes = FALSE, hjust = 0, size = 4
    ) +
    #scale_fill_manual(values = fill_vals) +
    coord_flip(ylim = c(0, upper_cap), clip = "off") +
    scale_y_continuous(
      breaks = seq(0, upper_cap, 5),
      expand = expansion(mult = c(0.02, 0.18))
    ) +
    # Dotted separator between Late phase (3) and Overall (4)
    geom_vline(xintercept = 3.5, linetype = "dashed", size = 0.5, color = "black") +
    labs(
      title    = panel_title,
      subtitle = panel_subtitle,
      x = "", y = "Incidence per 100 person-days"
    ) +
    theme_minimal() +
    theme(
      strip.text.y   = element_text(),
      plot.margin    = margin(5.5, 30, 5.5, 5.5),
      plot.title     = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle  = element_text(size = base_size + 2),
      axis.title     = element_text(size = base_size + 2),
      axis.text      = element_text(size = base_size + 2),
      plot.caption   = element_text(size = base_size - 3),
      legend.text    = element_text(size = 16),
      legend.title   = element_blank(),
      legend.position = "bottom",
      panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
      panel.grid.minor = element_line(linetype = "dotted", color = "grey85")
    )
}

  
# PLOT TWO SCENARIOS (Scenario 1 and scenario 2 with basecase scenario)
#---------------------------------------------------------------------------------------------

# Inputs
results_path_1 <- "./Output/Model_results/Observed_data/Two_step_sine_NonAdd_NonCol_seasonalityobserved.rds"
results_path_2 <- "./Output/Model_results/Observed_data/One_step_sine_NonAdd_NonCol_seasonalityobserved.rds"  # example
title_1 <- "Scenario 1"
subtitle_1 <- "Two-step intervention effect"
title_2 <- "Scenario 2"
subtitle_2 <- "One-step intervention effect"
  
upper_cap <- 25 #
  
# Build segs 
stan_data_fit <- readRDS("./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")
segs <- build_segs(stan_data_fit)

# Prepare data for both fits
data_s1 <- compute_from_fit(results_path_1, segs)
data_s2 <- compute_from_fit(results_path_2, segs)
  
# Save dataframes
data_s1$rates_draws; data_s1$rates_cri; data_s1$hr_tbl
data_s1$hr_tbl_HR
write.csv(data_s1$rates_draws, "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario1/rates_draws_s1.csv")
write.csv(data_s1$rates_cri,   "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario1/rates_summary_s1.csv")
write.csv(data_s1$hr_tbl,      "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario1/hr_tbl_s1.csv")

data_s2$rates_draws; data_s2$rates_cri; data_s2$hr_tbl
write.csv(data_s1$rates_draws, "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario2/rates_draws_s2.csv")
write.csv(data_s1$rates_cri,   "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario2/rates_summary_s2.csv")
write.csv(data_s1$hr_tbl,      "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario2/hr_tbl_s2.csv")


# Plot Manuscript Figure 2
#-------------------------------------------------------------------------------
# Scenario 1
data_s1$ann_box <- make_ann_box_from_hr_tbl(
    rates_draws = data_s1$rates_draws,
    hr_tbl      = data_s1$irr_tbl_DiD,   # your HR table (Phase 1 / Phase 2)
    IRR_overall_did_summary = data_s1$IRR_overall_did_summary,
    upper_cap   = 25,
    collapse_post = FALSE           # <- set TRUE if you want one combined post label
    )

data_s1$ann_box$hr_label[1] = ""
data_s1$ann_box$hr_label[2] = ""
# Scenario 2 (single HR label for both post phases)
data_s2$ann_box <- make_ann_box_from_hr_tbl(
    rates_draws = data_s2$rates_draws,
    hr_tbl      = data_s2$irr_tbl_DiD,
    IRR_overall_did_summary = data_s2$IRR_overall_did_summary,
    upper_cap   = 25,
    collapse_post = TRUE            # <- one label over Post-Intervention 2
  )
  
data_s2$ann_box$hr_label[1] = ""

# pick one key per group to show in the legend
legend_breaks <- c(
    "Intervention.Intervention phase 1\n~6 months",
    "Control.Intervention phase 1\n~6 months"
  )
  
p1 <- plot_incidence_panel(data_s1, title_1, subtitle_1, upper_cap = 25) +
    scale_fill_manual(
      values = c(
        "Intervention.Baseline" = "#FEE5D9",
        "Intervention.Early intervention phase" = "#CB181D",
        "Intervention.Late intervention phase" = "#CB181D",
        "Intervention.Intervention overall"    = "#CB181D",
        "Control.Baseline" = "#F0F0F0",
        "Control.Early intervention phase" = "#BDBDBD",
        "Control.Late intervention phase" = "#BDBDBD",
        "Control.Intervention overall"    = "#BDBDBD"
      ),
      breaks = legend_breaks,                                   
      labels = c("Intervention", "Control"), 
      name   = NULL
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "bottom")
p1  

p2 <- plot_incidence_panel(data_s2, title_2, subtitle_2, upper_cap = 25) +
    scale_fill_manual(
      values = c(
        "Intervention.Baseline" = "#FEE5D9",
        "Intervention.Early intervention phase" = "#CB181D",
        "Intervention.Late intervention phase" = "#CB181D",
        "Intervention.Intervention overall"    = "#CB181D",
        "Control.Baseline" = "#F0F0F0",
        "Control.Early intervention phase" = "#BDBDBD",
        "Control.Late intervention phase" = "#BDBDBD",
        "Control.Intervention overall"    = "#BDBDBD"
      ),
      breaks = legend_breaks,                                   # <-- show only 2 legend entries
      labels = c("Intervention", "Control"), 
      name   = NULL
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "bottom",  # legend only from p1
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  
p_both <-  p1 + p2 + patchwork::plot_layout(ncol = 2)
p_both

# STORE Figure 2
#----------------------------------------------------------------------------------
ggsave(
    filename = "./Output/Figures_and_tables/Paper/Figure2_Intervention_effect_baseline.png",
    plot     = p_both,
    width    = 15,   # in inches
    height   = 6,    # in inches
    dpi      = 300,   
    bg = "white"
  )

# PLOT SEASONAL PATTERN
#------------------------------------------------------------

# LOAD IN BASELINE SCENARIO
fit   <- readRDS(results_path_1)
D <- as_draws_df(fit)

# Baseline acquisition rate
median_acq = median(exp(D$q_1_2_base)); 
ci_acq = quantile(exp(D$q_1_2_base), probs=c(0.025, 0.975)); 

ratio = median_acq/0.0026
ratio

median(exp(D$q_2_1_base))

# Number of individuals with four observations
dfls0complete %>% group_by(intervention.text) %>%
  summarise(n = length(unique(menage_id_member))
  )

data_source = "observed"
global_interval_start <- ifelse(data_source=="simulated", as.Date("2022-10-01"), as.Date("2022-10-03"))
global_interval_end   <- as.Date("2024-02-19")
interval_length <- 28

# Convert to numeric
global_interval_start_numeric <- as.numeric(global_interval_start)
global_interval_end_numeric   <- as.numeric(global_interval_end)

# Number of intervals
num_intervals <- ceiling((global_interval_end_numeric - global_interval_start_numeric) / interval_length)

# Interval starts/ends
interval_starts <- seq(global_interval_start_numeric,
                       global_interval_start_numeric + (num_intervals - 1) * interval_length,
                       by = interval_length)
interval_ends <- interval_starts + (interval_length - 1)

# Midpoints (special rule for short last interval)
X_midpoints <- (interval_starts + pmin(interval_starts + interval_length - 1, global_interval_end_numeric)) / 2
if ((global_interval_end_numeric - interval_starts[num_intervals]) + 1 < interval_length) {
  X_midpoints[num_intervals] <- interval_starts[num_intervals]
}

phase_grid <- seq(0, 1, length.out = 400)  # smooth one-year grid

# One-year calendar helper (anchors plotting year to Jan 1 of year_anchor)
make_calendar <- function(start_date, end_date, interval_len) {
  gs <- as.numeric(start_date); ge <- as.numeric(end_date)
  n  <- ceiling((ge - gs) / interval_len)
  starts <- seq(gs, gs + (n - 1) * interval_len, by = interval_len)
  mids <- (starts + pmin(starts + interval_len - 1, ge)) / 2
  if ((ge - starts[n]) + 1 < interval_len) mids[n] <- starts[n]
  year_anchor <- 2023
  date0 <- as.Date(paste0(year_anchor, "-01-01"))
  list(midpoints = mids, date0 = date0)
}
cal <- make_calendar(global_interval_start, global_interval_end, interval_length)

summarise_draws <- function(mat, probs = c(.05, .5, .95)) {
  qs <- apply(mat, 2, quantile, probs = probs)
  tibble(q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,])
}

# Plot Figure 3 - season (HR vs yearly median = 1)
#------------------------------------------------------------
# Rainy season bands (adjust as needed)
rainy <- data.frame(start = as.Date("2023-06-01"), end = as.Date("2023-10-31"))

# Chose sine or spline
seasonality= "sine"

`%||%` <- function(a, b) if (is.null(a)) b else a

# Phase/date function
year_len <- 365  # must match Stan
phase_to_date <- function(phase, date0, year_len = year_len) date0 + round(phase * year_len)
date_to_phase <- function(d, origin, year_len = year_len) ((as.numeric(d) - as.numeric(origin)) / year_len) %% 1
plot_phase_curve <- function(
    df,                       # columns: Date, q_low, q50, q_hi (on LOG scale)
    title,
    subtitle,
    ylab = "Hazard ratio (reference = yearly median)",
    rainy = NULL,             # data.frame with columns start, end (Date)
    ref = c("mean","median"),
    show_peak_trough = TRUE,
    show_hr1 = TRUE,          # draw HR=1 reference
    year_median_qs = NULL,    # optional LOG-scale c(0.025,0.5,0.975)
    peak_qs = NULL,           # optional LOG-scale c(0.025,0.5,0.975)
    trough_qs = NULL,         # optional LOG-scale c(0.025,0.5,0.975)
    base_size = 18
) {
  ref <- match.arg(ref)
  fmt <- function(x) formatC(x, format = "f", digits = 2, drop0trailing = TRUE)
  
  # HR for plotting (normal scale)
  dfp <- df
  dfp$q_low <- exp(dfp$q_low); dfp$q50 <- exp(dfp$q50); dfp$q_hi <- exp(dfp$q_hi)
  if (!is.null(year_median_qs)) year_median_qs <- exp(year_median_qs)
  if (!is.null(peak_qs))        peak_qs        <- exp(peak_qs)
  if (!is.null(trough_qs))      trough_qs      <- exp(trough_qs)
  
  # Base layers
  p <- ggplot(dfp, aes(Date, q50)) +
    { if (!is.null(rainy))
      geom_rect(
        data = rainy,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = "forestgreen", alpha = 0.10
      )
    } +
    geom_ribbon(aes(ymin = q_low, ymax = q_hi), alpha = 0.25, fill = "cyan4") +
    geom_line(linewidth = 2, color = "cyan4") +
    { if (show_hr1) geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.9) }
  
  # Dynamic axis limits with extra space BELOW 0 so the label is visible
  ymax <- max(dfp$q_hi, na.rm = TRUE)
  ymin <- min(dfp$q_low, na.rm = TRUE)
  
  # how much space to show below zero (both relative and absolute cushions)
  bottom_frac <- 0.10                 # 10% of top scale as room below 0
  bottom_abs  <- 0.08                 # or at least 0.08 HR units
  bottom_pad  <- max(bottom_frac * ymax, bottom_abs)
  
  # lower limit goes BELOW 0 to create visual breathing room
  lower_limit <- min(0 - bottom_pad, ymin * 0.92)
  upper_limit <- ymax * 1.12
  
  # rainy label: bottom-left of (visible) band; never clipped
  if (!is.null(rainy) && nrow(rainy) > 0) {
    xmin <- min(dfp$Date, na.rm = TRUE)
    xmax <- max(dfp$Date, na.rm = TRUE)
    rainy2 <- subset(rainy, end >= xmin & start <= xmax)
    if (nrow(rainy2) > 0) {
      y_lab <- lower_limit + 0.02 * (upper_limit - lower_limit)  # 6% above bottom
      rainy_labs <- transform(rainy2, x_lab = pmax(start, xmin))
      p <- p + geom_text(
        data = rainy_labs,
        aes(x = x_lab, y = y_lab, label = "Rainy season"),
        inherit.aes = FALSE,
        hjust = 0, vjust = 0,
        color = "darkgreen", fontface = "bold", size = base_size/3.2
      )
    }
  }
  
  # peak & trough labels (never overlap the curve)
  if (show_peak_trough) {
    # Peak ABOVE ribbon
    pk_ix   <- which.max(dfp$q50)
    pk_date <- dfp$Date[pk_ix]
    pk_top  <- dfp$q_hi[pk_ix] * 1.04
    pk_txt  <- if (!is.null(peak_qs)) {
      sprintf("Peak month: %s\nHR = %s [%s–%s]",
              format(pk_date, "%b"), fmt(peak_qs[2]), fmt(peak_qs[1]), fmt(peak_qs[3]))
    } else {
      sprintf("Peak month: %s\nHR = %s", format(pk_date, "%b"), fmt(dfp$q50[pk_ix]))
    }
    p <- p +
      geom_vline(xintercept = pk_date, linetype = "dotdash",
                 col = "darkred", linewidth = 0.9, alpha = 0.7) +
      annotate("label", x = pk_date, y = pk_top,
               label = pk_txt, vjust = 0, hjust = 0.5,
               label.size = 0, alpha = 0.95, color = "black", size = base_size/3.2)
    
    # Trough BELOW ribbon but safely above axis
    tr_ix   <- which.min(dfp$q50)
    tr_date <- dfp$Date[tr_ix]
    
    # safe floor a bit above the bottom axis
    safe_floor <- lower_limit + 0.08 * (upper_limit - lower_limit)  # 8% above bottom
    tr_bot     <- max(dfp$q_low[tr_ix] * 0.96, safe_floor)
    
    tr_txt  <- if (!is.null(trough_qs)) {
      sprintf("Trough month: %s\nHR = %s [%s–%s]",
              format(tr_date, "%b"), fmt(trough_qs[2]), fmt(trough_qs[1]), fmt(trough_qs[3]))
    } else {
      sprintf("Trough month: %s\nHR = %s", format(tr_date, "%b"), fmt(dfp$q50[tr_ix]))
    }
    p <- p +
      geom_vline(xintercept = tr_date, linetype = "dotdash", alpha = 0.7,
                 col = "darkred", linewidth = 0.9) +
      annotate("label", x = tr_date, y = tr_bot,
               label = tr_txt, vjust = 1, hjust = 0.5,   # TOP of label at tr_bot
               label.size = 0, alpha = 0.95, color = "black", size = base_size/3.2)
  }
  
  # Axis & theme (allow negative bottom for space; prevent clipping)
  p +
    scale_x_date(date_labels = "%b", date_breaks = "1 month") +
    scale_y_continuous(limits = c(lower_limit, upper_limit),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Calendar month",
      y = ylab,
      title = title,
      subtitle = subtitle,
      caption = ""
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(size = base_size + 4, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 2),
      plot.caption  = element_text(size = base_size - 3),
      panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
      panel.grid.minor = element_line(linetype = "dotted", color = "grey85"),
      plot.margin = grid::unit(c(10, 10, 18, 10), "pt")
    )
}

#-----------------------------------------------------------------
# Now plot seasonal patterns, choose branche (spline / sine)
#-----------------------------------------------------------------

if (trimws(seasonality) == "spline_") {
  
  # Spline 
  K_fit   <- num_knots
  deg_fit <- spline_degree
  knots_vec <- knots
  stopifnot(length(knots_vec) == K_fit)
  nbasis_fit <- K_fit + deg_fit
  
  # Stan-matching periodic B-spline basis
  build_b_spline_vec <- function(t, ext_knots, ind, order) {
    if (order == 1) {
      as.numeric(ext_knots[ind] <= t & t < ext_knots[ind + 1])
    } else {
      d1 <- ext_knots[ind + order - 1] - ext_knots[ind]
      d2 <- ext_knots[ind + order]     - ext_knots[ind + 1]
      w1 <- if (d1 > 0) (t - ext_knots[ind]) / d1 else 0
      w2 <- if (d2 > 0) 1 - (t - ext_knots[ind + 1]) / d2 else 0
      w1 * build_b_spline_vec(t, ext_knots, ind, order - 1) +
        w2 * build_b_spline_vec(t, ext_knots, ind + 1, order - 1)
    }
  }
  
  make_B_matrix <- function(X_in, knots, spline_degree) {
    X_mod <- pmin(X_in %% 1, 1 - 1e-12)
    K  <- length(knots)
    nb <- K + spline_degree
    ext_knots <- c(knots[(K - spline_degree + 1):K] - 1,
                   knots,
                   knots[1:(spline_degree + 1)] + 1)
    order <- spline_degree + 1
    B <- sapply(1:nb, function(ind) build_b_spline_vec(X_mod, ext_knots, ind, order))
    as.matrix(B)
  }
  
  # Build basis on the phase grid
  B_phase <- make_B_matrix(phase_grid, knots_vec, deg_fit)
  
  vars_a  <- sprintf("a_1_2[%d]", seq_len(nbasis_fit))
  A_mat   <- as.matrix(D[, vars_a])  # draws × nbasis
  stopifnot(ncol(A_mat) == ncol(B_phase))
  
  # draws × n_phase (LOG scale)
  Y_log <- A_mat %*% t(B_phase)
  
  # Center per draw, yearly median HR = 1 after exp()
  row_med <- apply(Y_log, 1, median)
  Y_log   <- sweep(Y_log, 1, row_med, "-")
  
  # Timewise 95% CrIs (LOG)
  qs <- apply(Y_log, 2, quantile, probs = c(.05, .5, .95))
  spline_df <- tibble(
    phase = phase_grid,
    q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,]
  ) |>
    mutate(Date = cal$date0 + round(phase * year_len))
  
  # Peak & Trough across draws (LOG) for annotations
  peak_draws   <- apply(Y_log, 1, max)
  trough_draws <- apply(Y_log, 1, min)
  peak_qs      <- quantile(peak_draws,   c(.025, .5, .975), na.rm = TRUE)
  trough_qs    <- quantile(trough_draws, c(.025, .5, .975), na.rm = TRUE)
  
  p_season <- plot_phase_curve(
    df = spline_df,
    title = "Estimated seasonal variation in ESBL-E acquisition risk",
    subtitle = "Cubic spline seasonal function (median hazard ratio with 95% CrI)",
    ylab = NULL,
    rainy = rainy,
    ref = "median",
    show_peak_trough = TRUE,
    show_hr1 = TRUE,
    peak_qs = peak_qs,
    trough_qs = trough_qs
  )
  
} else {
  
  # Sine (evaluate on calendar; auto-align phi convention) 
  jan1 <- cal$date0
  date_grid <- seq(jan1, by = "1 day", length.out = year_len)
  
  # Model origin used to define phase in Stan (change to X_midpoints[1] if that’s what Stan used)
  model_origin <- as.Date(global_interval_start_numeric, origin = "1970-01-01")
  
  # Phase per calendar date
  phase_date <- ((as.numeric(date_grid) - as.numeric(model_origin)) / year_len) %% 1  # length = n_dates
  
  # Draws: amplitude (log scale) and phi
  s <- D[, c("a1","phi")]
  n_draws <- nrow(s); n_dates <- length(date_grid)
  
  # Try common phi conventions; pick the one peaking nearest August
  build_Y <- function(convention = c("sin_plus","sin_minus","cos_plus","sin_plus_cycles")) {
    convention <- match.arg(convention)
    angle_base <- matrix(2*pi*phase_date, nrow = n_dates, ncol = n_draws)
    phi_mat <- matrix(rep(s$phi, each = n_dates), nrow = n_dates, ncol = n_draws)
    a1_mat  <- matrix(rep(s$a1,  each = n_dates), nrow = n_dates, ncol = n_draws)
    Y <- switch(convention,
                sin_plus         = sin(angle_base +  phi_mat) * a1_mat,
                sin_minus        = sin(angle_base -  phi_mat) * a1_mat,
                cos_plus         = cos(angle_base +  phi_mat) * a1_mat,
                sin_plus_cycles  = { phi_rad <- 2*pi*phi_mat; sin(angle_base + phi_rad) * a1_mat }
    )
    t(Y)  # draws × dates
  }
  
  conventions <- c("sin_plus","sin_minus","cos_plus","sin_plus_cycles")
  Y_list <- lapply(conventions, build_Y)
  
  # Center per draw
  Y_list <- lapply(Y_list, function(Y) {
    rm <- apply(Y, 1, median)
    sweep(Y, 1, rm, "-")
  })
  
  # Choose convention whose median curve peaks closest to August
  target_month <- 8L
  med_curves <- lapply(Y_list, function(Y) apply(Y, 2, median))
  peak_dist <- function(dates, y_med) {
    pk <- dates[which.max(y_med)]
    pm <- as.integer(format(pk, "%m"))
    abs(((pm - target_month + 6) %% 12) - 6)
  }
  dists <- vapply(med_curves, function(m) peak_dist(date_grid, m), numeric(1))
  Y_log <- Y_list[[which.min(dists)]]
  
  # Timewise 95% CrIs (LOG)
  qs <- apply(Y_log, 2, quantile, probs = c(.05, .5, .95), na.rm = TRUE)
  sine_df <- tibble(Date = date_grid, q_low = qs[1,], q50 = qs[2,], q_hi = qs[3,])
  
  # Peak & Trough across draws (LOG) for annotations
  peak_draws   <- apply(Y_log, 1, max)
  trough_draws <- apply(Y_log, 1, min)
  peak_qs      <- quantile(peak_draws,   c(.025, .5, .975), na.rm = TRUE)
  trough_qs    <- quantile(trough_draws, c(.025, .5, .975), na.rm = TRUE)
  
  p_season <- plot_phase_curve(
    df = sine_df,
    title = "Scenario 1: Two-step intervention effect",
    subtitle = "Sinusoidal seasonal fit",
    ylab = NULL,
    rainy = rainy,
    ref = "median",
    show_peak_trough = TRUE,
    show_hr1 = TRUE,
    peak_qs = peak_qs,
    trough_qs = trough_qs
  )
}

p_season  

# STORE FIGURE 3
#-----------------------------------------------------------------------------

ggsave(
  filename = "./Output/Figures/Paper/Figure3_Seasonal_effect_baseline.png",
  plot     = p_season,
  width    = 12,   # in inches
  height   = 9,    # in inches
  dpi      = 300,   
  bg = "white"
)

# TABLES 
#-------------------------------------------------------------------------------

data <- read.csv(file = "./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")

villages <- data %>% select(menage_id, HouseID, village, village_name)
villages$VillageID = villages$village_name

data %>% 
  group_by(intervention.text) %>%
  summarise(
    n = length(unique(menage_id_member)),
    n_hh = length(unique(menage_id))
  )

#-------------------------------------------------------------------------------
# SUPPLEMENTARY MATERIAL
#-------------------------------------------------------------------------------

# individuals with at least 2 observations

# -------------------------------------------------------------------
# Figure S1: Calendar-time–adjusted per-village rates (and matched controls)
# -------------------------------------------------------------------

# Village vector
get_village_vec <- function(stan_data_fit) {
  # 1) Prefer a village/cluster field inside stan_data_fit
  cand <- c("VillageID","village_id","ClusterID","cluster_id","village","cluster","village_name")
  for (nm in cand) {
    if (!is.null(stan_data_fit[[nm]])) return(as.character(stan_data_fit[[nm]]))
  }
  
  # 2) If not present, try the 'villages' data frame (already built above)
  if (exists("villages", inherits = TRUE)) {
    vdf <- get("villages", inherits = TRUE)
    hcol <- intersect(c("HouseID","house_id","house"), names(vdf))[1]
    vcol <- intersect(c("VillageID","village_name","village","cluster"), names(vdf))[1]
    if (!is.na(hcol) && !is.na(vcol)) {
      vm <- vdf[[vcol]][match(stan_data_fit$HouseID, vdf[[hcol]])]
      if (any(is.na(vm))) {
        warning("Some HouseID values not found in 'villages' mapping; falling back to HouseID for those rows.")
        vm[is.na(vm)] <- as.character(stan_data_fit$HouseID[is.na(vm)])
      }
      return(as.character(vm))
    } else {
      warning("Could not find suitable columns in 'villages' (need HouseID and a village/cluster column).")
    }
  }
  
  # Last resort: use HouseID as a proxy
  message("No village field found; using HouseID as village proxy.")
  as.character(stan_data_fit$HouseID)
}


# Build segments exactly like Stan, now including village id
extract_segments_village <- function(stan_data_fit) {
  N  <- stan_data_fit$N
  L  <- stan_data_fit$interval_length
  g0 <- stan_data_fit$global_interval_start
  K  <- stan_data_fit$num_intervals
  village_vec <- get_village_vec(stan_data_fit)
  
  obs <- tibble(
    n        = seq_len(N),
    house    = stan_data_fit$HouseID,
    village  = village_vec,
    id       = stan_data_fit$menage_id_member,
    age      = stan_data_fit$age,
    sexe     = stan_data_fit$sexe,
    t        = stan_data_fit$date_use,
    int_flag = as.integer(stan_data_fit$intervention > 0.5),
    t_star1  = stan_data_fit$intervention_date,
    t_star2  = stan_data_fit$intervention_date2
  ) %>% arrange(id, n) %>%
    mutate(t0 = lag(t), first_obs = is.na(t0) | (id != lag(id)))
  
  clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
  overlap_len <- function(a, b, Lb, Ub) {
    if (is.na(a) || is.na(b)) return(0)
    as.numeric(max(0, min(b, Ub) - max(a, Lb) + 1))
  }
  make_one_obs_subints <- function(i) {
    if (obs$first_obs[i]) return(NULL)
    t0 <- obs$t0[i]; t1 <- obs$t[i]
    raw_start <- floor((t0 - g0) / L); raw_end <- ceiling((t1 - g0) / L)
    gi_start <- clamp(raw_start + 1, 1, K); gi_end <- clamp(raw_end, 1, K)
    a1 <- g0 + (gi_start - 1) * L; b1 <- a1 + L - 1
    first_len <- min(c(b1 - t0 + 1, t1 - t0 + 1, L))
    first <- NULL
    if (first_len > 0) {
      mid1 <- (t0 + b1)/2; idx_first <- clamp(floor((mid1 - g0)/L) + 1, 1, K)
      first <- list(list(a = t0, b = t0 + first_len - 1, idx = idx_first))
    }
    n_mid <- gi_end - gi_start - 1; mids <- NULL
    if (n_mid > 0) {
      mids <- vector("list", n_mid)
      for (m in 1:n_mid) {
        a <- g0 + (gi_start + m - 1) * L; b <- a + L - 1
        midm <- (a + b)/2; idxm <- clamp(floor((midm - g0)/L) + 1, 1, K)
        mids[[m]] <- list(a=a,b=b,idx=idxm)
      }
    }
    last_len <- min(c(t1 - (g0 + (gi_end - 1) * L) + 1, L))
    last <- NULL
    if (last_len > 0) {
      al <- t1 - last_len + 1; bl <- t1
      midl <- ((g0 + (gi_end - 1) * L) + t1)/2; idx_last <- clamp(floor((midl - g0)/L) + 1, 1, K)
      last <- list(list(a=al,b=bl,idx=idx_last))
    }
    c(first, mids, last)
  }
  
  seg_list <- vector("list", N)
  for (i in seq_len(N)) {
    if (obs$first_obs[i]) next
    subints <- make_one_obs_subints(i)
    if (is.null(subints) || !length(subints)) next
    t1 <- obs$t_star1[i]; t2 <- obs$t_star2[i]
    grp <- ifelse(obs$int_flag[i] == 1L, "Intervention", "Control")
    for (si in subints) {
      a <- si$a; b <- si$b; idx <- si$idx
      if (grp == "Intervention") {
        d0 <- overlap_len(a,b, -Inf, t1-1)
        d1 <- overlap_len(a,b, t1,   t2-1)
        d2 <- overlap_len(a,b, t2,   Inf)
        rows <- NULL
        if (d0>0) rows <- rbind(rows, c(idx, obs$house[i], obs$village[i], obs$age[i], obs$sexe[i], 0L, d0, grp))
        if (d1>0) rows <- rbind(rows, c(idx, obs$house[i], obs$village[i], obs$age[i], obs$sexe[i], 1L, d1, grp))
        if (d2>0) rows <- rbind(rows, c(idx, obs$house[i], obs$village[i], obs$age[i], obs$sexe[i], 2L, d2, grp))
        if (!is.null(rows)) seg_list[[i]] <- rbind(seg_list[[i]], rows)
      } else {
        seg_list[[i]] <- rbind(seg_list[[i]], c(idx, obs$house[i], obs$village[i], obs$age[i], obs$sexe[i], NA, b-a+1, grp))
      }
    }
  }
  segs <- do.call(rbind, seg_list) %>% as_tibble()
  names(segs) <- c("idx","house","village","age","sexe","iperiod","dur","group")
  segs %>%
    mutate(
      idx    = as.integer(idx),
      house  = as.integer(house),
      village= as.character(village),
      iperiod= suppressWarnings(as.integer(iperiod)),
      dur    = as.numeric(dur),
      group  = factor(group, levels = c("Intervention","Control"))
    )
}

# Rebuild per-row hazards from the fit (optionally drop season term)
compute_lambda_mats <- function(fit, segs, include_season = TRUE) {
  draws <- as_draws_df(fit)
  D <- nrow(draws); R <- nrow(segs)
  
  # seasonal term over global bins
  Y_cols <- grep("^Y_hat_1_2_out\\[\\d+\\]$", names(draws), value = TRUE)
  Y_idx  <- as.integer(stringr::str_match(Y_cols, "\\[(\\d+)\\]")[,2]); Y_cols <- Y_cols[order(Y_idx)]
  Y_mat  <- if (length(Y_cols)) as.matrix(draws[, Y_cols, drop = FALSE]) else matrix(0, D, max(segs$idx, na.rm=TRUE))
  
  # seasonal term over global bins
  Y_cols <- grep("^Y_hat_1_2_out\\[\\d+\\]$", names(draws), value = TRUE)
  Y_idx  <- as.integer(str_match(Y_cols, "\\[(\\d+)\\]")[,2]); Y_cols <- Y_cols[order(Y_idx)]
  Y_mat  <- if (length(Y_cols)) as.matrix(draws[, Y_cols, drop = FALSE]) else matrix(0, D, max(segs$idx, na.rm=TRUE))
  
  # index consistency check (prevents out-of-bounds) 
  if (include_season && length(Y_cols)) {
    K_fit <- ncol(Y_mat)
    if (max(segs$idx, na.rm = TRUE) > K_fit) {
      stop(
        "Calendar-bin index in 'segs' exceeds seasonal vector length from the fit.\n",
        "Rebuild segments with the *same* interval grid as the model used:\n",
        "- global_interval_start\n",
        "- interval_length\n",
        "- num_intervals"
      )
    }
  }
  
  # house REs
  u_cols <- grep("^u\\[\\d+\\]$", names(draws), value = TRUE)
  u_idx  <- as.integer(str_match(u_cols, "\\[(\\d+)\\]")[,2]); u_cols <- u_cols[order(u_idx)]
  U_mat  <- as.matrix(draws[, u_cols, drop = FALSE])
  
  q12_base <- draws[["q_1_2_base"]]
  b_age    <- draws[["beta_1_2_age"]]
  b_sex    <- draws[["beta_1_2_sexe"]]
  d_int1   <- draws[["beta_int1_1"]]
  d_int2   <- draws[["beta_int1_2"]]
  
  # assemble LP
  Y_by_row <- if (include_season) Y_mat[, segs$idx, drop = FALSE] else matrix(0, D, R)
  U_by_row <- U_mat[, segs$house, drop = FALSE]
  AGE_by_row <- matrix(as.numeric(segs$age),  nrow = D, ncol = R, byrow = TRUE)
  SEX_by_row <- matrix(as.numeric(segs$sexe), nrow = D, ncol = R, byrow = TRUE)
  BASE       <- matrix(q12_base, nrow = D, ncol = R)
  
  LP <- BASE + Y_by_row + U_by_row + (b_age * AGE_by_row) + (b_sex * SEX_by_row)
  
  # intervention effects on treated rows
  int_eff <- matrix(0, D, R)
  ix_I1 <- which(segs$group=="Intervention" & segs$iperiod==1L)
  ix_I2 <- which(segs$group=="Intervention" & segs$iperiod==2L)
  if (length(ix_I1)) int_eff[, ix_I1] <- matrix(d_int1, D, length(ix_I1))
  if (length(ix_I2)) int_eff[, ix_I2] <- matrix(d_int2, D, length(ix_I2))
  
  lambda <- exp(LP + int_eff)                   # per-day hazard per row, per draw
  Dur    <- matrix(segs$dur, nrow = D, ncol = R, byrow = TRUE)
  list(draws = draws, lambda = lambda, Dur = Dur)
}

# Per-village calendar-matched rates (pooled for control villages by calendar time)
# per_village_calendar_adjusted <- function(fit, stan_data_fit, include_season = TRUE, scale = 100) {
#   segs <- extract_segments_village(stan_data_fit)
#   mats <- compute_lambda_mats(fit, segs, include_season = include_season)
#   lambda <- mats$lambda; Dur <- mats$Dur; D <- nrow(lambda)
#   
#   sum_w <- function(mat, idx) if (length(idx)) rowSums(mat[, idx, drop=FALSE]) else rep(NA_real_, D)
#   rate_from_rows <- function(idx) {
#     if (!length(idx)) return(rep(NA_real_, D))
#     num <- sum_w(lambda * Dur, idx)
#     den <- sum(segs$dur[idx])
#     num / den
#   }
#   
#   # control hazards per calendar bin (across all control rows)
#   ix_C_all <- which(segs$group=="Control")
#   ctrl_bins <- sort(unique(segs$idx[ix_C_all]))
#   H_C <- sapply(ctrl_bins, function(j) {
#     jj <- ix_C_all[segs$idx[ix_C_all] == j]
#     num <- sum_w(lambda * Dur, jj)
#     den <- sum(segs$dur[jj])
#     num / den
#   })
#   if (is.null(dim(H_C))) H_C <- matrix(H_C, nrow = D)
#   colnames(H_C) <- as.character(ctrl_bins)
#   
#   weights_from_rows <- function(idx) {
#     tab <- tapply(segs$dur[idx], segs$idx[idx], sum)
#     idxs <- as.integer(names(tab)); w <- as.numeric(tab) / sum(tab)
#     list(idx = idxs, w = w)
#   }
#   std_rate <- function(wlist) {
#     common <- intersect(colnames(H_C), as.character(wlist$idx))
#     if (!length(common)) return(rep(NA_real_, D))
#     H <- H_C[, match(common, colnames(H_C)), drop = FALSE]
#     w <- wlist$w[match(as.integer(common), wlist$idx)]
#     w <- w / sum(w)
#     as.vector(H %*% w)
#   }
#   
#   qci <- function(x) c(med = median(x, na.rm=TRUE),
#                        "lo.2.5%" = quantile(x, 0.025, na.rm=TRUE),
#                        "hi.97.5%"= quantile(x, 0.975, na.rm=TRUE))
#   
#   # villages with intervention person-time
#   vill_int <- segs %>% filter(group=="Intervention") %>% distinct(village) %>% pull(village)
#   
#   # per-village outputs
#   out_rates_draws <- list()
#   out_rates_summ  <- list()
#   
#   for (v in vill_int) {
#     ix_v0 <- which(segs$group=="Intervention" & segs$village==v & segs$iperiod==0L)
#     ix_v1 <- which(segs$group=="Intervention" & segs$village==v & segs$iperiod==1L)
#     ix_v2 <- which(segs$group=="Intervention" & segs$village==v & segs$iperiod==2L)
#     
#     R0_I <- rate_from_rows(ix_v0); R1_I <- rate_from_rows(ix_v1); R2_I <- rate_from_rows(ix_v2)
#     
#     w0 <- weights_from_rows(ix_v0); w1 <- weights_from_rows(ix_v1); w2 <- weights_from_rows(ix_v2)
#     R0_C <- std_rate(w0); R1_C <- std_rate(w1); R2_C <- std_rate(w2)
#     
#     # per-draw long (for boxplots)
#     out_rates_draws[[v]] <- bind_rows(
#       tibble(village=v, Group="Intervention", phase_label="Baseline",      ID12_per100 = scale * R0_I),
#       tibble(village=v, Group="Intervention", phase_label="Phase 1",       ID12_per100 = scale * R1_I),
#       tibble(village=v, Group="Intervention", phase_label="Phase 2",       ID12_per100 = scale * R2_I),
#       tibble(village=v, Group="Control (matched)", phase_label="Baseline", ID12_per100 = scale * R0_C),
#       tibble(village=v, Group="Control (matched)", phase_label="Phase 1",  ID12_per100 = scale * R1_C),
#       tibble(village=v, Group="Control (matched)", phase_label="Phase 2",  ID12_per100 = scale * R2_C)
#     ) %>% filter(is.finite(ID12_per100))
#     
#     # summaries
#     out_rates_summ[[v]] <- bind_rows(
#       tibble(village=v, Group="Intervention",      phase_label="Baseline") %>% bind_cols(as.list(qci(scale*R0_I))),
#       tibble(village=v, Group="Intervention",      phase_label="Phase 1")  %>% bind_cols(as.list(qci(scale*R1_I))),
#       tibble(village=v, Group="Intervention",      phase_label="Phase 2")  %>% bind_cols(as.list(qci(scale*R2_I))),
#       tibble(village=v, Group="Control (matched)", phase_label="Baseline") %>% bind_cols(as.list(qci(scale*R0_C))),
#       tibble(village=v, Group="Control (matched)", phase_label="Phase 1")  %>% bind_cols(as.list(qci(scale*R1_C))),
#       tibble(village=v, Group="Control (matched)", phase_label="Phase 2")  %>% bind_cols(as.list(qci(scale*R2_C)))
#     )
#   }
#   
#   rates_draws_village <- bind_rows(out_rates_draws) %>%
#     mutate(
#       Group = factor(Group, levels = c("Intervention","Control (matched)")),
#       phase_label = factor(phase_label, levels = c("Baseline","Phase 1","Phase 2"))
#     )
#   
#   rates_summary_village <- bind_rows(out_rates_summ) %>%
#     #rename(Median = med, Lower = `lo.2.5%`, Upper = `hi.97.5%`) %>%
#     mutate(
#       Group = factor(Group, levels = c("Intervention","Control (matched)")),
#       phase_label = factor(phase_label, levels = c("Baseline","Phase 1","Phase 2"))
#     )
#   
#   list(segs = segs,
#        rates_draws_village = rates_draws_village,
#        rates_summary_village = rates_summary_village)
# }

# Village-explicit acquisition with or without pooled controls 
per_village_calendar_adjusted <- function(
    fit,
    stan_data_fit,
    include_season = TRUE,
    scale = 100,
    control_strategy = c("none","pooled_matched","per_village_self")
) {
  control_strategy <- match.arg(control_strategy)
  
  # Build Stan-like segments with village id
  segs <- extract_segments_village(stan_data_fit)
  
  # Map draws -> per-row hazards
  mats <- compute_lambda_mats(fit, segs, include_season = include_season)
  lambda <- mats$lambda
  Dur    <- mats$Dur
  D      <- nrow(lambda)
  
  sum_w <- function(mat, idx) if (length(idx)) rowSums(mat[, idx, drop = FALSE]) else rep(NA_real_, D)
  rate_from_rows <- function(idx) {
    if (!length(idx)) return(rep(NA_real_, D))
    num <- sum_w(lambda * Dur, idx)
    den <- sum(segs$dur[idx])
    num / den
  }
  
  # Function for CrIs
  qci <- function(x) c(med = stats::median(x, na.rm = TRUE),
                       "lo.2.5%" = quantile(x, 0.025, na.rm = TRUE),
                       "hi.97.5%"= quantile(x, 0.975, na.rm = TRUE))
  
 
  # INTERVENTION villages ONLY
  # -------------------------
  vill_int <- segs %>%
    dplyr::filter(group == "Intervention") %>%
    dplyr::distinct(village) %>%
    dplyr::pull(village)
  
  out_rates_draws <- vector("list", length(vill_int))
  names(out_rates_draws) <- vill_int
  out_rates_summ  <- vector("list", length(vill_int))
  names(out_rates_summ)  <- vill_int
  
  for (v in vill_int) {
    # rows for this village & phase (as defined by your Stan intervention split)
    ix_v0 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 0L)
    ix_v1 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 1L)
    ix_v2 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 2L)
    
    R0_I <- rate_from_rows(ix_v0)
    R1_I <- rate_from_rows(ix_v1)
    R2_I <- rate_from_rows(ix_v2)
    
    # Per-draw long (village-explicit, no controls)
    out_rates_draws[[v]] <- dplyr::bind_rows(
      tibble(village = v, Group = "Intervention", phase_label = "Baseline", ID12_per100 = scale * R0_I),
      tibble(village = v, Group = "Intervention", phase_label = "Phase 1",  ID12_per100 = scale * R1_I),
      tibble(village = v, Group = "Intervention", phase_label = "Phase 2",  ID12_per100 = scale * R2_I)
    ) %>% filter(is.finite(ID12_per100))
    
    # Summaries (village × phase)
    out_rates_summ[[v]] <- dplyr::bind_rows(
      tibble(village = v, Group = "Intervention", phase_label = "Baseline") %>%
        bind_cols(as.list(qci(scale * R0_I))),
      tibble(village = v, Group = "Intervention", phase_label = "Phase 1")  %>%
        bind_cols(as.list(qci(scale * R1_I))),
      tibble(village = v, Group = "Intervention", phase_label = "Phase 2")  %>%
        bind_cols(as.list(qci(scale * R2_I)))
    )
  }
  
  rates_draws_village <- dplyr::bind_rows(out_rates_draws) %>%
    mutate(
      Group = factor(Group, levels = c("Intervention")),
      phase_label = factor(phase_label, levels = c("Baseline","Phase 1","Phase 2"))
    )
  
  # add posterior variance directly here (one row per village × phase)
  rates_summary_village <- rates_draws_village %>%
    group_by(village, Group, phase_label) %>%
    dplyr::summarise(
      mean = mean(ID12_per100, na.rm = TRUE),
      median = stats::median(ID12_per100, na.rm = TRUE),
      lo = stats::quantile(ID12_per100, 0.025, na.rm = TRUE),
      hi = stats::quantile(ID12_per100, 0.975, na.rm = TRUE),
      var = stats::var(ID12_per100, na.rm = TRUE),
      n_draws = n(),
      .groups = "drop"
    )
  
  # -------------------------
  # Control handling
  # -------------------------
  if (control_strategy != "none") {
    
    if (control_strategy == "pooled_matched") {
      ix_C_all <- which(segs$group == "Control")
      ctrl_bins <- sort(unique(segs$idx[ix_C_all]))
      H_C <- sapply(ctrl_bins, function(j) {
        jj <- ix_C_all[segs$idx[ix_C_all] == j]
        num <- sum_w(lambda * Dur, jj)
        den <- sum(segs$dur[jj])
        num / den
      })
      if (is.null(dim(H_C))) H_C <- matrix(H_C, nrow = D)
      colnames(H_C) <- as.character(ctrl_bins)
      
      weights_from_rows <- function(idx) {
        tab <- tapply(segs$dur[idx], segs$idx[idx], sum)
        idxs <- as.integer(names(tab)); w <- as.numeric(tab) / sum(tab)
        list(idx = idxs, w = w)
      }
      std_rate <- function(wlist) {
        common <- intersect(colnames(H_C), as.character(wlist$idx))
        if (!length(common)) return(rep(NA_real_, D))
        H <- H_C[, match(common, colnames(H_C)), drop = FALSE]
        w <- wlist$w[match(as.integer(common), wlist$idx)]
        w <- w / sum(w)
        as.vector(H %*% w)
      }
      
      # Add matched pooled-controls per village/phase
      add_ctrl <- lapply(vill_int, function(v) {
        ix_v0 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 0L)
        ix_v1 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 1L)
        ix_v2 <- which(segs$group == "Intervention" & segs$village == v & segs$iperiod == 2L)
        w0 <- weights_from_rows(ix_v0); w1 <- weights_from_rows(ix_v1); w2 <- weights_from_rows(ix_v2)
        R0_C <- std_rate(w0); R1_C <- std_rate(w1); R2_C <- std_rate(w2)
        
        bind_rows(
          tibble(village = v, Group = "Control (matched)", phase_label = "Baseline", ID12_per100 = scale * R0_C),
          tibble(village = v, Group = "Control (matched)", phase_label = "Phase 1",  ID12_per100 = scale * R1_C),
          tibble(village = v, Group = "Control (matched)", phase_label = "Phase 2",  ID12_per100 = scale * R2_C)
        )
      })
      ctrl_draws <- dplyr::bind_rows(add_ctrl)
      
      rates_draws_village   <- dplyr::bind_rows(rates_draws_village, ctrl_draws)
      rates_summary_ctrl <- ctrl_draws %>%
        group_by(village, Group, phase_label) %>%
        dplyr::summarise(
          mean = mean(ID12_per100, na.rm = TRUE),
          median = median(ID12_per100, na.rm = TRUE),
          lo = quantile(ID12_per100, 0.025, na.rm = TRUE),
          hi = quantile(ID12_per100, 0.975, na.rm = TRUE),
          var = var(ID12_per100, na.rm = TRUE),
          n_draws = dplyr::n(),
          .groups = "drop"
        )
      rates_summary_village <- dplyr::bind_rows(rates_summary_village, rates_summary_ctrl)
    }
    
    if (control_strategy == "per_village_self") {
      .date_to_bin <- function(t, g0, L) as.integer(floor((t - g0) / L) + 1L)
      clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
      add_months <- function(t, m) t + round(30.4375 * m)
      
      g0 <- stan_data_fit$global_interval_start
      L  <- stan_data_fit$interval_length
      K  <- stan_data_fit$num_intervals
      
      # Map House -> (village, group)
      hg_map <- segs %>%
        distinct(house, village, group) %>%
        dplyr::rename(HouseID = house)
      
      # Pull observation dates
      date_vec <- if (!is.null(stan_data_fit$date_use)) {
        as.numeric(stan_data_fit$date_use)
      } else if (!is.null(stan_data_fit$date.use)) {
        as.numeric(stan_data_fit$date.use)
      } else {
        stop("stan_data_fit must contain 'date_use' (or 'date.use') to build control phases.")
      }
      
      obs_df <- tibble::tibble(
        HouseID = as.integer(stan_data_fit$HouseID),
        t = as.numeric(date_vec)
      ) %>%
        inner_join(hg_map, by = "HouseID") %>%
        filter(group == "Control", is.finite(t))
      
      # First & second sampling time per house (drop houses w/o a second)
      house_1st2nd <- obs_df %>%
        group_by(HouseID, village) %>%
        dplyr::summarise(
          t_first  = min(t, na.rm = TRUE),
          t_second = suppressWarnings(dplyr::nth(sort(unique(t)), 2)),
          .groups = "drop"
        ) %>%
        filter(is.finite(t_first), is.finite(t_second), t_second > t_first)
      
      # Median anchors per village
      piv <- house_1st2nd %>%
        group_by(village) %>%
        dplyr::summarise(
          t_first_med  = stats::median(t_first,  na.rm = TRUE),
          t_second_med = stats::median(t_second, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          # Baseline: [first, second)
          b_first  = clamp(.date_to_bin(t_first_med,  g0, L), 1L, K),
          b_second = clamp(.date_to_bin(t_second_med, g0, L), 1L, K),
          # Phase 1: [second, second + 6m)
          b_second_p6m = clamp(.date_to_bin(add_months(t_second_med, 6), g0, L), 1L, K),
          # Phase 2: [second + 6m, second + 9m)
          b_second_p9m = clamp(.date_to_bin(add_months(t_second_med, 9), g0, L), 1L, K)
        ) %>%
        filter(b_second > b_first)   # need a positive baseline window
      
      # compute calendar-weighted rate for one village & bin window
      rate_self_in_bins <- function(v, b_start, b_end) {
        if (!is.finite(b_start) || !is.finite(b_end) || b_end <= b_start) return(rep(NA_real_, D))
        idx_v <- which(segs$group == "Control" & segs$village == v & segs$idx >= b_start & segs$idx < b_end)
        if (!length(idx_v)) return(rep(NA_real_, D))
        
        # person-time weights by calendar bin within the window
        tab <- tapply(segs$dur[idx_v], segs$idx[idx_v], sum)
        bins <- as.integer(names(tab))
        w    <- as.numeric(tab) / sum(tab)
        
        # per-bin hazards for this village & window
        H <- sapply(bins, function(j) {
          jj <- idx_v[segs$idx[idx_v] == j]
          num <- sum_w(lambda * Dur, jj)
          den <- sum(segs$dur[jj])
          num / den
        })
        if (is.null(dim(H))) H <- matrix(H, nrow = D)
        colnames(H) <- as.character(bins)
        
        as.vector(H %*% w)
      }
      
      # Build draws for all control villages (Baseline, Phase 1, Phase 2)
      ctrl_phase_draws <- lapply(seq_len(nrow(piv)), function(i) {
        v  <- piv$village[i]
        R_B  <- rate_self_in_bins(v, piv$b_first[i],       piv$b_second[i])
        R_P1 <- rate_self_in_bins(v, piv$b_second[i],      piv$b_second_p6m[i])
        R_P2 <- rate_self_in_bins(v, piv$b_second_p6m[i],  piv$b_second_p9m[i])
        
        bind_rows(
          tibble(village = v, Group = "Control", phase_label = "Baseline", ID12_per100 = scale * R_B),
          tibble(village = v, Group = "Control", phase_label = "Phase 1",  ID12_per100 = scale * R_P1),
          tibble(village = v, Group = "Control", phase_label = "Phase 2",  ID12_per100 = scale * R_P2)
        ) %>% filter(is.finite(ID12_per100))
      })
      ctrl_draws <- bind_rows(ctrl_phase_draws)
      
      # bind to outputs and summarise
      rates_draws_village <- bind_rows(rates_draws_village, ctrl_draws) %>%
        mutate(
          Group = factor(Group, levels = c("Intervention","Control")),
          phase_label = factor(phase_label, levels = c("Baseline","Phase 1","Phase 2"))
        )
      
      rates_summary_ctrl <- ctrl_draws %>%
        group_by(village, Group, phase_label) %>%
        dplyr::summarise(
          mean   = mean(ID12_per100, na.rm = TRUE),
          median = median(ID12_per100, na.rm = TRUE),
          lo     = quantile(ID12_per100, 0.025, na.rm = TRUE),
          hi     = quantile(ID12_per100, 0.975, na.rm = TRUE),
          var    = var(ID12_per100, na.rm = TRUE),
          n_draws = n(),
          .groups = "drop"
        )
      
      rates_summary_village <- dplyr::bind_rows(rates_summary_village, rates_summary_ctrl)
    }
  }
  
  list(
    segs = segs,
    rates_draws_village = rates_draws_village,
    rates_summary_village = rates_summary_village
  )
}


# Function for a compact plot for selected villages (horizontal, with coral vs grey)
plot_village_rates <- function(rates_draws_village,lab, villages = NULL, upper_cap = 50) {
  df <- rates_draws_village
  if (!is.null(villages)) df <- df %>% filter(village %in% villages)
  
  fill_vals <- c(
    "Intervention.Baseline"      = "#FEE5D9",
    "Intervention.Phase 1"       = "#CB181D",
    "Intervention.Phase 2"       = "#CB181D",
    "Control (matched).Baseline" = "#F0F0F0",
    "Control (matched).Phase 1"  = "#BDBDBD",
    "Control (matched).Phase 2"  = "#BDBDBD"
  )
  
  df <- df %>%
    mutate(
      fill_key = interaction(Group, phase_label, drop = TRUE, lex.order = TRUE),
      village  = factor(village)  # keep your own order, or reorder by start time if you add it
    )
  
  ggplot(df, aes(x = phase_label, y = ID12_per100, fill = fill_key)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    coord_flip(ylim = c(0, upper_cap), clip = "off") +
    scale_y_continuous(breaks = seq(0, upper_cap, 5),
                       expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = NULL, y = "Incidence per 100 person-days",
         title = "Baseline scenario 1",
         subtitle = "Village-level intervention effect (calender-time matched control)",) +
    geom_text(data = lab,
              aes(x = phase_label, y = 30, label = label),
              inherit.aes = FALSE, hjust = 2, size = 3.5)+
    facet_wrap(~ village, ncol = 3) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"),
          plot.title    = element_text(size = base_size + 4, face = "bold"),
          plot.subtitle = element_text(size = base_size),
          axis.title    = element_text(size = base_size),
          axis.text     = element_text(size = base_size - 2),
          plot.caption  = element_text(size = base_size - 3),
          panel.grid.major = element_line(linetype = "dotted", color = "grey70"),
          panel.grid.minor = element_line(linetype = "dotted", color = "grey85"),
          plot.margin = grid::unit(c(10, 10, 18, 10), "pt"))
}

# DiD IRR = (I_Ph2 / I_Base) / (C_Ph2 / C_Base) per village
village_did_labels <- function(rates_draws_village, upper_cap = 25) {
  # 1) Make a per-phase per-group draw index so we can pivot wide
  rates_draws_village <- rates_draws_village %>%
    group_by(village, Group, phase_label) %>%
    mutate(draw = row_number()) %>%
    ungroup()
  
  # Map group + phase into compact codes and go wide
  rd <- rates_draws_village %>%
    mutate(
      G = ifelse(grepl("^Intervention", Group), "I", "C"),
      phase = dplyr::recode(
        phase_label,
        "Baseline" = "B",
        "Phase 1"  = "P1",
        "Phase 2"  = "P2",
        .default   = NA_character_
      )
    ) %>%
    filter(!is.na(phase)) %>%
    select(draw, village, G, phase, rate = ID12_per100) %>%
    pivot_wider(
      names_from  = c(G, phase),
      values_from = rate
    )
  
  # Ensure all needed columns exist (fill with NA if a village/phase is missing)
  needed <- c("I_B","I_P1","I_P2","C_B","C_P1","C_P2")
  missing <- setdiff(needed, names(rd))
  for (m in missing) rd[[m]] <- NA_real_
  
  # Per-draw DiD by phase (guard against zero/NA denominators)
  safe_div <- function(num, den) ifelse(is.finite(num) & is.finite(den) & den != 0, num/den, NA_real_)
  
  did_long <- rd %>%
    transmute(
      draw, village,
      `Phase 1` = safe_div( safe_div(I_P1, I_B), safe_div(C_P1, C_B) ),
      `Phase 2` = safe_div( safe_div(I_P2, I_B), safe_div(C_P2, C_B) )
    ) %>%
    tidyr::pivot_longer(
      cols = c(`Phase 1`, `Phase 2`),
      names_to = "phase_label",
      values_to = "did"
    )
  
  # Summaries + pretty label
  lab <- did_long %>%
    group_by(village, phase_label) %>%
    summarise(
      med = stats::median(did, na.rm = TRUE),
      lo  = stats::quantile(did, 0.025, na.rm = TRUE),
      hi  = stats::quantile(did, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("IRR %.2f (%.2f–%.2f)", med, lo, hi))
  
  lab
}

# -----------------------------------------------------------
# ESTIMATE VILLAGE LEVEL ACQUISITION RATES PER PHASE
# -----------------------------------------------------------
# Load in data used
stan_data_fit <- readRDS("./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")

# Load in model fit 
fit = readRDS(results_path_1)
results_path_1 <- "./Output/Model_results/Observed_data/Two_step_sine_NonAdd_NonCol_seasonalityobserved.rds"

# load in data with house and village ids
load("./Data/BF/clean/use_in_analyses/bf_esbl0123_long_all.rda")
data <- read.csv(file = "./Data/BF/clean/use_in_analyses/bf_df_model_all.csv")

villages <- data %>% select(menage_id, HouseID, village, village_name)
villages$VillageID = villages$village_name

# include_season = TRUE  -> season-adjusted (uses model’s seasonal term)
pv_cm <- per_village_calendar_adjusted(
  fit, stan_data_fit,
  include_season = TRUE,
  scale = 100,
  control_strategy = "pooled_matched"   # <- village-explicit, no pooled control
)


# ESTIMATES ACQUISITION RATES FOR INTERVENTION GROUP, and CALENDAR MATCHED CONTROLS
#------------------------------------------------------------------------------------
pv_cm$rates_draws_village$scen = paste0(pv_cm$rates_draws_village$village,"_", pv_cm$rates_draws_village$Group)
unique(pv_cm$rates_draws_village$scen)
draws = dim(pv_cm$rates_draws_village)[1]/22

pv_cm$rates_draws_village$draws = rep(c(1:draws), 22)

# Plot village-level IRR compared to calendar-matched
lab_DID <- village_did_labels(pv_cm$rates_draws_village, upper_cap = 25)
y = pv_cm$rates_summary_village%>%filter(phase_label!="Baseline", Group=="Intervention")%>%
  rename(hi_incidence = hi)

lab_DID = left_join(lab_DID, y%>%select(village, phase_label,hi_incidence))

p_vill <-plot_village_rates(pv_cm$rates_draws_village, lab = lab_DID, upper_cap = 25) 
print(p_vill)

# NOW GET VILLAGE EXPLICIT ESTIMATES FOR CONTROL GROUP PER VILLAGE TO LINK WITH WASH AND ABX USE
#---------------------------------------------------------------------------------------
pv <- per_village_calendar_adjusted(
  fit, stan_data_fit,
  include_season = TRUE,
  scale = 100,
  control_strategy = "per_village_self"   # <- alternatives are "none", "pooled_matched"
)

draws = dim(pv$rates_draws_village)[1]/22

pv$rates_draws_village$draws = rep(c(1:draws), 22)

vill_irr_tbl <- lab_DID %>%
  mutate(
    irr = sprintf("%.2f (%.2f–%.2f)", med, lo, hi)
  ) %>%
  select(village, phase_label, irr) %>%
  pivot_wider(
    id_cols   = village,
    names_from = c(phase_label),
    values_from = irr  # pretty column headers
  )

vill_irr_tbl 

# INTERVENTION pretty strings
# Intervention — keep median/lo/hi/var numeric; make pretty string without altering them
inc_base  <- pv_cm$rates_summary_village %>% filter(phase_label=="Baseline", Group=="Intervention")
inc_base$Var_base <- inc_base$var
inc_base$Inc_base <- with(inc_base, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

inc_phase1 <- pv_cm$rates_summary_village %>% filter(phase_label=="Phase 1", Group=="Intervention")
inc_phase1$Var_phase1 <- inc_phase1$var
inc_phase1$Inc_phase1 <- with(inc_phase1, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

inc_phase2 <- pv_cm$rates_summary_village %>% filter(phase_label=="Phase 2", Group=="Intervention")
inc_phase2$Var_phase2 <- inc_phase2$var
inc_phase2$Inc_phase2 <- with(inc_phase2, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

# Attach to IRR table (Intervention rows)
vill_irr_tbl <- left_join(vill_irr_tbl, inc_base %>% select(village, Inc_base)) %>%
  select(village, Inc_base, `Phase 1`, `Phase 2`)
vill_irr_tbl <- left_join(vill_irr_tbl, inc_phase1 %>% select(village, Inc_phase1)) %>%
  select(village, Inc_base, Inc_phase1, `Phase 1`, `Phase 2`)
vill_irr_tbl <- left_join(vill_irr_tbl, inc_phase2 %>% select(village, Inc_phase2)) %>%
  select(village, Inc_base, Inc_phase1, Inc_phase2, `Phase 1`, `Phase 2`)

# add Intervention variances
vill_irr_tbl <- left_join(vill_irr_tbl, var_base   %>% select(village, Var_base))
vill_irr_tbl <- left_join(vill_irr_tbl, var_phase1 %>% select(village, Var_phase1))
vill_irr_tbl <- left_join(vill_irr_tbl, var_phase2 %>% select(village, Var_phase2))

vill_irr_tbl$Group <- "Intervention"

#-----------------------------------------
# CONTROLS (same pattern)
#-----------------------------------------
inc_base_c  <- pv$rates_summary_village %>% filter(phase_label=="Baseline", Group!="Intervention")
inc_base_c$Var_base <- inc_base_c$var
inc_base_c$Inc_base <- with(inc_base_c, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

inc_phase1_c <- pv$rates_summary_village %>% filter(phase_label=="Phase 1", Group!="Intervention")
inc_phase1_c$Var_phase1 <- inc_phase1_c$var
inc_phase1_c$Inc_phase1 <- with(inc_phase1_c, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

inc_phase2_c <- pv$rates_summary_village %>% filter(phase_label=="Phase 2", Group!="Intervention")
inc_phase2_c$Var_phase2 <- inc_phase2_c$var
inc_phase2_c$Inc_phase2 <- with(inc_phase2_c, sprintf("%.2f (%.2f-%.2f)", median, lo, hi))

# assemble a Control table; IRR columns are NA for controls
vill_control <- cbind(
  village     = inc_base_c$village,
  Group       = "Control",
  Inc_base    = inc_base_c$Inc_base,
  Inc_phase1  = inc_phase1_c$Inc_phase1,
  Inc_phase2  = inc_phase2_c$Inc_phase2,
  Var_base    = var_base_c$Var_base,
  Var_phase1  = var_phase1_c$Var_phase1,
  Var_phase2  = var_phase2_c$Var_phase2,
  `Phase 1`   = rep(NA_character_, nrow(inc_base_c)),
  `Phase 2`   = rep(NA_character_, nrow(inc_base_c))
) %>% as.data.frame()

# stack, sort, title-case, select final cols
vill_irr_tbl <- rbind(vill_irr_tbl, vill_control)
vill_irr_tbl <- vill_irr_tbl[order(vill_irr_tbl$village), ]
vill_irr_tbl$village <- tools::toTitleCase(vill_irr_tbl$village)

# keep both the pretty incidence strings and the variances, plus the IRR columns
vill_irr_tbl_final <- vill_irr_tbl %>%
  select(village, Group,
         Inc_base, Inc_phase1, Inc_phase2,
         `Phase 1`, `Phase 2`)
vill_irr_tbl_final

# DATA FILE FOR ANALYSES WITH WASH AND ABX
Inc_final_i  <- pv_cm$rates_summary_village %>% filter(Group=="Intervention")
Inc_final_c  <- pv$rates_summary_village %>% filter(Group=="Control")
Inc_final <- rbind(Inc_final_i, Inc_final_c)

# Store plot and summary table
#------------------------------------------------

# Store Figure S1
ggsave(
  filename = "./Output/Figures_and_tables/Paper/FigureS2_Village-level_intervention_effect_basecase.png",
  plot     = p_vill,
  width    = 14,   # in inches
  height   = 12,    # in inches
  dpi      = 300,   
  bg = "white"
)

# Store table
writexl::write_xlsx(vill_irr_tbl , path = "./Output/Figures_and_tables/Paper/Table_S1_Village_level_IRR.xlsx")
write_csv(Inc_final, file = "./Output/Model_results/Model_summaries/Observed_data/Summary_tables/Scenario1/inc_s1.csv")
