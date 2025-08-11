###############################################################################
# MODEL WITH SEASONAL SPLINES
###############################################################################

# ------------------------------------------------------------------------------
# Seasonality Setup – Simulation vs. Fitting
# SIMULATED DATA: Seasonality
# MODEL FITTING: Seasonality with spline
# ------------------------------------------------------------------------------

# Date created: 1 July 2025
# Date last updated: 1 August 2025
# Author: Raneem Aizouk

rm(list=ls())

#------------------------------------------------------------------------------
# Cluster set up
#------------------------------------------------------------------------------

# Get command-line arguments or use defaults
args <- commandArgs(trailingOnly = TRUE)

#scenario    <- if (length(args) >= 1) args[[1]] else Sys.getenv("scenario", "Two_step_spline_seasonality")
#data_source <- if (length(args) >= 2) args[[2]] else Sys.getenv("data_source", "simulated")

scenario    <- Sys.getenv("scenario", "Two_step_spline_seasonality")
data_source <- Sys.getenv("data_source", "simulated")

cat("R sees scenario:", scenario, "\n")
cat("R sees data_source:", data_source, "\n")

# Set output directory based on data source
if (data_source == "simulated") {
  output_dir <- "./CABU_EICO/model-output/Simulated_data/"
} else if (data_source == "observed") {
  output_dir <- "./CABU_EICO/model-output/Observed_data/"
}

output_dir <- file.path(output_dir)

# Load libraries
pacman::p_load(rstan,dplyr,lubridate, tidyr)


# Simulated data:
#  Acquisition rate (1→2) includes seasonal effect: sine wave (sin(2πx)).
#  Recovery rate (2→1) is constant (no seasonality).
#  Rates passed to Stan via lambda12_vec and lambda21_vec.

# Fitting:
# This Stan model fits the acquisition rate (1 → 2) using a B-spline function:
# The recovery rate (2 → 1) is assumed constant (no seasonal variation).
# Spline-based seasonality has been removed and replaced with parametric sin/cos terms.
# Goal:
# Evaluate if spline-based fitting can recover sinusoidal acquisition pattern.

#-----------------------------------------------------
# Timing Variables Explained
# ------------------------------------------------------
# date1: The start date of the observation interval for a specific individual (usually previous visit date)
# date2: The end date of the observation interval for a specific individual (usually current visit date)

# interval_start_first: Start of the first global interval (e.g., 28-day chunk) that overlaps with date1
# interval_end_first: End of the first global interval that overlaps with date1

# first_subinterval: Number of days from date1 to interval_end_first
#  If the interval starts in the middle of a global bin, this is the length of that partial overlap
# If the interval starts exactly on the bin boundary, this will be 0

# interval_start_last: Start of the last global interval that overlaps with date2

# last_subinterval: Number of days from interval_start_last to date2
# If the interval ends before the end of a bin, this captures the tail chunk
# If the interval ends exactly at the bin boundary, this will be 0

# num_middle_subintervals: Number of full 28-day global intervals strictly between date1 and date2


# Set up spline basis
#------------------------------------------------------------------------------
num_knots <- 5
knots <- seq(0, 1, length.out = num_knots)

knots <- as.numeric(knots)
spline_degree <- 3
num_basis <- num_knots + spline_degree 


#------------------------------------------------------------------------------
# Load in data
#------------------------------------------------------------------------------

if (data_source == "simulated") {
  scen = ifelse(scenario=="Two_step_spline_seasonality",
                "Simulated_data_seasonality_stan_data",
                "Simulated_data_noseasonality_stan_data")
  scen_df = ifelse(scenario=="Two_step_spline_seasonality",
                   "Simulated_data_seasonality",
                   "Simulated_data_noseasonality")
  
  sim_stan_data <- readRDS(paste0("./CABU_EICO/data/Simulated_data/", scen, ".rds"))
  sim_df <- readRDS(paste0("./CABU_EICO/data/Simulated_data/", scen_df, ".rds"))
  
  # When run locally
  #sim_stan_data <- readRDS("./Data/Simulated_data/Simulated_data_seasonality_stan_data.rds")
  #sim_df <- readRDS("./Data/Simulated_data/Simulated_data_seasonality.rds")
  
  stan_data_fit <- list(
    N = nrow(sim_df),
    N_persons = length(unique(sim_df$MenagememberID)),
    H = length(unique(sim_df$HouseID)),
    HouseID = sim_df$HouseID,
    menage_id_member = sim_df$MenagememberID,
    age = sim_df$Age,
    round = sim_df$Round,
    sexe = sim_df$Sexe,
    observed_state = sim_df$Observed_State_Sim,
    date_use = sim_df$date_use_num,
    intervention = sim_df$Intervention,
    intervention_date = sim_stan_data$intervention_date,
    intervention_date2 = sim_stan_data$intervention_date2,
    global_interval_start = sim_stan_data$global_interval_start,
    global_interval_end = sim_stan_data$global_interval_end,
    interval_length = sim_stan_data$interval_length,
    num_data = length(sim_stan_data$X),
    X = sim_stan_data$X,
    max_middle = sim_stan_data$max_middle,
    num_intervals = sim_stan_data$num_data,
    idx_first_sim = sim_df$Idx_First_Sim,
    idx_last_sim = sim_df$Idx_Last_Sim,
    idx_middle_sim = sim_stan_data[["idx_middle_sim"]],
    first_subinterval_sim = sim_df$First_Subinterval_Sim,
    last_subinterval_sim = sim_df$Last_Subinterval_Sim,
    num_middle_subintervals_sim = sim_df$Num_Middle_Subintervals_Sim,
    global_interval_index_start = sim_df$Global_Interval_Index_Start,
    num_knots = num_knots,
    knots = knots,
    spline_degree = spline_degree,
    num_basis = num_basis
  )
  
} else if (data_source == "observed") {
  stan_data_fit <- readRDS("./CABU_EICO/data/Observed_data/bf_stan_data_all.rds")
  # If run locally
  #stan_data_fit <- readRDS("./Data/BF/clean/use_in_analyses/bf_stan_data_all.rds")
}

names(stan_data_fit)

#---------------------------------------------------------------------------
# Stan model code
#---------------------------------------------------------------------------

stan_code <- "functions {
 matrix transition_matrix(real t, real lambda_1_2, real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;
    matrix[2,2] P;
    
    //--- Fail safe, guard against division by zero if rates are very small
    //-----------------------------------------------------------------------
    if (total_lambda <= 0) {
      // No flow possible (or underflow): stay where you are
      P[1,1] = 1; P[1,2] = 0;
      P[2,1] = 0; P[2,2] = 1;
      return P;
    }
    {
      real exp_term = exp(-total_lambda * t);
      P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;
      P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;
      P[1,2] = (lambda_1_2 / total_lambda) * (1 - exp_term);
      P[2,1] = (lambda_2_1 / total_lambda) * (1 - exp_term);
    }
    return P;
  }
  // Recursive function to compute B-spline basis functions
  vector build_b_spline(real[] t,
                        real[] ext_knots,
                        int ind,
                        int order) {
    int M = size(t);
    vector[M] b_spline = rep_vector(0, M);
    if (order == 1) {
      for (i in 1:M)
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    } else {
      real denom1 = ext_knots[ind+order-1] - ext_knots[ind];
      real denom2 = ext_knots[ind+order] - ext_knots[ind+1];
      vector[M] w1 = denom1 > 0
        ? (to_vector(t) - rep_vector(ext_knots[ind], M)) / denom1
        : rep_vector(0, M);
      vector[M] w2 = denom2 > 0
        ? 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], M)) / denom2
        : rep_vector(0, M);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1)
               + w2 .* build_b_spline(t, ext_knots, ind + 1, order-1);
    }
    return b_spline;
  }

}

data {
  int<lower=1> N;
  array[N] int HouseID;
  int<lower=1> H;
  array[N] int<lower=1, upper=2> observed_state;
  array[N] int menage_id_member;
  array[N] int age;
  array[N] int round;
  array[N] int sexe;
  real global_interval_start;
  real global_interval_end;
  int<lower=1> interval_length;
  array[N] real date_use;
  // vector[N] intervention;   
  array[N] int<lower=0, upper=1> intervention; 
  int<lower=1> num_data;
  array[num_data] real X;
  array[N] real intervention_date;
  array[N] real intervention_date2;
  int<lower=1> max_middle;
  int<lower=1> num_intervals;
  //SPLINE
  int<lower=1> num_knots;
  vector[num_knots] knots;
  int<lower=0> spline_degree;
  //int<lower=1> num_basis;
}

transformed data {
  array[N] real first_subinterval;
  array[N] real middle_subinterval;
  array[N] real last_subinterval;
  array[N] int num_middle_subintervals;
  array[N] int global_interval_index_start;
  array[N] int global_interval_index_end;
  array[N] real interval_start_first;
  array[N] real interval_end_first;
  array[N, max_middle] real interval_start_middle;
  array[N, max_middle] real interval_end_middle;
  array[N] real interval_start_last;
  array[N] real interval_end_last;
  // Declare and initialize the index arrays to 1
  array[N] int idx_first;
  array[N, max_middle] int idx_middle;
  array[N] int idx_last;

  idx_first = rep_array(1, N);
  idx_middle = rep_array(1, N, max_middle);
  idx_last = rep_array(1, N);
  
  // SET UP SPLINE
  //-------------------------------------------------------------------------------------------
    int num_basis = num_knots + spline_degree ;
   matrix[num_data, num_basis] B;

  // this part is mapping X=[0, 1.38] to [0,1]
  array[num_data] real X_mod;  // Declare X_mod
  // Compute X_mod = fmod(X, 1.0)
  //X_mod = [0, 0.25, 0.5, 0.75, 0, 0.25, 0.38] This maps time points beyond 1 (e.g., 1.25 = 15 months) back to their equivalent position in the first year (e.g., 0.25 = 3 months).
  
  for (i in 1:num_data) {
    X_mod[i] = fmod(X[i], 1.0);
  }

  // Extend knots
  
  vector[2 * spline_degree + num_knots + 1] ext_knots;
  ext_knots = append_row(
              knots[(num_knots - spline_degree + 1):num_knots] - 1,
              append_row(knots, knots[1:(spline_degree+1)] + 1)
            );

  for (ind in 1:num_basis)
    B[:, ind] = build_b_spline(X_mod, to_array_1d(ext_knots), ind, spline_degree + 1);
  
  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n-1]) {
      real date1 = date_use[n-1];
      real date2 = date_use[n];
      
      // 1) Map date1/date2 into global-interval bins
      int raw_start = to_int(floor((date1 - global_interval_start) / interval_length));
      int raw_end = to_int(ceil((date2 - global_interval_start) / interval_length));
      
      // 1.1) Clamp into [1, num_intervals]
      global_interval_index_start[n] = max(1, min(num_intervals, raw_start + 1));
      global_interval_index_end[n] = max(1, min(num_intervals, raw_end));
      
      // 2) Number of full intervals strictly between
      num_middle_subintervals[n] = global_interval_index_end[n] - global_interval_index_start[n] - 1;
      
      // 3) First subinterval calendar start/end
      interval_start_first[n] = global_interval_start + (global_interval_index_start[n] - 1) * interval_length;
      interval_end_first[n] = interval_start_first[n] + interval_length - 1;
      
      // 4) Middle subintervals calendar start/end
      for (m in 1:num_middle_subintervals[n]) {
        int global_idx = global_interval_index_start[n] + m;
        interval_start_middle[n,m] = global_interval_start + (global_idx - 1) * interval_length;
        interval_end_middle[n,m] = interval_start_middle[n,m] + interval_length - 1;
      }
      
      // 5) Last subinterval calendar start/end
      interval_start_last[n] = global_interval_start + (raw_end - 1) * interval_length;
      interval_end_last[n] = date2;
      
      // 6) Durations
      // First subinterval duration
      first_subinterval[n] = fmin(fmin(interval_end_first[n] - date1 + 1, date2 - date1 + 1), interval_length * 1.0);
      
      // Last subinterval duration
      last_subinterval[n] = fmin(date2 - (global_interval_start + (global_interval_index_end[n] - 1) * interval_length) + 1, interval_length);
      
      middle_subinterval[n] = num_middle_subintervals[n] * interval_length;
      
      // 7) For each subinterval, compute the midpoint and map it to the corresponding global-interval index (1…num_data)
      // 7.1) First midpoint → idx_first
      {
        real mid1 = (date1 + interval_end_first[n]) / 2;
        int raw_mid1 = to_int(floor((mid1 - global_interval_start) / interval_length)) + 1;
        idx_first[n] = max(1, min(num_data, raw_mid1));
      }
      
      // 7.2) Middle midpoints → idx_middle
      for (m in 1:num_middle_subintervals[n]) {
        real ms = interval_start_middle[n,m];
        real me = interval_end_middle[n,m];
        real midm = (ms + me) / 2;
        int raw_midm = to_int(floor((midm - global_interval_start) / interval_length)) + 1;
        idx_middle[n,m] = max(1, min(num_data, raw_midm));
      }
      
      // 7.3) Last midpoint → idx_last
      {
        real ml = (interval_start_last[n] + date2) / 2;
        int raw_ml = to_int(floor((ml - global_interval_start) / interval_length)) + 1;
        idx_last[n] = max(1, min(num_data, raw_ml));
      }
    }
  }
}

parameters {
  real q_1_2_raw;
  real q_2_1_raw;
  real beta_1_2_age;
  real beta_2_1_age;
  real beta_1_2_sexe;
  real beta_2_1_sexe;
  vector[H] u_raw;
  real<lower=0> sigma_u;
  real<lower=0> sigma_q_1_2;
  real<lower=0> sigma_q_2_1;
  real<lower=0> sigma;
 
 
  real beta_int1_1;
  real beta_int1_2;
  real beta_int2_1;
  real beta_int2_2;
  
  // Infection spline (periodic)

  row_vector[num_basis - spline_degree] a_raw_1_2_free;
  real log_tau_raw_1_2;


}

transformed parameters {
  real q_1_2_base = -3.5 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
  vector[H] u = u_raw * sigma_u;
  
  // --- Infection spline coefficients 
  //-------------------------------------------------------------------------------------
  row_vector[num_basis] a_raw_1_2;

 // Complete the periodic spline by repeating the first 3 basis weights

  for (i in 1:(num_basis - spline_degree))
    a_raw_1_2[i] = a_raw_1_2_free[i];

  for (j in 1:spline_degree)
    a_raw_1_2[num_basis - spline_degree + j] = a_raw_1_2[j];

  // Scale by global variance
  real tau_1_2 = exp(log_tau_raw_1_2);
  row_vector[num_basis] a_1_2 = a_raw_1_2 * tau_1_2;

  // Linear trend term


  // Final seasonal effect for infection (with linear + periodic terms)
 vector[num_data] Y_hat_1_2 = to_vector(B * a_1_2');


}

model {
  q_1_2_raw ~ normal(0, 0.3);
  q_2_1_raw ~ normal(0, 0.38);
  beta_1_2_age ~ normal(0, 1);
  beta_2_1_age ~ normal(0, 1);
  beta_1_2_sexe ~ normal(0, 1);
  beta_2_1_sexe ~ normal(0, 1);
  beta_int1_1 ~ normal(0, 1);
  beta_int1_2 ~ normal(0, 1);
  beta_int2_1 ~ normal(0, 1);
  beta_int2_2 ~ normal(0, 1);
  sigma_q_1_2 ~ normal(0, 1);
  sigma_q_2_1 ~ normal(0, 1);
  u_raw ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
 //SPLINE
   a_raw_1_2_free ~ normal(0, 1);
  // a_raw_1_2  ~ normal(0, 1);
    log_tau_raw_1_2 ~ normal(0, 0.5);
     
     
  
  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n-1]) {
      real log12_base = q_1_2_base + u[HouseID[n]] + beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n];
      real log21_base = q_2_1_base + u[HouseID[n]] + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n];
       matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));
      real t_star1 = intervention_date[n];
      real t_star2 = intervention_date2[n];
 
      
    //----------------------------------------------------------------------------------------------------------
    // --- First subinterval
    //----------------------------------------------------------------------------------------------------------
    // The first subintervals between each of the observations, concerns the first interval between a transition, which is <28 days

     if (first_subinterval[n] > 0) {
        int i1 = idx_first[n];
        real s12 = Y_hat_1_2[i1];
        real t0 = date_use[n-1];
        real t1 = t0 + first_subinterval[n] - 1;
        real log_lambda_1_2;
        real log_lambda_2_1;
        
        //----------------------------------------------------------------------------------------------------------
        // Case 1: Intervention 1 starts within the interval → split in pre-intervention and post-intervention 1
        //----------------------------------------------------------------------------------------------------------
        if (t0 < t_star1 && t_star1 < t1 && intervention[n] == 1) {
          real d1a = t_star1 - t0;
          real d1b = t1 - t_star1 + 1;
          
          // --- Pre-intervention 1
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12;
          log_lambda_2_1 = log21_base;
          {
            matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
          // --- Post-intervention 1
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12 + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          //----------------------------------------------------------------------------------------------------------
          // Case 2: Intervention 2 starts within the interval (intervention 1 already active)
          //----------------------------------------------------------------------------------------------------------
        } else if (t0 < t_star2 && t_star2 < t1 && intervention[n] == 1) {
          real d1a = t_star2 - t0;
          real d1b = t1 - t_star2 + 1;
          
          // --- Pre-intervention 2 (intervention 1 only)
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12 + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
         // --- Post-intervention 2
         //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12 + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          {
            matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          //----------------------------------------------------------------------------------------------------------
          // Cases 3a–3c: No intervention starts inside → use midpoint logic
          //----------------------------------------------------------------------------------------------------------
        } else {
          real midpoint = (t0 + t1) / 2;
          
          // Case 3a: Interval lies fully within post-intervention 2
          //----------------------------------------------------------------------------------------------------------
          if (midpoint >= t_star2 && intervention[n] == 1) {
            log_lambda_1_2 = log12_base + s12 + beta_int1_1 + beta_int1_2;
            log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          // Case 3b: Interval lies fully within post-intervention 1
          //----------------------------------------------------------------------------------------------------------
          } else if (midpoint >= t_star1 && intervention[n] == 1) {
            log_lambda_1_2 = log12_base + s12 + beta_int1_1;
            log_lambda_2_1 = log21_base + beta_int2_1;
         // Case 3c: Interval is fully pre-intervention
         //---------------------------------------------------------------------------------------------------------
          } else {
            log_lambda_1_2 = log12_base + s12;
            log_lambda_2_1 = log21_base;
          }
          {
            matrix[2,2] P = transition_matrix(first_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
        }
      }
      
//----------------------------------------------------------------------------------------------------------
// --- Middle subintervals
//----------------------------------------------------------------------------------------------------------
// The middle subintervals between each of the transitions are all exactly 28 days

      for (m in 1:num_middle_subintervals[n]) {
        int im = idx_middle[n,m];
        real s12m = Y_hat_1_2[im];
        real t0m = global_interval_start + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;
        real log_lambda_1_2;
        real log_lambda_2_1;
        //----------------------------------------------------------------------------------------------------------
        // Case 1: Intervention 1 starts within the interval → split in pre-intervention and post-intervention 1
        //----------------------------------------------------------------------------------------------------------
        if (t0m < t_star1 && t_star1 < t1m && intervention[n] == 1) {
          real d2a = t_star1 - t0m;
          real d2b = t1m - t_star1 + 1;
          
          // Pre-intervention 1
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12m;
          log_lambda_2_1 = log21_base;
          {
            matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
          // Post-intervention 1 and pre-intervention 2
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12m + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
        //----------------------------------------------------------------------------------------------------------
        // Case 2: Intervention 2 starts within the interval (intervention 1 already active)
        //----------------------------------------------------------------------------------------------------------
        } else if (t0m < t_star2 && t_star2 < t1m && intervention[n] == 1) {
          real d2a = t_star2 - t0m;
          real d2b = t1m - t_star2 + 1;
          
          // Pre-intervention 2
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12m + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
          // Post-intervention 2
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12m + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          {
            matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
        //----------------------------------------------------------------------------------------------------------
        // Cases 3a–3c: No intervention starts inside → use midpoint logic
        //----------------------------------------------------------------------------------------------------------
        } else {
          real midpoint = (t0m + t1m) / 2;
          
          // Case 3a: Full interval during intervention 2
          //----------------------------------------------------------------------------------------------------------
          if (midpoint >= t_star2 && intervention[n] == 1) {
            log_lambda_1_2 = log12_base + s12m + beta_int1_1 + beta_int1_2;
            log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          } else if (midpoint >= t_star1 && intervention[n] == 1) {
            // Case 3b: Full interval during intervention 1
            //----------------------------------------------------------------------------------------------------------
            log_lambda_1_2 = log12_base + s12m + beta_int1_1;
            log_lambda_2_1 = log21_base + beta_int2_1;
          } else {
            // Case 3c: Full interval pre-intervention
            //----------------------------------------------------------------------------------------------------------
            log_lambda_1_2 = log12_base + s12m;
            log_lambda_2_1 = log21_base;
          }
          {
            matrix[2,2] P = transition_matrix(interval_length, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
        }
      }
      
//----------------------------------------------------------------------------------------------------------
// --- Last subinterval
//----------------------------------------------------------------------------------------------------------
// The last subintervals between each of the observations, concerns the last interval between a transition, which again is <28 days

      if (last_subinterval[n] > 0) {
        int il = idx_last[n];
        real s12l = Y_hat_1_2[il];
        real t1l = date_use[n];
        real t0l = t1l - last_subinterval[n] + 1;
        real log_lambda_1_2;
        real log_lambda_2_1;
    //----------------------------------------------------------------------------------------------------------
    // Case 1: Intervention 1 starts during interval → split in pre-intervention and post-intervention 1
    //----------------------------------------------------------------------------------------------------------
        if (t0l < t_star1 && t_star1 < t1l && intervention[n] == 1) {
          real d1a = t_star1 - t0l;
          real d1b = t1l - t_star1 + 1;
          
          // Pre-intervention 1
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12l;
          log_lambda_2_1 = log21_base;
          {
            matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
          // Post-intervention 1
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12l + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
      //----------------------------------------------------------------------------------------------------------
      // Case 2: Intervention 2 starts during interval → split in intervention 1 and post-intervention 2
      //----------------------------------------------------------------------------------------------------------
        } else if (t0l < t_star2 && t_star2 < t1l && intervention[n] == 1) {
          real d2a = t_star2 - t0l;
          real d2b = t1l - t_star2 + 1;
          
          // Pre-intervention 2 (intervention 1 only)
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12l + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          {
            matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
          
          // Post-intervention 2 (intervention 1 + 2)
          //----------------------------------------------------------------------------------------------------------
          log_lambda_1_2 = log12_base + s12l + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          {
            matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
      //----------------------------------------------------------------------------------------------------------
      // Case 3: No intervention date lies inside the interval
      //----------------------------------------------------------------------------------------------------------
        } else {
          real midpoint = (t0l + t1l) / 2;
          
          // Case 3a: Entire interval after intervention 2
          //----------------------------------------------------------------------------------------------------------
          if (midpoint >= t_star2 && intervention[n] == 1) {
            log_lambda_1_2 = log12_base + s12l + beta_int1_1 + beta_int1_2;
            log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
            
          // Case 3b: Entire interval after intervention 1 but before intervention 2
          //----------------------------------------------------------------------------------------------------------
          } else if (midpoint >= t_star1 && intervention[n] == 1) {
            log_lambda_1_2 = log12_base + s12l + beta_int1_1;
            log_lambda_2_1 = log21_base + beta_int2_1;
            
          // Case 3c: Entire interval before any intervention
          //----------------------------------------------------------------------------------------------------------
          } else {
            log_lambda_1_2 = log12_base + s12l;
            log_lambda_2_1 = log21_base;
          }
          {
            matrix[2,2] P = transition_matrix(last_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
            P_total *= P;
          }
        }
      }
      target += log(P_total[observed_state[n-1], observed_state[n]] + 1e-9);
    }
  }
}
generated quantities {
  vector[N] log_lambda_1_2_out;          // Log acquisition rate per observation
  vector[N] log_lambda_2_1_out;          // Log decolonisation rate per observation
  vector[num_data] Y_hat_1_2_out;        // Seasonal intervention effect over time
  array[N] int second_intervention_used;  // Whether second intervention was active
  array[N] real total_subinterval_duration; // Total duration of subintervals
  array[N] int y_rep;                    // Replicated data for posterior predictive checks

  // NEW: counters
  array[N] int acquisitions;
  array[N] int decolonisations;
  array[N] int at_risk_acquisition;
  array[N] int at_risk_decolonisation;

  // Initialize arrays
  second_intervention_used = rep_array(0, N);
  y_rep = rep_array(0, N);
  acquisitions = rep_array(0, N);
  decolonisations = rep_array(0, N);
  at_risk_acquisition = rep_array(0, N);
  at_risk_decolonisation = rep_array(0, N);

  // Copy seasonal curve
  for (i in 1:num_data) {
    Y_hat_1_2_out[i] = Y_hat_1_2[i];
  }

  for (n in 1:N) {
    real log12_base = q_1_2_base + u[HouseID[n]] + beta_1_2_age * age[n] + beta_1_2_sexe * sexe[n];
    real log21_base = q_2_1_base + u[HouseID[n]] + beta_2_1_age * age[n] + beta_2_1_sexe * sexe[n];
    real t_star1 = intervention_date[n];
    real t_star2 = intervention_date2[n];
    real log_lambda_1_2;
    real log_lambda_2_1;

    // Handle first observation in household
    if (n == 1 || menage_id_member[n] != menage_id_member[n-1]) {
      int i0 = idx_first[n];
      real s12 = Y_hat_1_2_out[i0];
      real midpoint = date_use[n];
      if (midpoint >= t_star2 && intervention[n] == 1) {
        log_lambda_1_2 = log12_base + s12 + beta_int1_1 + beta_int1_2;
        log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
        second_intervention_used[n] = 1;
      } else if (midpoint >= t_star1 && intervention[n] == 1) {
        log_lambda_1_2 = log12_base + s12 + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        second_intervention_used[n] = 0;
      } else {
        log_lambda_1_2 = log12_base + s12;
        log_lambda_2_1 = log21_base;
        second_intervention_used[n] = 0;
      }
      // Keep first observation as is
      y_rep[n] = observed_state[n];
      log_lambda_1_2_out[n] = log_lambda_1_2;
      log_lambda_2_1_out[n] = log_lambda_2_1;
      total_subinterval_duration[n] = 0;
      continue;
    }

    // Cumulative transition matrix
    matrix[2, 2] P_total = diag_matrix(rep_vector(1.0, 2));
    int current_state = observed_state[n-1]; // start state for this gap

//----------------------------------------------------------------------------------------------------------
// --- First subinterval
//----------------------------------------------------------------------------------------------------------
    if (first_subinterval[n] > 0) {
      int i1 = idx_first[n];
      real s12 = Y_hat_1_2_out[i1];
      real t0 = date_use[n-1];
      real t1 = t0 + first_subinterval[n] - 1;
      
      //----------------------------------------------------------------------------------------------------------
      // Case 1: Intervention 1 starts within the interval → split in pre-intervention and post-intervention 1
      //----------------------------------------------------------------------------------------------------------
      if (t0 < t_star1 && t_star1 < t1 && intervention[n] == 1) {
        real d1a = t_star1 - t0;
        real d1b = t1 - t_star1 + 1;
        
        // Pre-intervention 1
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12;
        log_lambda_2_1 = log21_base;
        {
          matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        // Post-intervention 1
        //----------------------------------------------------------------------------------------------------------

        log_lambda_1_2 = log12_base + s12 + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 0;
        
      //----------------------------------------------------------------------------------------------------------
      // Case 2: Intervention 2 starts during interval → split in intervention 1 and post-intervention 2
      //----------------------------------------------------------------------------------------------------------
      } else if (t0 < t_star2 && t_star2 < t1 && intervention[n] == 1) {
        real d1a = t_star2 - t0;
        real d1b = t1 - t_star2 + 1;
        
        // Pre-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12 + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        // Post-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12 + beta_int1_1 + beta_int1_2;
        log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
        {
          matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 1;
        
      //----------------------------------------------------------------------------------------------------------
      // Case 3: No intervention date lies inside the interval
      //----------------------------------------------------------------------------------------------------------
      } else {
        real midpoint = (t0 + t1) / 2;
        
        // Case 3a: Entire interval after intervention 2
        //----------------------------------------------------------------------------------------------------------
        if (midpoint >= t_star2 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12 + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          second_intervention_used[n] = 1;
          
        // Case 3b: Entire interval after intervention 1 but before intervention 2
        //----------------------------------------------------------------------------------------------------------
        } else if (midpoint >= t_star1 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12 + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          second_intervention_used[n] = 0;
          
        // Case 3c: Entire interval before any intervention
        //----------------------------------------------------------------------------------------------------------
        } else {
          log_lambda_1_2 = log12_base + s12;
          log_lambda_2_1 = log21_base;
          second_intervention_used[n] = 0;
        }
        {
          matrix[2,2] P = transition_matrix(first_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
      }
    }

//----------------------------------------------------------------------------------------------------------
// --- Middle subintervals
//----------------------------------------------------------------------------------------------------------
    for (m in 1:num_middle_subintervals[n]) {
      int im = idx_middle[n,m];
      real s12m = Y_hat_1_2_out[im];
      real t0m = global_interval_start + (global_interval_index_start[n] + m - 1) * interval_length;
      real t1m = t0m + interval_length - 1;

      //----------------------------------------------------------------------------------------------------------
      // Case 1: Intervention 1 starts within the interval → split in pre-intervention and post-intervention 1
      //----------------------------------------------------------------------------------------------------------
      if (t0m < t_star1 && t_star1 < t1m && intervention[n] == 1) {
        real d2a = t_star1 - t0m;
        real d2b = t1m - t_star1 + 1;
        
        // Pre-intervention 1
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12m;
        log_lambda_2_1 = log21_base;
        {
          matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        // Post-intervention 1
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12m + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 0;
      //----------------------------------------------------------------------------------------------------------
      // Case 2: Intervention 2 starts during interval → split in intervention 1 and post-intervention 2
      //----------------------------------------------------------------------------------------------------------
      } else if (t0m < t_star2 && t_star2 < t1m && intervention[n] == 1) {
        real d2a = t_star2 - t0m;
        real d2b = t1m - t_star2 + 1;
        
        // Pre-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12m + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        // Post-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12m + beta_int1_1 + beta_int1_2;
        log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
        {
          matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 1;
        
      //----------------------------------------------------------------------------------------------------------
      // Case 3: No intervention date lies inside the interval
      //----------------------------------------------------------------------------------------------------------
      } else {
        real midpoint = (t0m + t1m) / 2;
        
        // Case 3a: Entire interval after intervention 2
        //----------------------------------------------------------------------------------------------------------
        if (midpoint >= t_star2 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12m + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          second_intervention_used[n] = 1;
          
          // Case 3b: Entire interval after intervention 1 but before intervention 2
          //----------------------------------------------------------------------------------------------------------
        } else if (midpoint >= t_star1 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12m + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          second_intervention_used[n] = 0;
          
          // Case 3c: Entire interval before any intervention
          //----------------------------------------------------------------------------------------------------------
        } else {
          log_lambda_1_2 = log12_base + s12m;
          log_lambda_2_1 = log21_base;
          second_intervention_used[n] = 0;
        }
        {
          matrix[2,2] P = transition_matrix(interval_length, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
      }
    }

//----------------------------------------------------------------------------------------------------------
// --- Last subinterval
//----------------------------------------------------------------------------------------------------------
    if (last_subinterval[n] > 0) {
      int il = idx_last[n];
      real s12l = Y_hat_1_2_out[il];
      real t1l = date_use[n];
      real t0l = t1l - last_subinterval[n] + 1;
      
      //----------------------------------------------------------------------------------------------------------
      // Case 1: Intervention 1 starts within the interval → split in pre-intervention and post-intervention 1
      //----------------------------------------------------------------------------------------------------------    
      if (t0l < t_star1 && t_star1 < t1l && intervention[n] == 1) {
        real d1a = t_star1 - t0l;
        real d1b = t1l - t_star1 + 1;
        
        // Pre-intervention 1
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12l;
        log_lambda_2_1 = log21_base;
        {
          matrix[2,2] P = transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        
        // Post-intervention 1
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12l + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 0;
        
      //----------------------------------------------------------------------------------------------------------
      // Case 2: Intervention 2 starts during interval → split in intervention 1 and post-intervention 2
      //----------------------------------------------------------------------------------------------------------
      } else if (t0l < t_star2 && t_star2 < t1l && intervention[n] == 1) {
        real d2a = t_star2 - t0l;
        real d2b = t1l - t_star2 + 1;
        
        // Pre-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12l + beta_int1_1;
        log_lambda_2_1 = log21_base + beta_int2_1;
        {
          matrix[2,2] P = transition_matrix(d2a, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        // Post-intervention 2
        //----------------------------------------------------------------------------------------------------------
        log_lambda_1_2 = log12_base + s12l + beta_int1_1 + beta_int1_2;
        log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
        {
          matrix[2,2] P = transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
        second_intervention_used[n] = 1;
        
      //----------------------------------------------------------------------------------------------------------
      // Case 3: No intervention date lies inside the interval
      //----------------------------------------------------------------------------------------------------------
      } else {
        real midpoint = (t0l + t1l) / 2;
        
        // Case 3a: Entire interval after intervention 2
        //----------------------------------------------------------------------------------------------------------
        if (midpoint >= t_star2 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12l + beta_int1_1 + beta_int1_2;
          log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
          second_intervention_used[n] = 1;
          
        // Case 3b: Entire interval after intervention 1 and before intervention 2
        //----------------------------------------------------------------------------------------------------------
        } else if (midpoint >= t_star1 && intervention[n] == 1) {
          log_lambda_1_2 = log12_base + s12l + beta_int1_1;
          log_lambda_2_1 = log21_base + beta_int2_1;
          second_intervention_used[n] = 0;
          
        // Case 3c: Entire interval before any intervention
        //----------------------------------------------------------------------------------------------------------
        } else {
          log_lambda_1_2 = log12_base + s12l;
          log_lambda_2_1 = log21_base;
          second_intervention_used[n] = 0;
        }
        {
          matrix[2,2] P = transition_matrix(last_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
          int next_state = categorical_rng(to_vector(P[current_state]) / sum(P[current_state]));
          if (current_state == 1) at_risk_acquisition[n] += 1;
          if (current_state == 2) at_risk_decolonisation[n] += 1;
          if (current_state == 1 && next_state == 2) acquisitions[n] += 1;
          if (current_state == 2 && next_state == 1) decolonisations[n] += 1;
          current_state = next_state;
          P_total *= P;
        }
      }
    }

    // Final simulated state for this observation
    y_rep[n] = current_state;

    log_lambda_1_2_out[n] = log_lambda_1_2;
    log_lambda_2_1_out[n] = log_lambda_2_1;
    total_subinterval_duration[n] = first_subinterval[n] + middle_subinterval[n] + last_subinterval[n];
  }
}
 "  
#---------------------------------------------------------------------------
# Compile the Stan model
#---------------------------------------------------------------------------

compiled_model_fit <- stan_model(model_code = stan_code, verbose = TRUE)

#---------------------------------------------------------------------------
# Fit the model
#---------------------------------------------------------------------------

stan_fit <- sampling(
  object = compiled_model_fit,
  data = stan_data_fit,
  algorithm = "NUTS", # Use NUTS for fitting (not Fixed_param)
  iter = 5000,        # Number of iterations
  warmup = 2000,      # Warmup iterations
  chains = 4,         # Number of chains
  cores = 4,          # Parallel processing
  seed = 123,
  verbose = TRUE,
  control = list(adapt_delta = 0.95, max_treedepth = 12) # Adjust for convergence
)

#---------------------------------------------------------------------------
# Save output 
#---------------------------------------------------------------------------

saveRDS(stan_fit, file = paste0(output_dir,scenario,data_source, ".rds"))

print(summary(stan_fit)$summary)








