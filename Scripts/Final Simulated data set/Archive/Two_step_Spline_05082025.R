

####################
num_knots <- 5
knots <- seq(0, 1, length.out = num_knots)

knots <- as.numeric(knots)
spline_degree <- 3
num_basis <- num_knots + spline_degree 

library(rstan)


# Dates need to be numeric for Stan
simulated_dataset$date_use <- as.numeric(simulated_dataset$Date)
simulated_dataset$intervention_date <- as.numeric(simulated_dataset$Intervention_Date1)
simulated_dataset$intervention_date2 <- as.numeric(simulated_dataset$Intervention_Date2)


# Use these simulated variables from your simulation
stan_data <- list(
  menage_id_member = as.integer(simulated_dataset$MenageID),
  HouseID = as.integer(simulated_dataset$HouseID),
  H = length(unique(simulated_dataset$HouseID)),
  VillageID = as.integer(simulated_dataset$VillageID),
  intervention = as.integer(simulated_dataset$Intervention),
  round = as.integer(simulated_dataset$Round),
  age = as.integer(simulated_dataset$Age),
  sexe = as.integer(simulated_dataset$Sexe),
  intervention_date = as.numeric(simulated_dataset$Intervention_Date1),
  intervention_date2 = as.numeric(simulated_dataset$Intervention_Date2),
  N = N,
  observed_state = as.integer(observed_state_sim),
  global_interval_start = as.numeric(global_interval_start),
  global_interval_end = as.numeric(global_interval_end),
  interval_length = interval_length,
  date_use = as.numeric(simulated_dataset$date_use),
  X = X_seasonal, 
  num_data = length(X_seasonal),
  max_middle = max_middle,
  num_intervals = num_intervals,
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  num_basis = num_basis
)


str(stan_data)
str(simulated_dataset)
# Prepare data for Stan
stan_model_full <- "
functions {
  matrix transition_matrix(real t, real lambda_1_2, real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;                           // Total transition rate (acquisitions and decolonisations)
    real exp_term = total_lambda > 0 ? exp(-total_lambda * t) : 1.0;      // Exponential decay term
    matrix[2,2] P;
    if (total_lambda == 0) {
      P[1,1] = 1.0; P[1,2] = 0.0;                                          // No transition possible
      P[2,1] = 0.0; P[2,2] = 1.0;
    } else {
      P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;  // Prob staying in state 1
      P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;  // Prob staying in state 2
      P[1,2] = (lambda_1_2 / total_lambda) * (1 - exp_term);                           // Prob moving from 1 to 2
      P[2,1] = (lambda_2_1 / total_lambda) * (1 - exp_term);                           // Prob moving from 2 to 1
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
  int<lower=1> H;
  array[N] int<lower=1, upper=2> observed_state;
  array[N] int menage_id_member;
  array[N] int HouseID;
  array[N] int age;
  array[N] int round;
  array[N] int sexe;
  real global_interval_start;
  real global_interval_end;
  int<lower=1> interval_length;
  array[N] real date_use;
  array[N] int intervention;
  int<lower=1> num_knots;
  vector[num_knots] knots;
  int<lower=0> spline_degree;
  //int<lower=1> num_basis;
  int<lower=1> num_data;
  array[num_data] real X;

  //array[N] real Intervention_start_date;
  array[N] real intervention_date;  // or int if stored as day numbers
  array[N] real intervention_date2; 
  int<lower=1> max_middle;
  int<lower=1> num_intervals;


}

transformed data {
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

  array[N] real first_subinterval;
  array[N] real middle_subinterval;
  array[N] real last_subinterval;
  array[N] int  num_middle_subintervals;
  array[N] int  global_interval_index_start;
  array[N] int  global_interval_index_end;
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
  // Fill every position with “1” so that unused slots remain in range
  idx_first  = rep_array(1, N);
  idx_middle = rep_array(1, N, max_middle);
  idx_last   = rep_array(1, N);

  for (n in 2:N) {
    if (menage_id_member[n] == menage_id_member[n-1]) {
      real date1 = date_use[n-1];
      real date2 = date_use[n];

      // 1) map date1/date2 into global-interval bins
      int raw_start = to_int(floor((date1 - global_interval_start) / interval_length));
      int raw_end   = to_int(ceil((date2 - global_interval_start) / interval_length));

      // 1.1) clamp into [1, num_intervals]
      global_interval_index_start[n] = max(1, min(num_intervals, raw_start + 1));
      global_interval_index_end  [n] = max(1, min(num_intervals, raw_end));

      // 2) number of full intervals strictly between
      num_middle_subintervals[n] = global_interval_index_end[n]
                                - global_interval_index_start[n]
                                - 1;

      // 3) first subinterval calendar start/end
      interval_start_first[n] = global_interval_start
                              + (global_interval_index_start[n] - 1) * interval_length;
      interval_end_first  [n] = interval_start_first[n] + interval_length - 1;

      // 4) middle subintervals calendar start/end
      for (m in 1:num_middle_subintervals[n]) {
        int global_idx = global_interval_index_start[n] + m;
        interval_start_middle[n,m] = global_interval_start
                                   + (global_idx - 1) * interval_length;
        interval_end_middle  [n,m] = interval_start_middle[n,m] + interval_length - 1;
      }

      // 5) last subinterval calendar start/end
      interval_start_last[n] = global_interval_start + (raw_end - 1) * interval_length;
      interval_end_last  [n] = date2;

      // 6) durations
      first_subinterval[n]  = interval_end_first[n] - date1 + 1;
      last_subinterval [n]  = date2
                            - (global_interval_start
                               + (global_interval_index_end[n] - 1) * interval_length)
                            + 1;
      middle_subinterval[n] = num_middle_subintervals[n] * interval_length;
    // 7) For each subinterval, compute the midpoint and map it to the corresponding global‐interval index (1…num_data):
    // 7.1) first midpoint → idx_first
      {
        real mid1 = (date1 + interval_end_first[n]) / 2;
        int raw_mid1 = to_int(floor((mid1 - global_interval_start) / interval_length)) + 1;
        idx_first[n] = max(1, min(num_data, raw_mid1));
      }

      // 7.2) middle midpoints → idx_middle
      for (m in 1:num_middle_subintervals[n]) {
        real ms   = interval_start_middle[n,m];
        real me   = interval_end_middle  [n,m];
        real midm = (ms + me) / 2;
        int raw_midm = to_int(floor((midm - global_interval_start) / interval_length)) + 1;
        idx_middle[n,m] = max(1, min(num_data, raw_midm));
      }

      // 7.3) last midpoint → idx_last
      {
        real ml = (interval_start_last[n] + date2) / 2;
        int raw_ml = to_int(floor((ml - global_interval_start) / interval_length)) + 1;
        idx_last[n] = max(1, min(num_data, raw_ml));
      }
    }
  }
}


 
parameters {
  // Transition rate base parameters (non-centered)
  real q_1_2_raw;
  real q_2_1_raw;

  real<lower=0> sigma_q_1_2;
  real<lower=0> sigma_q_2_1;

  // Random effects
  vector[H] u_raw;
  real<lower=0> sigma_u;

  // Covariate effects
  real beta_1_2_age;
  real beta_2_1_age;
  real beta_1_2_sexe;
  real beta_2_1_sexe;

  // Intervention and seasonal covariate effects
  real beta_int1_1;
  real beta_int1_2;
  real beta_int2_1;
  real beta_int2_2;

  // Infection spline (periodic)


  row_vector[num_basis - spline_degree] a_raw_1_2_free;
  real log_tau_raw_1_2;
  real a0_raw_1_2;
  real<lower=0> sigma_a0_1_2;


  // Noise parameter for observation model
  real<lower=0> sigma;
}

transformed parameters {
  // ─── Base transition log-rates ───────────────────────────────────────────────
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;

  vector[H] u = u_raw * sigma_u;

  // ─── Infection spline coefficients ───────────────────────────────────────────
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
  // ─────────────────────────────────────────────────────────────────────────────
  // PRIORS
  // ─────────────────────────────────────────────────────────────────────────────
  // Base CTMC transition coefficients (raw scale → centered in transformed parameters)
  q_1_2_raw ~ normal(0, 0.3);
  q_2_1_raw ~ normal(0, 0.38);
  
  // Age and sex effects on infection/recovery hazards
  beta_1_2_age  ~ normal(0, 1);
  beta_2_1_age  ~ normal(0, 1);
  beta_1_2_sexe ~ normal(0, 1);
  beta_2_1_sexe ~ normal(0, 1);

  // Household‐level random effects
  u_raw   ~ normal(0, 0.5);
  sigma_u ~ normal(0, 0.5);

  // Spline parameters for seasonal variation (infection )
  a_raw_1_2  ~ normal(0, 1);

   beta_int1_1 ~ normal(0, 1);
  beta_int1_2 ~ normal(0, 1);
  beta_int2_1 ~ normal(0, 1);
  beta_int2_2 ~ normal(0, 1);

  a_raw_1_2_free ~ normal(0, 1);
 
  log_tau_raw_1_2 ~ normal(0, 0.5);


  // Noise on continuous‐time transitions
  sigma ~ normal(0, 1);

  // Intervention‐specific effects (first and second)

   // ─────────────────────────────────────────────────────────────────────────────
  // LIKELIHOOD: FOR EACH PAIR OF SUCCESSIVE OBSERVATIONS (n−1 → n) WITHIN A PERSON
  // ─────────────────────────────────────────────────────────────────────────────
  for (n in 2:N) {
    // Only compute transitions if this row (n) belongs to the same individual
    // as the previous row (n−1).  Otherwise, skip (no baseline → no transition).
    if (menage_id_member[n] == menage_id_member[n-1]) {
      // 1) Compute baseline log‐hazards for infection (1→2) and recovery (2→1):
      //    q_1_2_base and q_2_1_base are the “population‐average” log‐rates,
      //    u[HouseID[n]] is the household‐specific random intercept,
      //    beta_*_age and beta_*_sexe are individual‐level covariates.
      real log12_base = q_1_2_base
                     + u[HouseID[n]]
                     + beta_1_2_age  * age[n]
                     + beta_1_2_sexe * sexe[n];
      real log21_base = q_2_1_base
                     + u[HouseID[n]]
                     + beta_2_1_age  * age[n]
                     + beta_2_1_sexe * sexe[n];

      // Start with an “identity” 2×2 matrix P_total; we will multiply in pieces.
      matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));

      // Intervention dates (first and second) in numeric form (days since epoch)
      real t_star1 = intervention_date[n];
      real t_star2 = intervention_date2[n];

      // ─────────────────────────────────────────────────────────────────────────
      // (A) FIRST SUBINTERVAL: from (date_use[n−1]) up to FIRST interval boundary
      // ─────────────────────────────────────────────────────────────────────────
      if (first_subinterval[n] > 0) {
        // idx_first[n] gives the spline index (1…num_data) for this piece’s midpoint.
        
        int    i1  = idx_first[n];
        real   s12 = Y_hat_1_2[i1];  // seasonal adjustment at that midpoint for infection
      

        // t0 = the actual calendar date of observation (n−1)
        
          real t0  = date_use[n-1];
          
          
        // t1 = end of “first subinterval,” measured in days.
        // (Because first_subinterval[n] = length of that piece in days.)
        
        
        real t1  = t0 + first_subinterval[n] - 1;

        // (A.1) If t_star1 lies strictly inside [t0, t1], we split this subinterval
        
        real log_lambda_1_2;
        real log_lambda_2_1;
        
        if (t0 < t_star1 && t_star1 < t1   && intervention[n] == 1 ){
        real d1a = t_star1 - t0;        // duration before first intervention
        real d1b = t1 - t_star1 + 1;    // duration after first intervention

         // --- Pre-intervention 1: baseline (no interventions)
           log_lambda_1_2 = log12_base + s12;
            log_lambda_2_1 = log21_base;
         P_total *= transition_matrix(d1a, exp(log_lambda_1_2), exp(log_lambda_2_1));


           // Post-intervention 1: add first intervention effect
          log_lambda_1_2 = log12_base + s12 +  beta_int1_1;
          log_lambda_2_1 = log21_base +  beta_int2_1;
           P_total *= transition_matrix(d1b, exp(log_lambda_1_2), exp(log_lambda_2_1));

        } else {
      // No first intervention effect in this subinterval:
    // apply only baseline + seasonality as if the first intervention doesnt belong to the first subinterval the 
          log_lambda_1_2 = log12_base + s12;
           log_lambda_2_1 = log21_base;
    P_total *= transition_matrix(first_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}
      // ─────────────────────────────────────────────────────────────────────────
      // (B) MIDDLE SUBINTERVALS: intervals fully between 1st and last global bins
      // ─────────────────────────────────────────────────────────────────────────
      for (m in 1:num_middle_subintervals[n]) {
        int    im = idx_middle[n,m];
        real   s12m = Y_hat_1_2[im];

        // t0m = start of the m-th middle bin, t1m = end of that bin
        real t0m = global_interval_start
                   + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;

        // (B.1) If the first intervention date t_star1 falls inside this middle bin
        // Declare log_lambda variables at the start of the loop
         real log_lambda_1_2;
         real log_lambda_2_1;
        
        if (t0m < t_star1 && t_star1 < t1m && intervention[n] == 1) {
          real d2a = t_star1 - t0m;          // time before first‐intervention in this bin
          real d2b = t1m - t_star1 + 1;  // time after first‐intervention (but before second)
          
       
          // Pre‐first in this middle bin (no interventions yet/ baseline):
          
           log_lambda_1_2 = log12_base + s12m;
              log_lambda_2_1 = log21_base;
          
          P_total *= transition_matrix(
            d2a,
            exp( log_lambda_1_2 ),
            exp(log_lambda_2_1 )
          );

          // Post-intervention part (only first intervention effect)
         log_lambda_1_2 = log12_base + s12m + beta_int1_1;
         log_lambda_2_1 = log21_base +  beta_int2_1;
          P_total *= transition_matrix(d2b, exp(log_lambda_1_2), exp(log_lambda_2_1));
     
       } else {
       // CASE 2: First intervention does NOT lie inside -> use midpoint logic
    real midpoint = (t0m + t1m) / 2;

    if (midpoint >= t_star2) {
      // Both interventions cumulative
    log_lambda_1_2 = log12_base + s12m +  beta_int1_1 +  beta_int1_2;
       log_lambda_2_1 = log21_base +  beta_int2_1 +  beta_int2_2;
    } else if (midpoint >= t_star1) {
      // Only first intervention
      log_lambda_1_2 = log12_base + s12m +  beta_int1_1;
    log_lambda_2_1 = log21_base +  beta_int2_1;
    } else {
      // Baseline
     log_lambda_1_2 = log12_base + s12m;
     log_lambda_2_1 = log21_base;
    }

    P_total *= transition_matrix(interval_length, exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}
      // ─────────────────────────────────────────────────────────────────────────
      // (C) LAST SUBINTERVAL: from the last global‐bin boundary up to date_use[n]
      // ─────────────────────────────────────────────────────────────────────────
      // --- Last subinterval
      // If t_star2 (second intervention) lies within the subinterval:
// Split into two parts:
// Pre-second intervention: Apply the first intervention only .
// Post-second intervention: Apply both interventions.

// --- Last subinterval
if (last_subinterval[n] > 0) {
  int il = idx_last[n];
  real s12l = Y_hat_1_2[il];
  real t1l = date_use[n];
  real t0l = t1l - last_subinterval[n] + 1;
  real log_lambda_1_2;
  real log_lambda_2_1;

  if (t0l < t_star2 && t_star2 < t1l  && intervention[n] == 1) {
    // Split at second intervention
    real d3a = t_star2 - t0l;       // before second intervention
    real d3b = t1l - t_star2 + 1;   // after second intervention

    // --- Pre-second intervention: only first intervention effect
    log_lambda_1_2 = log12_base + s12l +  beta_int1_1;
    log_lambda_2_1 = log21_base + beta_int2_1;
    P_total *= transition_matrix(d3a, exp(log_lambda_1_2), exp(log_lambda_2_1));

    // --- Post-second intervention: both interventions
  log_lambda_1_2 = log12_base + s12l +  beta_int1_1 +  beta_int1_2;
  log_lambda_2_1 = log21_base +  beta_int2_1 +  beta_int2_2;
    P_total *= transition_matrix(d3b, exp(log_lambda_1_2), exp(log_lambda_2_1));

  } else {
    // No split: use midpoint check
    real midpoint = (t0l + t1l) / 2;

    if (midpoint >= t_star2) {
      // Both interventions cumulative
     log_lambda_1_2 = log12_base + s12l +  beta_int1_1 + beta_int1_2;
      log_lambda_2_1 = log21_base + beta_int2_1 + beta_int2_2;
    } else if (midpoint >= t_star1) {
      // Only first intervention
     log_lambda_1_2 = log12_base + s12l +  beta_int1_1;
     log_lambda_2_1 = log21_base +  beta_int2_1;
    } else {
      // Baseline (no interventions)
    log_lambda_1_2 = log12_base + s12l;
     log_lambda_2_1 = log21_base;
    }

    P_total *= transition_matrix(last_subinterval[n], exp(log_lambda_1_2), exp(log_lambda_2_1));
  }
}

      // ─────────────────────────────────────────────────────────────────────────
      // (D) ADD LOG‐LIKELIHOOD CONTRIBUTION FOR THIS TRANSITION
      // ─────────────────────────────────────────────────────────────────────────
      // P_total is now the 2×2 transition probability matrix for going from
      // state at (n−1) to state at n, after multiplying through every “piece”
      // (first/middle/last) with the appropriate interventions flags.
      target += log(
        P_total[ observed_state[n-1], observed_state[n] ] + 1e-9
      );
    }
  }
}

generated quantities {
  vector[num_data] seasonality;
for (i in 1:num_data) {
vector[num_data] Y_hat_1_2 =  to_vector(a_1_2 * B);

}
  }

"
compiled_model <- stan_model(model_code = stan_model_full)

fit <- sampling(
  compiled_model,
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 123
)






