
functions {
  // 2×2 CTMC transition matrix for a 2-state chain
  matrix transition_matrix(real t,
                           real lambda_1_2,
                           real lambda_2_1) {
    real total_lambda = lambda_1_2 + lambda_2_1;
    real exp_term     = exp(-total_lambda * t);

    matrix[2,2] P;
    P[1,1] = (lambda_2_1 / total_lambda) + (lambda_1_2 / total_lambda) * exp_term;
    P[2,2] = (lambda_1_2 / total_lambda) + (lambda_2_1 / total_lambda) * exp_term;
    P[1,2] = (lambda_2_1 / total_lambda) * (1 - exp_term);
    P[2,1] = (lambda_1_2 / total_lambda) * (1 - exp_term);
    return P;
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
  //array[N] int intervention_village;
 


  int<lower=1> num_data;
  array[num_data] real X;

  //array[N] real Intervention_start_date;
  array[N] real intervention_date;  // or int if stored as day numbers
  array[N] real intervention_date2; 
  int<lower=1> max_middle;
  int<lower=1> num_intervals;


}

transformed data {
  
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
  real beta_int1;   // first‐intervention effect
  real beta_int2;   // second‐intervention effect

  // Replace spline parameters with sinusoid coefficients:
  real a1;  // Amplitude for sin(2piX)
  real b1;  // Amplitude for cos(2piX)


}

transformed parameters {
  real q_1_2_base = -4.79 + q_1_2_raw * sigma_q_1_2;
  real q_2_1_base = -4.71 + q_2_1_raw * sigma_q_2_1;
  vector[H] u = u_raw * sigma_u;
  
  // Define seasonality as sinusoidal:
  vector[num_data] Y_hat_1_2;
  for (i in 1:num_data) {
    Y_hat_1_2[i] = a1 * sin(2 * pi() * X[i]) + b1 * cos(2 * pi() * X[i]);
  }


 
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

  
  // Noise on continuous‐time transitions
  sigma ~ normal(0, 1);

  // Intervention‐specific effects (first and second)
  beta_int1 ~ normal(0, 1);
  beta_int2 ~ normal(0, 1);
  // New priors for sinusoidal coefficients:
   a1 ~ normal(0, 0.5);
   b1 ~ normal(0, 0.5);

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
      real lambda12 = q_1_2_base
                     + u[HouseID[n]]
                     + beta_1_2_age  * age[n]
                     + beta_1_2_sexe * sexe[n];
      real lambda21 = q_2_1_base
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
       // real   s21 = Y_hat_2_1[i1];  // seasonal adjustment at that midpoint for recovery

        // t0 = the actual calendar date of observation (n−1)
        real t0    = date_use[n-1];
        // t1 = end of “first subinterval,” measured in days.
        // (Because first_subinterval[n] = length of that piece in days.)
        real t1    = t0 + first_subinterval[n] - 1;

        // (A.1) If t_star1 lies strictly inside [t0, t1], we split this subinterval
        if (t0 < t_star1 && t_star1 < t1) {
          real d1a = t_star1 - t0;          // # days before first‐intervention
          real d1b = t1      - t_star1 + 1; // # days after first‐intervention (but before second)

          // — Pre‐first piece (no intervention active):
          P_total *= transition_matrix(
            d1a,
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 + s12),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2 )
          );

          // — Post‐first piece (first intervention active, but second still not yet active):
          // Compute midpoint of [t_star1, t1]
          real midpoint_post1 = (t_star1 + t1) / 2;
          // Because we are after t_star1, flag1_post1 = 1
          int  flag1_post1    = 1;
          // flag2_post1 = 1 only if that midpoint ≥ t_star2 (rare, since t_star2 >> t_star1)
          int  flag2_post1    = (midpoint_post1 >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            d1b,
            exp(lambda12 + flag1_post1*beta_int1 + flag2_post1*beta_int2 + s12),
            exp(lambda21 + flag1_post1*beta_int1 + flag2_post1*beta_int2 )
          );

        } else {
          // (A.2) No split by t_star1 inside this “first” piece
          real midpoint1 = (t0 + t1) / 2;
          int  flag1_1   = (midpoint1 >= t_star1) ? 1 : 0;
          int  flag2_1   = (midpoint1 >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            first_subinterval[n],
            exp(lambda12 + flag1_1*beta_int1 + flag2_1*beta_int2 + s12),
            exp(lambda21 + flag1_1*beta_int1 + flag2_1*beta_int2 )
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (B) MIDDLE SUBINTERVALS: intervals fully between 1st and last global bins
      // ─────────────────────────────────────────────────────────────────────────
      for (m in 1:num_middle_subintervals[n]) {
        int    im = idx_middle[n,m];
        real   s12_m = Y_hat_1_2[im];
        //real   s21_m = Y_hat_2_1[im];

        // t0m = start of the m-th middle bin, t1m = end of that bin
        real t0m = global_interval_start
                   + (global_interval_index_start[n] + m - 1) * interval_length;
        real t1m = t0m + interval_length - 1;

        // (B.1) If the first intervention date t_star1 falls inside this middle bin
        if (t0m < t_star1 && t_star1 < t1m) {
          real d2a = t_star1 - t0m;          // time before first‐intervention in this bin
          real d2b = t1m     - t_star1 + 1;  // time after first‐intervention (but before second)

          // — Pre‐first in this middle bin (no interventions yet):
          P_total *= transition_matrix(
            d2a,
            exp(lambda12 + 0*beta_int1 + 0*beta_int2 + s12_m),
            exp(lambda21 + 0*beta_int1 + 0*beta_int2)
          );

          // — Post‐first in this same bin (first active, second not yet active):
          real midpoint_post1m = (t_star1 + t1m) / 2;
          int  flag1_post1m    = 1;
          int  flag2_post1m    = (midpoint_post1m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            d2b,
            exp(lambda12 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 + s12_m),
            exp(lambda21 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 )
          );

        } else {
          // (B.2) No split by t_star1 inside this middle bin
          real midpoint2m = (t0m + t1m) / 2;
          int  flag1_2m   = (midpoint2m >= t_star1) ? 1 : 0;
          int  flag2_2m   = (midpoint2m >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            interval_length,
            exp(lambda12 + flag1_2m*beta_int1 + flag2_2m*beta_int2 + s12_m),
            exp(lambda21 + flag1_2m*beta_int1 + flag2_2m*beta_int2 )
          );
        }
      }

      // ─────────────────────────────────────────────────────────────────────────
      // (C) LAST SUBINTERVAL: from the last global‐bin boundary up to date_use[n]
      // ─────────────────────────────────────────────────────────────────────────
      if (last_subinterval[n] > 0) {
        int    il = idx_last[n];
        real   s12_l = Y_hat_1_2[il];
       // real   s21_l = Y_hat_2_1[il];

        real t1l = date_use[n];
        real t0l = t1l - last_subinterval[n] + 1;

        // By design, t0l > t_star1, since the “last subinterval” always starts
        // after the end of the first‐intervention window.  Therefore, we skip
        // any check for t_star1 here.  The only possible split in this last piece
        // is around t_star2 (September 2023), which may or may not lie inside [t0l, t1l].

        // (C.1) If second intervention falls inside this last subinterval:
        if (t0l < t_star2 && t_star2 < t1l) {
          real d3a = t_star2 - t0l;          // time from t0l until just before t_star2
          real d3b = t1l     - t_star2 + 1;  // time from t_star2 through t1l

          // — Pre‐second piece:  first intervention is already active on entire last piece —
          //   (Because t0l > t_star1 implies we are already post‐first for all days in [t0l, t_star2).)
          real flag1_pre2 = 1;
          real flag2_pre2 = 0;
          P_total *= transition_matrix(
            d3a,
            exp(lambda12 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 + s12_l),
            exp(lambda21 + flag1_pre2*beta_int1 + flag2_pre2*beta_int2 )
          );

          // — Post‐second piece: both interventions active —
          real flag1_post2 = 1;
          real flag2_post2 = 1;
          P_total *= transition_matrix(
            d3b,
            exp(lambda12 + flag1_post2*beta_int1 + flag2_post2*beta_int2 + s12_l),
            exp(lambda21 + flag1_post2*beta_int1 + flag2_post2*beta_int2 )
          );

        } else {
          // (C.2) No split by t_star2 inside this last piece:
          //    first intervention is already active (flag1 = 1).
          //    second intervention may or may not yet be active, depending on midpoint.
          real midpoint_l = (t0l + t1l) / 2;
          int  flag1_l    = 1;  // always post‐first
          int  flag2_l    = (midpoint_l >= t_star2) ? 1 : 0;

          P_total *= transition_matrix(
            last_subinterval[n],
            exp(lambda12 + flag1_l*beta_int1 + flag2_l*beta_int2 + s12_l),
            exp(lambda21 + flag1_l*beta_int1 + flag2_l*beta_int2 )
          );
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
  vector[N]             log_lambda_1_2_out;
  vector[N]             log_lambda_2_1_out;
  vector[num_data]      Y_hat_1_2_out;
 // vector[num_data]      Y_hat_2_1_out;
  int                   y_rep[N];
  

  
  // Recompute seasonal effect with sin/cos
  for (i in 1:num_data) {
  Y_hat_1_2_out[i] = a1 * sin(2 * pi() * X[i]) + b1 * cos(2 * pi() * X[i]);
}

  for (n in 1:N) {
    // 1) Baseline log‐rates (same as in transformed parameters)
real log12 = q_1_2_base
+ u[HouseID[n]]
+ beta_1_2_age  * age[n]
+ beta_1_2_sexe * sexe[n];
real log21 = q_2_1_base
+ u[HouseID[n]]
+ beta_2_1_age  * age[n]
+ beta_2_1_sexe * sexe[n];

// First observation of each individual: no transition
if (n == 1 || menage_id_member[n] != menage_id_member[n-1]) {
  y_rep[n] = -1;  // no simulated transition
  // Record the “initial” rates (use idx_first to pick spline)
  int i0 = idx_first[n];
  log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[i0];
  log_lambda_2_1_out[n] = log21 ;
  continue;
}

// Otherwise, build P_total exactly as in the model block:
  matrix[2,2] P_total = diag_matrix(rep_vector(1.0, 2));

real t_star1 = intervention_date[n];   // first intervention date
real t_star2 = intervention_date2[n];  // second intervention date

// ─────────────────────────────────────────────────────────────────────
// (A) FIRST SUBINTERVAL: [date_use[n−1], date_use[n−1] + first_subinterval[n] − 1]
// ─────────────────────────────────────────────────────────────────────
if (first_subinterval[n] > 0) {
  int    i1  = idx_first[n];
  real   s12 = Y_hat_1_2_out[i1];
  //real   s21 = Y_hat_2_1_out[i1];
  
  real t0 = date_use[n-1];
  real t1 = t0 + first_subinterval[n] - 1;
  
  // (A.1) split if t_star1 lies inside
  if (t0 < t_star1 && t_star1 < t1) {
    real d1a = t_star1 - t0;
    real d1b = t1      - t_star1 + 1;
    
    // Pre‐first: no interventions
    P_total *= transition_matrix(
      d1a,
      exp(log12 + 0*beta_int1 + 0*beta_int2 + s12),
      exp(log21 + 0*beta_int1 + 0*beta_int2 )
    );
    
    // Post‐first (but still pre‐second)
    real midpoint_post1 = (t_star1 + t1) / 2;
    int  flag1_post1    = 1;
    int  flag2_post1    = (midpoint_post1 >= t_star2) ? 1 : 0;
    
    P_total *= transition_matrix(
      d1b,
      exp(log12 + flag1_post1*beta_int1 + flag2_post1*beta_int2 + s12),
      exp(log21 + flag1_post1*beta_int1 + flag2_post1*beta_int2 )
    );
    
  } else {
    // (A.2) no split by t_star1
    real midpoint1 = (t0 + t1) / 2;
    int  flag1_1   = (midpoint1 >= t_star1) ? 1 : 0;
    int  flag2_1   = (midpoint1 >= t_star2) ? 1 : 0;
    
    P_total *= transition_matrix(
      first_subinterval[n],
      exp(log12 + flag1_1*beta_int1 + flag2_1*beta_int2 + s12),
      exp(log21 + flag1_1*beta_int1 + flag2_1*beta_int2 )
    );
  }
}

// ─────────────────────────────────────────────────────────────────────
// (B) MIDDLE SUBINTERVALS
// ─────────────────────────────────────────────────────────────────────
for (m in 1:num_middle_subintervals[n]) {
  int    im   = idx_middle[n,m];
  real   s12m = Y_hat_1_2_out[im];
  //real   s21m = Y_hat_2_1_out[im];
  
  real t0m = global_interval_start
  + (global_interval_index_start[n] + m - 1) * interval_length;
  real t1m = t0m + interval_length - 1;
  
  // (B.1) if t_star1 falls inside this bin
  if (t0m < t_star1 && t_star1 < t1m) {
    real d2a = t_star1 - t0m;
    real d2b = t1m     - t_star1 + 1;
    
    // Pre‐first in this middle bin
    P_total *= transition_matrix(
      d2a,
      exp(log12 + 0*beta_int1 + 0*beta_int2 + s12m),
      exp(log21 + 0*beta_int1 + 0*beta_int2 )
    );
    
    // Post‐first (but pre‐second)
    real midpoint_post1m = (t_star1 + t1m) / 2;
    int  flag1_post1m    = 1;
    int  flag2_post1m    = (midpoint_post1m >= t_star2) ? 1 : 0;
    
    P_total *= transition_matrix(
      d2b,
      exp(log12 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 + s12m),
      exp(log21 + flag1_post1m*beta_int1 + flag2_post1m*beta_int2 )
    );
    
  } else {
    // (B.2) no split by t_star1 in this bin
    real midpoint2m = (t0m + t1m) / 2;
    int  flag1_2m   = (midpoint2m >= t_star1) ? 1 : 0;
    int  flag2_2m   = (midpoint2m >= t_star2) ? 1 : 0;
    
    P_total *= transition_matrix(
      interval_length,
      exp(log12 + flag1_2m*beta_int1 + flag2_2m*beta_int2 + s12m),
      exp(log21 + flag1_2m*beta_int1 + flag2_2m*beta_int2 )
    );
  }
}

// ─────────────────────────────────────────────────────────────────────
// (C) LAST SUBINTERVAL: entirely after t_star1, may split around t_star2
// ─────────────────────────────────────────────────────────────────────
if (last_subinterval[n] > 0) {
  int    il   = idx_last[n];
  real   s12l = Y_hat_1_2_out[il];
 // real   s21l = Y_hat_2_1_out[il];
  
  real t1l = date_use[n];
  real t0l = t1l - last_subinterval[n] + 1;
  
  // (C.1) if t_star2 lies inside this last piece
  if (t0l < t_star2 && t_star2 < t1l) {
    real d3a = t_star2 - t0l;
    real d3b = t1l     - t_star2 + 1;
    
    // Pre‐second: first intervention is already active (t0l > t_star1)
    P_total *= transition_matrix(
      d3a,
      exp(log12 + 1*beta_int1 + 0*beta_int2 + s12l),
      exp(log21 + 1*beta_int1 + 0*beta_int2 )
    );
    
    // Post‐second: both interventions active
    P_total *= transition_matrix(
      d3b,
      exp(log12 + 1*beta_int1 + 1*beta_int2 + s12l),
      exp(log21 + 1*beta_int1 + 1*beta_int2 )
    );
    
  } else {
    // (C.2) no split by t_star2; entire last piece has first active,
    // second may or may not be active depending on midpoint
    real midpoint_l = (t0l + t1l) / 2;
    int  flag1_l    = 1;  // always post‐first here
    int  flag2_l    = (midpoint_l >= t_star2) ? 1 : 0;
    
    P_total *= transition_matrix(
      last_subinterval[n],
      exp(log12 + flag1_l*beta_int1 + flag2_l*beta_int2 + s12l),
      exp(log21 + flag1_l*beta_int1 + flag2_l*beta_int2 )
    );
  }
}

// Draw the next state for observation n
y_rep[n] = categorical_rng(
  to_vector(P_total[ observed_state[n-1] ]) 
  / sum(P_total[ observed_state[n-1] ])
);

// Record the log‐hazards at the very end of this transition (using idx_last):
  log_lambda_1_2_out[n] = log12 + Y_hat_1_2_out[idx_last[n]];
log_lambda_2_1_out[n] = log21 ;
  }
}

