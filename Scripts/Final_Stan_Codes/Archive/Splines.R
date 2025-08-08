# Example data
num_knots <- 5
knots <- c(0, 0.25, 0.5, 0.75, 1)  # Knots in [0, 1]
spline_degree <- 3
X <- seq(0, 1.38, length.out = 100)  # Time points for ~15 months
num_data <- length(X)

stan_data <- list(
  num_knots = num_knots,
  knots = knots,
  spline_degree = spline_degree,
  num_data = num_data,
  X = X
)

stan_model <- "
functions {
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
  int<lower=1> num_knots;
  vector[num_knots] knots;
  int<lower=0> spline_degree;
  int<lower=1> num_data;
  array[num_data] real X;
}

transformed data {
  int num_basis = num_knots + spline_degree;
  matrix[num_data, num_basis] B;
  array[num_data] real X_mod;

  for (i in 1:num_data)
    X_mod[i] = fmod(X[i], 1.0); // wrap around (periodic)

  vector[2 * spline_degree + num_knots + 1] ext_knots;
  ext_knots = append_row(
              knots[(num_knots - spline_degree + 1):num_knots] - 1,
              append_row(knots, knots[1:(spline_degree+1)] + 1)
            );

  for (ind in 1:num_basis)
    B[:, ind] = build_b_spline(X_mod, to_array_1d(ext_knots), ind, spline_degree + 1);
}

parameters {
  row_vector[num_basis - spline_degree] a_raw_1_2_free;
  real log_tau_raw_1_2;
  real a0_raw_1_2;
  real<lower=0> sigma_a0_1_2;

  row_vector[num_basis - spline_degree] a_raw_2_1_free;
  real log_tau_raw_2_1;
  real a0_raw_2_1;
  real<lower=0> sigma_a0_2_1;
}

transformed parameters {
  row_vector[num_basis] a_raw_1_2;
  row_vector[num_basis] a_raw_2_1;

  real tau_1_2 = exp(log_tau_raw_1_2);
  real tau_2_1 = exp(log_tau_raw_2_1);

  for (i in 1:(num_basis - spline_degree))
    a_raw_1_2[i] = a_raw_1_2_free[i];
  for (k in 1:spline_degree)
    a_raw_1_2[num_basis - spline_degree + k] = a_raw_1_2[k];

  for (i in 1:(num_basis - spline_degree))
    a_raw_2_1[i] = a_raw_2_1_free[i];
  for (k in 1:spline_degree)
    a_raw_2_1[num_basis - spline_degree + k] = a_raw_2_1[k];
}

model {
  // Priors for infection spline
  a_raw_1_2_free ~ normal(0, 1);
  log_tau_raw_1_2 ~ normal(0, 0.5);
  a0_raw_1_2 ~ normal(0, 0.5);
  sigma_a0_1_2 ~ normal(0, 0.05);

  // Priors for recovery spline
  a_raw_2_1_free ~ normal(0, 2);
  log_tau_raw_2_1 ~ normal(0, 0.5);
  a0_raw_2_1 ~ normal(0, 0.5);
  sigma_a0_2_1 ~ normal(0, 0.05);

  // Smoothness prior
  for (k in 3:(num_basis - spline_degree)) {
    target += normal_lpdf(a_raw_1_2_free[k] - 2 * a_raw_1_2_free[k-1] + a_raw_1_2_free[k-2] | 0, 0.5);
    target += normal_lpdf(a_raw_2_1_free[k] - 2 * a_raw_2_1_free[k-1] + a_raw_2_1_free[k-2] | 0, 0.5);
  }
}

generated quantities {
  row_vector[num_basis] a_1_2 = a_raw_1_2 * exp(log_tau_raw_1_2);
  row_vector[num_basis] a_2_1 = a_raw_2_1 * exp(log_tau_raw_2_1);

  real a0_1_2 = a0_raw_1_2 * sigma_a0_1_2;
  real a0_2_1 = a0_raw_2_1 * sigma_a0_2_1;

  vector[num_data] Y_hat_1_2 = a0_1_2 * to_vector(X) + to_vector(a_1_2 * B');
  vector[num_data] Y_hat_2_1 = a0_2_1 * to_vector(X) + to_vector(a_2_1 * B');

  vector[num_data] Y_hat_1_2_out = Y_hat_1_2;
  vector[num_data] Y_hat_2_1_out = Y_hat_2_1;
}

"

library(rstan)
writeLines(stan_model, "periodic_splines.stan")
compiled_model <- stan_model(file = "periodic_splines.stan")

fit <- sampling(compiled_model, data = stan_data, iter = 2000, warmup = 1000, chains = 4, cores = 4, seed = 123)

library(ggplot2)

# ---- Extract posterior samples ----
post <- rstan::extract(fit)

# Sample spline curves for plotting
num_samples <- 50
sample_idx <- sample(1:nrow(post$Y_hat_2_1_out), num_samples)

plot_data <- data.frame()
for (i in sample_idx) {
  temp_df <- data.frame(
    X = stan_data$X,
    Y = post$Y_hat_2_1_out[i, ],
    sample = paste0("sample_", i)
  )
  plot_data <- rbind(plot_data, temp_df)
}

# Plot
ggplot(plot_data, aes(x = X, y = Y, group = sample)) +
  geom_line(alpha = 0.2, color = "red") +
  labs(
    title = "Recovery Spline (Y_hat_2_1) - Prior Predictive Samples",
    x = "Time",
    y = "Spline Value"
  ) +
  theme_minimal()
#########################################winter effect
