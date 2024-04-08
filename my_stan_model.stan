data {
  int<lower=1> num_individuals;
  int<lower=1> num_days;
  matrix[num_individuals, num_days] initial_infection_matrix;
  real lamda1;  // Background acquisition rate
  real lamda2;  // Additional acquisition rate
}

parameters {
  matrix[num_individuals, num_days] acquisition_status;  // Binary matrix for acquisition
}

model {
  for (day in 2:num_days) {
    for (individual_number in 1:num_individuals) {
      if (acquisition_status[individual_number, day - 1] != 1) {
        real force_of infection;
        force_of_infection = lamda1 + lamda2 * sum(acquisition_status[, day - 1]);
        force_of_infection = fmax(0, fmin(1, force_of_infection));
        real p;
        p = 1 - exp(-force_of_infection);
        
        // Generate a uniform random number [0, 1]
        real u_r_n = uniform_rng(0, 1);
        
        // Check if the random value (u_r_n) is less than the probability p
        // If yes, update the infection status to 1 (infected)
        if (u_r_n < p) {
          acquisition_status[individual_number, day] = 1;
        }
      }
    }
  }
}
