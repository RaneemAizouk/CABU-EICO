# Install the rstan package
install.packages("rstan", dependencies = TRUE)
install.packages("V8", dependencies = TRUE)
install.packages("rstan")

# Load the required library
library(rstan)

# Define the data (initial_infection_matrix) and the number of individuals and days
num_individuals <- 5
num_days <- 10
lamda1 <- 0.1 
lamda2 <-0.2
# Create an empty matrix filled with zeros
initial_infection_matrix <- matrix(0, nrow = num_individuals, ncol = num_days)
# Create an empty matrix filled with zeros
initial_infection_matrix <- matrix(0, nrow = num_individuals, ncol = num_days)

# Introduce random infections
for (i in 1:num_individuals) {
  for (j in 1:num_days) {
    # Generate a random number between 0 and 1
    random_number <- runif(1)
    
    # Define a threshold for introducing infections (e.g., 0.2 for a 20% chance of infection)
    threshold <- 0.2
    
    # If the random number is below the threshold, mark the individual as infected
    if (random_number < threshold) {
      initial_infection_matrix[i, j] <- 1
    }
  }
}

# Compile the Stan model
library(rstan)
model <- stan_model(file = 'code_2.stan')

# Define data
data_list <- list(
  num_individuals = num_individuals,
  num_days = num_days,
  initial_infection_matrix = matrix(0, nrow = num_individuals, ncol = num_days),
  lamda1 = lamda1,
  lamda2 = lamda2
)

# Run the MCMC simulation, pass data to stan 
stan_samples <- sampling(model, data = data_list, iter = 1000, chains = 4)

# Extract and analyse the results
stan_summary <- summary(stan_samples)
print(stan_summary)

