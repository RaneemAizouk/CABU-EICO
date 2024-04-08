# Define parameters
num_individuals <- 10
num_days <- 1000
lamda1 <- 0.2
lamda2 <- 0.3
lamda3 <- 0.2
recovery_Prob <- 0.15
num_simulations <- 1000

# Create an empty matrix to store the status of individuals
initial_infection_matrix <- matrix(0, nrow = num_individuals, ncol = num_days)

# Simulate infections
for (individual in 1:num_individuals) {
  random_number <- runif(1)
  if (random_number < lamda3) {
    initial_infection_matrix[individual, 1] <- 1
  }
}
# Create a function to simulate infection spread for a single simulation
simulate_infection <- function() {
  initial_infection <- sample(c(0, 1), num_individuals, replace = TRUE, prob = c(1 - lamda3, lamda3))
  
  for (day in 2:num_days) {
    for (individual in 1:num_individuals) {
      was_infected <- initial_infection[individual]
      if (was_infected == 1 && runif(1) < recovery_Prob) {
        initial_infection[individual] <- 0
      } else {
        if (was_infected == 0) {
          force_of_infection <- lamda1 + lamda2 * sum(initial_infection)
          p <- 1 - exp(-force_of_infection)
          if (runif(1) < p) {
            initial_infection[individual] <- 1
          }
        }
      }
    }
  }
  sum(initial_infection)
}

# Simulate infection spread over multiple simulations and count infected individuals in households
infected_counts <- replicate(num_simulations, simulate_infection())

# Create a histogram of the number of infected individuals in households
hist(infected_counts, breaks = seq(-0.5, num_individuals + 0.5, by = 1), xlab = "Number of Infected Individuals",
     ylab = "Frequency", main = "Distribution of Infected Individuals in Households")

