# Define the dimensions of the matrix
num_individuals <- 5  # Number of rows
num_days <- 10  # Number of columns
lamda1 <- 0.2
lamda2 <- 0.3
recovery_rate <- 0.15 # Rate at which individuals recover
num_simulations <- 1000

# Create an empty list to store the results of multiple simulations
sim_results <- vector("list", length = num_simulations)

# Create an empty matrix filled with zeros to represent uninfected individuals
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

# Loop through each day
for (day in 2:num_days) {
  # Loop through each individual
  for (individual_number in 1:num_individuals) {
    # Check if the individual is already infected
    if (initial_infection_matrix[individual_number, day - 1] == 1) {
      # If already infected, keep them infected
      initial_infection_matrix[individual_number, day] <- 1
    } else {
      # Calculate the force of infection for the specified individual at day
      force_of_infection <- lamda1 + lamda2 * sum(initial_infection_matrix[, day - 1])
      
      # Limit the force_of_infection to be between 0 and 1
      force_of_infection <- pmin(1, pmax(0, force_of_infection))
      
      # Calculate the probability of getting infected
      p <- 1 - exp(-force_of_infection)
      
      # Generate a uniform random number [0, 1]
      u.r.n <- runif(1)
      
      # Check if the random value (u.r.n) is less than the probability p
      # If yes, update the infection status to 1 (infected)
      if (u.r.n < p) {
        initial_infection_matrix[individual_number, day] <- 1
      }
    }
  }
}
# Apply the recovery process for each individual
for (individual_number in 1:num_individuals) {
  if (initial_infection_matrix[individual_number, day] == 1 && runif(1) < recovery_rate) {
    initial_infection_matrix[individual_number, day] <- 0  # Change the status to recovered (0)
  }
}

# Print the final infection status matrix
print(initial_infection_matrix)

# Store the simulation result in sim_results list
for (sim in 1:num_simulations) {
  sim_results[[sim]] <- initial_infection_matrix
}

# Create an empty matrix to store simulation results
simulation_matrix <- array(dim = c(num_individuals, num_days, num_simulations))

# Store simulation results in the matrix
for (sim in 1:num_simulations) {
  simulation_matrix[, , sim] <- sim_results[[sim]]
}

# Calculate the percentage of infected individuals for each day and simulation
infected_percentage <- apply(simulation_matrix, c(2, 3), function(x) mean(x == 1))

# Plot the results
matplot(1:num_days, infected_percentage, type = "l", xlab = "Days", ylab = "Percentage Infected",
        main = "Infection Spread Over 10 Days for 1000 Simulations")

# Prepare the scatter plots for each day
par(mfrow = c(2, 5))  # Divides the plotting area into a 2x5 grid for 10 plots

for (day_to_plot in 1:num_days) {
  # Calculate the percentage of infected individuals at the chosen day for each simulation
  infected_percentage_at_day <- infected_percentage[day_to_plot, ]
  
  # Create scatter plots for each day
  plot(rep(day_to_plot, num_simulations), infected_percentage_at_day, 
       xlab = paste("Day", day_to_plot), ylab = "Percentage Infected",
       main = paste("Percentage of Infected Individuals at Day", day_to_plot, "for 1000 Simulations"),
       pch = 19, col = "blue")  # Using blue color for dots
}


# Calculate the recovery status for each individual in each simulation
recovery_matrix <- array(dim = c(num_individuals, num_days, num_simulations))

for (sim in 1:num_simulations) {
  for (individual in 1:num_individuals) {
    for (day in 1:num_days) {
      recovery_matrix[individual, day, sim] <- ifelse(sim_results[[sim]][individual, day] == 0, 1, 0)
    }
  }
}








