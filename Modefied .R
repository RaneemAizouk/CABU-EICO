#  The code simulates the infection and recovery process over several days for a defined number of individuals.
# Define the dimensions of the matrix
num_individuals <- 3  # Number of rows
num_days <- 10  # Number of columns
lamda1 <- 0.2  # Initialise lamda1 value, background rate
lamda2 <- 0.3  # Set lamda2, the rate associated with infectious individuals within the household
lamda3 <- 0.2  # Set lamda3, the likelihood of individuals being initially infected on the first day
recovery_Prob <- 0.15  # Set the likelihood of an infected individual transitioning from the infected state to the recovered state within a specific timeframe
num_simulations <- 1000

# Create an empty list to store the results of multiple simulations and after running the simulations sim_results stores the outcomes of each simulation in a list[[1]]
#[1] 0 1 1 0 1 0 0 1 1 0

#[[2]]
#[[1] 1 0 0 1 1 0 1 1 1 0 simulation day 1 then days 2 for the individuals in that family



#sim_results <- vector("list", length = num_simulations)

# Create an empty matrix filled with zeros to represent uninfected individuals
initial_infection_matrix <- matrix(0, nrow = num_individuals, ncol = num_days)

# Introduce initial infections for day 1
for (individual in 1:num_individuals) {
  
    random_number <- runif(1)  # Generate a random number between 0 and 1
   
    if (random_number < lamda3)
      # put 1 inside the matrix to refer to first day
      {
      initial_infection_matrix[individual, 1] <- 1  # Set the individual as infected
    }
    
}

# simulate_infection performs the infection spread simulation for a single household. It uses the parameters set earlier to simulate the infection process for the defined number of days.
simulate_infection <- function(){
for (day in 2:num_days) {
  for (individual in 1:num_individuals) {

    was_infected <- initial_infection_matrix[individual, day - 1]  # Store the previous infection status
    
    # Apply the recovery process to previously infected individuals
    if (was_infected == 1 && runif(1) < recovery_Prob) {
      initial_infection_matrix[individual, day] <- 0  # Change the status to recovered
    } else {
      # Check if the individual was not infected at the previous time step
      if (was_infected == 0) {
        force_of_infection <- lamda1 + lamda2 * sum(initial_infection_matrix[, day - 1])
        
        # Calculate the probability of getting infected
        p <- 1 - exp(-force_of_infection)
        
        # Generate a uniform random number [0, 1]
        u_r_n <- runif(1)
        
        # Check if the random value is less than the probability p
        if (u_r_n < p) {
          initial_infection_matrix[individual, day] <- 1  # Change the status to infected
      
        }
      }
    }
  }
}
  return(sum(initial_infection_matrix[, num_days]))  # Return count of infected individuals in the household
}

# New code for multiple simulations and displaying results
simulate_multiple_infections <- function() {
  result <- vector("numeric", length = num_simulations)
  for (sim in 1:num_simulations) {
    result[sim] <- simulate_infection()  # Execute the simulation function
  }
  return(result)
}

# Execute multiple simulations and gather the results
simulation_results <- simulate_multiple_infections()

# Show the final results
print(simulation_results)

# Calculate mean infected individuals per day
mean_infected_per_day <- colMeans(matrix(simulation_results, nrow = num_simulations))
print(mean_infected_per_day)

# Create a density plot
density_plot <- density(simulation_results)

# Plot the density
plot(density_plot, main = "Density Plot of Infected Individuals", xlab = "Number of Infected People", ylab = "Density")

# Create a histogram of the number of infected individuals in households
hist(simulation_results, breaks = seq(0, max(simulation_results) + 0.5, by = 1), 
     xlab = "Number of Infected People", ylab = "Frequency of Infections", 
     main = "Distribution of Infected Individuals in Households")

# Plot the mean count of infected individuals over time
plot(1:num_days, mean_infected_per_day, type = "l", xlab = "Days", ylab = "Mean Infected Count",
     main = "Mean Count of Infected Individuals Over Time")
abline(h = mean(mean_infected_per_day), col = "red")  # Add a horizontal line for the mean value
