
# Define the dimensions of the matrix
setwd("/Users/raizouk/Desktop/R files")

num_individuals <- 5  # Number of rows
num_days <- 10  # Number of columns

# Create an empty matrix filled with zeros to represent uninfected individuals
initial_infection_matrix <- matrix(0, nrow = num_individuals, ncol = num_days)

# Set the infection status for specific individuals at time step 1
# initial_infection_matrix[1, 2] <- 1  # Row 1, Column 2 is infected
# initial_infection_matrix[2, 3] <- 1  # Row 2, Column 3 is infected
# print(initial_infection_matrix)

# Define the background infection rate from outside the household
lamda1 <- 0.2
lamda1<- 0.3

# Loop through each day
for (day in 2:num_days) {
  # Loop through each individual
  for (individual_number in 1:num_individuals) {
    # Check if the individual is already infected
    if (initial_infection_matrix[individual_number, day - 1] == 1) {
      # Calculate the probability of recovery
      p_recovery <- 0.3
      recovery_probability <- p_recovery
      
      # Generate a uniform random number [0, 1]
      u.r.n <- runif(1)
      
      # Check if the random value (u.r.n) is less than the probability of recovery
      # If yes, set the entry to 0 (recovered)
      if (u.r.n < recovery_probability) {
        initial_infection_matrix[individual_number, day] <- 0
      } else {
        # Individual remains infected, proceed with calculating the force of infection
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
}


# Print the final infection status matrix
print(initial_infection_matrix)
# Define the transition matrix
transition_matrix <- matrix(c(0.7, 0.3,  # Probability of staying infected
                              0.6, 0.4),  # Probability of recovery
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("Infected", "Not Infected"), c("Infected", "Not Infected")))

# Initial state
# initial_state <- "Infected"
# Number of simulations
num_simulations <- 1000
# Number of time steps
num_steps <- 10
# Create a list to store results
simulation_results <- vector("list", num_simulations)
# Loop through simulations
for (sim in 1:num_simulations) {
  # Initialize the current state as random (Infected or Not Infected)
  current_state <- sample(c("Infected", "Not Infected"), size = 1)
  
  # Simulate one Markov chain
  state_sequence <- character(num_steps)
  
  
  for (i in 1:num_steps) {
    # Record the current state
    state_sequence[i] <- current_state
    
    # Transition to the next state based on probabilities
    next_state <- sample(c("Infected", "Not Infected"), size = 1, prob = transition_matrix[current_state, ])
    
    # Update the current state
    current_state <- next_state
  }
  # Store the results for this simulation
  simulation_results[[sim]] <- state_sequence
  print(simulation_results)
}

library(ggplot2)

# Example simulation results (replace with your simulation data)
# simulation_results is a list containing the state sequences
# For instance, simulation_results[[1]] is the state sequence of the first simulation
# Replace this with your actual data
simulation_results <- list(
  c("Infected", "Infected", "Not Infected", "Not Infected", "Infected", "Infected", "Not Infected", "Infected", "Not Infected", "Infected"),
  c("Infected", "Infected", "Infected", "Infected", "NotInfected", " Infected", "Not Infected", "Not Infected", "Not Infected", "Not Infected")
)

# Create a data frame for plotting
simulation_df <- data.frame(
  Simulation = rep(1:length(simulation_results), each = num_steps),
  Time = rep(1:num_steps, times = length(simulation_results)),
  State = unlist(simulation_results)
)

# Create a step plot
ggplot(simulation_df, aes(x = Time, y = Simulation, color = State)) +
  geom_line() +
  labs(x = "Time Step", y = "Simulation", color = "State") +
  scale_color_manual(values = c("Infected" = "red", "Not Infected" = "blue")) +
  theme_minimal()

print(simulation_results[[1]])

