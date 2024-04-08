# Define the dimensions of the matrix
num_individuals <- 5  # Number of individuals
num_days <- 10  # Number of days
lamda1 <- 0.1
lamda2 <- 0.2
recovery_rate <- 0.15  # Assumed recovery rate
num_simulations <- 1000  # Number of simulations

# Create an empty list to store the results of multiple simulations
sim_results <- vector("list", length = num_simulations)

# Simulation loop to model disease progression using a two-state Markov chain
for (sim in 1:num_simulations) {
  # Create an empty matrix to represent infection status for each simulation
  initial_infection_matrix <- matrix("Susceptible", nrow = num_individuals, ncol = num_days)
  
  for (individual in 1:num_individuals) {
    for (day in 2:num_days) {
      current_state <- initial_infection_matrix[individual, day - 1]
      
      # Calculate the force of infection for the specified individual at day
      force_of_infection <- lamda1 + lamda2 * sum(initial_infection_matrix[, day - 1] == "Infected")
      force_of_infection <- pmin(1, pmax(0, force_of_infection))
      
      # Define transition probabilities for Susceptible to Infected based on force of infection
      transition_matrix <- matrix(c(1 - force_of_infection, force_of_infection,  # Transition probability for Susceptible to Infected
                                    0, 1),  # Transition probability for Infected to stay Infected
                                  nrow = 2, byrow = TRUE,
                                  dimnames = list(c("Susceptible", "Infected"), c("Susceptible", "Infected")))
      
      # Transition to the next state based on the adjusted transition probabilities
      next_state <- sample(colnames(transition_matrix), size = 1, prob = transition_matrix[current_state, ])
      
      # Apply recovery
      if (current_state == "Infected" && runif(1) < recovery_rate) {
        initial_infection_matrix[individual, day] <- "Recovered"
      } else {
        initial_infection_matrix[individual, day] <- next_state
      }
    }
  }
  # Store the simulation result in sim_results list
  sim_results[[sim]] <- initial_infection_matrix
}

# Create an empty matrix to store simulation results
simulation_matrix <- array(dim = c(num_individuals, num_days, num_simulations))

# Store simulation results in the matrix
for (sim in 1:num_simulations) {
  simulation_matrix[, , sim] <- ifelse(sim_results[[sim]] == "Susceptible", 0, 
                                       ifelse(sim_results[[sim]] == "Infected", 1, 2))
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
#EXTRA FROM PREVIOUS ONE
# Define the transition matrix
# Define transition probabilities based on infection model parameters
transition_matrix <- matrix(c(0.5, 0.5,  # Probability of staying infected
                              0.6, 0.4),  # Probability of recovery
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("Infected", "Not Infected"), c("Infected", "Not Infected")))
# Number of simulations
num_simulations <- 1000
# Number of time steps
num_steps <- 10
# Simulation of a Markov chain reflecting infection progression
# Set the initial state of an individual based on the initial infection matrix
initial_state <- ifelse(initial_infection_matrix[individual_number, 1] == 1, "Infected", "Not Infected")

# Simulate one Markov chain for an individual
state_sequence <- character(num_days)

for (day in 1:num_days) {
  # Record the current state
  state_sequence[day] <- initial_state
  
  # Transition to the next state based on probabilities
  next_state <- sample(c("Infected", "Not Infected"), size = 1, prob = transition_matrix[initial_state, ])
  
  # Update the current state
  initial_state <- next_state
}

# Now state_sequence represents the individual's status over the simulated days
print(state_sequence)



















transition_matrix <- matrix(c(0.5, 0.5,  # Probability of staying infected
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


