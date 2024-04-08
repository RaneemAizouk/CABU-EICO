# Function to simulate sampling with imperfect sensitivity
sample_household <- function(sim_result, sensitivity) {
  # Simulate sampling with imperfect sensitivity
  sampled_S <- rbinom(1, sim_result$S, sensitivity)
  sampled_I <- pmin(sim_result$S, sampled_S)
  
  # sampled_I <- sim_result$S - sampled_S
  return(data.frame(t = sim_result$t, S = sampled_S, I = sampled_I))
}

# Example usage
# Assuming sim_result is a data frame with columns t, S, and I representing simulation results
sim_result <- data.frame(t = c(1, 2, 3), S = c(5, 4, 3), I = c(1, 2, 3))

# Specify sensitivity (e.g., 90%)
sensitivity <- 0.9

# Simulate sampling
sampled_data <- sample_household(sim_result, sensitivity)

# Display sampled data
print(sampled_data)

