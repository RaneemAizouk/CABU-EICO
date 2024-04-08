gillespe1 <- function(lambda3, lambda4, rho1, intital_states, tmax) {
  # Initialise variables
  time <- c(0)
  current_states <- intital_states
  states_m <- matrix( current_states,nrow=1,ncol=3)
  S_m <- numeric(0)  # Initialize vector to store the number of susceptibles
  I_m <- numeric(0)  # Initialize vector to store the number of infected individuals
  
  # Main loop
  t <- 0
  
  
  while (t < tmax) {
    # Calculate number of susceptibles (S) and infected individuals (I)
    S_m <- sum(current_states == 0)
    I_m <- sum(current_states == 1)
    
    # Calculate rates
    pf3 <- (lambda3 + lambda4 * I_m) * S_m  # Acquisition rate
    pf4 <- rho1 * I_m  # Recovery rate
    pf_m <- pf3 + pf4  # Total rate
    
    # ru is the event could happen next.
    ru <- runif(1)  # Sample for one event
    # Generate a random time step from an exponential distribution
    dt <- rexp(1, rate =(pf_m))
    t <- t + dt
    # Check if the new_time exceeds tmax before proceeding
    if (t > tmax) break
    
    # Determine the type of event based on ru
    if (ru < pf3 / pf_m) {
      # Infection event
      susceptibles <- which(current_states == 0)
      if (length(susceptibles) > 1) {
        # If there are susceptibles, randomly select one for infection
        individual_to_update <- sample(susceptibles, size = 1)
        current_states[individual_to_update] <- 1
      } else {
        current_states[susceptibles] <- 1
      }
    } else {
      # Recovery event
      infected <- which(current_states == 1)
      if (length(infected) > 1) {
        # If there are infected individuals, randomly select one for recovery
        individual_to_update <- sample(infected, size = 1)
        current_states[individual_to_update] <- 0
      } else {
        current_states[infected] <- 0
      }
    }
    
    # Record states
    # if (all(states[nrow(states),] == initial_states)) print(paste(ru, pf1/pf))
    states_m <- rbind(states_m, current_states)
    time <- c(time, t)
  }
  
  # Update the number of susceptibles (S) and infected individuals (I) after the loop
  #S <- sum(initial_states == 0)
  # I <- sum(initial_states == 1)
  # Combine time and states into a data frame
  result_df1 <- data.frame(time = time, X1 = states_m[, 1], X2 = states_m[, 2], X3 = states_m[, 3], state_m = apply(states_m, 1, paste, collapse = ""))
  
  # Set column names
  colnames(result_df1) <- c("time", "X1", "X2", "X3", "state")
  
  return(result_df1)
  
}

# Example usage
lambda3 <- 0.2
lambda4 <- 0.3
rho1 <- 0.1
current_states <- c(0, 1, 1)  # Example initial states
tmax <- 14

result2 <- gillespe1(lambda3, lambda4, rho1, current_states, tmax)

print(result2)

