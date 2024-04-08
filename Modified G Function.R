gillespe <- function(lambda1, lambda2, rho,   intital_states, tmax) {
  # Initialise variables
  time <- c(0)
  current_states <- intital_states
  states <- matrix(  current_states,nrow=1,ncol=3)
  S <- numeric(0)  # Initialize vector to store the number of susceptibles
  I <- numeric(0)  # Initialize vector to store the number of infected individuals
  
  # Main loop
  t <- 0
  
 
  while (t < tmax) {
    # Calculate number of susceptibles (S) and infected individuals (I)
    S <- sum(current_states == 0)
    I <- sum(current_states == 1)
    
    # Calculate rates
    pf1 <- (lambda1 + lambda2 * I) * S  # Acquisition rate
    pf2 <- rho * I  # Recovery rate
    pf <- pf1 + pf2  # Total rate
    
    # ru is the event could happen next.
    ru <- runif(1)  # Sample for one event
    # Generate a random time step from an exponential distribution
    dt <- rexp(1, rate =(pf))
    t <- t + dt
    # Check if the new_time exceeds tmax before proceeding
    if (t > tmax) break
   
    # Determine the type of event based on ru
    if (ru < pf1 / pf) {
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
    states <- rbind(states, current_states)
    time <- c(time, t)
  }
 
  # Update the number of susceptibles (S) and infected individuals (I) after the loop
  #S <- sum(initial_states == 0)
 # I <- sum(initial_states == 1)
  # Combine time and states into a data frame
  result_df <- data.frame(time = time, X1 = states[, 1], X2 = states[, 2], X3 = states[, 3], state = apply(states, 1, paste, collapse = ""))

  # Set column names
  colnames(result_df) <- c("time", "X1", "X2", "X3", "state")
  
  return(result_df)

}

# Example usage
lambda1 <- 0.2
lambda2 <- 0.3
rho <- 0.1
current_states <- c(0, 1, 1)  # Example initial states
tmax <- 14

result <- gillespe(lambda1, lambda2, rho, current_states, tmax)

print(result)

