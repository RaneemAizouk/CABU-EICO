gillespe <- function(tvec, parameters, initial_state) {
  # parameters
  lambda1 <- parameters["lambda1"]  # Background infections from outside the household
  lambda2 <- parameters["lambda2"]  # Infection rate within the household
  rho <- parameters["rho"] # recovery rate
  I <- initial_state["I"] # initial infection state 
  
  # Initialization of time and state vectors
  t0 <- min(tvec)
  X <- initial_state
  
  states <- matrix( nrow =length(tvec) , ncol = length(initial_state), dimnames = list(tvec, names(initial_state)))
  
  # Transitions matrix
  transitions <- matrix(c(1, -1, -1, 1), ncol = 2, byrow = TRUE, dimnames = list(c("Infection", "Recovery"), c("S", "I")))
  
  # Rates function
  ratefun <- function(X1, parameters) {
    with(as.list(c(X1, parameters)), {
      return(c((lambda1 + lambda2 * X1["I"]) * X1["S"],rho * X1["I"]))
    })
  }
  states[1,]<- initial_state # assigning the values of the initial_state vector to the first row of the states matrix.
  for (ctr in 1:(length(tvec)-1)){
    condition = t0<tvec[ctr+1]
    while(condition) {       
      # Calculate rates
      rates <- ratefun(X, parameters)
      
      # Check if all rates are zero, then break
      if (all(rates == 0)) break
      
      # Sample time until the next event
      totrate <- sum(rates)
      elapsed <- rexp(1, totrate)
      
      # Pick transition
      which.trans <- sample(1:nrow(transitions), size = 1, prob = rates / totrate)
      # Update time and state based on the chosen event
      t0 <- t0 + elapsed
      X <- X + transitions[which.trans, ]
      condition = t0<tvec[ctr+1]
    }
    states[ctr+1, ] <- X
  } 
  # Return the results
  return(data.frame(time = tvec, states))
}

# Example usage
tvec <- seq(0, 14, 1)
parameters <- c(lambda1 = 0.2, lambda2 = 0.3, rho = 0.2)
initial_state <- c(I = 2, S = 7)

result <- gillespe(tvec, parameters, initial_state)
print(result)
