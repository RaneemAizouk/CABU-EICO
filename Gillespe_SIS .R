library(dplyr)

gillespe <- function(lambda1, lambda2, rho, S0, I0, tmax,N) {
  t <- 0 # initial time
  S <- S0 #  initial S state 
  I <- I0 #  initial I state 
  
  ta <- numeric() # Initialise empty numeric vector to store time t
  Sa <- numeric() # Initialise empty numeric vector to store state S
  Ia <- numeric() # Initialise empty numeric vector to store state I
  
  # In Gillespie algorithm, typically update the states within the loop that represents the progression of time.
  # (S and I) are updated within the while loop at each iteration.
  while (t < tmax) {
    ta <- c(ta, t) #append current time to vector ta 
    Sa <- c(Sa, S) #append current state to vector Sa
    Ia <- c(Ia, I) #append current state to vector Ia
    # Calculate rates
    pf1 <- (lambda1 + lambda2 * I) * S # Acquisition rate
    pf2 <- rho * I # Recovery rate
    pf <- pf1 + pf2 # Total rate
    
    dt <- rexp(1, rate = pf) #  Generate a random time step (dt) from an exponential distribution with the totalrate parameter 
    t <- t + dt
  
    
    if (t > tmax) {
      break
    }
   
    #  sample next event 
    ru <- runif(1)   #  a random number that is used to determine which event occurs next.

    if (ru < pf1/pf) {
      S <- S - 1
      I <- I + 1
    } else {
      I <- I - 1
      S <- S + 1
    }
   
  }
  # Create results data frame
  results <- data.frame(t = ta, S = Sa, I = Ia)
  return(results)
}

# Example usage
#lambda1 <- 0.1
#lambda2 <- 0.2
#rho<-0.2
#N <- 9
#S0 <- 7
#I0 <- 2
#tmax <- 14

#result <- gillespe(lambda1, lambda2, rho, S0, I0, tmax)
#print(result)

