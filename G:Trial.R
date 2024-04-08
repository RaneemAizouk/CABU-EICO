gillespe <- function(lambda1, lambda2, rho, S0, I0, tmax, N) {
  t <- 0
  S <- S0
  I <- I0
  
  ta <- numeric()
  Sa <- numeric()
  Ia <- numeric()
  HouseIDa <- numeric()  # New vector to store individual IDs
  
  while (t < tmax) {
    ta <- c(ta, t)
    Sa <- c(Sa, S)
    Ia <- c(Ia, I)
    HouseIDa <- c(HouseIDa, rep(seq_along(S), each = 1))  # Assign individual IDs
    
    pf1 <- (lambda1 + lambda2 * I) * S
    pf2 <- rho * I
    pf <- pf1 + pf2
    dt <- rexp(1, rate = pf)
    t <- t + dt
    
    if (t > tmax) {
      break
    }
    
    ru <- runif(1)
    if (ru < pf1/pf) {
      S <- S - 1
      I <- I + 1
    } else {
      I <- I - 1
      S <- S + 1
    }
  }
  
  results <- data.frame(t = ta, S = Sa, I = Ia, HouseID = HouseIDa)

  return(results)
}
# Example usage
lambda1 <- 0.1
lambda2 <- 0.2
rho <- 0.2
N <- 9
S0 <- 7
I0 <- 2
tmax <- 14

result <- gillespe(lambda1, lambda2, rho, S0, I0, tmax, N)

# Extract HouseIDa vector from the result
HouseIDa <- result$HouseID

# Print the HouseIDa vector
print(HouseIDa)

