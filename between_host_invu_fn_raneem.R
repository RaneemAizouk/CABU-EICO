
betweenhost.stoch = function(t, state, parameters)
{
  #parameters["beta.R"] <- parameters["beta.S"] * (1 - parameters["c"])
  
  lambda1 <- parameters["lambda1"]  # Background infections from outside the household
  lambda2 <- parameters["lambda2"]  # Infection rate within the household
  rho <- parameters["rho"] 
  ## Transitions matrix
  transitions <- list(
    c(S=-1,I=1), # Infection
    c(S=1,I=-1) # Recovery
  )
  ## Rates function
  rates <- function(state, parameters, t)
  {
    with(as.list(c(state, parameters)),{ ### Need source state vars
      return(c((lambda1 + lambda2 * I) * S,
               rho * I))
    })
  }
  out.temp <- as.data.frame(ssa.adaptivetau(state, transitions, rates, parameters,
                                            tf=diff(range(t))))
  out.temp$time <- out.temp$time + min(t)
  out <- cbind(time = t, apply(out.temp[, -1], 2, function(col) {
    approx(x = out.temp[, 1], y = col, xout = t, method = "constant")$y
  }))
  # approx function with the method set to "constant," it implies that if you request the value of the function at a time point between two observed time points, the function will return the value of the nearest observed time point.
  return(as.data.frame(out))
}
getwd()
# Set the working directory to Desktop/Rfiles

