# Setup function to assign unique IDs to individuals in each household
setup_households <- function(num_households, individuals_per_household) {
  households <- list()
  id_counter <- 1
  
  for (house in 1:num_households) {
    for (indiv in 1:individuals_per_household) {
      households[[house]][indiv] <- list(
        id = id_counter,
        initial_state = sample(c(0, 1), 1) # Random initial state for example
      )
      id_counter <- id_counter + 1
    }
  }
  
  return(households)
}


# Create households
households <- setup_households(2, 3)  # For example, 2 households with 3 individuals each

# Print the structure of the first household to check
print(households)

