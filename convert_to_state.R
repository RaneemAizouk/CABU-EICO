binary_vector <- c("000", "001", "010","001","110","101","011","111")
decimal_vector <- sapply(binary_vector, function(x) strtoi(x, base = 2))

print("Binary numbers:")
print(binary_vector)
print("Decimal numbers:")
print(decimal_vector)

convert_to_state <- function(states_df) {
  state_numeric <- apply(states_df, 1, function(states) {
    state_str <- paste(states, collapse = "")
    return(ifelse(state_str == "111", 1,
                  ifelse(state_str == "000", 2,
                         ifelse(state_str == "110", 3,
                                ifelse(state_str == "101", 4,
                                       ifelse(state_str == "011", 5,
                                              ifelse(state_str == "001", 6,
                                                     ifelse(state_str == "010", 7, 8))))))))
  })
  return(state_numeric)
}

# Example usage
initial_states <- c(0, 1, 1)
result_states <- matrix(c(0, 0, 1, 1, 0, 1, 1, 1, 1), ncol = 3)
result_df <- data.frame(X1 = result_states[, 1], X2 = result_states[, 2], X3 = result_states[, 3])

# Convert states to decimal
result_df$state_numeric <- convert_to_state(result_df)

# Print the result
print(result_df)


# Define a function to generate individual IDs based on household and individual numbers
generate_individual_id <- function(household_id, individual_id) {
  return(paste(household_id, individual_id, sep = "-"))
}

# Define a class for Individual
Individual <- setClass(
  "Individual",
  slots = list(
    household_id = "numeric",
    individual_id = "character"
  )
)

# Create a constructor function for Individual
Individual$new <- function(household_id, individual_id) {
  return(new("Individual", household_id = household_id, individual_id = individual_id))
}

# Example usage
household_id <- 1
individual_id <- 1

# Generate individual ID
individual_id <- generate_individual_id(household_id, individual_id)

# Create individual object
individual <- Individual$new(household_id, individual_id)

# Print individual ID
print(paste("Individual ID:", individual@individual_id))

