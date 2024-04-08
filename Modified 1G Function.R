# Load packages, checks if the "pacman" package is installed. If not, it installs the package and loads it.
setwd ("/Users/raizouk/Desktop/R files")
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(viridis, adaptivetau)
if (!require(msm)) {
  install.packages("msm")
  library(msm)
}
p_load(viridis, adaptivetau)
if (!require(msm)) {
  install.packages("tidyverse")
  library(msm)
}

# Clear workspace, remove all objects from the current environment.
rm(list = ls())

# Source model functions, sources (loads) functions from the file named "between_host_invu_fn_raneem.R"
source("./Modified G Function.R")

convert_to_state <- function(states_df) {
  state_numeric <- apply(states_df, 1, function(states) {
    state_str <- paste(states, collapse = "")
    return(ifelse(state_str == "111", 7,
                  ifelse(state_str == "110", 6,
                         ifelse(state_str == "101", 5,
                                ifelse(state_str == "100", 4,
                                       ifelse(state_str == "011", 3,
                                              ifelse(state_str == "010", 2,
                                                     ifelse(state_str == "001",1, 0))))))))
  })
  return(state_numeric)
}


# Initialize an empty list to store results
results_list <- list()

# Simulation Params
# a sequence of time points from 0 to 14 with a step of 1.
tvec <- seq(0, 14, 1)
tmax <- max(tvec)
lambda1 <- 0.2
lambda2 <- 0.3
rho <- 0.2

# Color palette for plotting
cols <- viridis(3)

# Initialize an empty list to store results
results_list <- list()
n.sims <- 1
n.clusters <- 22
houses_per_cluster <- 12

for (sim in 1:n.sims) {
  # Run simulations for each cluster
  for (cluster in 1:n.clusters) {
    # Run simulations for each house in the cluster
    for (house in 1:houses_per_cluster) {
      # Call the gillespe function
      # N <- sample(1:5, 1)  # Random number of individuals in the house with a maximum of 5 individuals per house.
      N <- 3
      initial_states <- c(0, 0, 1)  # Example initial states
      
      result <- gillespe(lambda1, lambda2, rho, initial_states, tmax)
      
      # Calculate HouseID
      HouseID <- rep((cluster - 1) * houses_per_cluster + house, length(result$t))
      
      # Convert S and I to states
      # result$state_numeric <- convert_to_state(result$S, result$I)
      
      # Combine results and identifiers
      result <- cbind(result, HouseID)
      
      # Convert states to decimal
      result$state_numeric <- convert_to_state(result[, c("X1", "X2", "X3")])
      
      
      
      # Store the result in the list
      results_list[[length(results_list) + 1]] <- result
    }}}
library(dplyr)
library(ggplot2)
library(tidyr)



# # Plotting state frequencies for all simulations
# for (i in seq_along(results_list)) {
#   # Extract result
#   result <- results_list[[i]]
#   
#   # Calculate state frequencies
#   state_frequencies <- table(result$state_numeric)
#   
#   # Convert to data frame for ggplot
#   state_freq_df <- data.frame(state_numeric = as.numeric(names(state_frequencies)),
#                               frequency = as.numeric(state_frequencies))
#   
#   # Plotting the state frequencies
#   plot <- ggplot(state_freq_df, aes(x = factor(state_numeric), y = frequency, fill = factor(state_numeric))) +
#     geom_bar(stat = "identity") +
#     labs(title = "State Frequencies",
#          x = "State",
#          y = "Frequency") +
#     scale_fill_viridis_d() +  # Use viridis color scale
#     theme_minimal()
#   
#   # Save or print the plot as desired
#   print(plot)
# }



# 3) Drewing samples by choosing the nearst time point to the required time points c(3,6,9,12)

# Assuming results_list is your list of results
desired_time_points <- c(3, 6, 9, 12)

nearest_times_list <- lapply(results_list, function(item) {
  sapply(desired_time_points, function(time) {
    closest_points <- item$t[item$t <= time]
    if (length(closest_points) > 0) {
      return(tail(closest_points, 1))
    } else {
      return(NA)
    }
  })
})

# Create Time_list
Time_list <- list()

for (i in seq_along(nearest_times_list)) {
  time_points <- unique(na.omit(unlist(nearest_times_list[[i]])))
  if (length(time_points) > 0) {
    X1_values <- results_list[[i]]$X1[results_list[[i]]$time %in% time_points]
    X2_values <- results_list[[i]]$X2[results_list[[i]]$time %in% time_points]
    X3_values <- results_list[[i]]$X3[results_list[[i]]$time %in% time_points]
    STATES_values <- results_list[[i]]$state[results_list[[i]]$time %in% time_points]
    state_numeric_values <- results_list[[i]]$state_numeric[results_list[[i]]$time %in% time_points]
    HouseID_values <- results_list[[i]]$HouseID[results_list[[i]]$time %in% time_points]
    
    # Create a data frame with time, X1, X2, X3, STATES, state_numeric, and HouseID
    common_data <- data.frame(
      time = time_points,
      X1 = X1_values,
      X2 = X2_values,
      X3 = X3_values,
      STATES = STATES_values,
      state_numeric = state_numeric_values,
      HouseID = HouseID_values
    )
    # Filter out rows where time is zero
    common_data <- common_data[common_data$time != 0, ]
    
    Time_list[[i]] <- common_data
  }
}

# Combine data from Time_list into a single data frame
all_data <- do.call(rbind, Time_list)

# Sort the data frame by time and IndividualID
all_data <- all_data[order(all_data$time, all_data$  HouseID), ]

# Print the sorted data frame
print(all_data)

# Assuming 'all_data' is your data frame
library(reshape2)
library(tidyr)
# Assuming 'long_data' is your data frame

library(reshape2)

# Assuming you have already created the IndividualID and HouseID variables

# Melt the data frame
long_data <- melt(all_data, id.vars = c("HouseID",  "time"), measure.vars = c("X1", "X2", "X3", "state_numeric"))

# Separate data for X1, X2, X3, STATES, state_numeric
long_data <- spread(long_data, key = "variable", value = "value")

# Rename columns for clarity
colnames(long_data) <- c("HouseID",  "time", "X1", "X2", "X3", "state_numeric")
# Add 1 to state_numeric
long_data$state_numeric <- long_data$state_numeric + 1

long_data <- long_data[order(long_data$HouseID, long_data$time), ]

#Remove subjects that have only one observation, as they cannot contribute to the estimation of transition probabilities.
#Eg:HouseID = c(1, 1, 2, 3, 3, 3, 4),house hold 2 has one observation here so we remove it.
sufficient_data <- long_data %>%
  group_by(HouseID) %>%
  filter(n() > 1) %>%
  ungroup()

result <- statetable.msm(state = sufficient_data$state_numeric, subject =  sufficient_data$HouseID, data = sufficient_data)
print(result)
# Define the Q-matrix with constraints
Q <- matrix(0, nrow = 8, ncol = 8) # Initialize Q-matrix with zeros
r1 <- lambda1
r2 <- lambda1+lambda2
r3<- rho
r4 <- lambda1+2*lambda2
 
# Fill in the rates based on the constraints
# Assuming direct transitions are allowed as per the defined rules
Q[1, c(2,3,5)] <- r1 # From 000 to states with one infected
Q[2, c(6,4)] <- r2  # From states with one infected to two
Q[3, c(7,4)] <- r2 # From states with one infected to two
Q[5, c(6,7)] <- r2
Q[c(2,3,5), 1] <- r3
Q[7, c( 5,3)] <- r3
Q[6, c( 5,2)] <- r3
Q[4, c( 3,2)] <- r3
Q[c(4,6,7), 8] <- r4 # From states with two infected to three
Q[8, c(4,6,7)] <- r3# From 111 to states with two infected 
# Ensur e rows sum to 0
diag(Q) <- -rowSums(Q)
print(Q)

# Fit the msm model
result <- msm(state_numeric ~ time, subject = HouseID, data = sufficient_data, qmatrix = Q) 
# method = "Nelder-Mead",  control = list(maxit = 100000, reltol = 1e-12))
summary(result)

