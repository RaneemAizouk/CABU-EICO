# Load required packages using pacman for simplicity
if (!require(pacman)) {
  install.packages("pacman")
}
library(pacman)
p_load(viridis, adaptivetau, msm, tidyverse, ggplot2)

# Clear workspace, remove all objects from the current environment
rm(list = ls())

# Source model functions, sources (loads) functions from the file named "Modified G Function.R"
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
                                                     ifelse(state_str == "001", 1, 0))))))))
  })
  return(state_numeric) 
}

# Define simulation parameters
tvec <- seq(0, 14, 1)
tmax <- max(tvec)
lambda1 <- 0.02
lambda2 <- 0.03
rho <- 0.02
#current_states <- c(0, 1, 1)  
# Initialize an empty list to store results
results_list <- list()
n.sims <- 1
n.clusters <- 12
houses_per_cluster <- 22
individuals_per_house <- 3


# Adjust the loop to generate result_df with the desired structure
for (sim in 1:n.sims) {
  for (cluster in 1:n.clusters) {
    for (house in 1:houses_per_cluster) {
      
      HouseID <- paste((cluster - 1) * houses_per_cluster + house,"_",sim)
      #HouseID <- paste((cluster - 1) * houses_per_cluster + house)
      
      # Assuming gillespe function or similar returns states for the house
      # This should return a data frame with columns for time and states (X1, X2, X3)
      result_df <- gillespe(lambda1, lambda2, rho, current_states, tmax)
      
      # Initialize columns for individual IDs in result_df
      result_df$X1_individual_id <- NA
      result_df$X2_individual_id <- NA
      result_df$X3_individual_id <- NA
      
      # Assign individual IDs for each state
      for (individual in 1:individuals_per_house) {
        individual_id <- paste(HouseID, individual, sep="-")
        
        # Assuming each individual corresponds to a state (X1, X2, X3)
        # Adjust the assignment of individual IDs accordingly
        if (individual == 1) {
          result_df$X1_individual_id <- individual_id
        } else if (individual == 2) {
          result_df$X2_individual_id <- individual_id
        } else if (individual == 3) {
          result_df$X3_individual_id <- individual_id
        }
      }
      # Add HouseID to result_df before storing it
      result_df$HouseID <- HouseID  # This line adds HouseID to your dataframe
      
      # Convert states to decimal
      
      result_df$state_numeric <- convert_to_state(result_df[, c("X1", "X2", "X3")])
      result_df$simulation <- sim
      
      # Store the modified result_df in the results_list
      results_list[[length(results_list) + 1]] <- result_df
    }
  }
}


# Assuming you later process results_list to create all_data or similar
# This part would be greatly dependent on the structure of 'result' from gillespe

# Example of processing results_list (Adjust according to your actual data structure)
# This is a placeholder and needs actual implementation details

# Now all_data would contain Individual_IDs if the above loop correctly assigns them

# Further processing and plotting can now utilize Individual_ID


#     library(ggplot2)
#     # Plotting the results with a gradient color scale
#     ggplot(result, aes(x = time, y = X1 + X2 + X3, color = state_numeric)) +
#       geom_step() +
#       labs(title = "Simulation Results",
#            x = "Time",
#            y = "Population",
#            color = "State") +
#       scale_color_viridis_c() +  # Use viridis color scale
#       theme_minimal()
#     
#     # Print the plot
#     print(ggplot(result, aes(x = time, y = X1 + X2 + X3, color = state_numeric)) +
#             geom_step() +
#             labs(title = "Simulation Results",
#                  x = "Time",
#                  y = "Population",
#                  color = "State") +
#             scale_color_viridis_c() +  # Use viridis color scale
#             theme_minimal())
#   }
# }
# Print the result
#print(result)



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
    #Individual_IDs <- paste(HouseID_values, "_Individual", seq_along(X1_values), sep="")
    #X1_individual_ids <- results_list[[i]]$X1_individual_id[results_list[[i]]$time %in% time_points]
    #X2_individual_ids <- results_list[[i]]$X2_individual_id[results_list[[i]]$time %in% time_points]
    #X3_individual_ids <- results_list[[i]]$X3_individual_id[results_list[[i]]$time %in% time_points]
    simulation_values <-results_list[[i]]$simulation[results_list[[i]]$time %in% time_points]
    
    # Create a data frame with time, X1, X2, X3, STATES, state_numeric, and HouseID
    common_data <- data.frame(
      time = time_points,
      X1 = X1_values,
      X2 = X2_values,
      X3 = X3_values,
      STATES = STATES_values,
      state_numeric = state_numeric_values,
      HouseID = HouseID_values,
      simulation = simulation_values
      #I1_individual_id = X1_individual_ids,  # Correct column name
      #I2_individual_id = X2_individual_ids,  # Correct column name
      #I3_individual_id = X3_individual_ids   # Correct column name
    )
    
    
    # Filter out rows where time is zero
    common_data <- common_data[common_data$time != 0, ]
    
    Time_list[[i]] <- common_data
  }
}

# Combine data from Time_list into a single data frame
all_data <- do.call(rbind, Time_list)

# Print the unique time points for each simulation
print(all_data)


# Assuming 'all_data' is your data frame
library(reshape2)
library(tidyr)

# Melting the data frame to long format
long_data <- melt(all_data, id.vars = c("HouseID", "time","simulation"), 
                  measure.vars = c("X1", "X2", "X3", "state_numeric"))


library(tidyr)
# Assuming 'all_data' contains 'HouseID', 'time', and 'state_numeric' among other columns
# Prepare long_data by selecting only the relevant columns for MSM analysis
long_data <- all_data[, c("HouseID", "time","state_numeric","simulation")]

# Ensure 'state_numeric' is numeric, as required by the msm package
long_data$state_numeric <- as.numeric(long_data$state_numeric)

# Optionally, adjust 'state_numeric' values if necessary for your analysis
# For example, incrementing state_numeric by 1 if your state numbering starts from 1
long_data$state_numeric <- long_data$state_numeric + 1
# # Melt the data frame
# long_data <- melt(all_data, id.vars = c("HouseID", "time"), measure.vars = c("X1", "X2", "X3", "state_numeric"))
# 
# # Separate data for X1, X2, X3, STATES, state_numeric
# long_data <- spread(long_data, key = "variable", value = "value")
# 
# # Rename columns for clarity
# colnames(long_data) <- c("HouseID", "time", "X1", "X2", "X3", "state_numeric")
# # Add 1 to state_numeric
# long_data$state_numeric <- long_data$state_numeric + 1
# 
# 
# # Arrange the data by HouseID and time
# long_data <- long_data[order(long_data$HouseID, long_data$t), ]
# Add IndividualID
#long_data$IndividualID <- with(long_data, ave(1:nrow(long_data), HouseID, FUN = seq_along))

sufficient_data <- long_data %>%
  arrange(HouseID,time) %>%
  group_by(HouseID) %>%
  filter(n() > 1) %>%
  #filter(simulation == 3) %>%
  ungroup()

result <- statetable.msm(state = sufficient_data$state_numeric, subject = sufficient_data$HouseID, data = sufficient_data)
print(result)
# Define the Q matrix
#Q <- matrix(0, nrow = 8, ncol = 8) # Initialize Q-matrix with zeros
r1 <- lambda1
r2 <- lambda1+lambda2
r3<- rho
r4 <- lambda1+2*lambda2

Q <- matrix(c(
  0,  r1, r1, 0,  r1, 0,  0,  0,
  r3, 0,  0,  r2, 0,  r2, 0,  0,
  r3, 0,  0,  r2, 0,  0, r2,  0,
  0,  r3, r3, 0,  0,  0,  0,  r4,
  r3, 0,  0,  0,  0,  r2, r2, 0,
  0,  r3, 0,  0,  r3, 0,  0,  r4,
  0,  0,  r3, 0,  r3, 0,  0,  r4,
  0,  0,  0,  r3, 0,  r3, r3, 0
), nrow=8, byrow =TRUE)

#qconstraint = c(1,1,1,3,2,2,3,2,2,3,3,4,3,2,2,3,2,4,3,3,4,3,3,3)
# Ensure rows sum to 0
diag(Q) <- -rowSums(Q)
print(Q)
# Fit the msm model
result <- msm(state_numeric ~ time, subject = HouseID, data = sufficient_data, qmatrix = Q, qconstraint = c(1,1,1,3,2,2,3,2,2,3,2,4,3,2,2,3,3,4,3,3,4,3,3,3) )
summary(result)


# # Fill in the rates based on the constraints
# # Assuming direct transitions are allowed as per the defined rules
# Q[1,c(2,3,5)] <- r1 # From 000 to states with one infected
# Q[2, c(6,4)] <- r2  # From states with one infected to two
# Q[3, c(7,4)] <- r2 # From states with one infected to two
# Q[5, c(6,7)] <- r2
# Q[c(2,3,5), 1] <- r3
# Q[7, c( 5,3)] <- r3
# Q[6, c( 5,2)] <- r3
# Q[4, c( 3,2)] <- r3
# Q[c(4,6,7), 8] <- r4 # From states with two infected to three
# Q[8, c(4,6,7)] <- r3# From 111 to states with two infected 
# # Ensur e rows sum to 0
# diag(Q) <- -rowSums(Q)
# print(Q)
# 
# 
# summary(result)

