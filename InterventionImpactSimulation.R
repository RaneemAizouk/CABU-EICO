#The code uses simulation techniques and statistical
#modelling to evaluate the effects of public health 
#interventions on disease dynamics in a population.
#It simulates transmission events before and after interventions,
#providing insights into potential shifts in disease spread patterns. 
#The Gillespie algorithm is used in the pre-intervention phase to model random
#and variable probabilities of events where it uses parameters 
#like transmission rates (lambda1 and lambda2) and recovery rate (rho) 
#to mirror the baseline dynamics of disease spread within the population. 
#In the post-intervention phase, a modified version of the algorithm is used to 
#simulate the impact of intervention strategies. This adaptation adjusts the parameters to 
#reflect anticipated changes due to interventions—namely, a 50% reduction in transmission rates (lambda3 and lambda4)
#and a 20% increase in the recovery rate (rho1).


# Load required packages using pacman for simplicity
if (!require(pacman)) {
  install.packages("pacman")
}
# Ensure msmtools is installed
if (!requireNamespace("msmtools", quietly = TRUE)) {
  install.packages("msmtools")
}

library(msmtools)

library(pacman)
p_load(viridis, adaptivetau, msm, tidyverse, ggplot2)

# Clear workspace, remove all objects from the current environment
rm(list = ls())
setwd("/Users/raizouk/Desktop/R files")
# Source model functions, sources (loads) functions from the file named "Modified G Function.R"
source("./Modified G Function.R") # Control group function
source("./Modified G Function 2.R") # Intervention group function




#Pre-intervention,  the function "convert_to_state" takes the data frame "states_df" where each row in it represents a binary state composed of three binary digits (0 or 1)
#and convert  each row into a numeric state.
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

# Post-intervention, the function "convert_to_state_m" takes the data frame "states_df_m" where each row in it represents a binary state composed of three binary digits (0 or 1)
# and convert  each row into a numeric state.
convert_to_state_m <- function(states_df_m) {
  state_numeric_m <- apply(states_df_m, 1, function(states_m) {
    state_str_m <- paste(states_m, collapse = "")
    return(ifelse(state_str_m == "111", 15,
                  ifelse(state_str_m == "110", 14,
                         ifelse(state_str_m == "101", 13,
                                ifelse(state_str_m == "100", 12,
                                       ifelse(state_str_m == "011", 11,
                                              ifelse(state_str_m == "010", 10,
                                                     ifelse(state_str_m == "001", 9, 8))))))))
  })
  return(state_numeric_m)
}
# Define simulation parameters
tvec <- seq(0, 20, 1)
tmax <- max(tvec)
lambda1 <- 0.02 #before interventions.
lambda2 <- 0.03 #before interventions.
lambda3 <- 0.01 #after interventions (50% reduction). 
lambda4 <- 0.015 # after interventions (50% reduction).
rho <- 0.05 # before interventions.
rho1 <- 0.06 #after interventions (20% increase).
current_states <- c(0, 1, 1)  
# Initialize an empty list to store results
results_list <- list()
results_list_m <- list()
n.sims <- 1
n.clusters <- 12
houses_per_cluster <- 22
individuals_per_house <- 3
intervention_threshold <- 0.6
# Initialize named lists for better indexing
results_list <- list()
results_list_m <- list()

# Loop for simulations
for (sim in 1:n.sims) {
  for (cluster in 1:n.clusters) {
    for (house in 1:houses_per_cluster) {
      # ASSIGN HOUSE ID  
      HouseID <- paste((cluster - 1) * houses_per_cluster + house, "_", sim)
      
      # Randomly decide on the intervention status based on threshold
      intervention_status <- ifelse(runif(1) <= intervention_threshold, 1, 0)
      
      # Choose simulation function based on intervention status
      if (intervention_status == 1) {
        result_df <- gillespe1(lambda3, lambda4, rho1, current_states, tmax)
        result_df <- result_df %>% mutate(state_numeric = convert_to_state_m(select(., starts_with("X"))))
          #mutate(state_numeric = convert_to_state_m(.))
      } else {
        result_df <- gillespe(lambda1, lambda2, rho, current_states, tmax)
        result_df <- result_df %>% mutate(state_numeric = convert_to_state(select(., starts_with("X"))))
          #mutate(state_numeric = convert_to_state(.))
      }
      
      
      # Common processing for both scenarios
      result_df$HouseID <- HouseID
      result_df$simulation <- sim
      result_df$intervention_status <- intervention_status
     # result_df <- result_df %>%
      #   mutate(state_numeric = convert_to_state(select(., starts_with("X"))),
      #          X1_individual_id = paste(HouseID, "1", sep="-"),
      #          X2_individual_id = paste(HouseID, "2", sep="-"),
      #          X3_individual_id = paste(HouseID, "3", sep="-"))
      # 
      # # Further processing, if needed
      # result_df <- result_df %>%
      #   mutate(
      #     X1_individual_id = paste(HouseID, "1", sep="-"),
      #     X2_individual_id = paste(HouseID, "2", sep="-"),
      #     X3_individual_id = paste(HouseID, "3", sep="-")
      #   )
      results_list[[length(results_list) + 1]] <- result_df
      # Use HouseID as the key for direct access
     # results_list[HouseID] <- result_df
    }
  }
}

# # Loop for simulations
# for (sim in 1:n.sims) {
#   for (cluster in 1:n.clusters) {
#     for (house in 1:houses_per_cluster) {
#       # ASSIGN HOUSE ID  
#        HouseID <- paste((cluster - 1) * houses_per_cluster + house, "_", sim)
#      
# # This condition is to define those who are committed to the interventions after imposing it, so 60% of them will follow the interventions. 
#       if (runif(1) <= intervention_threshold) {
#         # Post-intervention scenario, follow the gillespe1 algorithm to simulate the data and estimate the assigned parameters(lambda3, lambda4, rho1) post interventions. 
#         result_df_m <- gillespe1(lambda3, lambda4, rho1, current_states, tmax)
#         result_df_m$HouseID <- HouseID
#         result_df_m$simulation <- sim
#         result_df_m <- result_df_m %>%
#           mutate(state_numeric_m = convert_to_state_m(select(., X11, X22, X33)),
#                  X11_individual_id = paste(HouseID, "1", sep="-"),
#                  X22_individual_id = paste(HouseID, "2", sep="-"),
#                  X33_individual_id = paste(HouseID, "3", sep="-"),
#                  intervention_status = 1) # Indicating post-intervention
#         results_list_m[[length(results_list_m) + 1]] <- result_df_m
#       } else {
#         # Pre-intervention scenario,follow the gillespe algorithm to simulate the data and estimate the assigned parameters(lambda1, lambda2, rho) post interventions. 
#         result_df <- gillespe(lambda1, lambda2, rho, current_states, tmax)
#         result_df$HouseID <- HouseID
#         result_df$simulation <- sim
#         result_df <- result_df %>%
#           mutate(state_numeric = convert_to_state(select(., X1, X2, X3)),
#                  X1_individual_id = paste(HouseID, "1", sep="-"),
#                  X2_individual_id = paste(HouseID, "2", sep="-"),
#                  X3_individual_id = paste(HouseID, "3", sep="-"), 
#                  intervention_status = 0)# Indicating pre-intervention
#         results_list[[length(results_list) + 1]] <- result_df
#       }
#     }
#   }
# }



# Count data frames with only one row in post intervention part.
#one_row_counts <- sapply(results_list_m, function(df) nrow(df) == 1)
# Sum up true values to get the total count
#total_one_row_df = sum(one_row_counts)




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





# Pre-intervention,  drewing samples by choosing the nearst time point to the required time points c(3,6,9,12)

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

#Pre-intervention, Create Time_list
Time_list <- list()

for (i in seq_along(nearest_times_list)) {
  time_points <- unique(na.omit(unlist(nearest_times_list[[i]])))
  #zero_time_count <- sum(common_data$time == 0)
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
    intervention_status <- results_list[[i]]$intervention_status[results_list[[i]]$time %in% time_points]

    # Create a data frame with time, X1, X2, X3, STATES, state_numeric, and HouseID
    # Creating common_data for the current set of results
    common_data <- data.frame(time = time_points, X1 = X1_values, X2 = X2_values, X3 = X3_values,
                                   STATES = STATES_values, state_numeric = state_numeric_values,
                                  HouseID = HouseID_values, simulation = simulation_values, intervention_status =  intervention_status
 
     #I1_individual_id = X1_individual_ids,  # Correct column name
     #I2_individual_id = X2_individual_ids,  # Correct column name
     #I3_individual_id = X3_individual_ids   # Correct column name
    )

    # Filter out rows where time is zero
    common_data <- common_data[common_data$time != 0, ]

    Time_list[[i]] <- common_data
  }
}
# Time_list <- list()
# 
# for (i in seq_along(nearest_times_list)) {
#   # Extract unique and non-missing time points for the current dataset
#   time_points <- unique(na.omit(unlist(nearest_times_list[[i]])))
#   
#   # Only proceed if there are valid time points to process
#   if (length(time_points) > 0) {
#     # Extract data for the current time points
#     X1_values <- results_list[[i]]$X1[results_list[[i]]$time %in% time_points]
#     X2_values <- results_list[[i]]$X2[results_list[[i]]$time %in% time_points]
#     X3_values <- results_list[[i]]$X3[results_list[[i]]$time %in% time_points]
#     STATES_values <- results_list[[i]]$state[results_list[[i]]$time %in% time_points]
#     state_numeric_values <- results_list[[i]]$state_numeric[results_list[[i]]$time %in% time_points]
#     HouseID_values <- results_list[[i]]$HouseID[results_list[[i]]$time %in% time_points]
#     simulation_values <- results_list[[i]]$simulation[results_list[[i]]$time %in% time_points]
#     
#     # Creating common_data for the current set of results
#     common_data <- data.frame(time = time_points, X1 = X1_values, X2 = X2_values, X3 = X3_values,
#                               STATES = STATES_values, state_numeric = state_numeric_values,
#                               HouseID = HouseID_values, simulation = simulation_values)
#     
#     # Filter out rows where time is zero, if needed
#     common_data <- common_data[common_data$time != 0, ]
#     
#     # Save the processed data for the current iteration
#     Time_list[[i]] <- common_data
#   }
# }

# Combine data from Time_list into a single data frame
all_data <- do.call(rbind, Time_list)

# Print the unique time points for each simulation
print(all_data)

# 
# desired_time_points_m <- c(3, 6, 9, 12)
# 
# nearest_times_list_m <- lapply(results_list_m, function(item1) {
#   sapply(desired_time_points_m, function(time) {
#     closest_points_m <- item1$t[item1$t <= time]
#     if (length(closest_points_m) > 0) {
#       return(tail(closest_points_m, 1))
#     } else {
#       return(NA)
#     }
#   })
# })
# 
# #Post-intervention, Create Time_list_m
# Time_list_m <- list()
# 
# for (i in seq_along(nearest_times_list_m)) {
#   # Extract unique and non-missing time points for the current dataset
#   time_points_m <- unique(na.omit(unlist(nearest_times_list_m[[i]])))
#   
#   # Only proceed if there are valid time points to process
#   if (length(time_points_m) > 0) {
#     # Extract data for the current time points
#     X11_values <- results_list_m[[i]]$X11[results_list_m[[i]]$time %in% time_points_m]
#     X22_values <- results_list_m[[i]]$X22[results_list_m[[i]]$time %in% time_points_m]
#     X33_values <- results_list_m[[i]]$X33[results_list_m[[i]]$time %in% time_points_m]
#     STATES_values_m <- results_list_m[[i]]$state_m[results_list_m[[i]]$time %in% time_points_m]
#     state_numeric_values_m <- results_list_m[[i]]$state_numeric_m[results_list_m[[i]]$time %in% time_points_m]
#     HouseID_values_m <- results_list_m[[i]]$HouseID[results_list_m[[i]]$time %in% time_points_m]
#     simulation_values_m <- results_list_m[[i]]$simulation[results_list_m[[i]]$time %in% time_points_m]
#     
#     # Creating common_data_m for the current set of results
#     common_data_m <- data.frame(time = time_points_m, X11 = X11_values, X22 = X22_values, X33 = X33_values,
#                                 STATES = STATES_values_m, state_numeric_m = state_numeric_values_m,
#                                 HouseID = HouseID_values_m, simulation = simulation_values_m)
#     
#     # Save the processed data for the current iteration
#     Time_list_m[[i]] <- common_data_m
#   }
# }
# 
# # Combine data from Time_list_m into a single data frame
# all_data_m <- do.call(rbind, Time_list_m)
# 
# # Print the unique time points for each simulation
# print(all_data_m)




# #Post-intervention, drawing samples by choosing the nearest time point to the required time points c(3,6,9,12)
# 
# # Assuming results_list_m is your list of results
# desired_time_points_m <- c(3, 6, 9, 12)
# 
# nearest_times_list_m <- lapply(results_list_m, function(item1) {
#   sapply(desired_time_points_m, function(time) {
#     closest_points_m <- item1$t[item1$t <= time]
#     if (length(closest_points_m) > 0) {
#       return(tail(closest_points_m, 1))
#     } else {
#       return(NA)
#     }
#   })
# })
# 
# #Post-intervention, Create Time_list
# Time_list_m <- list()
# 
# for (i in seq_along(nearest_times_list_m)) {
#   time_points_m <- unique(na.omit(unlist(nearest_times_list_m[[i]])))
#   #zero_time_count_m <- sum(common_data_m$time == 0)
#   if (length(time_points_m) > 0) {
#     X11_values <- results_list_m[[i]]$X11[results_list_m[[i]]$time %in% time_points_m]
# 
#     X22_values <- results_list_m[[i]]$X22[results_list_m[[i]]$time %in% time_points_m]
# 
#     X33_values <- results_list_m[[i]]$X33[results_list_m[[i]]$time %in% time_points_m]
# 
#     STATES_values_m <- results_list_m[[i]]$state_m[results_list_m[[i]]$time %in% time_points_m]
#     state_numeric_values_m <- results_list_m[[i]]$state_numeric_m[results_list_m[[i]]$time %in% time_points_m]
#     HouseID_values_m <- results_list_m[[i]]$HouseID[results_list_m[[i]]$time %in% time_points_m]
#     #Individual_IDs <- paste(HouseID_values, "_Individual", seq_along(X1_values), sep="")
#     #X1_individual_ids <- results_list[[i]]$X1_individual_id[results_list[[i]]$time %in% time_points]
#     #X2_individual_ids <- results_list[[i]]$X2_individual_id[results_list[[i]]$time %in% time_points]
#     #X3_individual_ids <- results_list[[i]]$X3_individual_id[results_list[[i]]$time %in% time_points]
#     simulation_values_m <-results_list_m[[i]]$simulation[results_list_m[[i]]$time %in% time_points_m]
#     intervention_status <- results_list_m[[i]]$intervention_status[results_list_m[[i]]$time %in% time_points_m]
#     # Create a data frame with time, X1, X2, X3, STATES, state_numeric, and HouseID
#     common_data_m <- data.frame(
#       time = time_points_m,
#       X11 = X11_values,
#       X22 = X22_values,
#       X33 = X33_values,
#       STATES = STATES_values_m,
#       state_numeric_m = state_numeric_values_m,
#       HouseID = HouseID_values_m,
#       simulation = simulation_values_m,
#       intervention_status =   intervention_status
#       #I1_individual_id = X1_individual_ids,  # Correct column name
#       #I2_individual_id = X2_individual_ids,  # Correct column name
#       #I3_individual_id = X3_individual_ids   # Correct column name
#     )
#     # Filter out rows where time is zero
#     common_data_m <- common_data_m[common_data_m$time != 0, ]
# 
#     Time_list_m[[i]] <- common_data_m
#   }
# }
# all_data_m <- do.call(rbind, Time_list_m)
# 
# # Print the unique time points for each simulation
# print(all_data_m)

library(reshape2)
library(tidyr)

# #1) (all_data_m) Converting data to long format data to fit msm package 
# # Post intervention, combine data from "Time_list_m" into a single data frame
# 
# 
# long_data_m <-melt(all_data_m, id.vars=c("HouseID", "time","intervention_status"),
#                    measure.vars = c("X11","X22","X33","state_numeric_m"))
# long_data_m <- all_data_m[,c("HouseID", "time","state_numeric_m", "intervention_status")]
# long_data_m$state_numeric_m <- as.numeric(long_data_m$state_numeric_m)
# long_data_m$state_numeric_m <- long_data_m$state_numeric_m + 1

#2) (all_data) Converting data to long format data  to fit msm package 

 #Melting the data frame to long format
long_data <- melt(all_data, id.vars = c("HouseID", "time","intervention_status"), 
                  measure.vars = c("X1", "X2", "X3", "state_numeric"))

# Pre intervention, assuming 'all_data' contains 'HouseID', 'time', and 'state_numeric' among other columns
# Prepare long_data by selecting only the relevant columns for MSM analysis
long_data <- all_data[, c("HouseID", "time","state_numeric", "intervention_status")]

# Ensure 'state_numeric' is numeric, as required by the msm package
long_data$state_numeric <- as.numeric(long_data$state_numeric)


# Incrementing state_numeric by 1 if  state numbering starts from 0 so it will start from 1.
long_data$state_numeric <- long_data$state_numeric + 1

# # Post intervention, remove all observations from "long_data_m" that only have one observation at a certain time point for each HouseID.
# sufficient_data_m <- long_data_m %>%
#   arrange(HouseID, time) %>%
#   group_by(HouseID) %>%
#   filter(n() > 1) %>%
#   #filter(simulation == 1) %>%
#   ungroup()


# Pre intervention, remove all observations from "long_data" that only have one observation at a certain time point for each HouseID.

sufficient_data <- long_data %>%
 arrange(HouseID, time) %>%
  group_by(HouseID) %>%
   filter(n() > 1) %>%
#filter(simulation == 1) %>%
    ungroup()

# Create a table showing the observed numbers of transitions between each pair of states in a multi-state model. 
# # Post intervention
# result_m <- statetable.msm(state = sufficient_data_m$state_numeric_m, subject = sufficient_data_m$HouseID, data = sufficient_data_m)
# print(result_m)
# Pre intervention
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




#qconstraint = c(1,1,1,3,2,2,3,2,2,3,3,4,3,2,2,3,3,4,3,3,4,3,3,3)
# Ensure rows sum to 0
diag(Q) <- -rowSums(Q)
print(Q)

r1_m <- lambda3
r2_m <- lambda3+lambda4
r3_m<- rho1
r4_m <- lambda3+2*lambda4


Q_m <- matrix(c(
  0,  r1_m, r1_m, 0,  r1_m, 0,  0,  0,
  r3_m, 0,  0,  r2_m , 0,  r2_m , 0,  0,
  r3_m, 0,  0,  r2_m , 0,  0, r2_m ,  0,
  0,  r3_m, r3_m, 0,  0,  0,  0,  r4_m,
  r3_m, 0,  0,  0,  0, r2_m , r2_m , 0,
  0,  r3_m, 0,  0,  r3_m, 0,  0,  r4_m,
  0,  0,  r3_m, 0,  r3_m, 0,  0,  r4_m,
  0,  0,  0,  r3_m, 0, r3_m, r3_m, 0
), nrow=8, byrow =TRUE)
diag(Q_m) <- -rowSums(Q_m)
print(Q_m)

# Assuming Q and Q_m are already defined as shown in your snippet

# Combine Q and Q_m into a larger matrix, Q_combined
# This step assumes you're not considering transitions between pre- and post- intervention states
Q_combined <- rbind(cbind(Q, matrix(0, nrow=8, ncol=8)),
                    cbind(matrix(0, nrow=8, ncol=8), Q_m))
print(Q_combined)
# Fit the msm model for Post intervention data. 
result <- msm(state_numeric ~ time, subject = HouseID, data = sufficient_data, qmatrix =  Q_combined,covariates = ~ intervention_status , qconstraint = c(1,1,1,3,2,2,3,2,2,3,3,4,3,2,2,3,3,4,3,3,4,3,3,3,5,5,5,7,6,6,7,6,6,7,7,8,7,6,6,7,7,8,7,7,8,7,7,7))

summary(result)



# # Fit the msm model for Pre intervention data. 
# result_m <- msm(state_numeric_m ~ time, subject = HouseID, data = sufficient_data, qmatrix = Q_combined , qconstraint = c(1,1,1,3,2,2,3,2,2,3,3,4,3,2,2,3,3,4,3,3,4,3,3,3))
# summary(result_m)
# # To plot observed vs. expected state prevalences
# prevplot(result)

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

