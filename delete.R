# Set the working directory to Desktop/Rfiles
setwd("/Users/raizouk/Desktop/R files")

# Load packages, checks if the "pacman" package is installed. If not, it installs the package and loads it.
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(viridis, adaptivetau)
if (!require(msm)){
  install.packages("msm")
  library(msm)
}
p_load(viridis, adaptivetau)
if (!require(msm)){
  install.packages("tidyverse")
  library(msm)
}

# Clear workspace, remove all objects from the current environment.
rm(list = ls())

# Source model functions, sources (loads) functions from the file named "between_host_invu_fn_raneem.R"
source("./between_host_invu_fn_raneem3.R")

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

# Results list and result are the lists for the general simulations.
# Stochastic simulation
# we consider the simulations over the clusters and houses.
# Repeat simulations for (n.sims) times
# Initialize an empty list to store results
results_list <- list()
n.sims <- 1
n.clusters <- 22
houses_per_cluster <- 12
# Create a function to convert S and I values to states
#convert_to_state <- function(S, I) {
 # state_str <- paste(S, I, sep = "")
 # state_numeric <- ifelse(state_str == "111", 1,
  #                        ifelse(state_str == "000", 2,
  #                               ifelse(state_str == "110", 3,
  #                                      ifelse(state_str == "101", 4,
   #                                            ifelse(state_str == "011", 5,
   #                                                   ifelse(state_str == "001", 6,
    #                                                         ifelse(state_str == "010", 7, 8)))))))
  #return(state_numeric)
#}

for (sim in 1:n.sims) {
  # Run simulations for each cluster
  for (cluster in 1:n.clusters) {
    # Run simulations for each house in the cluster
    for (house in 1:houses_per_cluster) {
      # Call the gillespe function
      #N <- sample(1:5, 1)  # Random number of individuals in the house with a maximum of 5 individuals per house.
      N<-3
      S1 <- sample(0:N, 1)  # Random number of susceptible individuals in the house
      I1 <- N - S1  # Ensure that S0 + I0 equals N
      state.orig <- c(I = I1, S = S1)
      
      result <- gillespe(lambda1, lambda2, rho, S1, I1, tmax)
       
      # Calculate HouseID
      HouseID <- rep((cluster - 1) * houses_per_cluster + house, length(result$t))
      
       # Convert S and I to states
     # result$state_numeric <- convert_to_state(result$S, result$I)
    
      # Combine results and identifiers
      result <- cbind(result, HouseID)
       # Store the result in the list
      results_list[[length(results_list) + 1]] <- result
    
      }
      
    }
  }



print(result)
# 2) Find the mea
#n which include a) Rounding the time points using the floor function. This step acknowledges that the mean cannot be calculated for continuous time points, so rounding them down to the nearest whole number provides discrete time points for aggregation
# b) put them in a new list with the corresponding S and I values. c) find the mean 
# To find the mean, we took the floor function on the time values to round them down to the nearest whole number and store the corresponding values for S and I in a data frame.
# The resulting data frame (approx_list) provides a summarized view of the simulation results at discrete time points, making it easier to calculate and visualize the mean behavior over time.
# approx_list: is the list for the aproximated values to draw the mean.
approx_list <- list() # create an empty list()
for (item in results_list){
  approx_list[[length(approx_list)+1]] <- data.frame(t=floor(item$t), S=item$S, I=item$I)
}

# we will have repeated time values so we wanted to ensures to aggregate or summarise data at each unique time point. why we use group_by(t)
grouped_approx_list <- do.call(rbind, approx_list) %>%
  group_by(t) %>%
  add_count(name="group_count") %>%
  summarize(S = sum(S), I = sum(I), group_count=unique(group_count))
# to combine the data frames in the list approx_list into a single data frame and then perform further operations using the %>% (pipe) operator.
out.mean <- do.call(rbind,approx_list) %>%
  # to combine similar time values. 
  group_by(t) %>%
  #group_count: provides the count of observations (group_count) for each unique time point in the data.
  add_count(name="group_count") %>%
  summarize(mean_S = sum(S)/ n(), mean_I = sum(I)/n(), group_count=unique(group_count))

print(out.mean)


#SHALL WE APROOXIMATE THE MEANS TO UNIQIE TIMES POINTS OR LEAVE IT AS DECIMALS ?
# Plot output

names(out.mean) <- c("time", "S", "I")
alphaval <- 0.3

par(mfrow=c(1,1))
plot(NULL,
     las=1,xaxs="i",yaxs="i",ylim=c(0,6),bty="n",xlim=range(tvec),
     xlab="time",ylab="Frequency")
# Plot legend
legend("topleft", bty="n", col=cols[1:2], legend=c("Susceptible", "Infected"), lwd=3)


# Plot individual simulations
for (ii in 1:n.sims) {
  lines(results_list[[ii]]$t, results_list[[ii]]$S, col=adjustcolor(cols[1], alphaval))
  lines(results_list[[ii]]$t, results_list[[ii]]$I, col=adjustcolor(cols[2], alphaval))
}

# Corrected lines for plotting the mean
lines(out.mean[["time"]],out.mean[["S"]],col=adjustcolor(cols[1]),lwd=3)
lines(out.mean[["time"]],out.mean[["I"]],col=adjustcolor(cols[2]),lwd=3)

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
  #na.omit.:Removes NA values 
  if (length(time_points) > 0) {
    S_values <- results_list[[i]]$S[results_list[[i]]$t %in% time_points]
    I_values <- results_list[[i]]$I[results_list[[i]]$t %in% time_points]
    HouseID_values <- results_list[[i]]$HouseID[results_list[[i]]$t %in% time_points]
    # Create a data frame with time, S, I, and HouseID
    common_data <- data.frame(t = time_points, S = S_values, I = I_values,HouseID = HouseID_values)
    # Filter out rows where time is zero
    common_data <- common_data[common_data$t != 0, ]
    
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
# Assuming 'long_data' is your data frame
# Assuming 'all_data' is your data frame
library(reshape2)


# Melt the data frame
long_data <- melt(all_data, id.vars = c("HouseID", "t"), measure.vars = c( "S","I"))

# Separate data for S and I
long_data <- spread(long_data, key = "variable", value = "value")

# Rename columns for clarity
colnames(long_data) <- c("HouseID", "t", "count_S")

# Arrange the data by HouseID and time
long_data <- long_data[order(long_data$HouseID, long_data$t), ]

print(long_data)


# Define the Q matrix
Q <- matrix(c(-0.5, 0.5, 0.2, -0.2), nrow = 2, ncol = 2, byrow = TRUE)

# Q.init  <- init.msm(state ~ t, id, data=long_data, qmatrix=Q)
                                                
# Fit the msm model
result <- msm(
  #state  ~ c(count_S, count_I),
  state ~ t, 
  subject = HouseID,
  data = long_data,
  qmatrix = Q,
  #gen.inits = True
)

