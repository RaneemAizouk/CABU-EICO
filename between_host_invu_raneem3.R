# Set the working directory to Desktop/Rfiles
setwd("/Users/raizouk/Desktop/R files")

# Load packages, checks if the "pacman" package is installed. If not, it installs the package and loads it.
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}
p_load(viridis, adaptivetau)

# Clear workspace, remove all objects from the current environment.
rm(list=ls())

# Source model functions, sources (loads) functions from the file named "between_host_originalG_fn_raneem.R"
#source("./between_host_originalG_fn_raneem.R")

source("./Gillespe_SIS .R")

# Initialize an empty list to store results
results_list <- list()

# Set simulation parameters
lambda1 <- 0.2
lambda2 <- 0.3
rho <- 0.2
n.clusters <- 22

 houses_per_cluster <- 12

# Set initial conditions

S0 <- 3
I0 <-2 
tmax <- 365
n.sims <- 10
cols <- viridis(3)  # Assuming you want three colors
#here
# Run simulations
#for (i in 1:n.sims) {
  # Call the gillespe function
  # N <- sample(1:5, 1)
  # result <- gillespe(lambda1, lambda2, rho, S0, I0, tmax, N = N)
  
  # Store the result in the list
  # results_list[[i]] <- result
  #} here
# Repeat simulations for 100 times
for (sim in 1:n.sims) {
  # Run simulations for each cluster
  for (cluster in 1:n.clusters) {
    # Run simulations for each house in the cluster
    for (house in 1:houses_per_cluster) {
      # Call the gillespe function
      N <- sample(1:5, 1)  # Random number of individuals in the house
 
      S <- sample(0:N, 1)  # Random number of susceptible individuals in the house  
      I <- N - S  # Ensure that S0 + I0 equals N
      result <- gillespe(lambda1, lambda2, rho, S0, I0, tmax,N)
      sampled_result <- sample_household(result, sensitivity = 0.9)  # Adjust sensitivity as needed
      
      # Store the result in the list
      results_list[[length(results_list) + 1]] <- result
    }
  }
}
print(results_list)
# Plot output
alphaval <- 0.3
out.mean <- apply(do.call(rbind, results_list), c(1, 2), mean)
# Plot legend
print(do.call(rbind, results_list))
par(mfrow=c(1,1))
plot(NULL,
     las=1,xaxs="i",yaxs="i",ylim=c(0,6),bty="n",xlim=c(0,tmax),
     xlab="Day",ylab="Number")
legend("topleft", bty="n", col=cols[1:2], legend=c("Susceptible", "Infected"), lwd=3)


# Plot individual simulations
for (ii in 1:n.sims) {
  lines(results_list[[ii]]$t, results_list[[ii]]$S, col=adjustcolor(cols[1],alphaval))
  lines(results_list[[ii]]$t, results_list[[ii]]$I, col=adjustcolor(cols[2], alphaval))
}

# Plot average

#lines(out.mean["t"], out.mean["S"], col=adjustcolor(cols[1]), lwd=3)
#lines(out.mean["t"], out.mean["I"], col=adjustcolor(cols[2]), lwd=3)

