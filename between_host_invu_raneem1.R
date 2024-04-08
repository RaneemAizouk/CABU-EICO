# Set the working directory to Desktop/Rfiles
setwd("/Users/raizouk/Desktop/R files")
# Load packages,checks if the "pacman" package is installed. If not, it installs the package and loads it.
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(viridis,
       adaptivetau)
# Clear workspace, remove all objects from the current environment.
rm(list=ls())
# Source model functions,sources (loads) functions from the file named"between_host_invu_fn_raneem.R"

# Initialize an empty list to store results

source("./between_host_invu_fn_raneem.R")
# Params
# a sequence of time points from 0 to 14 with a step of 1.
results_list <- list()
tvec <- seq(0,365,1)
# Parameters for the model
parms <- c(lambda1 = 0.2, lambda2 = 0.3, rho = 0.2)
# Initial state of the system
state.orig <- c(I=2,S=7)
# Color palette for plotting
cols <- viridis(3)
# Stochastic simulation
n.sims <- 100
N <- 9 # number of people
# as.vector: convert matrix into a vector.
# rmultinom :generates random multinomial variates, size=N represents the number of trials 
# prob parameter is a vector of probabilities for each category.
# generates the initial state vector for the simulation
state <- as.vector(rmultinom(1,size=N,prob=c(0.2,0.8)))
# Assigns names to the elements of the state vector,
names(state) <- names(state.orig)
# conduct the stochastic simulation.
# It uses the replicate function to run the (betweenhost.stoch) function n.sims times, passing the initial state and parameters parms. The t function transposes the resulting matrix for each simulation, and the output is stored in the out variable.
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
alphaval <- 0.15
par(mfrow=c(1,1))
plot(NULL,
     las=1,xaxs="i",yaxs="i",ylim=c(0,10),bty="n",xlim=range(tvec),
     xlab="Day",ylab="Number")
for (ii in 1:n.sims){
  lines(out["time",,ii],out["S",,ii],col=adjustcolor(cols[1],alphaval))
  lines(out["time",,ii],out["I",,ii],col=adjustcolor(cols[2],alphaval))
}
legend("topleft",bty="n",col=cols[1:2],legend=c("Susceptible","Infected"),
       lwd=3)
out.mean <- apply(out,c(1,2),mean)
names(out.mean)=names(state)
lines(out.mean["time",],out.mean["S",],col=adjustcolor(cols[1]),lwd=3)
lines(out.mean["time",],out.mean["I",],col=adjustcolor(cols[2]),lwd=3)

