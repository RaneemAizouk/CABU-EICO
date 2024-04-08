# ... (your existing code)

# 1) Create a list for the n.sim * n.clusters * houses_per_cluster
total_houses <- n.sims * n.clusters * houses_per_cluster
time_points <- c(3, 9, 12)  # Time points to check

# Initialize an empty list to store results
chosen_samples <- list()

# 2) Loop over each time point
for (time_point in time_points) {
  # Find the closest time point in your simulation
  closest_time <- sapply(results_list, function(result) result$t[which.min(abs(result$t - time_point))])
  
  # Choose the sample from the simulation which is closer to the defined time point
  closest_sample <- results_list[[which.min(abs(closest_time - time_point))]]
  
  # 3) Take the sample from the susceptible and infected individuals for the chosen time points
  chosen_samples[[length(chosen_samples) + 1]] <- data.frame(
    t = floor(closest_sample$t),
    S = closest_sample$S,
    I = closest_sample$I
  )
}

# 4) Run the Gillespie algorithm on the chosen sample
# (assuming the Gillespie function is available, replace with your actual Gillespie function)
gillespie_result <- gillespe(lambda1, lambda2, rho, chosen_samples[[1]]$S[1], chosen_samples[[1]]$I[1], tmax, N)

# 5) Draw a graph for the sample
# (assuming you want to plot the Gillespie result, replace with your actual plotting code)
par(mfrow = c(1, 1))
plot(gillespie_result$t, gillespie_result$S, type = "l", col = "blue", xlab = "Time", ylab = "Population", main = "Gillespie Algorithm Result")
lines(gillespie_result$t, gillespie_result$I, col = "red")

# ... (your existing code)

