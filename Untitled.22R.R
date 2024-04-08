# Given data
infected <- c(37, 10, 73, 130, 55)
deaths <- c(10, 1, 21, 30, 12)

# Function to calculate negative log likelihood
calculate_nll <- function(p, infected, deaths) {
  # Calculate individual likelihoods for each outbreak
  individual_likelihoods <- dbinom(deaths, size = infected, prob = p, log = TRUE)
  
  # Sum the log likelihoods
  total_nll <- -sum(individual_likelihoods)
  
  return(total_nll)
}

# Values of p between 0 and 1 with steps of 0.01
p_values <- seq(0, 1, 0.01)

# Calculate NLL for each value of p
nll_values <- sapply(p_values, function(p) calculate_nll(p, infected, deaths))

# Plot NLL
plot(p_values, nll_values, type = "l", col = "blue", lwd = 2,
     xlab = "Probability of Death given Infected (p)",
     ylab = "Negative Log Likelihood",
     main = "Negative Log Likelihood vs. Probability of Death")

# Find the minimum NLL and corresponding p
min_nll <- min(nll_values)
optimal_p <- p_values[which.min(nll_values)]

# Highlight the minimum NLL point on the plot
points(optimal_p, min_nll, col = "red", pch = 16)
text(optimal_p, min_nll, sprintf("Min NLL at p = %.2f", optimal_p), pos = 4, col = "red")

N<-c(37, 10, 73, 130, 55 )
d<-c(10, 1,  21, 30, 12)
NLL<-NULL
p<-seq(0, 1, .01)
for(p_i in p){
  NLL <- c(NLL,-sum(dbinom(d,  N, p_i, log = TRUE)))
}
# Assuming params[2] is the log of the reciprocal of the infectious period
log_reciprocal_infectious_period <- -2  # Replace with your actual value

# Taking the exponential to get the reciprocal of the infectious period
reciprocal_infectious_period <- exp(log_reciprocal_infectious_period)

# Print the results
cat("Log of Reciprocal Infectious Period:", log_reciprocal_infectious_period, "\n")
cat("Reciprocal Infectious Period:", reciprocal_infectious_period, "\n")
# Assuming the observed number of infections on Tuesday is 10
fakedataforTuesday <- data.frame(cases = 10)

# Let's say the model predicts the following number of new infections for each time point
# This is just an illustrative example; you would replace this with your actual predictions
new.infections <- c(5, 8, 12, 10, 15)

# Assuming the probability of observing an infection is 0.2 (20%)
probability.of.observing.infection <- 0.2

# Calculate the negative log-likelihood using the binomial distribution
NLL <- -dbinom(fakedataforTuesday$cases, new.infections, probability.of.observing.infection, log=TRUE)

# Print the result
print(NLL)


