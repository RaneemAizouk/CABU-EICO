# Load required library
if (!require(pacman)) install.packages("pacman")
pacman::p_load(viridis)

# Number of rolls
num_rolls <- 10

# Probability of each outcome on the die
probabilities <- rep(1/6, 6)

# Simulate 1000 sets of 10 die rolls
simulations <- rmultinom(1000, size = num_rolls, prob = probabilities)

# Plot the results with default colors and labels
barplot(simulations, beside = TRUE,
        main = "Multinomial Distribution of Die Rolls",
        xlab = "Die Face", ylab = "Frequency",
        col = viridis(10),
        legend.text = TRUE, args.legend = list(x = "topright", bty = "n"))

