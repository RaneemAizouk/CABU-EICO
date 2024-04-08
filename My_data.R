# Creating the N_E_Cframe
N_E_C<- data.frame(
  Egg_complement = c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23), 
  Number_1 = c(0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
  Number_2 = c(2,5,11,5,2,1,4,3,1,2,0,3,2,2,0,0,0,1,0,0),
  Number_3 = c(1,1,3,1,1,0,3,4,6,4,3,4,6,4,2,6,2,1,0,1),
  Number_4 = c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
)
library(dplyr)

# Assuming N_E_C is your data frame with the first column for Egg Complement (E)
# and the next four columns for clutch sizes C1 to C4 counts.
#cf: Clutch factor, theoretical clutch size for comparison.
# this quantifies how well different fixed clutch sizes (1, 2, 3, and 4) match the observed data, with lower SSQ 
SSQ <- function(cf) {
  sum2 <- 0
  for(i in 1:nrow(N_E_C)) {
    sum1 <- 0
    for(j in 1:4) {
      eq <- ((j - cf)^2) * (N_E_C[i, j + 1]) # Adjust column index for clutch sizes
      sum1 <- sum1 + eq
    }
    sum2 <- sum2 + sum1
  }
  return(sum2)
}

# Calculate SSQ for each clutch factor
ssq1 <- SSQ(1)
ssq2 <- SSQ(2)
ssq3 <- SSQ(3)
ssq4 <- SSQ(4)

# Data frame to hold the SSQ values and their comparison
ssqDF <- data.frame(
  cf = c(1, 2, 3, 4),
  ssq = c(ssq1, ssq2, ssq3, ssq4)
) %>%
  mutate(compare = ssq / sum(N_E_C[, 2:ncol(N_E_C)])) # Adjust to sum all clutch sizes

# Display the results
print(ssqDF)
# Assuming N_E_C is your data frame with columns: EggComplement, and four more columns for clutch sizes 1 to 4

#6.2

# Assuming N_E_C is structured with EggComplement as the first column and counts for clutch sizes 1 to 4 in the subsequent columns.

# Initialize a large value for comparison
min_ssq <- 10^6
optimal_params <- list(C1=NULL, C2=NULL, e1=NULL)

# Function to calculate weighted SSQ based on clutch size switch
#the sum of squared differences between observed clutch sizes and a theoretical model 
#threshold e1,You set a threshold e1 to 5, meaning you'll switch from using C1 to C2 when the egg complement exceeds 5.
calculate_ssq_switch <- function(C1, C2, e1, N_E_C) {
  ssq <- 0
  for (i in 1:nrow(N_E_C)) {
    E <- N_E_C$Egg_complement[i]
    for (C in 1:4) {
      actual_clutch_size <- ifelse(E <= e1, C1, C2)
      count <- N_E_C[i, C+1]  # Adjusted for the layout of N_E_C, corrected from 'data' to 'N_E_C'
      ssq <- ssq + ((C - actual_clutch_size)^2) * count
    }
  }
  return(ssq)
}

# Search over all combinations of C1, C2 (where C2 > C1), and e1
for (C1 in 1:3) {
  for (C2 in (C1+1):4) {
    for (e1 in unique(N_E_C$Egg_complement)) {
      current_ssq <- calculate_ssq_switch(C1, C2, e1, N_E_C)
      if (current_ssq < min_ssq) {
        min_ssq <- current_ssq
        optimal_params$C1 <- C1
        optimal_params$C2 <- C2
        optimal_params$e1 <- e1
      }
    }
  }
}

# Display the optimal parameters and their corresponding SSQ
print(optimal_params)
print(min_ssq)
library(ggplot2)

# Example: Plot SSQ for varying e1, with fixed C1 and C2.
C1_fixed <- 2  # Example fixed values
C2_fixed <- 3  # Example fixed values

ssq_values <- numeric(length = length(unique(N_E_C$Egg_complement)))
e1_values <- unique(N_E_C$Egg_complement)

for (e1 in e1_values) {
  ssq_values[which(e1_values == e1)] <- calculate_ssq_switch(C1_fixed, C2_fixed, e1, N_E_C)
}

# Create a data frame for plotting
plot_data <- data.frame(e1 = e1_values, SSQ = ssq_values)

# Plot
ggplot(plot_data, aes(x = e1, y = SSQ)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = paste("SSQ for C1 =", C1_fixed, "and C2 =", C2_fixed),
       x = "Egg Complement Threshold (e1)",
       y = "Sum of Squared Differences (SSQ)")

