# R code for modelling challenge Monday 17th Jan 2021, 1-4pm

# Install and load the deSolve package

data <- readRDS("fakedataforTuesday.RDS")
if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}
library(deSolve)
# View the structure of the loaded data
str(fakedataforTuesday)

# View the first few rows of the data
head(fakedataforTuesday.RDS)
# The model function
sir <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    beta <- gamma * R0
    dSdt <- -beta * S * I / N
    dIdt <- beta * S * I / N - gamma * I
    dcumIdt <- beta * S * I / N  # cumulative new infections
    
    list(c(dSdt, dIdt, dcumIdt))
  })
}

# Set initial parameters and conditions
time <- 0
N0 <- 10000000
proportion.initially.susceptible <- 1
initial.SI <- c(S = proportion.initially.susceptible * N0, I = 20.5, cumI = 0)
values <- c(R0 = 3, gamma = 1/8, N = N0)

# Specify times to output model results
time.out <- seq(0, 365, by = 1)

# Run the ODEs
sir.output <- data.frame(lsoda(
  y = initial.SI,
  times = time.out,
  func = sir,
  parms = values
))

# Function to generate fake data
gen.fake.data <- function(probability.of.observing.infection = 0.01, siroutput.df) {
  num.timepoints <- dim(siroutput.df)[1]
  new.infections <- round(c(siroutput.df$cumI[1], siroutput.df$cumI[2:num.timepoints] - siroutput.df$cumI[1:(num.timepoints - 1)]))
  fake.data <- data.frame(t = siroutput.df$time, cases = rbinom(rep(1, num.timepoints), new.infections, probability.of.observing.infection))
  return(fake.data)
}

# Generate fake data
fake.data <- gen.fake.data(0.01, sir.output)

# Function to calculate SSQ for SIR model
getSSQforSIRmodel <- function(params) {
  proportion.initially.susceptible <- 1
  initial.SI <- c(S = proportion.initially.susceptible * N0, I = 20.5, cumI = 0)
  probability.of.observing.infection <- 0.01
  values <- c(R0 = params[1], gamma= params[2], N = N0)
  time.out <- seq(0, 365, by = 1)
  siroutput.df <- data.frame(lsoda(
    y = initial.SI,
    times = time.out,
    func = sir,
    parms = values
  ))
  num.timepoints <- dim(siroutput.df)[1]
  new.infections <- round(c(siroutput.df$cumI[1], siroutput.df$cumI[2:num.timepoints] - siroutput.df$cumI[1:(num.timepoints - 1)]))
  predicted.observed.infections <- new.infections * probability.of.observing.infection
  SSQ <- sum((predicted.observed.infections - fake.data$cases)^2)
  return(SSQ)
}

# Grid search for R0
R0.values <- seq(1, 10, 0.1)
SSQvec <- rep(NA, length(R0.values))
for (i in 1:length(R0.values)) {
  ssq <- getSSQforSIRmodel(c(R0.values[i]))
  SSQvec[i] <- ssq
}

# Plot results
plot(R0.values, SSQvec, type = "l", xlab = "R0", ylab = "SSQ")

# Function to calculate SSQ for SIR model with two parameters
getSSQforSIRmodel2 <- function(params) {
  proportion.initially.susceptible <- 1
  initial.SI <- c(S = proportion.initially.susceptible * N0, I = 20.5, cumI = 0)
  probability.of.observing.infection <- 0.01
  values <- c(R0 = params[1], gamma = params[2], N = N0)
  time.out <- seq(0, 365, by = 1)
  siroutput.df <- data.frame(lsoda(
    y = initial.SI,
    times = time.out,
    func = sir,
    parms = values
  ))
  num.timepoints <- dim(siroutput.df)[1]
  new.infections <- round(c(siroutput.df$cumI[1], siroutput.df$cumI[2:num.timepoints] - siroutput.df$cumI[1:(num.timepoints - 1)]))
  predicted.observed.infections <- new.infections * probability.of.observing.infection
  SSQ <- sum((predicted.observed.infections - fake.data$cases)^2)
  return(SSQ)
}

# Grid search for R0 and gamma
R0.values <- seq(1, 10, 0.1)
gamma.values <- 1/seq(0.5, 16, 0.5)
SSQmat <- matrix(NA, nrow = length(R0.values), ncol = length(gamma.values))
min.ssq <- Inf
for (i in 1:length(R0.values)) {
  for (j in 1:length(gamma.values)) {
    ssq <- getSSQforSIRmodel2(c(R0.values[i], gamma.values[j]))
    SSQmat[i, j] <- ssq
    if (ssq < min.ssq) {
      min.ssq <- ssq
      best.fit.R0 <- R0.values[i]
      best.fit.gamma <- gamma.values[j]
    }
  }
}

# Plot results
image(R0.values, 1/gamma.values, log(SSQmat), col = heat.colors(100), xlab = "R0", ylab = "1/gamma.values")
points(best.fit.R0, 1/best.fit.gamma, pch = 19, col = "green")

# Explore optim starting points
opt.out1 <- optim(c(2, 1/8), getSSQforSIRmodel)
opt.out2 <- optim(c(6, 1/16), getSSQforSIRmodel)
opt.out3 <- optim(c(7, 1/4), getSSQforSIRmodel)

opt.out1
opt.out2
opt.out3

# Explore optim starting points with reduced tolerance
opt.out1_low_tol <- optim(c(3, 1/8), getSSQforSIRmodel, control = list(reltol = 1e-12))
opt.out2_low_tol <- optim(c(6, 1/16), getSSQforSIRmodel, control = list(reltol = 1e-12))
opt.out3_low_tol <- optim(c(1.5, 1/4), getSSQforSIRmodel, control = list(reltol = 1e-12))

opt.out1_low_tol
opt.out2_low_tol
opt.out3_low_tol
p <-seq(0,1,0.1)
L<-dbinom(30,100,p)
print (L)
 
                    