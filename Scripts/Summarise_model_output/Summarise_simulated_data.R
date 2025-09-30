####################################################################
# SUMMARISE SIMULATED DATA OUTPUT
####################################################################

rm(list=ls())

os = readRDS("./Output/Model_results/Validation_with_simulated_data/Two_step_sine_seasonality.rds")
ons = readRDS("./Output/Model_results/Validation_with_simulated_data/Two_step_sine_noseasonality.rds")

