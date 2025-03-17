##############################################################
# DESCRIPTIVES DRIVERS VS COLONISATION DATA
##############################################################

# Author: Esther van Kleef
# Date: 17 March 2025

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, msm, skimr,scales)

# SET DIRECTORY
DirectoryData <- "./Data/BF/clean"

source("./Scripts/functions/multiplot.R")

# Load household data
hh = read.csv(paste0(DirectoryData,"/individual_datasets_final/Household_WASH_BF.csv"))
hh_ind = read.csv(paste0(DirectoryData,"/individual_datasets_final/Individual_WASH_BF.csv"))

# Load stool data
col = read.csv(paste0(DirectoryData,"/individual_datasets_final/Individual_Stool_BF.csv"))

# Link stool data with household data
d_r1 = left_join(col%>%filter(redcap_event_name=="round_0_arm_1"), hh_ind%>%filter(redcap_event_name=="round_0_arm_1"))
d_r4 = left_join(col%>%filter(redcap_event_name=="round_3_arm_1"), hh_ind%>%filter(redcap_event_name=="round_3_arm_1"))

# PREAMBLE
#-----------------------------------------------------------

# Define custom order for bars
custom_order <- c(
  "main.drinking.water.dry.binary",
  "main.drinking.water.rainy.binary",
  "cleaning.water.storage.binary",
  "improved.sanitation.binary",
  "correct.handwashing.binary",
  "livestock.access.house.binary",
  "animal.excrement.floor.binary"
)
# Make long format for WASH INDICATORS
# Convert variables to factor with custom order
hh_l <- hh %>%
  select(menage_id, village_name, intervention.text, redcap_event_name, date.enquete,
         main.drinking.water.dry.binary, main.drinking.water.rainy.binary, 
         cleaning.water.storage.binary, correct.handwashing.binary, 
         improved.sanitation.binary, livestock.access.house.binary, 
         animal.excrement.floor.binary) %>%
  mutate(menage_id = row_number()) %>%  # Create a unique ID if none exists
  pivot_longer(
    cols = all_of(custom_order), 
    names_to = "variable", 
    values_to = "value"
  ) %>%
  mutate(variable = factor(variable, levels = custom_order))  # Apply custom order

# VISUALISE TWO ROUNDS
#-------------------------------------------------------------------

# Summarize proportions per round
df_summary <- hh_l %>%
  group_by(redcap_event_name, variable, intervention.text) %>%
  summarise(prop_yes = mean(value == "Yes" | value == "Improved", na.rm = TRUE), .groups = "drop")

color_palette <- c("round_0_arm_1" = "steelblue", "round_3_arm_1" = "tomato")

# Plot with updated theme
p = ggplot(df_summary, aes(x = variable, y = prop_yes, fill = redcap_event_name)) +
  geom_col(position = "dodge") +
  facet_wrap(~intervention.text)+
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "WASH Indicator", 
    y = "Proportion of Yes or Improved responses", 
    title = "WASH Indicators Baseline vs Post-intervention",
    fill = "Survey round"  # Change legend title
  ) +
  scale_fill_manual(
    values = color_palette,    # Apply custom colors
    labels = c("Baseline", "Post-intervention")      # Change legend labels
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  # Striped background grid
    panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
    axis.text.y = element_text(size=12),  # Improve readability of x-axis labels
    axis.text.x = element_text(size=12),
    legend.position = "bottom"
  ) +
  coord_flip()  # Keep flipped coordinates
p

df_diff <- df_summary %>%
  pivot_wider(names_from = redcap_event_name, values_from = prop_yes) %>%
  mutate(Difference = round_0_arm_1 - round_3_arm_1)

p1 = ggplot(df_diff, aes(x = variable, y = Difference, fill = intervention.text)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "WASH Indicator", 
    y = "Difference (Baseline - Post)",
    title = "Change in WASH Indicators",
    subtitle ="Intervention vs Control",
    fill = "Group"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c("darkgreen", "coral"),    # Apply custom colors
    labels = c("Controle", "Intervention")      # Change legend labels
  ) +
  theme(
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  # Striped background grid
    panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
    axis.text.y = element_text(size=12),  # Improve readability of x-axis labels
    axis.text.x = element_text(size=12),
    legend.position = "bottom"
  )+
  coord_flip()
p1

pdf("./Output/Figures/BF/WASH_pre_post_diff.pdf", width=12,height=5)
multiplot(p, cols=1)
dev.off()

pdf("./Output/Figures/BF/WASH_pre_post_abs_diff.pdf", width=7,height=5)
multiplot(p1, cols=1)
dev.off()

# BY VILLAGE CLUSTER
#----------------------------------------------------------
# Summarize proportions per round per village
df_summary_c <- hh_l %>%
  group_by(redcap_event_name, variable, intervention.text, village_name) %>%
  summarise(prop_yes = mean(value == "Yes" | value == "Improved", na.rm = TRUE), .groups = "drop")

df_summary_c <- df_summary_c %>%
  mutate(x_axis_var = ifelse(intervention.text == "Intervention", variable, village_name))

# Plot
p2 = ggplot(df_summary_c, aes(x = x_axis_var, y = prop_yes, fill = redcap_event_name)) +
  geom_col(position = "dodge") +
  facet_wrap(intervention.text~variable, scales = "free_y") +  # Free x-axis per intervention/control
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "",  # Remove x-axis label to avoid confusion
    y = "Proportion of Yes or Improved responses", 
    title = "WASH Indicators: Intervention vs Control",
    fill = "Survey Round"
  ) +
  scale_fill_manual(
    values = color_palette,    
    labels = c("Baseline", "Post-intervention")  
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  
    panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
    axis.text.x = element_text(size=12),  # Rotate for readability
    axis.text.y = element_text(size=12),
    legend.position = "bottom"
  )+coord_flip()

p2

output_dir <- paste0("./Output/Figures/BF/", output_dir)

output_file <- paste0(output_dir, "/WASH_Multiplot.png")  # Output file path

# Initialize an empty list to store plots
plot_list <- list()

# Loop over each variable and create a plot
for (var in unique(df_summary_c$variable)) {
  
  # Subset data for the current variable
  df_subset <- df_summary_c %>% filter(variable == var)
  
  # Create plot
  p <- ggplot(df_subset, aes(x = x_axis_var, y = prop_yes, fill = redcap_event_name)) +
    geom_col(position = "dodge") +
    facet_wrap(~intervention.text, scales = "free_y") +  
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "",  
      y = "Proportion of Yes or Improved responses", 
      title = var,  # Use variable name as title
      fill = "Survey Round"
    ) +
    scale_fill_manual(
      values = color_palette,    
      labels = c("Baseline", "Post-intervention")  
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  
      panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.position = "bottom"
    ) +
    coord_flip()
  
  # Store the plot in the list
  plot_list[[var]] <- p
}
# Save the combined figure using the multiplot function
png(output_file, width = 12, height = 10, units = "in", res = 300, bg = "white")  # Start PNG device with white background
multiplot(plotlist = plot_list, cols = 2)  # Arrange in 2 columns
dev.off() 

output_file <- paste0(output_dir, "/WASH_Changes_Multiplot.png")  # Output file path


df_diff_c <- df_summary_c %>%
  select(village_name, intervention.text, variable, redcap_event_name, prop_yes) %>%  # Ensure village_name is included
  pivot_wider(names_from = redcap_event_name, values_from = prop_yes) %>%  
  mutate(Difference = `round_0_arm_1` - `round_3_arm_1`) 

# Initialize an empty list to store plots
plot_list <- list()

# Loop over each variable (WASH Indicator)
for (var in unique(df_diff_c$variable)) {
  
  # Subset data for the current variable
  df_subset <- df_diff_c %>% filter(variable == var)
  
  # Create plot for difference per village
  p <- ggplot(df_subset, aes(x = village_name, y = Difference, fill = intervention.text)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = percent_format()) +
    labs(
      x = "Village",  
      y = "Difference (Baseline - Post)", 
      title = paste("Change in:", var),
      subtitle = "Intervention vs Control",
      fill = "Group"
    ) +
    theme_minimal() +
    scale_fill_manual(
      values = c("darkgreen", "coral"),  # Custom colors
      labels = c("Control", "Intervention")  
    ) +
    theme(
      panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  
      panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
      axis.text.x = element_text(size=12),  # Rotate village names for readability
      axis.text.y = element_text(size=12),
      legend.position = "bottom"
    ) +
    coord_flip()  # Flip for better visualization
  
  # Store the plot in the list
  plot_list[[var]] <- p
}


# Save the combined figure using the multiplot function
png(output_file, width = 10, height = 19, units = "in", res = 300)  # Start PNG device
multiplot(plotlist = plot_list, cols = 2)  # Adjust 'cols' to control layout
dev.off()  # Close PNG device




#-----------------------------------------------------------
# DESCRIBE HOUSEHOLD DATA PER CLUSTER
#-----------------------------------------------------------
