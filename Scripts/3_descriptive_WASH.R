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
hh = read.csv(paste0(DirectoryData,"/FINAL_FOR_SHARING/Household_WASH_BF.csv"))
hh_ind = read.csv(paste0(DirectoryData,"/FINAL_FOR_SHARING/Household_stool_WASH_BF.csv"))

load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_all.rda"))
load(paste0(DirectoryData,"/use_in_analyses/bf_esbl0123_long_completecases_UPDATE.rda"))


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

dfls0 = dfls0 %>% mutate(
  agegr = factor(agegr, levels=c("0-4", "5-17", "18-49", "50+"))
)

dfls0complete = dfls0complete %>% mutate(
  agegr = factor(agegr, levels=c("0-4", "5-17", "18-49", "50+"))
)

#---------------------------------------------------------
# DESCRIPTIVES
#----------------------------------------------------------
# R0 characteristics 
table_wash = table1(~ n.householdmember +
                      n.child.0to5 +
                      n.households.concession +
                      #n.households +
                      ses.quintile +
                      main.drinking.water.dry.binary+
                      main.drinking.water.rainy.binary+ 
                      cleaning.water.storage.binary+
                      correct.handwashing.binary+
                      improved.sanitation.binary+
                      livestock.access.house.binary+
                      animal.excrement.floor.binary+
                      factor(q1.diar.prev.water.pot.covered) +
                      #q1_diarrhee_prevenu +
                      factor(q1.diar.prev.no.finger.in.waterglass) +  
                      factor(q1.diar.prev.utensil.to.take.water.from.pot) + 
                      #factor(q1_diarrhee_prevenu.filtrer_l_eau_de_boisson) +
                      factor(q1.diar.prev.cover.food) +
                      factor(q1.diar.prev.boil.water) + 
                      factor(q1.diar.prev.filter.water) + 
                      factor(q1.diar.prev.other)+
                      factor(q1.diar.prev.cant.be.avoided) + 
                      factor(q1.diar.prev.dont.know) + 
                      factor(q2.main.water.source.dry) + 
                      factor(q3.main.water.source.rainy) + 
                      factor(q4.cans.storage.water) + 
                      factor(q5a.cans.storage.water.closed.when.filled) +
                      factor(q5b.cans.storage.water.closed.wen.empty) + 
                      factor(q5c.cans.cleaned.before.reuse) +
                      factor(q6.treatment.water) + 
                      factor(q7.principle.defication) + 
                      factor(q8.other.defecation.flush.toiled.septic) + 
                      factor(q8.other.defecation.pit.latrine.ventilation) + 
                      factor(q8.other.defecation.pit.latrine.slab) + 
                      factor(q8.other.defecation.pit.latrine.no.slab) +
                      factor(q8.other.defecation.open.defecation) + 
                      factor(q8.other.defecation.other) +                                
                      #q8.other.defecation.none + 
                      #q9.shared.toilet.c + 
                      as.numeric(q10.n.shared.toilet) +
                      factor(q11.toilet.last.cleaned) + 
                      factor(q12.disposal.child.stool) + 
                      factor(q13.disposal.latrine.pit) +
                      factor(q14.handwashing.product) + 
                      factor(q15.handwashing.defecation) + 
                      factor(q16.handwashing.meals) +  
                      factor(q17.animals.around.household) + 
                      factor(q18.animal.inside.cow) +  
                      factor(q18.animal.inside.sheep.goat) + 
                      factor(q18.animal.inside.pig) + 
                      factor(q18.animal.inside.donkey.horse) + 
                      factor(q18.animal.inside.chicken.goose.duck) +  
                      factor(q18.animal.inside.other) +  
                      factor(q19.animal.outside.cow) + 
                      factor(q19.animal.outside.sheep.goat) +
                      factor(q19.animal.outside.pig) +
                      factor(q19.animal.outside.donkey.horse) +
                      factor(q19.animal.outside.chicken.goose.duck) + 
                      factor(q19.animal.outside.other) +   
                      factor(q20.animal.excrement.floor) +    
                      factor(q21.when.animal.ill.treatment.with.vet) + 
                      factor(q21.when.animal.ill.treatment.without.vet) +
                      factor(q21.when.animal.ill.sell.bucher) +
                      factor(q21.when.animal.ill.slaugter.eat.at.home) +     
                      factor(q21.when.animal.ill.dies.burie.dispose) +
                      factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text), data=dfls0%>%filter(redcap_event_name=="round_0_arm_1"))
table_wash

# R2 characteristics 
table_wash_r2 = table1(~ n.householdmember +
                         n.child.0to5 +
                         n.households.concession +
                         #n.households +
                         ses.quintile +
                         main.drinking.water.dry.binary+
                         main.drinking.water.rainy.binary+ 
                         cleaning.water.storage.binary+
                         correct.handwashing.binary+
                         improved.sanitation.binary+
                         livestock.access.house.binary+
                         animal.excrement.floor.binary+
                         factor(q1.diar.prev.water.pot.covered) +
                         #q1_diarrhee_prevenu +
                         factor(q1.diar.prev.no.finger.in.waterglass) +  
                         factor(q1.diar.prev.utensil.to.take.water.from.pot) + 
                         #factor(q1_diarrhee_prevenu.filtrer_l_eau_de_boisson) +
                         factor(q1.diar.prev.cover.food) +
                         factor(q1.diar.prev.boil.water) + 
                         factor(q1.diar.prev.filter.water) + 
                         factor(q1.diar.prev.other)+
                         factor(q1.diar.prev.cant.be.avoided) + 
                         factor(q1.diar.prev.dont.know) + 
                         factor(q2.main.water.source.dry) + 
                         factor(q3.main.water.source.rainy) + 
                         factor(q4.cans.storage.water) + 
                         factor(q5a.cans.storage.water.closed.when.filled) +
                         factor(q5b.cans.storage.water.closed.wen.empty) + 
                         factor(q5c.cans.cleaned.before.reuse) +
                         factor(q6.treatment.water) + 
                         factor(q7.principle.defication) + 
                         factor(q8.other.defecation.flush.toiled.septic) + 
                         factor(q8.other.defecation.pit.latrine.ventilation) + 
                         factor(q8.other.defecation.pit.latrine.slab) + 
                         factor(q8.other.defecation.pit.latrine.no.slab) +
                         factor(q8.other.defecation.open.defecation) + 
                         factor(q8.other.defecation.other) +                                
                         #q8.other.defecation.none + 
                         #q9.shared.toilet.c + 
                         as.numeric(q10.n.shared.toilet) +
                         factor(q11.toilet.last.cleaned) + 
                         factor(q12.disposal.child.stool) + 
                         factor(q13.disposal.latrine.pit) +
                         factor(q14.handwashing.product) + 
                         factor(q15.handwashing.defecation) + 
                         factor(q16.handwashing.meals) +  
                         factor(q17.animals.around.household) + 
                         factor(q18.animal.inside.cow) +  
                         factor(q18.animal.inside.sheep.goat) + 
                         factor(q18.animal.inside.pig) + 
                         factor(q18.animal.inside.donkey.horse) + 
                         factor(q18.animal.inside.chicken.goose.duck) +  
                         factor(q18.animal.inside.other) +  
                         factor(q19.animal.outside.cow) + 
                         factor(q19.animal.outside.sheep.goat) +
                         factor(q19.animal.outside.pig) +
                         factor(q19.animal.outside.donkey.horse) +
                         factor(q19.animal.outside.chicken.goose.duck) + 
                         factor(q19.animal.outside.other) +   
                         factor(q20.animal.excrement.floor) +    
                         factor(q21.when.animal.ill.treatment.with.vet) + 
                         factor(q21.when.animal.ill.treatment.without.vet) +
                         factor(q21.when.animal.ill.sell.bucher) +
                         factor(q21.when.animal.ill.slaugter.eat.at.home) +     
                         factor(q21.when.animal.ill.dies.burie.dispose) +
                         factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text),na.rm=T, data=dfls0%>%filter(redcap_event_name=="round_3_arm_1"))
table_wash_r2

# FOR PAPER - BASELINE
table_wash_sel = table1(~ n.householdmember +
                          n.child.0to5 +
                          n.households + 
                          n.population +
                          ses.quintile +
                          age +  
                          agegr +
                          sexe +
                          main.drinking.water.dry.binary+
                          main.drinking.water.rainy.binary+ 
                          cleaning.water.storage.binary+
                          correct.handwashing.binary+
                          improved.sanitation.binary+
                          livestock.access.house.binary+
                          animal.excrement.floor.binary
                        | factor(intervention.text), data=dfls0%>%filter(redcap_event_name=="round_0_arm_1"))
table_wash_sel

# TABLE 1 - Post intervention
table_wash_r2_sel = table1(~ n.householdmember +
                          n.child.0to5 +
                            n.households + 
                            n.population +
                            ses.quintile +
                          age +  
                          sexe +
                          main.drinking.water.dry.binary+
                          main.drinking.water.rainy.binary+ 
                          cleaning.water.storage.binary+
                          correct.handwashing.binary+
                          improved.sanitation.binary+
                          livestock.access.house.binary+
                          animal.excrement.floor.binary
                        | factor(intervention.text), data=dfls0%>%filter(redcap_event_name=="round_3_arm_1"))
table_wash_r2_sel

table_wash_sel_complete = table1(~ n.householdmember +
                          n.child.0to5 +
                          n.households + 
                          n.population +
                          ses.quintile +
                          age +  
                          agegr +
                          sexe +
                          main.drinking.water.dry.binary+
                          main.drinking.water.rainy.binary+ 
                          cleaning.water.storage.binary+
                          correct.handwashing.binary+
                          improved.sanitation.binary+
                          livestock.access.house.binary+
                          animal.excrement.floor.binary
                        | factor(intervention.text), data=dfls0complete%>%filter(redcap_event_name=="round_0_arm_1"))

table_wash_sel_r2_complete = table1(~ n.householdmember +
                                   n.child.0to5 +
                                   n.households + 
                                   n.population +
                                   ses.quintile +
                                   age +  
                                   agegr +
                                   sexe +
                                   main.drinking.water.dry.binary+
                                   main.drinking.water.rainy.binary+ 
                                   cleaning.water.storage.binary+
                                   correct.handwashing.binary+
                                   improved.sanitation.binary+
                                   livestock.access.house.binary+
                                   animal.excrement.floor.binary
                                 | factor(intervention.text), data=dfls0complete%>%filter(redcap_event_name=="round_3_arm_1"))

table_wash_sel_complete
table_wash_sel_r2_complete

table_wash_sel

t1flex(table_wash) %>% 
  save_as_docx(path="./Output/Tables/BF/wash_baseline.docx")
t1flex(table_wash_sel) %>% 
  save_as_docx(path="./Output/Tables/BF/wash_baseline_table1.docx")

write.table(table_wash_sel, "./Output/Tables/BF/wash_baseline_table1.csv", col.names = T, row.names=F, append= F, sep=';')
write.table(table_wash_sel_complete, "./Output/Tables/BF/wash_baseline_table1_completecases.csv", col.names = T, row.names=F, append= F, sep=';')
write.table(table_wash_r2_sel, "./Output/Tables/BF/wash_post-intervention.csv", col.names = T, row.names=F, append= F, sep=';')
write.table(table_wash_sel_r2_complete, "./Output/Tables/BF/wash_post-intervention_completecases.csv", col.names = T, row.names=F, append= F, sep=';')

# NEED TO CHECK SES LINKAGE AS I HAVE SOME MISSING in ROUND 1!

hh %>% filter(redcap_event_name=="round_0_arm_1") %>%
  group_by(intervention.text) %>%
  summarise(n = length(unique(menage_id))
)

# Number of households per cluster
hh %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(intervention.text, village_name) %>%
  summarise(n = n_distinct(menage_id)) %>%
  group_by(intervention.text) %>%
  summarise(
    median = median(n),
    q1 = quantile(n, probs = 0.25),
    q3 = quantile(n, probs = 0.75)
  )

# Number of individuals per intervention arm
dfls0 %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(intervention.text) %>%
  summarise(n = n_distinct(menage_id_member)) #%>%
  #group_by(intervention.text) %>%
  #summarise(
  #  median = median(n),
  #  q1 = quantile(n, probs = 0.25),
  #  q3 = quantile(n, probs = 0.75)
  #)

dfls0complete %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(intervention.text) %>%
  summarise(n = n_distinct(menage_id_member)) #%>%
#group_by(intervention.text) %>%
#summarise(
#  median = median(n),
#  q1 = quantile(n, probs = 0.25),
#  q3 = quantile(n, probs = 0.75)
#)

hh %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(village_name) %>%
  summarise(n = n_distinct(menage_id)) %>%
  summarise(
    median = median(n),
    q1 = quantile(n, probs = 0.25),
    q3 = quantile(n, probs = 0.75)
  )

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
# DESCRIBE PREVALENCE PER HOUSEHOLD AND VILLAGE
#-----------------------------------------------------------

# Compute ESBL prevalence per household per round
hh_prev <- dfls0 %>%
  group_by(redcap_event_name, menage_id) %>%
  summarise(
    n.positive = sum(esble, na.rm = TRUE),
    n.tested = unique(n.tested),
    prevalence = n.positive / unique(n.tested)
  )

hh_prev = left_join(hh_prev, hh %>%select(c(menage_id, redcap_event_name,main.drinking.water.dry.binary,main.drinking.water.rainy.binary,
                                            cleaning.water.storage.binary,correct.handwashing.binary,improved.sanitation.binary,
                                            livestock.access.house.binary,animal.excrement.floor.binary)))

hh_prev_long <- hh_prev %>% filter(redcap_event_name%in%c("round_0_arm_1","round_3_arm_1"))%>%
  pivot_longer(
    cols = c(main.drinking.water.dry.binary, main.drinking.water.rainy.binary,
             cleaning.water.storage.binary, correct.handwashing.binary,
             improved.sanitation.binary, livestock.access.house.binary,
             animal.excrement.floor.binary),
    names_to = "WASH_indicator",
    values_to = "WASH_value"
  ) %>%
  mutate(WASH_value = ifelse(is.na(WASH_value), "Missing", as.character(WASH_value)),
    WASH_value = factor(WASH_value, levels=c("Improved", "Unimproved", "Yes", "No", "Missing")))

p5 = ggplot(data = hh_prev_long %>% filter(redcap_event_name=="round_0_arm_1"), aes(x = factor(WASH_indicator), y = prevalence, fill=WASH_value)) +
  theme_minimal()+
  geom_boxplot() +
  #geom_jitter(width = 0.3, alpha = 0.5, aes(col=WASH_value)) +
  labs(x = "WASH Indicator", y = "Prevalence per Household",
       title ="Baseline") +
  theme(
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  
    panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
    axis.text.x = element_text(angle = 0, size = 10),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
    legend.position = "bottom"
  ) +
  #facet_wrap(~redcap_event_name, ncol=2) +
  coord_flip()
p5

p6 = ggplot(data = hh_prev_long %>% filter(redcap_event_name=="round_3_arm_1"), aes(x = factor(WASH_indicator), y = prevalence, fill=WASH_value)) +
  theme_minimal()+
  geom_boxplot() +
  #geom_jitter(width = 0.3, alpha = 0.5, aes(col=WASH_value)) +
  labs(x = "WASH Indicator", y = "Prevalence per Household",
       title ="Post-intervention") +
  theme(
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed"),  
    panel.grid.minor = element_line(colour = "gray95", linetype = "dotted"),
    axis.text.x = element_text(angle = 0, size = 10),  # Rotate x-axis labels
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    strip.text = element_text(size = 14, face = "bold"),  # Increase facet label size
    legend.position = "bottom"
  ) +
  #facet_wrap(~redcap_event_name, ncol=2) +
  coord_flip()
p6

png("./Output/Figures/BF/WASH_prevalence_pre.png",width = 6, height = 4, units = "in", res = 300)  # Start PNG device
multiplot(p5, cols = 1)  # Adjust 'cols' to control layout
dev.off()  # Close PNG device

png("./Output/Figures/BF/WASH_prevalence_post.png",width = 6, height = 4, units = "in", res = 300)  # Start PNG device
multiplot(p6, cols = 1)  # Adjust 'cols' to control layout
dev.off()  # Close PNG device
