################################################################
# DESCRIPTIVES; RANDOM-EFFECTS MODEL and FREQUENTIST MSM MODEL
################################################################

# 20 April 2024
# Author: Esther van Kleef

rm(list=ls())

# Load functions
source("./Scripts/functions/multiplot.R")

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, summarytools, msm,
               survey,tsModel, mgcv, wesanderson, patchwork)

# Colours
palette <- wes_palette("Darjeeling1", n = 5)
palette2 <- wes_palette("BottleRocket2", n = 1)
palette3 <- wes_palette("GrandBudapest1", n = 2)[2]
palette4 <- wes_palette("BottleRocket2", n = 2)[2]
palette5 = c(palette3, palette[2],palette2,palette[5], palette[4],palette4)

# Load datasets

# ALL INDIVIDUALS - Wide format
df = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_esble_r0123_wide.csv")
length(unique(df$menage_id_member)) # 1237

# ALL INDIVIDUALS - long format
dfl = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_esble_r0123.csv")
length(unique(dfl$menage_id_member)) # 1237

# SELECTED INDIVIDUALS WITH FOUR OBSERVATIONS - long format
load("./Data/BF/clean/use_in_analyses/bf_esbl0123_long_completecases.rda")
length(unique(dfls0$menage_id_member)) # 785

# SELECTED INDIVIDUALS WITH FOUR OBSERVATIONS - wide format
load("./Data/BF/clean/use_in_analyses/bf_esbl_wide_completecases.rda") # 785 individuals
dfc = df[complete.cases(df[, c("r0.esble", "r1.esble", "r2.esble", "r3.esble")]),] # 803 individuals; Difference between dfa and dfc comes from 
# those with a missing date of consent. They are filtered out in the dfa dataset, as that one is the same as the data in long format, where only 
# those with a date of consent present were included (as otherwise the time since intervention could not be calculated for the msm model)

# ADD TIME VARIABLE TO FULL DATASET
dfl = dfl %>% mutate(
  date = as.Date(date, format="%Y-%m-%d"),
  date.stool.collection = as.Date(date.stool.collection, format="%Y-%m-%d"),
  date.consent = as.Date(date.consent, format = "%Y-%m-%d"),
  date.stool.collection = as.Date(date.stool.collection, format="%Y-%m-%d"),
  date.check = date - date.consent,
  date.check.collect = date - date.stool.collection,
  date.use = as.character(date.stool.collection),
  date.use = ifelse(is.na(as.character(date.stool.collection)), 
                    as.character(date.consent),as.character(date.stool.collection)), # This replaced the empty dates with consent date as proxy
  date.use = as.Date(date.use, format = "%Y-%m-%d"),
  month = format(date.use,"%m"),
  date.check.use = date - date.use,
  time = ifelse(redcap_event_name =="round_0_arm_1", 0,
                ifelse(redcap_event_name == "round_1_arm_1", 1,
                       ifelse(redcap_event_name == "round_2_arm_1", 2,3))),
  rainy = ifelse(month%in%c("06","07","08","09"),"yes", "no"))

table(dfl$date.use, useNA="always")
dfl$date.use[is.na(dfl$date.use)] = dfl$date[is.na(dfl$date.use)] 

# ECCMID
contr = dfc %>% filter(intervention.text=="control")
table(contr$r0.esble)

# INTERVENTION DATES
# Round 1: 13 Feb - 26 April 2023
# Round 2: 8 May - 13 June 2023
# Round 3: 21 Aug - 3 Nov 2023

village_order <- dfls0 %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(village_name) %>%
  dplyr::summarize(earliest_date = min(date.use, na.rm = TRUE), .groups = "drop") %>%
  arrange(earliest_date)
village_order
print(n=22,village_order)
unique(dfls0$village_name)

# Step 2: Reorder village_name in the original dataframe
dfls0 <- dfls0 %>%
  mutate(village_name = factor(village_name,
                               levels = c("KOKOLO", "BOULPON","POESSI","SOALA","Zimidin","NAZOANGA","SEGUEDIN",
                                          "DACISSE","KOURIA","GOULOURE","SIGLE",
                                          "BOLOGHO","SOUM","BALOGHO","LALLE","NANORO","MOGODIN","KALWAKA",
                                          "RAKALO","POESSE","SOAW","PELLA"))
  )


# Convert the dates to Date type if not already
vline_dates <- as.Date(c("2023-02-13", "2023-04-26", "2023-05-08", 
                         "2023-06-13", "2023-08-21", "2023-11-03"))

# Define colors for each pair of shaded areas
shaded_colors <- c("red", "blue", "black")

# Plot
p = ggplot(dfls0 %>% filter(intervention.text=="intervention"), aes(x = date.use, fill = redcap_event_name)) +
  geom_bar() +
  facet_grid(rows = vars(village_name), scales = "free_x") +  # Use facet_grid with free x-scales
  labs(x = "Date", y = "Count", title = "Village sampling times",
       subtitle="Intervention group", fill="Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m-%d") +  # Set weekly breaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = vline_dates[c(1, 3, 5)], 
                              xmax = vline_dates[c(2, 4, 6)], 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("round_0_arm_1" = palette5[1], 
                               "round_1_arm_1" = palette5[2], 
                               "round_2_arm_1" =palette5[3], 
                               "round_3_arm_1" = palette5[4])) +
  # Add vertical lines in the specified colors
  geom_vline(xintercept = vline_dates[1], col = "red", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[2], col = "red", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[3], col = "blue", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[4], col = "blue", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[5], col = "black", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[6], col = "black", linetype = 2, size = 0.5)
p   

p2 = ggplot(dfls0 %>% filter(intervention.text == "control"), aes(x = date.use, fill = redcap_event_name)) +
  geom_bar() +
  facet_grid(rows = vars(village_name), scales = "free_x") +  # Use facet_grid with free x-scales
  labs(x = "Date", y = "Count", title = "Village sampling times",
       subtitle="Control group", fill="Sampling round") +
  theme_minimal() +
  scale_x_date(date_breaks = "1 week", date_labels = "%Y-%m-%d") +  # Set weekly breaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  # Add shaded areas between pairs of dates with the same colors as the vertical lines
  geom_rect(data = data.frame(xmin = vline_dates[c(1, 3, 5)], 
                              xmax = vline_dates[c(2, 4, 6)], 
                              ymin = -Inf, ymax = Inf,
                              fill_color = shaded_colors),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.2) +
  scale_fill_manual(values = c("round_0_arm_1" = palette5[1], 
                               "round_1_arm_1" = palette5[2], 
                               "round_2_arm_1" =palette5[3], 
                               "round_3_arm_1" = palette5[4])) +
  # Add vertical lines in the specified colors
  geom_vline(xintercept = vline_dates[1], col = "red", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[2], col = "red", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[3], col = "blue", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[4], col = "blue", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[5], col = "black", linetype = 2, size = 0.5) +
  geom_vline(xintercept = vline_dates[6], col = "black", linetype = 2, size = 0.5)
p2 

combined_plot <- p + p2 + plot_layout(ncol = 2)
combined_plot

ggsave(file = "./Output/Figures/BF/samplingtimes_villages.pdf", 
       plot = combined_plot, width = 32, height = 15)          

# CREATE DATASETS FOR REGRESSION ANALYSES
###############################################################
dfls0$intervention.start= ifelse(dfls0$time%in%c(0,1)&dfls0$intervention.text=="intervention", "control", dfls0$intervention.text) # Co-variate that indicates the start time of the intervention
table(dfls0$time, dfls0$intervention.start)

# Dataset with individuals that have observation at R0, R1 and R2 and where at risk at R2
df_r2 = df %>% filter((!is.na(r1.esble) & !is.na(r2.esble)) & r1.esble == 0 & (r2.esble == 1| r2.esble == 0)) # 395 Include only those at risk at time 2, regardless of whether they have an observation at R3
# Dataset with individuals that have observation at R2 and R3 where at risk at R3
df_r3 = df %>% filter((!is.na(r2.esble)) & r2.esble == 0 & (r3.esble == 1| r3.esble == 0))# 291

# Dataset with individuals that have observation at R0, R1, R2 and R3 where at risk at R2
df_r2_c = dfa %>% filter(r1.esble == 0 & (r2.esble == 1| r2.esble == 0)) # 340; Include only those at risk at time 2, only among those with complete observations
# Dataset with individuals that have observation at R0, R1, R2 and R3 where at risk at R2
df_r3_c = dfa %>% filter(r2.esble == 0 & (r3.esble == 1| r3.esble == 0)) # 267; Include only those at risk at time 2, only among those with complete observations

# Create the same dataset as df_r2_c and df_r3_c with the longformat dataset for for checking
dfls2 = dfls0 %>% filter(time==2 & !is.na(acquisition)) # 338 individuals, why the difference of 2 individuals?
ids = df_r2_c$menage_id_member[which(!df_r2_c$menage_id_member%in%dfls2$menage_id_member)]
#View(dfl[dfl$menage_id_member%in%ids,]) # The difference is coming from those with the date of consent missing for the round of concern. They
# are thrown out in the long format, as that is the dataset for the MSM model which needs a date, to calculate the 'day since intervention'.

# Create the same dataset as df_r2_c with the longformat dataset for for checking
dfls3 = dfls0 %>% filter(time==3 & !is.na(acquisition)) # 264 individuals

# To keep comparison between the long and wide format, and thus logistic regression and MSM model, use the dfls datasets
rm(df_r2_c, df_r3_c, dfc)


#############################################################
# DESCRIPTIVES
#############################################################
# Number of individuals per round
length(unique(df$menage_id_member[!is.na(df$r0.esble)])) # 1209
length(unique(df$menage_id_member[!is.na(df$r1.esble)])) # 1029
length(unique(df$menage_id_member[!is.na(df$r2.esble)])) # 1034
length(unique(df$menage_id_member[!is.na(df$r3.esble)])) # 1006

# Number of interventions per arm
table(df$intervention.text[!is.na(df$r0.esble)]) # Control: 606
table(df$intervention.text[!is.na(df$r1.esble)]) # Control: 522
table(df$intervention.text[!is.na(df$r2.esble)]) # Control: 520
table(df$intervention.text[!is.na(df$r3.esble)]) # Control: 532

# Number of individuals with R0 and R1
length(unique(df$menage_id_member[!is.na(df$r0.esble)&!is.na(df$r1.esble)])) # 1018
table(df$intervention.text[!is.na(df$r0.esble)&!is.na(df$r1.esble)]) # Control: 518

# Number of individuals with R0 and R1 and R2
length(unique(df$menage_id_member[!is.na(df$r0.esble)&!is.na(df$r1.esble)&!is.na(df$r2.esble)])) # 903
table(df$intervention.text[!is.na(df$r0.esble)&!is.na(df$r1.esble)&!is.na(df$r2.esble)]) # Control: 461

# Number of individuals with R0 and R1 and R2 and R3
length(unique(df$menage_id_member[!is.na(df$r0.esble)&!is.na(df$r1.esble)&!is.na(df$r2.esble)&!is.na(df$r3.esble)])) # 803
table(df$intervention.text[!is.na(df$r0.esble)&!is.na(df$r1.esble)&!is.na(df$r2.esble)&!is.na(df$r3.esble)]) # Control: 425

# Number of individuals with R0 and R1 and R2 and R3 AND NO MISSING CONSENT DATA
length(unique(dfa$menage_id_member)) # 785
table(dfa$intervention.text) # Control: 425



# Table 1
dfa0 = dfls0%>%filter(time==0 & !is.na(esble) & !is.na(intervention.text)) # check why these NA's are here
dfa1 = dfls0%>%filter(time==1)
dfa2 = dfls0%>%filter(time==2)
dfa3 = dfls0%>%filter(time==3& !is.na(esble) & !is.na(intervention.text))

# R0 characteristics among those with four observation rounds
wash_r0_table1 = table1(~ factor(esble) + age + sexe + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
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
                          factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text), data=dfa0)
wash_r0_table1

# R0 characteristics among those with four observation rounds
wash_r3_table1 = table1(~ factor(esble) + age + sexe + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
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
                          factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text), data=dfa3)
wash_r3_table1

t1flex(wash_r0_table1) %>% 
  save_as_docx(path="./Output/Tables/BF/wash_r0_table1.docx")
write.table(wash_r0_table1 , "./Output/Tables/BF/wash_r0_table1.csv", col.names = T, row.names=F, append= F, sep=';')


t1flex(wash_r3_table1) %>% 
  save_as_docx(path="./Output/Tables/BF/wash_r3_table1.docx")
write.table(wash_r3_table1 , "./Output/Tables/BF/wash_r3_table1.csv", col.names = T, row.names=F, append= F, sep=';')

# R0 characteristics selected
wash_r0_table1 = table1(~ age + sexe + 
                          n.householdmember +
                          n.tested + 
                          factor(q2.main.water.source.dry) + 
                          factor(q3.main.water.source.rainy) + 
                          factor(q6.treatment.water) + 
                          factor(q7.principle.defication) +
                          factor(q17.animals.around.household) | factor(intervention.text), data=dfa0)
wash_r0_table1

# Baseline characteristics  - acquisitions?
wash_r1_table1 = table1(~ factor(intervention.text) + age + 
                          factor(sexe) + 
                          antibiotic + 
                          watch  + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          factor(q2.main.water.source.dry) + 
                          factor(q3.main.water.source.rainy) + 
                          factor(q6.treatment.water) + 
                          factor(q7.principle.defication) +
                          factor(q17.animals.around.household) +
                          factor(village_name) | factor(acquisition), data=dfa1%>%filter(!is.na(acquisition)))
wash_r1_table1

# R2 characteristics  - acquisitions?
wash_r2_table1 = table1(~ factor(intervention.text) + age + 
                          factor(sexe) + 
                          antibiotic + 
                          watch  + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          factor(q2.main.water.source.dry) + 
                          factor(q3.main.water.source.rainy) + 
                          factor(q6.treatment.water) + 
                          factor(q7.principle.defication) +
                          factor(q17.animals.around.household) +
                          factor(village_name) | factor(acquisition), data=dfa2%>%filter(!is.na(acquisition)))
wash_r2_table1

# Effect of the intervention M9?
# R3 characteristics  - acquisitions?
wash_r3_table1 = table1(~ factor(intervention.text) + age + 
                          factor(sexe) + 
                          antibiotic + 
                          watch  + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          factor(q2.main.water.source.dry) + 
                          factor(q3.main.water.source.rainy) + 
                          factor(q6.treatment.water) + 
                          factor(q7.principle.defication) +
                          factor(q17.animals.around.household) +
                          factor(village_name) | factor(acquisition), data=dfa3%>%filter(!is.na(acquisition)))
wash_r3_table1

###############################################################
# TABLES
###############################################################

da = dfls0 # Change this subject to whether wanting to plot the whole dataset or the selected dataset

# positive per round - full dataset
dfl %>% group_by(time) %>%
  summarise(
    esblsum = sum(esble),
    n = length(unique(menage_id_member)),
    prev = esblsum/n
  )

# positive per round - selected dataset
dfls0 %>% filter(!is.na(intervention.text)) %>%
  group_by(intervention.text, time) %>% 
  summarise(
    esblsum = sum(esble-1),
    n = length(unique(menage_id_member)),
    prev = esblsum/n
  )

# Acquisitions per round
# Of note:
# variable 'acquisition' =  Counts individuals at risk that were for R2, negative at R1; for R3, negative at R2
# variable 'acquisition2' is for a sensitivity analyses =  Counts individuals at risk that were for R2, negative at R1; for R3, either negative or postive at R2 (assuming a colonisation duration of 4M (paper Bootsma))
# variable 'acquisition3' =  Counts individuals at risk that were for R2, negative at R0 and R1; for R3, negative at R2

# ACQUISITIONS AMONG THOSE WITH FOUR OBSERVATION ROUNDS
pa = dfls0 %>% filter(!is.na(intervention.text) &time!=0) %>% group_by(time, intervention.text) %>%
  summarise(
    acqsum = sum(acquisition, na.rm=T),
    acqsum2 = sum(acquisition2, na.rm=T),
    acqsum3 = sum(acquisition3, na.rm=T),
    n = length(unique(menage_id_member)[!is.na(acquisition)]),
    n2 = length(unique(menage_id_member)),
    n3 = length(unique(menage_id_member)[!is.na(acquisition3)])
  )%>%
  mutate(
    n2 = ifelse(time==0, 0,
                ifelse(time%in%c(1,2),n,n2)),
    prev = acqsum/n,
    prev2= acqsum2/n2,
    prev3= acqsum3/n3
  )
pa

# DECOLONISATIONS AMONG THOSE WITH FOUR OBSERVATION ROUNDS
pad = dfls0 %>% filter(!is.na(intervention.text) &time!=0) %>% group_by(time, intervention.text) %>%
  summarise(
    decolsum = sum(decolonisation, na.rm=T),
    decolsum2 = sum(decolonisation2, na.rm=T),
    decolsum3 = sum(decolonisation3, na.rm=T),
    n = length(unique(menage_id_member)[!is.na(decolonisation)]),
    n2 = length(unique(menage_id_member)),
    n3 = length(unique(menage_id_member)[!is.na(decolonisation3)])
  )%>%
  mutate(
    n2 = ifelse(time==0, 0,
                ifelse(time%in%c(1,2),n,n2)),
    prev = decolsum/n,
    prev2= decolsum2/n2,
    prev3= decolsum3/n3
  )
pad


# acquisitions per round - selected dataset among those at risk
da$esbllag = lag(da$esble)
da %>% filter(!is.na(intervention.text)) %>% group_by(time,intervention.text) %>%
  summarise(
    acqsum = sum(acquisition, na.rm=T),
    n = length(unique(menage_id_member)[esbllag!=2]),
    prev = acqsum/n
  )


# decolonisations per round - selected dataset
da %>% filter(!is.na(intervention.text))%>% group_by(time, intervention.text) %>%
  summarise(
    decolsum = sum(decolonisation, na.rm=T),
    n = length(unique(menage_id_member)[esbllag!=1]),
    prev = decolsum/n
  )


# Per household, how many positive and tested
dh = da %>% filter(!is.na(intervention.text)) %>% 
  group_by(time, village_name,intervention.text,menage_id) %>%
  summarise(n.tested = unique(n.tested)[1],
            nesbl = sum(esble-1),
            n.atrisk = n.tested-unique(nesbl_tminus)[1],
            eprev = nesbl/n.tested,
            nacquisition =sum(acquisition,na.rm=T),
            aprev = nacquisition/n.atrisk
            )
dh = dh[order(dh$menage_id,dh$time),] 

dh%>%group_by(time) %>%
  summarise(sum(nesbl),
            sum(nacquisition, na.rm=T),
            sum(n.atrisk),
            sum(n.tested))

dv = dh %>% filter(!is.na(intervention.text)) %>% 
  group_by(time, village_name,intervention.text) %>%
  summarise(n.tested = sum(n.tested),
            n.atrisk = sum(n.atrisk,na.rm=T),
            nesbl = sum(nesbl),
            eprev = nesbl/n.tested,
            nacquisition = sum(nacquisition, na.rm=T),
            aprev = nacquisition/n.atrisk
  )

dv[order(dv$village_name,dv$time),]
dv%>%group_by(time) %>%
  summarise(sum(nesbl),
            sum(nacquisition, na.rm=T),
            sum(n.atrisk),
            sum(n.tested))


# mean and median colonisation prevalence over village cluster
dv %>% filter(!is.na(intervention.text)) %>% 
  group_by(time, intervention.text) %>%
  summarise(
    meanp = mean(eprev),
    lower_ci = ggplot2::mean_cl_normal(eprev, conf.int = 0.95)[2],
    upper_ci = ggplot2::mean_cl_normal(eprev, conf.int = 0.95)[3],
    median = median(eprev),
    q1 = quantile(eprev,probs = 0.25),
    q3 = quantile(eprev,probs=0.75)
  )

# mean and median acquisition prevalence over village cluster
dv %>% filter(!is.na(intervention.text)) %>% 
  group_by(time, intervention.text) %>%
  summarise(
    meanp = mean(aprev,na.rm=T),
    lower_ci = ggplot2::mean_cl_normal(aprev, conf.int = 0.95)[2],
    upper_ci = ggplot2::mean_cl_normal(aprev, conf.int = 0.95)[3],
    median = median(aprev,na.rm=T),
    q1 = quantile(aprev,probs = 0.25,na.rm=T),
    q3 = quantile(aprev,probs=0.75,na.rm=T)
  )

# month vs round
start.collection = dfls0 %>% filter(!is.na(intervention.text)) %>% 
  select(c(time,intervention.text,village_name,menage_id, date.use,month)) %>%
  group_by(time,village_name) %>%
  summarise(start.collection.village = min(date.use,na.rm = T),
            start.collection.month.village = month(start.collection.village)
  )

start.collection

dh = left_join(dh,start.collection, by=c("time", "village_name"))
dv = left_join(dv,start.collection, by=c("time", "village_name"))

table(dv$time,dv$start.collection.month.village)

###############################################################
# FIGURES
###############################################################


# Plot prevalence per household per round
p1=ggplot(dh, aes(x=start.collection.village, y=eprev, fill=menage_id, group=menage_id)) + geom_point(aes(colour=menage_id)) +
  geom_line(aes(colour=menage_id)) + 
  theme_bw() + 
  theme(legend.position="none",
          axis.text.x=element_text(size=15, angle=90),
          axis.title=element_text(size=16,face="bold"),
          axis.text=element_text(size=14,face="bold"),
          strip.text=element_text(size=14,face="bold")) + 
  facet_wrap(.~village_name)

# Plot prevalence per village per round
p2 = ggplot(dv, aes(x=start.collection.village, y=eprev, colour=village_name, 
               group=village_name)) + geom_point() +
  geom_line() + 
  theme_bw() + 
  theme(legend.position="none",
        axis.text.x=element_text(size=15),
        axis.title=element_text(size=16,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        strip.text=element_text(size=14,face="bold")) + 
  facet_wrap(.~village_name,ncol=4)

# Plot acquisitions per village per round
ggplot(dv, aes(x=start.collection.village, y=aprev, colour=village_name, 
               group=village_name)) + geom_point(size=3) +
  geom_line(size=1) + 
  theme_bw() + 
  theme(legend.position="none",
        axis.text.x=element_text(size=15),
        axis.title=element_text(size=16,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        strip.text=element_text(size=14,face="bold")) + 
  facet_wrap(.~intervention.text,ncol=4)

# Number of individuals included per household
hist(table(da$menage_id[da$time==1]))

# Trajectories of individuals
pdf(file = "./Output/Figures/BF/trajectories_all_villages.pdf", width = 10, height = 10)

# Loop through each village and create plots
for(i in unique(da$village_name)){
  d = da %>% filter(village_name == i)
  
  tid = ggplot(d, aes(x = time, y = esble - 1, col = menage_id_member, group = menage_id_member)) + 
    geom_point(size = 2) + 
    geom_line() +
    facet_wrap(~ menage_id_member) +
    theme(legend.position = "none") + 
    ggtitle(paste0("Village ", i))
  
  # Print the plot for each village
  print(tid)
}
dev.off()

# number of esbls by date (as villages have had data collection consequatively. Just to see if some seasonality can be found)
# This to check seasonal effect
dat = dfls0 %>% filter(!is.na(intervention.text)) %>% 
  group_by(date.use,intervention.text) %>%
  summarise(nesbl = sum(esble-1))

ggplot(dat, aes(x=date.use, y=nesbl, col=intervention.text)) + geom_point() + geom_smooth()

dat = dfls0 %>% filter(!is.na(intervention.text)) %>% 
  group_by(month,intervention.text) %>%
  summarise(nesbl = sum(esble-1),
            ntest = length(menage_id_member),
            prev = nesbl/ntest)

dat = dat %>% filter(!is.na(month))

ggplot(dat, aes(x=month, y=nesbl, col=intervention.text,group=intervention.text)) + 
  geom_point(size=2) + geom_smooth() + geom_line() 

bs = ggplot(dat, aes(x=month, y=prev, col=intervention.text,group=intervention.text)) + 
  geom_vline(aes(xintercept =1), lty="dashed", size=1, col=c("grey"))+
  geom_vline(aes(xintercept =4), lty="dashed", size=1, col=c("grey"))+
  geom_vline(aes(xintercept =10), lty="dashed", size=1, col=c("grey"))+
  geom_vline(aes(xintercept =6), lty=1, size=1, col=c("darkgreen"))+
  
  geom_point(size=4) + geom_smooth(size=1) + 
  theme_bw() + 
  theme(
        axis.text.x=element_text(size=15),
        axis.title=element_text(size=16,face="bold"),
        axis.text=element_text(size=14,face="bold"),
        strip.text=element_text(size=14,face="bold"))+
  labs(title = "Smoothed fit of ESBL-E positive (%) per month",
       x = "Month",
       y = "% positive")# There appears to be a seasonal effect
bs


median_t = dv %>% filter(time==0) %>% group_by(intervention.text) %>%
  summarise(
    median = median(eprev),
    mean=mean(eprev)
  )

median_t

# Plot boxplot of fraction positive per village per intervention group based on household prevalence
bp = dh %>% filter(time==1& !n.tested %in%c(0,1)) %>% 
  ggplot(., aes(x = reorder(village_name, eprev),village_name, y = eprev, fill = village_name)) +
  geom_jitter(alpha=0.5) + theme_bw() +
  # geom_hline(yintercept=median, 
  #            color = "red", size=2) +
  geom_boxplot() + theme(legend.position="none",
                         axis.text.x=element_text(size=15, angle=90),
                         axis.title=element_text(size=16,face="bold"),
                         axis.text=element_text(size=14,face="bold"),
                         strip.text=element_text(size=14,face="bold")) +
  ylim(0,1) + 
  facet_wrap(~intervention.text, scales=("free_x")) + 
  geom_abline(data = median_t, aes(intercept = as.numeric(median), slope = 0),lty="dashed", size=1, col="red")+ 
  labs(title = "Boxplot of ESBL-E positive (%) per village cluster (Baseline)",
       x = "Village",
       y = "% positive")
print(bp)

# Number of individuals tested vs esbl positive - is there a pattern? If not 
# with prevalence then good
ggplot(dh, aes(x=n.tested, y=nesbl, col=factor(time))) + geom_point() + geom_jitter()
ggplot(dh, aes(x=n.tested, y=eprev, col=factor(time))) + geom_point() + geom_jitter()

# Boxplot of eprev by intervention per round
bpi = dv %>% 
  ggplot(., aes(x = factor(time), y = eprev, fill = intervention.text)) +
  geom_abline(data = median_t, aes(intercept = as.numeric(median), slope = 0),lty="dashed", size=1, col=c("red", "turquoise"))+ 
  geom_vline(data = dv%>% filter(time==1), aes(xintercept =time+0.5), lty="dashed", size=1, col=c("grey"))+
  geom_jitter(alpha=0.5,size=3) + theme_bw() +
  # geom_hline(yintercept=median, 
  #            color = "red", size=2) +
  geom_boxplot() + theme(
                         axis.text.x=element_text(size=15),
                         axis.title=element_text(size=16,face="bold"),
                         axis.text=element_text(size=14,face="bold"),
                         strip.text=element_text(size=14,face="bold")) +
  ylim(0,1) + 
  labs(title = "Boxplot of fraction ESBL-E positive across village clusters",
       x = "Round",
       y = "Positive/n tested")
print(bpi)

median_ta = dv %>% filter(time==1) %>% group_by(intervention.text) %>%
  summarise(
    median = median(aprev),
    mean=mean(aprev)
  )

median_ta

bai =  dv %>% filter(time!=0) %>%
  ggplot(., aes(x = factor(time), y = aprev, fill = intervention.text)) +
  geom_abline(data = median_ta, aes(intercept = as.numeric(median), slope = 0),lty="dashed", size=1, col=c("red", "turquoise"))+ 
  geom_jitter(alpha=0.5,size=3) + theme_bw() +
  # geom_hline(yintercept=median, 
  #            color = "red", size=2) +
  geom_boxplot() + theme(
    axis.text.x=element_text(size=15),
    axis.title=element_text(size=16,face="bold"),
    axis.text=element_text(size=14,face="bold"),
    strip.text=element_text(size=14,face="bold")) +
  ylim(0,1) + 
  labs(title = "Boxplot of ESBL-E acquisitions across village clusters",
       x = "Round",
       y = "Acquired/at risk")
print(bai)

pdf(file="./Output/Figures/BF/prevalence_per_village.pdf", width=10,height=7)
print(p1)
dev.off()

pdf(file="./Output/Figures/BF/prevalence_per_village_agg.pdf", width=10,height=7)
print(p2)
dev.off()


pdf(file="./Output/Figures/BF/prevalence_per_month.pdf", width=10,height=7)
print(bs)
dev.off()

pdf(file="./Output/Figures/BF/prevalence_per_village_meanv.pdf", width=10,height=7)
print(bp)
dev.off()

pdf(file="./Output/Figures/BF/prevalence_per_round_meanv.pdf", width=10,height=7)
print(bpi)
dev.off()

pdf(file="./Output/Figures/BF/prevalence_acq_per_round_meanv.pdf", width=10,height=7)
print(bai)
dev.off()

pdf(file="./Output/Figures/BF/prevalence_esbl_acq_per_village_meanv.pdf", width=15,height=5)
print(multiplot(bpi, bai, cols=2))
dev.off()

pdf(file="./Output/Figures/BF/trajectories_indv_v.pdf", width=15,height=15)
print(tid)
dev.off()

#####################################################################
# REGRESSION MODEL
#####################################################################

# Note that if acquisition occurs in between screening round 0 (baseline) and round 1 (at time of intervention)
# then an acquisition cannot occur between round 1 and round 2 (3 month later) because individual will 
# already be positive at round 1. And if already positive at round 1 individual will not be at risk 
# of a post-intervention acquisition

# Proposed analysi 1 s -  include in analysis only if a negative swab at round one
# i.e r1_esbl_pos =0   and a swab (whether pos or neg) at round 2 (i.e r2_esbl_pos is 1 or 0)
#  outcome is 1 if r2_esbl_pos  
# adjust for clustering at household level (menage_id) and include covariates for age (age), sex (sexe) and 
# intervention (as coded in field intervention_text)
# then as a sensitivity analysis consider only those with two initial negative swab

################################################################
# RANDOM-EFFECTS MODEL - Acquisition
################################################################

# first with no random effects or other covariates
# first with random effects and no other covariates
# first with random effects and other covariates
# table(dfls2$acquisition)
# 
# # Fit odds ratios for R2 acquisitions only
# ma2_r2 <- glmer(acquisition ~ intervention.text + age + sexe + (1|menage_id), 
#              data=dfls2,  family = binomial(link = "logit"))
# summary(ma2_r2)
# data.frame(exp(fixef(ma2_r2)))
# exp(confint(ma2_r2,method="Wald"))  # Effect intervention: OR = 1.13 (0.64 - 1.79);
# 
# # Fit odds ratios for R3 acquisitions only
# ma2_r3 <- glmer(acquisition ~ intervention.text + age + sexe + (1|menage_id), 
#                 data=dfls3,  family = binomial(link = "logit"))
# summary(ma2_r3)
# data.frame(exp(fixef(ma2_r3)))
# exp(confint(ma2_r3,method="Wald"))  # Effect intervention: OR = 0.90 (0.46 - 1.78); 
# 
# # Fit odds for R3 acquisition - sensitivity analyses
# ma2_r3s <- glmer(acquisition2 ~ intervention.text + age + sexe + (1|menage_id), 
#                 data=dfls3,  family = binomial(link = "logit"))
# summary(ma2_r3s)
# data.frame(exp(fixef(ma2_r3s)))
# exp(confint(ma2_r3s,method="Wald"))  # Effect intervention: OR = 0.82 (0.53 - 1.26); 
# 
# 
# # Now fit both together
# ma2 <- glmer(acquisition ~ intervention.start + age + sexe + (1|menage_id_member), 
#              data=dfls0,  family = binomial(link = "logit"))
# summary(ma2)
# data.frame(exp(fixef(ma2)))
# exp(confint(ma2,method="Wald"))  
# # exp(confint(ma2)) # Failed to converge with bootstrap method
# 
# 
# # MODEL WITH 'SEASONALITY' (CYCLIC SPLINES, R PICK) + covariates
# d= dfls0 %>% mutate(menage_id_member = as.factor(menage_id_member),
#                     village_name = as.factor(village_name))
# 
# # Without seasonality
# ma2gam = gamm(acquisition ~ age + sexe + intervention.start,
#             random=list(village_name=~1,menage_id_member=~1), 
#             family=binomial(link = "logit"), data=d)
# summary(ma2gam$gam)
# 
# # Without seasonality but with time
# ma2tgam = gamm(acquisition ~ age + sexe + intervention.start +
#               factor(time),
#               random=list(village_name=~1, menage_id_member=~1), 
#               family=binomial(link = "logit"), data=d)
# summary(ma2tgam$gam)
# exp(coef(ma2tgam$gam))
# # Without seasonality and different effect at time 2 and 3 by fitting interaction term with time gives error; with days3 not
# # mas2inter_gam = gamm(acquisition ~ age + sexe + intervention.start*factor(time) + 
# #                       s(as.numeric(month), bs="cc"),
# #                     random=list(village_name=~1,menage_id_member=~1), 
# #                     family=binomial(link = "logit"), data=d)
# # 
# # summary(mas2inter_gam$gam)
# d$days3 = d$days
# d$days3[d$time==0] = 0
# 
# mas2inter_gam = gamm(acquisition ~ age + sexe + intervention.start*days3 + 
#                        s(as.numeric(month), bs="cc"),
#                      random=list(village_name=~1,menage_id_member=~1), 
#                      family=binomial(link = "logit"), data=d)
# 
# summary(mas2inter_gam$gam)
# 
# 
# # With seasonality
# mas2gam = gamm(acquisition ~ age + sexe + intervention.start + 
#               s(as.numeric(month), bs="cc"),
#             random=list(village_name=~1,menage_id_member=~1), 
#             family=binomial(link = "logit"), data=d)
# 
# summary(mas2gam$gam)
# 
# 
# 
# # With seasonality - Sensitivity analyses (acquisition R3 = ESBL-E is 0 or 1 at R2)
# mas2sgam = gamm(acquisition2 ~ age + sexe + intervention.start +
#               s(as.numeric(month), bs="cc"),
#             random=list(village_name=~1,menage_id_member=~1), 
#             family=binomial(link = "logit"), data=d)
# 
# # Model fitted during ECCMID
# mas2gam_eccmid = gamm(acquisition ~ age + sexe + intervention.text*factor(time) + 
#                  s(as.numeric(month), bs="cc"),
#                random=list(village_name=~1,menage_id_member=~1), 
#                family=binomial(link = "logit"), data=d)
# 
# summary(mas2gam_eccmid$gam)
# 
# # Model fitted during ECCMID - sensitivity
# mas2sgam_eccmid = gamm(acquisition2 ~ age + sexe + intervention.text*factor(time) + 
#                         s(as.numeric(month), bs="cc"),
#                       random=list(village_name=~1,menage_id_member=~1), 
#                       family=binomial(link = "logit"), data=d)
# 
# summary(mas2sgam_eccmid$gam)
# 
# # PLOT SEASONAL TERM
# # The seasonal term seems to make sense
# plot(mas2gam$gam,
#      trans=plogis,
#      pages=1,
#      seWithMean = T)
# gam.check(mas2$gam)
# 
# # The seasonal term seems to make sense - sensitivity analyses
# plot(mas2sgam$gam,
#      trans=plogis,
#      pages=1,
#      seWithMean = T)
# 
# # The seasonal term seems to make sense - sensitivity analyses
# plot(mas2gam_eccmid$gam,
#      trans=plogis,
#      pages=1,
#      seWithMean = T)
# 
# plot(mas2sgam_eccmid$gam,
#      trans=plogis,
#      pages=1,
#      seWithMean = T)
# 
# # PLOT ODDS RATIOS - MODEL WITHOUT SEASONALITY
# # Extract the model coefficients and standard errors
# coef_summary <- summary(ma2gam$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios
# # May want to bootstrap those
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fitnos = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fitnos$variable =row.names(fitnos)
# fitnos$variable = factor(fitnos$variable, levels=unique(fitnos$variable))
# fitnos
# 
# # PLOT ODDS RATIOS - MODEL WITH SEASONALITY and INTERACTION TERMS DAYS
# # Extract the model coefficients and standard errors
# coef_summary <- summary(mas2inter_gam$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios
# # May want to bootstrap those
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fitinter = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fitinter$variable =row.names(fitinter)
# fitinter$variable = factor(fitinter$variable, levels=unique(fitinter$variable))
# fitinter
# 
# # PLOT ODDS RATIOS - MODEL WITH SEASONALITY
# # Extract the model coefficients and standard errors
# coef_summary <- summary(mas2gam$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios
# # May want to bootstrap those
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fit = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fit$variable =row.names(fit)
# fit$variable = factor(fit$variable, levels=unique(fit$variable))
# fit
# 
# 
# 
# 
# # PLOT ODDS RATIOS - SENSITIVITY analyses
# # Extract the model coefficients and standard errors
# coef_summary <- summary(mas2sgam$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios 
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fits = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fits$variable =row.names(fits)
# fits$variable = factor(fits$variable, levels=unique(fits$variable))
# fits
# 
# # PLOT ODDS RATIOS - ECCMID
# # Extract the model coefficients and standard errors
# coef_summary <- summary(mas2gam_eccmid$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios 
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fit_eccmid = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fit_eccmid$variable =row.names(fit_eccmid)
# fit_eccmid$variable = factor(fit_eccmid$variable, levels=unique(fit_eccmid$variable))
# fit_eccmid
# 
# # PLOT ODDS RATIOS - ECCMID SENSITIVITY ANALYSES
# # Extract the model coefficients and standard errors
# coef_summary <- summary(mas2sgam_eccmid$gam)
# 
# # Extract the coefficients and standard errors
# coefficients <- coef_summary$p.coeff
# standard_errors <- coef_summary$se[1:length(coefficients)]
# 
# # Compute the odds ratios
# odds_ratios <- exp(coefficients)
# 
# # Compute the confidence intervals for the odds ratios 
# lower_ci <- exp(coefficients - 1.96 * standard_errors)
# upper_ci <- exp(coefficients + 1.96 * standard_errors)
# 
# fit_eccmid_s = data.frame(cbind(odds_ratios,lower_ci, upper_ci))
# fit_eccmid_s$variable =row.names(fit_eccmid_s)
# fit_eccmid_s$variable = factor(fit_eccmid_s$variable, levels=unique(fit_eccmid_s$variable))
# fit_eccmid_s
# 
# # Plot the odds ratios with confidence intervals
# ms = ggplot(fits, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition (9M at risk: 3M ESBL-E = 0|1)",
#        x = "Odds ratio",
#        y = "")+
#   xlim(0,3)# +
#   # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
#   #                             "factor(time)2","factor(time)3",
#   #                             "intervention.textintervention:factor(time)2",
#   #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
#   #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
#   #                  "R2",
#   #                             "R3","Intervention: R2",
#   #                             "Intervention: R3"))
# 
# 
# ms
# 
# 
# 
# #fit = melt(fit,value.name="" )
# # Plot the odds ratios with confidence intervals
# mnos = ggplot(fitnos, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition - No seasonality",
#        x = "Odds ratio",
#        y = "")+
#   # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
#   #                             "factor(time)2","factor(time)3",
#   #                             "intervention.textintervention:factor(time)2",
#   #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
#   #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
#   #                             "R2",
#   #                             "R3","Intervention: R2",
#   #                             "Intervention: R3")) +
#   xlim(0,3)
# mnos
# 
# 
# 
# # Plot the odds ratios with confidence intervals
# ms = ggplot(fits, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition - Seasonality (9M at risk: 3M ESBL-E = 0|1)",
#        x = "Odds ratio",
#        y = "")+
#   xlim(0,3)# +
# # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
# #                             "factor(time)2","factor(time)3",
# #                             "intervention.textintervention:factor(time)2",
# #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
# #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
# #                  "R2",
# #                             "R3","Intervention: R2",
# #                             "Intervention: R3"))
# 
# 
# ms
# 
# 
# m = ggplot(fit, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition - Seasonality (9M at risk: 3M ESBL-E = 0)",
#        x = "Odds ratio",
#        y = "")+
#   # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
#   #                             "factor(time)2","factor(time)3",
#   #                             "intervention.textintervention:factor(time)2",
#   #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
#   #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
#   #                             "R2",
#   #                             "R3","Intervention: R2",
#   #                             "Intervention: R3")) +
#   xlim(0,3)
# m
# 
# 
# m_eccmid = ggplot(fit_eccmid, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition - ECCMID (seasonality + interaction time (R1,R2,R3)*intervention.group)",
#        x = "Odds ratio",
#        y = "")
#   # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
#   #                             "factor(time)2","factor(time)3",
#   #                             "intervention.textintervention:factor(time)2",
#   #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
#   #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
#   #                             "R2",
#   #                             "R3","Intervention: R2",
#   #                             "Intervention: R3")) +
#   
# m_eccmid
# 
# m_eccmid_s = ggplot(fit_eccmid_s, aes(x=odds_ratios, y=variable)) + geom_point(size=3)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower_ci, xmax = upper_ci), 
#                  position = position_dodge(width = 0.2),size=1)+
#   theme_bw()+ theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   labs(title = "Odds ratios for ESBL-acquisition - ECCMID (sensitivity, i.e. 9M at risk: 3M ESBL-E = 0|1)",
#        x = "Odds ratio",
#        y = "")
# # scale_y_discrete(breaks = c("age","sexeMale","intervention.textintervention",
# #                             "factor(time)2","factor(time)3",
# #                             "intervention.textintervention:factor(time)2",
# #                             "intervention.textintervention:factor(time)3"),  # Set custom breaks
# #                  labels = c("Age", "Sex: Male", "Intervention: Yes",
# #                             "R2",
# #                             "R3","Intervention: R2",
# #                             "Intervention: R3")) +
# 
# m_eccmid_s
# 
# pdf("./Output/Figures/fit_plotOR.pdf", width=25,height = 20)
# print(multiplot(m,m_eccmid,mnos,ms,m_eccmid_s,cols=2))
# dev.off()
# # Get probability of acquistion at T2 and T3
# ##############################################
# # Assuming 'intervention_group_value' is 1 for intervention group and 0 for control group
# # and 'time_value' is the specific time point for which you want to calculate the probability
# set.seed(40)
# village = sample(unique(d$village_name[d$intervention.text=="intervention"]),1)
# id = sample(unique(d$menage_id_member[d$village_name%in%village]), 1)
# new_data = d %>% filter(menage_id_member==id,time!=0) %>% select(time, month,village_name,
#                                  menage_id,menage_id_member,
#                                  intervention.text,age,sexe)
# 
# 
# new_data$age = median(d$age,na.rm=T)
# new_data$sexe = "Female"
# new_data$time = c(1:3)
# #new_data$month = c("03","07", "11")
# new_data$intervention.text = "intervention"
# #new_data$village_name = "BOLOGHO"
# new_data
# 
# # Get probabilities for time 1 - 3
# predict(mas2$gam,newdata=new_data,type="response")
# 
# # Bootstrap confidence intervals
# #############################################
# 
# # Define a function to obtain predicted values for a given combination of predictor values
# predict_values <- function(model, newdata) {
#   predict(model, newdata = newdata, type = "response")
# }
# 
# # Set the number of bootstrap samples
# num_bootstraps <- 500
# 
# # Initialize a vector to store predicted values
# predicted_probabilities <- matrix(NA, nrow = num_bootstraps, ncol = 6)
# predicted_probabilities_sens <- matrix(NA, nrow = num_bootstraps, ncol = 6)
# 
# set.seed(100)
# id = sample(unique(dfls0$menage_id_member[dfls0$village_name%in%village]), 1)
# village = sample(unique(d$village_name[d$intervention.text=="intervention"]),1)
# id;village
# 
# # Perform bootstrapping
# for (i in 1:num_bootstraps) {
#   # Generate a bootstrap sample by resampling with replacement
#   ids = sample(d$menage_id_member,size=length(unique(d$menage_id_member)),
#                replace=T)
#   ids = data.frame(ids) %>% rename(menage_id_member="ids")
#   bootstrap_data <- left_join(ids,d, by="menage_id_member")
#   bootstrap_data = rbind(bootstrap_data, d%>%filter(menage_id_member==id))
#   new_data = bootstrap_data %>% filter(time!=0)
#   new_data$intervention.text = "intervention"
#   
#   new_data2 = bootstrap_data %>% filter(time!=0)
#   new_data2$intervention.text = "control"
#   new_data = rbind(new_data,new_data2)
#   # new_data = bootstrap_data %>% filter(menage_id_member==id,time!=0) %>% select(time, month,village_name,
#   #                                                                  menage_id,menage_id_member,
#   #                                                                  intervention.text,age,sexe)
#   # 
#   # new_data = new_data %>% filter(!duplicated(time))
#   # new_data$age = median(d$age,na.rm=T)
#   # new_data$sexe = "Female"
#   # new_data$time = c(1:3)
#   # #new_data$month = c("03","07", "11")
#   # new_data$intervention.text = "intervention"
#   # new_data = rbind(new_data,new_data)
#   # new_data$intervention.text = c(rep("intervention",3),  c(rep("control",3)))
#   # #new_data$village_name = "BOLOGHO"
#   # 
#   # Fit the GAMM model to the bootstrap sample
#   bootstrap_model <- gamm(acquisition ~ age + sexe + intervention.text*factor(time) + 
#                             s(as.numeric(month), bs="cc"),
#                           random=list(village_name=~1, menage_id_member=~1), 
#                           family=binomial(link = "logit"), data=bootstrap_data)
#   
#   bootstrap_model_sens <- gamm(acquisition2 ~ age + sexe + intervention.text*factor(time) + 
#                             s(as.numeric(month), bs="cc"),
#                           random=list(village_name=~1, menage_id_member=~1), 
#                           family=binomial(link = "logit"), data=bootstrap_data)
#   
#   # Calculate predicted values using the fitted model
#   
#   p <- predict(bootstrap_model$gam, newdata = new_data, type = "response")
#   p = data.frame(cbind(p, new_data$time, new_data$intervention.text))
#   names(p) = c("prob","time", "intervention.text")
#   p = p %>% group_by(time,intervention.text) %>% summarise(
#     median = median(as.numeric(prob), na.rm=T))
#   
#   predicted_probabilities[i, ] = p$median
#   p <- predict(bootstrap_model_sens$gam, newdata = new_data, type = "response")
#   p = data.frame(cbind(p, new_data$time, new_data$intervention.text))
#   names(p) = c("prob","time", "intervention.text")
#   p = p %>% group_by(time,intervention.text) %>% summarise(
#     median = median(as.numeric(prob), na.rm=T))
#   predicted_probabilities_sens[i, ] = p$median
#   
#   
#   print(i)
# }
# 
# # Calculate confidence intervals using empirical quantiles
# estimate <- apply(predicted_probabilities, 2, quantile, probs = 0.50,na.rm=T)
# lower_bound <- apply(predicted_probabilities, 2, quantile, probs = 0.025,na.rm=T)
# upper_bound <- apply(predicted_probabilities, 2, quantile, probs = 0.975,na.rm=T)
# 
# acq_probs = data.frame(cbind(estimate,lower_bound,upper_bound)) 
# acq_probs
# 
# # Calculate confidence intervals using empirical quantiles
# estimate_s <- apply(predicted_probabilities_sens, 2, quantile, probs = 0.50,na.rm=T)
# lower_bound_s <- apply(predicted_probabilities_sens, 2, quantile, probs = 0.025,na.rm=T)
# upper_bound_s <- apply(predicted_probabilities_sens, 2, quantile, probs = 0.975,na.rm=T)
# 
# acq_probs_s = data.frame(cbind(estimate_s,lower_bound_s,upper_bound_s)) 
# acq_probs_s
# 
# # Calculate the OR of the prediction probabilities
# data_pred = data.frame(predicted_probabilities)
# names(data_pred) = c("R1_control","R1_int","R2_control","R2_int","R3_control","R3_int")
# data_pred_s = data.frame(predicted_probabilities_sens)
# names(data_pred_s) = c("R1_control","R1_int","R2_control","R2_int","R3_control","R3_int")
# 
# OR_pred1 = data_pred$R1_int/data_pred$R1_control
# OR_pred2 = (data_pred$R2_int/data_pred$R2_control)
# OR_pred3 = (data_pred$R3_int/data_pred$R3_control)
# OR_pred1_s = data_pred_s$R1_int/data_pred_s$R1_contol
# OR_pred2_s = (data_pred_s$R2_int/data_pred_s$R2_control)
# OR_pred3_s = (data_pred_s$R3_int/data_pred_s$R3_control)
# 
# OR_estimate1 <-quantile(OR_pred1, probs = 0.50,na.rm=T)
# OR_lower_bound1 <- quantile(OR_pred1, probs = 0.025,na.rm=T)
# OR_upper_bound1 <- quantile(OR_pred1, probs = 0.975,na.rm=T)
# 
# OR_estimate1 <-quantile(OR_pred1, probs = 0.50,na.rm=T)
# OR_lower_bound1 <- quantile(OR_pred1, probs = 0.01,na.rm=T)
# OR_upper_bound1 <- quantile(OR_pred1, probs = 0.99,na.rm=T)
# 
# OR_estimate2 <-quantile(OR_pred2, probs = 0.50,na.rm=T)
# OR_lower_bound2 <- quantile(OR_pred2, probs = 0.01,na.rm=T)
# OR_upper_bound2 <- quantile(OR_pred2, probs = 0.99,na.rm=T)
# 
# OR_estimate3 <-quantile(OR_pred3, probs = 0.50,na.rm=T)
# OR_lower_bound3 <- quantile(OR_pred3, probs = 0.01,na.rm=T)
# OR_upper_bound3 <- quantile(OR_pred3, probs = 0.99,na.rm=T)
# 
# OR_probs1 = data.frame(cbind(OR_estimate1,OR_lower_bound1,OR_upper_bound1)) 
# OR_probs2 = data.frame(cbind(OR_estimate2,OR_lower_bound2,OR_upper_bound2)) 
# OR_probs3 = data.frame(cbind(OR_estimate3,OR_lower_bound3,OR_upper_bound3)) 
# names(OR_probs1) = c("estimate", "lower", "upper")
# names(OR_probs2) = c("estimate", "lower", "upper")
# names(OR_probs3) = c("estimate", "lower", "upper")
# 
# # Sensitivity analyses
# OR_estimate1_s <-quantile(OR_pred1_s, probs = 0.50,na.rm=T)
# OR_lower_bound1_s <- quantile(OR_pred1_s, probs = 0.025,na.rm=T)
# OR_upper_bound1_s <- quantile(OR_pred1_s, probs = 0.975,na.rm=T)
# 
# OR_estimate2_s <-quantile(OR_pred2_s, probs = 0.50,na.rm=T)
# OR_lower_bound2_s <- quantile(OR_pred2_s, probs = 0.025,na.rm=T)
# OR_upper_bound2_s <- quantile(OR_pred2_s, probs = 0.975,na.rm=T)
# 
# OR_estimate3_s <-quantile(OR_pred3_s, probs = 0.50,na.rm=T)
# OR_lower_bound3_s <- quantile(OR_pred3_s, probs = 0.025,na.rm=T)
# OR_upper_bound3_s <- quantile(OR_pred3_s, probs = 0.975,na.rm=T)
# 
# OR_probs1_s = data.frame(cbind(OR_estimate1_s,OR_lower_bound1_s,OR_upper_bound1_s)) 
# OR_probs2_s = data.frame(cbind(OR_estimate2_s,OR_lower_bound2_s,OR_upper_bound2_s)) 
# OR_probs3_s = data.frame(cbind(OR_estimate3_s,OR_lower_bound3_s,OR_upper_bound3_s)) 
# names(OR_probs1_s) = c("estimate", "lower", "upper")
# names(OR_probs2_s) = c("estimate", "lower", "upper")
# names(OR_probs3_s) = c("estimate", "lower", "upper")
# 
# 
# OR_probs = rbind(OR_probs2,OR_probs2_s,OR_probs3,OR_probs3_s)
# 
# row.names(OR_probs) = c("R2","R2 (sensitivity)", "R3", "R3 (sensitivity)")
# 
# 
# #acq_l = melt(acq_probs, values=row.names(acq_probs))
# 
# # Plot the acquisition probabilities with confidence intervals
# p_or_pred = ggplot(OR_probs, aes(y=row.names(OR_probs), x=estimate)) + geom_point(size=5)+
#   geom_vline(xintercept = 1, linetype=2, col="red")+
#   geom_linerange(aes(xmin = lower, xmax = upper), 
#                  position = position_dodge(width = 0.2), size=1, linetype=c(1,2,1,2))+
#   theme_bw()+theme(
#     axis.text.x=element_text(size=12),
#     axis.title=element_text(size=12,face="bold"),
#     axis.text=element_text(size=12,face="bold"),
#     strip.text=element_text(size=12,face="bold")) +
#   xlim(0,1.5) + 
#   labs(title = "Odds Ratio ESBL-acquisition probability",
#        x = "Ratio intervention(t)/control(t)",
#        y = "") +
#   scale_y_discrete(breaks = c("R2","R2 (sensitivity)", "R3", "R3 (sensitivity)"),  # Set custom breaks
#                    labels = c(">3M",">3M (sensitivity)", 
#                               ">9M",">9M (sensitivity)" ))  # Set custom labels
# p_or_pred
# 
# 
# pdf("./Output/Figures/fit_plotOR_ratio.pdf", width=7,height = 5)
# print(p_or_pred)
# dev.off()
# 
# ################################################################
# # RANDOM-EFFECTS MODEL - Decolonisation
# ################################################################
# 
# # Fit odds ratios for R2 acquisitions only
# md2_r2 <- glmer(decolonisation ~ intervention.text + age + sexe + (1|menage_id), 
#                 data=dfls2,  family = binomial(link = "logit"))
# summary(md2_r2)
# data.frame(exp(fixef(md2_r2)))
# exp(confint(md2_r2,method="Wald"))  # Effect intervention: OR = 0.99 (0.57 - 1.75); 
# # Fit odds ratios for R3 acquisitions only
# md2_r3 <- glmer(decolonisation ~ intervention.text + age + sexe + (1|menage_id), 
#                 data=dfls3,  family = binomial(link = "logit"))
# summary(md2_r3)
# data.frame(exp(fixef(md2_r3)))
# exp(confint(md2_r3,method="Wald"))  # Effect intervention: OR = 1.10 (0.70 - 1.74); 
# 
# 
# # Now fit both together
# mds2 = gamm(decolonisation ~ age + sexe + intervention.text*factor(time) + 
#               s(as.numeric(month), bs="cc"),
#             random=list(village_name=~1, menage_id_member=~1), 
#             family=binomial(link = "logit"), data=d)
# 
# summary(mds2$gam)
# exp(coef(mds2$gam))
# exp(confint(mds2$gam))
# 
# ################################################################
# # RANDOM-EFFECTS MODEL - Colonisation
# ################################################################
# dfls = dfls %>% mutate(
#   intervention.inplace= ifelse(time%in%c(2,3)&intervention.text=="intervention","yes","no")
# )
# 
# dfls0 = dfls0 %>% mutate(
#   intervention.inplace = ifelse(time%in%c(2,3)&intervention.text=="intervention","yes","no")
# )
# 
# table(dfls$time,dfls$intervention.inplace)
# 
# # Fit R0, R1, R2 and R3 together --> Then can calculate probability of infection at each time point
# m2 <- glmer(esble-1 ~ intervention.text*factor(time) + harmonic(month,1,12)  + age + sexe + (1|village_name/menage_id_member), data=dfls0,  family = binomial(link = "logit"))
# summary(m2)
# data.frame(exp(fixef(m2)))
# exp(confint(m2,method="Wald"))  # Effect intervention: OR = 
# 
# 
# ##############################################
# # MSM model
# ##############################################
# 
# # Define the transition intensity parameters
# lambda <- 0.02  # Transition intensity from state 1 to state 2 (acquisition)
# rho <- 1/(5*30)  # Transition intensity from state 2 to state 1 (decolonisation)
# 
# # Create the Q matrix
# Q <- matrix(0, nrow = 2, ncol = 2)
# Q[1, 1] <- 0
# Q[2, 2] <- 0
# Q[1, 2] <- lambda
# Q[2, 1] <- rho
# diag(Q) <- -rowSums(Q)
# 
# rownames(Q) <- colnames(Q) <- c("no_esbl", "esbl")
# Q
# crudeinits.msm(esble ~ days, menage_id_member, data=dfls0, qmatrix=Q)
# 
# incl = unique(dfls0$menage_id_member[dfls0$time==0])
# dfls0$after = ifelse(dfls0$time%inc(0,1), "no", "yes")
# du = dfls0 %>% filter(menage_id_member%in%incl) %>%
#   select(time,days,days2,month,rainy, intervention.text,intervention.start,village_name, menage_id, menage_id_member, age,sexe,esble, nesbl_tminus)
# du$intervention.start = ifelse(du$time%in%c(0), "control", du$intervention.text)
# table(du$time, du$intervention.start)
# 
# du = du %>% filter(!menage_id_member=="EHL00002504-1")
# 
# # add seasonal term
# harmonics <- data.frame(harmonic(du$month, 2, period = 12))
# names(harmonics) = c("harmonic_1","harmonic_2", "harmonic_3", "harmonic_4")
# du = cbind(du,harmonics)
# 
# 
# # Fit the multi-state model
# msm_fit <- msm(
#   formula = esble ~ days2,          # Model formula
#   covariates = list("1-2" = ~factor(intervention.start) + nesbl_tminus+age+sexe, 
#                     "2-1" = ~factor(intervention.start)),
#   subject = menage_id_member,    # ID variable
#   data = du,                     # Dataset
#   qmatrix = Q,                     # Transition intensity matrix
#   gen.inits = T,                   # initial values generated by the model
#   method = "BFGS"                  # Optimization method
# )
# 
# # Summarize the fitted model
# summary(msm_fit)
# 
# qmatrix.msm(msm_fit)
# pmatrix.msm(msm_fit, t=max(dfls0$days2))
# sojourn.msm(msm_fit)
# 
