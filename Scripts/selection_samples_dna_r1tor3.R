#####################################################
# DATA EXPLORATION
#####################################################
# This code is selecting households per cluster to follow up over time for sequencing.
# The code is selecting from households  which were included in round 0 (see selection_samples_dna.R)

# Considering we tested up to 6 individuals in the household, and with the identifiability of the model in mind, 
# we might, for r1 – r4 reconsider which households we include for sequencing, i.e. exclude household larger 
# than 10 household members and take and additional number of households compared to the ones we selected 
# in the first round (in total we can select 850 household members over 4 rounds).



# Date: 7 June 2024
# Author: Esther van Kleef
rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

# Load in data from individuals selected in round 0
############################################################

iddnaR0 = read.csv("./Data/BF/clean/bf_sample_ids_dna_update.csv", sep=";")
length(unique(iddnaR0$id_ecantillon))

# Load in data from all rounds
############################################################
# ALL INDIVIDUALS - Wide format
dfe = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_esble_r0123_wide.csv")
length(unique(df$menage_id_member)) # 1237

# ALL INDIVIDUALS - long format
dfl = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_all_r0123.csv")
length(unique(dfl$menage_id_member)) # 1237

# SELECTED INDIVIDUALS WITH FOUR OBSERVATIONS - long format
load("./Data/BF/clean/use_in_analyses/bf_esbl0123_long_completecases.rda")
length(unique(dfls0$menage_id_member)) # 785


# Select only those that were sampled in R0
###############################################################
dfdnaR0 = dfl %>% filter(household%in%iddnaR0$household & esble==1 & redcap_event_name=="round_0_arm_1") # 188
dfdna = dfl %>% filter(household%in%iddnaR0$household & esble==1) # 675

# how many households? 
length(unique(dfdna$household)) # 

# which missing?
iddnaR0$id_ecantillon[!iddnaR0$id_ecantillon%in%dfdna$id_ecantillon]
iddnaR0$household[!iddnaR0$household%in%dfdna$household] # two households missing

# Now we have 850 samples, so another 175 left for four rounds. As some might fail, we could add 5 % = 184 samples
# What is the household size of these households?
table(dfdna$n.householdmember)
nperhh = (dfdna %>% group_by(village, household) %>%
            summarise(n =n()))

summary(nperhh$n) # mean of 12 samples per household. 184/12 = 15 more households
length(unique(dfl$household)) # 261 households
length(unique(dfl$village)) # 22 villages

# household sizes
table(dfdna$n.householdmember[dfdna$redcap_event_name=="round_0_arm_1"])

hh.size = dfdna %>% group_by(household, n.householdmember) %>%
  summarise(n = n())
table(hh.size$n.householdmember) # 17/57 households that are larger than 10. So in these households, quite some individuals with unobserved status, as sampled up to 6 individuals

hh.size.r0 = dfdna %>% filter(redcap_event_name=="round_0_arm_1") %>% group_by(household, n.householdmember) %>%
  summarise(n = n())
table(hh.size.r0$n.householdmember) # 17/57 households that are larger than 10. So in these households, quite some individuals with unobserved status, as sampled up to 6 individuals

# Suggest to exclude the households of larger than 20 individuals for the next rounds
hh_larg = hh.size.r0$household[hh.size.r0$n.householdmember>19]

# Remove from selected dataset
dfdna_lrm = dfdna %>% filter(!household%in%hh_larg & !redcap_event_name == "round_0_arm_1")
dfdna_lrm = rbind(dfdna_lrm, dfdna %>% filter(redcap_event_name=="round_0_arm_1"))

# Now add new households
# 850 - 602 = 248 more samples. 248*1.05 = 260 more samples / 12 on average per hh = 22 additional households, i.e. one per village

# Create datasets with households where at least one individual was positive in R0 and not yet included
dfl_eligible = dfl %>% filter(!household %in% iddnaR0$household & n.householdmember<15 & redcap_event_name=="round_0_arm_1" & esble==1)
length(unique(dfl_eligible$household)) # 168 households

length(unique(dfl_eligible$village)) # all villages included

new_hh = NULL
for(i in unique(dfl_eligible$village)){
  hh = unique(dfl_eligible$household[dfl_eligible$village==i])
  hh_s = sample(hh,2)
  new_hh = c(new_hh,hh_s)
  print(i)
}
new_hh

# Update list
hh_ids_dna = c(iddnaR0$household, new_hh)
length(unique(hh_ids_dna))

dfdna_update = dfl %>% filter(household%in%hh_ids_dna & esble==1 & !redcap_event_name=="round_0_arm_1"& n.householdmember<15)
dfdna_update = rbind(dfdnaR0,dfdna_update)

# check hh sizes
hh.size.rall = dfdna_update %>% group_by(household, n.householdmember) %>%
  summarise(n = n())
table(hh.size.rall$n.householdmember) # 17/57 households that are larger than 10. So in these households, quite some individuals with unobserved status, as sampled up to 6 individuals

nperround = dfdna_update %>% group_by(menage_id_member) %>%
  summarise(n = n())

nperroundhh = dfdna_update %>% group_by(household, n.householdmember) %>%
  summarise(n = n())
table(nperroundhh$n)

# Export list of IDs
dfdna_update_sel = dfdna_update %>% select(redcap_event_name,household, id_ecantillon, menage_id_member)
write.xlsx(dfdna_update_sel, "./Data/BF/clean/bf_sample_ids_dna_update_R0123.xlsx")
