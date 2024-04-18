#####################################################
# DATA CLEANING AND LINKAGE BURKINA FASO
#####################################################
# This code is cleaning the data and linking the different datasets

# 22 January 2024
# Author: Esther van Kleef
# Last update: 10 February 2024

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

# Lab data
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_manualchange_no_password.xlsx"))
#car_r0_old = read_xlsx(paste0(DirectoryData,"/householdsurvey/CABUBPortageAsymptom_DATA_2023-10-17.xlsx"))
#car_r1_old = read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_2024-02-05_1720.xlsx"), sheet=1)
#car_r2_old = read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_2024-02-05_1720.xlsx"), sheet=2)

car_r0123 = read.csv(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_2024-04-17_1527.csv"), sep=";")

car_r0 = car_r0123 %>% filter(redcap_event_name %in% "round_0_arm_1")
car_r1 = car_r0123 %>% filter(redcap_event_name %in% "round_1_arm_1")
car_r2 = car_r0123 %>% filter(redcap_event_name %in% "round_2_arm_1")
car_r3 = car_r0123 %>% filter(redcap_event_name %in% "round_3_arm_1")

#View(car_r0[!car_r0$record_id %in% car_r0_old$record_id,])


# Lab ids vs household ids
hh_lab_ids =  readxl::read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))
names(hh_lab_ids)
names(car_r0)
names(hh_lab_ids) = c("household", "menage_id", "bras")

# Household data
#wash_r0_old = read_xlsx(paste0(DirectoryData, "/householdsurvey/WP4_WASH_07_09_2023.xlsx"))
#wash_r1_old = read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24.xlsx"), sheet=1) # These wash_r1 and wash_r2 are not householdsurvey data; but do provide denominator data for lab R1 and R2
#wash_r2_old = read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24.xlsx"), sheet=2)
wash_r0123 = readxl::read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBWP4_DATA_2024-04-17_1528.xlsx"))

wash_r0 = wash_r0123 
wash_r1 = wash_r0123 %>% filter(redcap_event_name %in% "round_1_arm_1")
wash_r2 = wash_r0123 %>% filter(redcap_event_name %in% "round_2_arm_1")
wash_r3 = wash_r0123 %>% filter(redcap_event_name %in% "round_3_arm_1")


# Villages (that are the clusters) of CABU-EICO
villages = readxl::read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","ajouter")


# Check variables between r0 and r1 and r2 - seems we have received the household variables
# of those with a lab sample, not the lab data
# which(names(car_r1) %in% names(car_r0))
# which(names(car_r1) %in% names(wash_r0))
# length(which(names(car_r1) %in% names(wash_r0)))
# which(!names(car_r1) %in% names(wash_r0))
# names(car_r1[which(names(car_r1) %in% names(car_r0))])
# names(car_r1[which(!names(car_r1) %in% names(wash_r0))])
# names(car_r1[which(!names(car_r1) %in% names(wash_r0))])

# Antibiotic use data
abx = read.csv("./Data/BF/clean/watch_acute.csv")


# Add variables village and household
# ROUND 0
####################
car_r0$village = substr(car_r0$record_id, start = 1, stop = 2)
car_r0$household = str_extract(car_r0$record_id, "[^-]+")
car_r0 = merge(car_r0, villages, by="village")

# ROUND 1
####################
car_r1$village = substr(car_r1$record_id, start = 1, stop = 2)
car_r1$household = str_extract(car_r1$record_id, "[^-]+")
car_r1 = merge(car_r1, villages, by="village")

# ROUND 2
####################
car_r2$village = substr(car_r2$record_id, start = 1, stop = 2)
car_r2$household = str_extract(car_r2$record_id, "[^-]+")
car_r2 = merge(car_r2, villages, by="village")


############################################################
# CLEAN LAB DATA
############################################################

# Clean germe; diameters are CLSI guidelines, 
# Jan jacobs comm: ESBL positive is defined as cetriax/cefo <=22 (including I and R, based on CLSI guidelines,
# See in WP6 folder "word file: Interpretation antibiogramme des isolats PORTAGE asymptomatique_ESBL_E. coliKlebsielle.docx), 
# So I understood to use this, instead of esbltest == 1 
# This is then interpreted as, cetriax/cefo resistant following ESBL selective medium)

# Evk 20 September 2023:
# Would be good to do a final check with Yougbare if this is correct or whether
# we can just use all the observations in the car_r0 file (i.e. esbltest == 0 and esbltest == 1). 

# ROUND 0
####################
car_r0 = car_r0 %>%
  mutate(germe_c = ifelse(germe %in% c("E.COLI", "E-COLI", "ECOLI", "E.COLI 2", "E.COLI 1", 
                                       "E-COLI 2","E-CLI","E-CLOI"),"e.coli", 
                          ifelse(germe %in% c("SALMONELLA SP","SALMONELLA SPP","SALMONELLA SSP",
                                              "SALMONELLA", "SALMO"),"salmonella",NA)),
         germe_c = ifelse(morphotyp%in%c(1, NA),germe_c, paste0(germe_c, "_", 2)),
         esbl_pos = ifelse(diametr_cetriax_or_cefota <= 22, 1, 0))

# Number of cases positive
table(car_r0$germe_c, car_r0$esbl_pos)
table(car_r0$esbl_pos) # These are the ESBL positive patients based on cetriax_or_cefota, 769
table(car_r0$testesbl) # These are the ESBL positive patients based on esbl_pos, 772
table(car_r0$esbl_pos==1 & car_r0$testesbl==1) # difference of 7; we decided to ignore these differences

names(car_r0)

# ROUND 1
####################
table(car_r1$germe)

car_r1 = car_r1 %>%
  mutate(germe_c = ifelse(germe %in% c("E.COLI", "E-COLI", "ECOLI", "E.COLI 2", "E.COLI 1", 
                                       "E-COLI 2","E-CLI","E-CLOI"),"e.coli", 
                          ifelse(germe %in% c("SALMONELLA SP","SALMONELLA SPP","SALMONELLA SSP", "SALMO SPP", "SALMONELLE SPP", "SELMO",
                                              "SALMONELLA", "SALMO"),"salmonella",NA)),
         germe_c = ifelse(morphotyp%in%c(1, NA),germe_c, paste0(germe_c, "_", 2)),
         esbl_pos = ifelse(diametr_cetriax_or_cefota <= 22, 1, 0))

# Number of cases positive
table(car_r1$germe_c, car_r1$germe)
table(car_r1$germe, car_r1$esbl_pos)

table(car_r1$germe_c, car_r1$esbl_pos)
table(car_r1$esbl_pos) # These are the ESBL positive patients based on cetriax_or_cefota, 653
table(car_r1$testesbl) # These are the ESBL positive patients based on esbl_pos, 653
table(car_r1$esbl_pos==1 & car_r1$testesbl==1) # difference of 2. We decided to ignore these differences


# ROUND 2
####################
table(car_r2$germe, useNA = "always")

car_r2 = car_r2 %>%
  mutate(germe_c = ifelse(germe %in% c("E.COLI", "E-COLI", "ECOLI", "E.COLI 2", "E.COLI 1", "E.COLIE",
                                       "E-COLI 2","E-CLI","E-CLOI"),"e.coli", 
                          ifelse(germe %in% c("SALMONELLA SP","SALMONELLA SPP","SALMONELLA SSP",
                                              "SALMONELLA", "SALMO"),"salmonella",NA)),
         germe_c = ifelse(morphotyp%in%c(1, NA),germe_c, paste0(germe_c, "_", 2)),
         esbl_pos = ifelse(diametr_cetriax_or_cefota <= 22, 1, 0))

# Number of cases positive
table(car_r2$germe_c, car_r2$germe)

table(car_r2$germe_c, car_r2$esbl_pos)
table(car_r2$esbl_pos) # These are the ESBL positive patients based on cetriax_or_cefota, 653
table(car_r2$testesbl) # These are the ESBL positive patients based on esbl_pos, 653

table(car_r2$esbl_pos==1 & car_r2$testesbl==1) # difference of 17. We decided to ignore these differences

#################################################################
# CLEAN HOUSEHOLD DATA
#################################################################
# Good to note is that the household data is collected at household head level
# i.e. there is one survey taken by the household head. 

# Then there is a healthcare utilisation survey attached to the household survey
# Here for each household member, the household head is asked about the number of 
# healthcare visits in the last 30 days
# Each of these visits (per provider type) are one line of data
# Therefore there are multiple observations within a household

# ROUND 0
################################
wash_r0 = wash_r0 %>% mutate(
  dob = as.Date(dob, format = "%Y-%m-%d"),
  age = tolower(age),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d")
 # ageyears = gsub(" ", "", age), # remove spaces
#  ageyears = as.numeric(gsub("ans", "", age, fixed = TRUE))
) %>% select(-dob)

#wash_r0$ageyears[grepl("mois", wash_r0$age)==T] <- as.numeric(gsub("mois", "", wash_r0$age[grepl("mois", wash_r0$age)==T], fixed = TRUE))
#wash_r0$ageyears[grepl("mois", wash_r0$age)] <- wash_r0$ageyears[grepl("mois", wash_r0$age)] / 12
#wash_r0$ageyears[is.na(wash_r0$age) & !is.na(wash_r0$dob)] <- round((as.Date("2023-01-30") - wash_r0$dob[is.na(wash_r0$age) & !is.na(wash_r0$dob)])/365.25, 0)
  
# ROUND 1
################################
  
wash_r1 = wash_r1 %>% mutate(
  dob = as.Date(dob, format = "%Y-%m-%d"),
  age = tolower(age),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d")
  # ageyears = gsub(" ", "", age), # remove spaces
  #  ageyears = as.numeric(gsub("ans", "", age, fixed = TRUE))
) %>% select(-dob)


# ROUND 2
################################
wash_r2 = wash_r2 %>% mutate(
  dob = as.Date(dob, format = "%Y-%m-%d"),
  age = tolower(age),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d")
  # ageyears = gsub(" ", "", age), # remove spaces
  #  ageyears = as.numeric(gsub("ans", "", age, fixed = TRUE))
) %>% select(-dob)


#################################################################
# LINK LAB AND HOUSEHOLD DATA
#################################################################

# Create household variable in lab data (car_r0) to link to WASH survey by adding zero's.
###################################################################################
# Don't need the below section as we have now the Correspondande-Code_Lab-ID_Menage.xlsx file
###################################################################################
# df = data.frame(household = car_r0$household) 
# df = df %>%  separate(household,
#                       into = c("text", "num"),
#                       sep = "(?<=[A-Za-z])(?=[0-9])")
# 
# df$num_ad = NULL
# df$household = car_r0$household
# 
# for(i in 1:length(df$num)){
#   if(nchar(df$num)[i]==3){
#     p = "00000"
#     df$num_ad[i] = paste0(p,df$num[i])
#   }else if(nchar(df$num)[i]==4){
#     p = "0000"
#     df$num_ad[i] = paste0(p,df$num[i])
#   }else if(nchar(df$num)[i]==5){
#     p = "000"
#     df$num_ad[i] = paste0(p,df$num[i])
#   }
#   else if(nchar(df$num)[i]==6){
#     p = "00"
#     df$num_ad[i] = paste0(p,df$num[i])
#   }
# }
# 
# #car_r0 = left_join(car_r0, hh_lab_ids)
# 
# #This is not yet doing the trick fully as some menage_id have no zero's (see nchar == 8 for some id's)
# # Also some nchar == 12 for some in the new df$menage_id which should be 11
# df$menage_id = paste0(df$text,df$num_ad)
# nchar(df$menage_id)
# nchar(wash_r0$menage_id)

#car_r0$menage_id = df$menage_id

# Need to adjust the 28 IDs in the household dataset that have 0's lacking
# ids_m = as.table(cbind(unique(car_r0$household[car_r0$found_in_wash==0]),
#              unique(car_r0$menage_id[car_r0$found_in_wash==0])))
# # 11 will be resolved through this method
# which(ids_m[,1] %in% unique(wash_r0$menage_id))
# ids_m[which(ids_m[,1] %in% unique(wash_r0$menage_id))] # All "CRAS"
# unique(df$text[car_r0$found_in_wash==0]) 

# Need to check this (26 January 2024)
# p = 1
# for(i in ids_m[,1]){
#   repl =  wash_r0$menage_id[wash_r0$menage_id==i] 
#   wash_r0$menage_id[wash_r0$menage_id==i] = rep(ids_m[p,2], length(repl))
#   p = p+1
# }
# 
# # which left?
# car_r0$found_in_wash[which(car_r0$menage_id %in% wash_r0$menage_id)] = 1
# car_r0$found_in_wash[is.na(car_r0$found_in_wash)] = 0
# length(unique(car_r0$household[car_r0$found_in_wash==0])) # 17 households can not be found in WASH survey database; could be due to error's in ID or due to zero's that need to be removed, needs checking
# unique(car_r0$menage_id[car_r0$found_in_wash==0])
# unique(car_r0$household[car_r0$found_in_wash==0])
###################################################################################
# Code can be disgarded up till here
###################################################################################

# ROUND 0 (pre-intervention round)
###################################################################################
# Link village and cluster data to household survey
wash_r0$village = substr(wash_r0$menage_id, start = 1, stop = 2)
wash_r0 = merge(wash_r0, villages, by="village")


# Link lab data with lab vs hh survey ids (as IDs are in different format)
car_r0 = left_join(car_r0, hh_lab_ids, by="household")
# Check if all linked
table(car_r0$bras, useNA= "always")

# Link lab data with hh survey ids
length(unique(car_r0$menage_id))
length(unique(wash_r0$menage_id))

# Can all IDs be traced back from lab to wash survey?
car_r0$found_in_wash[which(car_r0$menage_id %in% wash_r0$menage_id)] = 1
car_r0$found_in_wash[is.na(car_r0$found_in_wash)] = 0
length(unique(car_r0$household[car_r0$found_in_wash==0])) # after cleaning, all households can be found in WASH survey database; 
unique(car_r0$menage_id[car_r0$found_in_wash==0])
unique(car_r0$household[car_r0$found_in_wash==0])

# HH ids "ETE00000201" "SBA00001601" "SJG00002501" not found in hh survey
# "ETE00201" "SDC04101" "SJG02501" --> 
# In WASH survey "ETE00101" can be found, typo?
# In hh_labs_ids, SDC04101 corresponds to SBA00001601, should be SDC00004101?
# In WASH survey "SJG02401" can be found, typo?


wash_r0$data_row = c(1:nrow(wash_r0)) # This variable we can use for identifying back those individuals that had a sample taken

# Make dataset with only those household individuals that had a stool sample taken and their individual variables
wash_r0_lab = wash_r0 %>% filter(!is.na(cs_id_individu)) %>% # ensure just 1 observation per person of whom individual is esbl positive
  select(data_row, cs_id_individu,num_echantillon, menage_id, village, age, sexe, date_consentement, date_recuperation_selle)
  


# Merge wash_r0 patient characteristics with lab data - NEED TO GET IDs IIN THE SAME FORMAT
which(car_r0$record_id %in% unique(wash_r0_lab$cs_id_individu)) # Have to change the format of both to make sure matching can be done
head(car_r0$record_id ); head(wash_r0_lab$cs_id_individu)
unique(car_r0$record_id)
unique(wash_r0$cs_id_individu)

# SEE HOW TO DO THE LINKAGE
# Wash database does no have "-" and puts and "M" before each household member. Also sometimes 01 and sometimes 1 as method of writing
# PROBABLY IF WE TAKE FROM wash_r0 cs_id_individue all "MX" (so M + number after ID), then put that in a seperate column.
# Then we take the household number from the WASH survey (so menage_id), combine these to with a hyphen

# THEN for car_r0 we take also menage_Id and digits after the "-", remove the 0's then two can be combined

# Change IDs to the same format

# Create a variable in which we will safe the new formatted and cleaned ID
wash_r0_lab$menage_id_member = NA

# Create a variable that will just store the household member number
wash_r0_lab$member = NA

# Extract the number after M to get household number
wash_r0_lab$member =  gsub(".*M", "", wash_r0_lab$cs_id_individu)
# Remove leading '0's
wash_r0_lab$member = as.character(as.numeric(wash_r0_lab$member))
# Check the one's which are now NA
table(is.na(wash_r0_lab$member))
wash_r0_lab$cs_id_individu[is.na(wash_r0_lab$member)] # 22 individuals have a household member number missing
#View(wash_r0_lab[is.na(wash_r0_lab$member),])
wash_r0_lab$num_echantillon[is.na(wash_r0_lab$member)] # of 0 we can get them from num_echantillon
#wash_r0_lab$member[wash_r0_lab$cs_id_individu=="CRAS1003"] = "1" not needed anymore after cleaned datasets 17/4/2024
#wash_r0_lab$member[wash_r0_lab$cs_id_individu=="SCA00101"] = "3" not needed anymore after cleaned datasets 17/4/2024
table(is.na(wash_r0_lab$member))


# Now create new variable for linking with lab dataset
wash_r0_lab$menage_id_member = paste0(wash_r0_lab$menage_id, "-", wash_r0_lab$member)
wash_r0_lab$menage_id_member
# Make NAs for the one's that still need checking
wash_r0_lab$menage_id_member[is.na(wash_r0_lab$member)] = NA
table(is.na(wash_r0_lab$menage_id_member))


# Check if dubplicates from wash_r0
table(duplicated(wash_r0_lab$menage_id_member)) # 21 duplicates
dups = unique(wash_r0_lab$menage_id_member[duplicated(wash_r0_lab$menage_id_member)])
#View(wash_r0_lab[which(wash_r0_lab$menage_id_member%in%dups),])
# For   "SCB00002402-7", "SCB00002402-8"  "SCB00002402-9"  "SCB00002402-10" "SCB00002402-11" "SCB00002402-12"
# THese are entered double, so keep the ones without NA for age
#nkeep = which(wash_r0_lab$menage_id_member%in%c("SCB00002402-7", "SCB00002402-8",  "SCB00002402-9" ,
#                                                "SCB00002402-10" ,"SCB00002402-11", "SCB00002402-12") & is.na(wash_r0_lab$age)) 

#wash_r0_lab_d = wash_r0_lab[-c(nkeep),]

# Then for remainder keep first record for now
wash_r0_lab_de = wash_r0_lab %>% filter(!duplicated(menage_id_member)) # STILL NEED TO DO CHECKING FOR THOSE THAT ARE DUPLICATES BUT KEEP LIKE THIS FOR NOW
table(duplicated(wash_r0_lab_de$menage_id_member)) # 




# Now create the same variable for the lab dataset
# Create a variable in which we will safe the new formatted and cleaned ID
car_r0$menage_id_member = NA

# Create a variable that will just store the household member number
car_r0$member = NA

# Extract the number after M to get household number
car_r0$member =  gsub(".*-", "", car_r0$record_id)
# Remove leading '0's
car_r0$member = as.character(as.numeric(car_r0$member))
# Check the one's which are now NA
table(is.na(car_r0$member))
car_r0$record_id[is.na(car_r0$member)] # 16 individuals have a household member number after M
car_r0$member[is.na(car_r0$member)] =  gsub(".*M", "", car_r0$record_id[is.na(car_r0$member)])
table(is.na(car_r0$member))

# Now create new variable for linking with lab dataset
car_r0$menage_id_member = paste0(car_r0$menage_id, "-", car_r0$member)
car_r0$menage_id_member


# Remove salmonellas
 humanR0 = car_r0 %>% filter(germe_c != "salmonella")

# Keep for those with two e. coli's only ESBL positive
# First make dataset with only those with more than one e coli
dups = humanR0$record_id[duplicated(humanR0$record_id)] # 97
dups_members=humanR0$menage_id_member[duplicated(humanR0$menage_id_member)] # 98

d = humanR0[humanR0$menage_id_member%in% dups_members,]
# Keep only ESBL positives
d_dup = d %>% filter(esbl_pos==1)

length(d_dup$menage_id_member[duplicated(d_dup$menage_id_member)]) # 97 individuals left
# Keep only one record of those as we only need to know IF they were ESBL positive, not how many strains
d_dup = d_dup%>%filter(!duplicated(menage_id_member))

# Now merge with those that only had 1 E. coli
d = humanR0[!humanR0$menage_id_member%in% dups_members,]
humanR0e = rbind(d_dup,d)
length(unique(humanR0e$record_id)) # one missing needs checking still
length(unique(humanR0$record_id)) # all household samples are included
table(!duplicated(humanR0e$record_id))

# Select relevant variables from lab data
names(humanR0e)

humanR0e_sel = humanR0e %>% select(village,village_name,menage_id, intervention_text, household,menage_id_member,record_id,
                                   id_ecantillon,germe_c, date,esbl_pos)

# MERGE wash_r0_lab with car_r0
HR0e = left_join(humanR0e_sel,wash_r0_lab_de, by= c("menage_id_member", "menage_id", "village"),  suffix = c("", ""))
table(HR0e$esbl_pos)
table(!duplicated(HR0e$record_id))


# MERGE OTHER WAY AROUND FOR DENOMINATOR DATA
HR0 = left_join(wash_r0_lab_de,humanR0e_sel, by= c("menage_id_member","menage_id","village"),  suffix = c("", ""))
table(HR0$esbl_pos, useNA="always") # We have 25 less indiivduals with an ESBL which need checking
names(HR0)

HR0 = HR0 %>%
  mutate(esbl_pos=ifelse(is.na(esbl_pos),0,esbl_pos)) # FOR NOW ASSUME ALSO 30 INDIVIDUALS WITH TYPE IN ID_MENAGE_MEMBER  ARE NEGATIVE

# ARE ALL MERGED? 
no_m = which(!car_r0$menage_id_member %in%HR0$menage_id_member) # THESE ONES CAN NOT BE MATCHED AND NEED CHECKING
length(no_m)
car_r0$menage_id_member[no_m] # WE LOSE 11 individuals and need checking


# Select relevant variables from the household survey
################################################################################################

# variables excluding healthcare seeking behaviour survey questions (and related medicine use); as these
# are not 1 observation per household
wash_r0_sel = wash_r0 %>% select(data_row,menage_id,village, village_name, intervention_text,   
                                 redcap_event_name,
                                 date_enquete,groupe,nmbre_personne_menage, nbre_enf_0_5ans,
                                 nbre_menage_conc,
                                 informations_gnrales_complete, q1_diarrhee_prevenu___1,         
                                 q1_diarrhee_prevenu___2, q1_diarrhee_prevenu___3,         
                                 q1_diarrhee_prevenu___4,q1_diarrhee_prevenu___5,         
                                 q1_diarrhee_prevenu___6,q1_diarrhee_prevenu___7,         
                                 q1_diarrhee_prevenu___8,q1_diarrhee_prevenu___9,         
                                 autr_mesur_prev_diarrhe,q2_source_princ_saison_seche,  
                                 q3_source_princ_saison_pluv,q4_bidon_stock,                 
                                 q5a_bidon_ferme_rempli,q5b_bidon_ferme_vide,           
                                 q5c_bidon_nettoye,q6_traite_eau,               
                                 q6_autre_traitmen_eau,q7_type_inst_sanitaire,          
                                 q7_autr_typ_ins_sanitair,q8_autr_lieu_defecation___1,     
                                 q8_autr_lieu_defecation___2,q8_autr_lieu_defecation___3,     
                                 q8_autr_lieu_defecation___4,q8_autr_lieu_defecation___5,     
                                 q8_autr_lieu_defecation___6,q8_autr_lieu_defecation___7,     
                                 q8_autre_preciser,  q9_toilette_partagee,            
                                 q10_combien_partag,q11_dernier_nettoyage, 
                                 q12_elimine_selle_enf, q12_autre_preciser,
                                 q13_vidange_toilette,q13_autre_preciser,            
                                 q14_produit_lavag_main,q15_lave_apr_defec,              
                                 q16_lave_apr_repas,q17_animaux_menage,              
                                 q18_animaux_interieur___1,q18_animaux_interieur___2,       
                                 q18_animaux_interieur___3,q18_animaux_interieur___4,       
                                 q18_animaux_interieur___5,q18_animaux_interieur___6,       
                                 q18_autre_specifie, q19_animaux_dehors___1,         
                                 q19_animaux_dehors___2,q19_animaux_dehors___3,          
                                 q19_animaux_dehors___4,q19_animaux_dehors___5,          
                                 q19_animaux_dehors___6,q19_autre_specifie,            
                                 q20_excrement_animaux,q21_animal_malade___1,           
                                 q21_animal_malade___2,q21_animal_malade___3,           
                                 q21_animal_malade___4,q21_animal_malade___5,           
                                 q21_animal_malade___6, eau_assainissement_hygine_complete) %>%
  filter(!is.na(date_enquete)) # Denominator data (i.e. people tested for esbl) for R0


# FINAL DATASET FOR ANALYSES WITH ALL HOUSEHOLD VARIABLES
HR0_final = left_join(HR0,wash_r0_sel, by="menage_id", suffix = c("", ""))
names(HR0_final)
table(is.na(HR0_final$menage_id))


# ROUND 1
#############################################################################
wash_r1$village = substr(wash_r1$menage_id, start = 1, stop = 2)
wash_r1 = merge(wash_r1, villages, by="village")


# Link lab data with lab vs hh survey ids (as IDs are in different format)
car_r1 = left_join(car_r1, hh_lab_ids, by="household")
#View(car_r1[is.na(car_r1$menage_id),]) # These are not part of the hh_lab_ids ref check of r0, thus new households with an ESBL Or 
                                       # For ECCMID, if we want to still submit something, we may want to ignore these IDs for now


# Check if all linked
table(car_r1$bras, useNA= "always") # 0 can not be linked

# Link lab data with hh survey ids
length(unique(car_r1$menage_id))
length(unique(wash_r1$menage_id))

wash_r1$data_row = c(1:nrow(wash_r1)) # This variable we can use for identifying back those individuals that had a sample taken

# Make dataset with only those household individuals that had a stool sample taken and their individual variables
wash_r1_lab = wash_r1 %>% filter(!is.na(cs_id_individu)) %>% # ensure just 1 observation per person of whom individual is esbl positive
  select(data_row, cs_id_individu,num_echantillon, menage_id, village, age, sexe, date_consentement, date_recuperation_selle)


# Create a variable in which we will safe the new formatted and cleaned ID
wash_r1_lab$menage_id_member = NA

# Create a variable that will just store the household member number
wash_r1_lab$member = NA

# Extract the number after M to get household number
wash_r1_lab$member =  gsub(".*M", "", wash_r1_lab$cs_id_individu)
# Remove leading '0's
wash_r1_lab$member = as.character(as.numeric(wash_r1_lab$member))
# Check the one's which are now NA
table(is.na(wash_r1_lab$member))
wash_r1_lab$cs_id_individu[is.na(wash_r1_lab$member)] # 0 missing
wash_r1_lab$num_echantillon[is.na(wash_r1_lab$member)] # of 0 we can get them from num_echantillon
#wash_r1_lab$member[wash_r1_lab$num_echantillon=="SEA06902M2"] = "2"
table(is.na(wash_r1_lab$member))


# Now create new variable for linking with lab dataset
wash_r1_lab$menage_id_member = paste0(wash_r1_lab$menage_id, "-", wash_r1_lab$member)
wash_r1_lab$menage_id_member
# Make NAs for the one's that still need checking
wash_r1_lab$menage_id_member[is.na(wash_r1_lab$member)] = NA
table(is.na(wash_r1_lab$menage_id_member))


# Check if dubplicates from wash_r0
table(duplicated(wash_r1_lab$menage_id_member)) # 4 duplicates
dups = unique(wash_r1_lab$menage_id_member[duplicated(wash_r1_lab$menage_id_member)])
#View(wash_r1_lab[which(wash_r1_lab$menage_id_member%in%dups),])

# FOR NOW KEEP first record BUT NEEDS CHECKING!!!
wash_r1_lab_de = wash_r1_lab %>% filter(!duplicated(menage_id_member)) # STILL NEED TO DO CHECKING FOR THOSE THAT ARE DUPLICATES BUT KEEP LIKE THIS FOR NOW
table(duplicated(wash_r1_lab_de$menage_id_member)) # 

######################


# Now create the same variable for the lab dataset
# Create a variable in which we will safe the new formatted and cleaned ID
car_r1$menage_id_member = NA

# Create a variable that will just store the household member number
car_r1$member = NA

# Extract the number after M to get household number
car_r1$member =  gsub(".*-", "", car_r1$record_id)

# Remove leading '0's
car_r1$member = as.character(as.numeric(car_r1$member))
# Check the one's which are now NA
table(is.na(car_r1$member))
car_r1$record_id[is.na(car_r1$member)] # 9 individuals have a household member number after M
car_r1$member[is.na(car_r1$member)] =  gsub(".*M", "", car_r1$record_id[is.na(car_r1$member)])
table(is.na(car_r1$member))

# Now create new variable for linking with lab dataset
car_r1$menage_id_member = paste0(car_r1$menage_id, "-", car_r1$member)
car_r1$menage_id_member

# Remove salmonellas
humanR1 = car_r1 %>% filter(germe_c != "salmonella")

# Keep for those with two e. coli's only ESBL positive
# First make dataset with only those with more than one e coli
dups = humanR1$record_id[duplicated(humanR1$record_id)] # 71
length(dups)
dups_members=humanR1$menage_id_member[duplicated(humanR1$menage_id_member)] # 71

d = humanR1[humanR1$menage_id_member%in% dups_members,]
# Keep only ESBL positives
d_dup = d %>% filter(esbl_pos==1)

length(d_dup$menage_id_member[duplicated(d_dup$menage_id_member)]) # 70 individuals left
# Keep only one record of those as we only need to know IF they were ESBL positive, not how many strains
d_dup = d_dup%>%filter(!duplicated(menage_id_member))

# Now merge with those that only had 1 E. coli
d = humanR1[!humanR1$menage_id_member%in% dups_members,]
humanR1e = rbind(d_dup,d)
length(unique(humanR1e$record_id)) # one missing needs checking still
length(unique(humanR1$record_id)) # all household samples are included
table(!duplicated(humanR1e$record_id))


# Select relevant variables from lab data
humanR1e_sel = humanR1e %>% select(village,village_name,menage_id, intervention_text, household,menage_id_member,record_id,
                                   id_ecantillon,germe_c, date,esbl_pos)

# ROUND 2
#############################################################################

wash_r2$village = substr(wash_r2$menage_id, start = 1, stop = 2)
wash_r2 = merge(wash_r2, villages, by="village")


# Link lab data with lab vs hh survey ids (as IDs are in different format)
car_r2 = left_join(car_r2, hh_lab_ids, by="household")
#View(car_r2[is.na(car_r2$menage_id),]) # These are not part of the hh_lab_ids ref check of r0, thus new households with an ESBL Or 
                                       # For ECCMID, if we want to still submit something, we may want to ignore these IDs for now

# Check if all linked
table(car_r2$bras, useNA= "always") # 0 can not be linked

# Link lab data with hh survey ids
length(unique(car_r2$menage_id))
length(unique(wash_r2$menage_id))

wash_r2$data_row = c(1:nrow(wash_r2)) # This variable we can use for identifying back those individuals that had a sample taken

# Make dataset with only those household individuals that had a stool sample taken and their individual variables
wash_r2_lab = wash_r2 %>% filter(!is.na(cs_id_individu)) %>% # ensure just 1 observation per person of whom individual is esbl positive
  select(data_row, cs_id_individu,num_echantillon, menage_id, village, age, sexe, date_consentement, date_recuperation_selle)


# Create a variable in which we will safe the new formatted and cleaned ID
wash_r2_lab$menage_id_member = NA

# Create a variable that will just store the household member number
wash_r2_lab$member = NA

# Extract the number after M to get household number
wash_r2_lab$member =  gsub(".*M", "", wash_r2_lab$cs_id_individu)
# Remove leading '0's
wash_r2_lab$member = as.character(as.numeric(wash_r2_lab$member))
# Check the one's which are now NA
table(is.na(wash_r2_lab$member))
wash_r2_lab$cs_id_individu[is.na(wash_r2_lab$member)] # 2 missing
wash_r2_lab$num_echantillon[is.na(wash_r2_lab$member)] # can not get it from num_echanillon
table(is.na(wash_r2_lab$member))


# Now create new variable for linking with lab dataset
wash_r2_lab$menage_id_member = paste0(wash_r2_lab$menage_id, "-", wash_r2_lab$member)
wash_r2_lab$menage_id_member
# Make NAs for the one's that still need checking
wash_r2_lab$menage_id_member[is.na(wash_r2_lab$member)] = NA
table(is.na(wash_r2_lab$menage_id_member))


# Check if dubplicates from wash_r2
table(duplicated(wash_r2_lab$menage_id_member)) # 4 duplicates
dups = unique(wash_r2_lab$menage_id_member[duplicated(wash_r2_lab$menage_id_member)])
#View(wash_r2_lab[which(wash_r2_lab$menage_id_member%in%dups),])

# FOR NOW KEEP first record BUT NEEDS CHECKING!!!
wash_r2_lab_de = wash_r2_lab %>% filter(!duplicated(menage_id_member)) # STILL NEED TO DO CHECKING FOR THOSE THAT ARE DUPLICATES BUT KEEP LIKE THIS FOR NOW
table(duplicated(wash_r2_lab_de$menage_id_member)) # 


# Now create the same variable for the lab dataset
# Create a variable in which we will safe the new formatted and cleaned ID
car_r2$menage_id_member = NA

# Create a variable that will just store the household member number
car_r2$member = NA

# Extract the number after M to get household number
car_r2$member =  gsub(".*-", "", car_r2$record_id)

# Remove leading '0's
car_r2$member = as.character(as.numeric(car_r2$member))

# Check the one's which are now NA
table(is.na(car_r2$member))
car_r2$record_id[is.na(car_r2$member)] # 13 individuals have a household member number after M
car_r2$member[is.na(car_r2$member)] =  gsub(".*M", "", car_r2$record_id[is.na(car_r2$member)])
table(is.na(car_r2$member))

# Now create new variable for linking with lab dataset
car_r2$menage_id_member = paste0(car_r2$menage_id, "-", car_r2$member)
car_r2$menage_id_member


# Remove salmonellas
humanR2 = car_r2 %>% filter(germe_c != "salmonella")

# Keep for those with two e. coli's only ESBL positive
# First make dataset with only those with more than one e coli
dups = humanR2$record_id[duplicated(humanR2$record_id)] # 78
dups_members=humanR2$menage_id_member[duplicated(humanR2$menage_id_member)] # 80

d = humanR2[humanR2$menage_id_member%in% dups_members,]
# Keep only ESBL positives
d_dup = d %>% filter(esbl_pos==1)

length(d_dup$menage_id_member[duplicated(d_dup$menage_id_member)]) # 79 individuals left
# Keep only one record of those as we only need to know IF they were ESBL positive, not how many strains
d_dup = d_dup%>%filter(!duplicated(menage_id_member))

# Now merge with those that only had 1 E. coli
d = humanR2[!humanR2$menage_id_member%in% dups_members,]
humanR2e = rbind(d_dup,d)
length(unique(humanR2e$record_id)) # one missing needs checking still
length(unique(humanR2$record_id)) # 2 households are not inlcuded and need checking
table(!duplicated(humanR2e$record_id))

# Select relevant variables from lab data
humanR2e_sel = humanR2e %>% select(village,village_name,menage_id, intervention_text, household,menage_id_member,record_id,
                                   id_ecantillon,germe_c, date,esbl_pos)



# LINK ROUNDS TOGETHER
##########################################################################

# DO A CHECK IF ALL VARIABLES CAN BE LINKED 
not_in_dataset0 <- anti_join(car_r0, wash_r0_lab, by = "menage_id_member")# 11 can not be found back in wash_r0; these are the ECC1801 which should be ECC01101 and are duplicates
not_in_dataset1 <- anti_join(car_r1, wash_r0_lab, by = "menage_id_member")# 1 can not be found back in wash_r0
not_in_dataset2 <- anti_join(car_r2, wash_r0_lab, by = "menage_id_member")# 0 can not be found back in wash_r0

length(unique(not_in_dataset0$menage_id)) # belonging to 1 different households
length(unique(not_in_dataset1$menage_id)) # belonging to 1 different households
length(unique(not_in_dataset2$menage_id)) # belonging to 0 different households

# export those that can not be found back
not_in_dataset0 = not_in_dataset0 %>% select(record_id, id_ecantillon,household,menage_id,menage_id_member,member)
not_in_dataset1 = not_in_dataset1 %>% select(record_id, id_ecantillon,household,menage_id,menage_id_member,member)
not_in_dataset2 = not_in_dataset2 %>% select(record_id, id_ecantillon,household,menage_id,menage_id_member,member)

write_xlsx(not_in_dataset0, paste0(DirectoryDataOut, "/need_checking/unlinked_car0_wash0.xlsx"))
write_xlsx(not_in_dataset1, paste0(DirectoryDataOut, "/need_checking/unlinked_car1_wash0.xlsx"))
write_xlsx(not_in_dataset2, paste0(DirectoryDataOut, "/need_checking/unlinked_car2_wash0.xlsx"))

# Check which ones can not be found back and why
print(cbind(not_in_dataset0$record_id, not_in_dataset0$menage_id_member)) # The translation to new variable worked well, so it is not that
print(cbind(not_in_dataset1$record_id, not_in_dataset1$menage_id_member)) # There are a few that generate NA's which need updating of lab IDs
# "CAG00701-01" --> in wash_r0 no household with that number, should be CAG00702?
# "CRAS01002-01" -->  in wash_r0 no household with that number, should be CRAS0101 OR CRAS0201?
# "CRAS01002-02" -->  in wash_r0 no household with that number, should be CRAS0101 OR CRAS0201?
# "CRAS01002-03" -->  in wash_r0 no household with that number, should be CRAS0101 OR CRAS0201?
# "ECC01801-5"   -->  in wash_r0 no household with that number, should be ECC08801?
print(cbind(not_in_dataset2$record_id, not_in_dataset2$menage_id_member))
# CRAS00701-06"  --> 
# "CRAS00901-05"           
# "CRAS00901-3"            
# "CRAS00901-3"          
# "EKE07110-01"           
# "EKE07110-02"           
# "EKE07110-03"           
# "EKE07110-05"   

not_in_dataset0 <- anti_join(car_r0, wash_r0_lab, by = "menage_id")# 0 household IDs can not be found back in wash_r0
not_in_dataset1 <- anti_join(car_r1, wash_r0_lab, by = "menage_id")# 0 can not be found back in wash_r0
not_in_dataset2 <- anti_join(car_r2, wash_r0_lab, by = "menage_id")# 0 can not be found back in wash_r0

# CLEAN UP AT LEAST FOR car_r0 ids to see if can be linked
#wash_r0_ch = wash_r0_lab%>% select(data_row,menage_id, num_echantillon, cs_id_individu, menage_id_member)
#wash_r0_lab$menage_id_member[wash_r0_lab$data_row==177] = "CRAS0602-3"

# HOUSEHOLD MEMBER NUMBER NOT IN wash_r0
# car_r0$record_id %in% c("EAA04801-08", "ECA00502-3", "ECC05301-09","ECC05301-8","ECC05301-10","ECC05301-11",
#                         "EDI02001-7","EEA01901-02",) # Not finalised


# FOR NOW (ECCMID ANALYSES); CREATE A DATASET WITH ONLY HOUSEHOLDS THAT CAN BE TRACED BACK IN ALL DATABASES AND JUST THE E. COLI's
##############################################################################################################

# Create datasets that just include the individuals that have three observations
# Create a variable that says whether individuals was tested in following rounds or not

w2 = anti_join(wash_r2_lab_de, wash_r0_lab_de, by = "menage_id_member")# 0 household IDs can not be found back in wash_r0
w1 = anti_join(wash_r1_lab_de, wash_r0_lab_de, by = "menage_id_member")# 0 can not be found back in wash_r0

#w21 = anti_join(wash_r2_lab_de, wash_r1_lab_de, by = "menage_id_member")# 134 household IDs can not be found back in wash_r1
#w12 = anti_join(wash_r1_lab_de, wash_r2_lab_de, by = "menage_id_member")# 135 household IDs can not be found back in wash_r2

e2 = anti_join(humanR1e_sel, humanR0, by = "menage_id_member")# 292 household IDs can not be found back in car_r0
e1 = anti_join(humanR2e_sel, humanR0, by = "menage_id_member")# 235 household IDs can not be found back in car_r0

h2 = anti_join(wash_r2_lab_de, HR0, by = "menage_id_member")# 0 household IDs can not be found back in HR0 should actually be the same as wash_r0
h1 = anti_join(wash_r1_lab_de, HR0, by = "menage_id_member")# 0 can not be found back in HR0


ids_r1 = unique(wash_r1_lab_de$menage_id_member)
ids_r2 = unique(wash_r2_lab_de$menage_id_member)

data_eccmid = HR0 %>% 
  mutate(r1_test = ifelse(menage_id_member %in% ids_r1, 1,0),
         r2_test = ifelse(menage_id_member %in% ids_r2, 1,0))

humanR1e_sel$r1_esbl_pos = humanR1e_sel$esbl_pos 
table(humanR1e_sel$r1_esbl_pos)
humanR1e_sel = humanR1e_sel %>% select(-c(esbl_pos))

humanR2e_sel$r2_esbl_pos = humanR2e_sel$esbl_pos
table(humanR2e_sel$r2_esbl_pos)
humanR2e_sel = humanR2e_sel %>% select(-c(esbl_pos))

# For linkage
R1e = humanR1e_sel %>% select(menage_id_member,r1_esbl_pos)
R2e = humanR2e_sel %>% select(menage_id_member,r2_esbl_pos)

#test = left_join(data_eccmid,humanR1e_sel)

data_eccmid = left_join(data_eccmid,R1e, by="menage_id_member")
data_eccmid = left_join(data_eccmid,R2e, by="menage_id_member")

data_eccmid$r1_esbl_pos[which(is.na(data_eccmid$r1_esbl_pos) & data_eccmid$r1_test==1)] = 0 
data_eccmid$r2_esbl_pos[which(is.na(data_eccmid$r2_esbl_pos) & data_eccmid$r2_test==1)] = 0 

# data_eccmid = data_eccmid %>%
#   mutate(esbl_pos_r1 = ifelse(is.na(r1_esbl_pos) & r1_test==1, 0, r1_esbl_pos),
#          esbl_pos_r2 = ifelse(is.na(r2_esbl_pos) & r2_test==1, 0, r2_esbl_pos))
data_eccmid = left_join(data_eccmid,villages, by="village")

# Some checks
table(data_eccmid$r1_esbl_pos, useNA="always")
table(car_r1$germe_c,car_r1$esbl_pos, useNA="always")

table(data_eccmid$r2_esbl_pos, useNA="always")
table(car_r2$germe_c,car_r2$esbl_pos, useNA="always")

# Further checks
table(data_eccmid$r1_test) # 19 could not be found back in R0 so seems to have gone right
table(data_eccmid$r2_test) # 28 could not be found back in R0 so seems to have gone right

names(data_eccmid)

data_eccmid = data_eccmid %>% select(-c(village_name.x, intervention_text.x, household, id_ecantillon,date,germe_c,ajouter)) %>%
  rename(village_name = "village_name.y",
         intervention_text = "intervention_text.y") %>%
  # ADD Acquisitions and decolonisations 
  mutate(r1_acquisition = ifelse(esbl_pos == 0 & r1_esbl_pos==1,1,
                                 ifelse(is.na(r1_esbl_pos), NA,0)),
         r2_acquisition = ifelse(r1_esbl_pos == 0 & r2_esbl_pos==1,1,
                                 ifelse(is.na(r2_esbl_pos), NA,0)),
         r1_decolonisation = ifelse(esbl_pos == 1 & r1_esbl_pos==0,1,
                                 ifelse(is.na(esbl_pos)|is.na(r1_esbl_pos), NA,0)),
         r2_decolonisation = ifelse(r1_esbl_pos == 1 & r2_esbl_pos==0,1,
                                    ifelse(is.na(r1_esbl_pos)|is.na(r2_esbl_pos), NA,0)),
         age = as.numeric(age),
         agegr10 = cut(age, 
                       breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                       labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                       include.lowest = TRUE)
         )

table(data_eccmid$agegr10)

# Number of individuals we have complete records of
sum(complete.cases(data_eccmid[, c("esbl_pos", "r1_esbl_pos", "r2_esbl_pos")])) # 892

# Add household variables
data_eccmid_hh = left_join(data_eccmid,HR0_final, by="menage_id_member", suffix=c("","")) %>% 
  select(-c(household,id_ecantillon,germe_c,date)) %>%
  mutate(age = as.numeric(age),
         sexe = as.character(sexe))
sapply(data_eccmid_hh, function(x) class(x))

names(data_eccmid_hh)
class(data_eccmid_hh$age)

# Acquisitions
table(data_eccmid$r1_acquisition)
table(data_eccmid$r2_acquisition)

# by intervention
table(data_eccmid$intervention_text,data_eccmid$r1_acquisition)
table(data_eccmid$intervention_text,data_eccmid$r2_acquisition)

prop.table(table(data_eccmid$intervention_text,data_eccmid$r1_acquisition),1)
prop.table(table(data_eccmid$intervention_text,data_eccmid$r2_acquisition),1)

# By age
table(data_eccmid$agegr10,data_eccmid$r1_acquisition)
prop.table(table(data_eccmid$agegr10,data_eccmid$r1_acquisition),1)



# Decolonisations
table(data_eccmid$r1_decolonisation)
table(data_eccmid$r2_decolonisation)

# by intervention
prop.table(table(data_eccmid$intervention_text,data_eccmid$r1_decolonisation),1)
prop.table(table(data_eccmid$intervention_text,data_eccmid$r2_decolonisation),1)

# Add antibiotic use by village using weigthed averages
test = abx %>% group_by(village.cluster, round, providertype) %>%
  summarise(n = n(),
            watch = sum(watch==1),
            antibiotic = sum(antibiotic==1),
            hcu = unique(hcu))%>%
  mutate(watch_p = watch/n,
         abx_p = antibiotic/n) %>%
  filter( round == "baseline") 

# Should take only data from BF
# test_w = test %>% group_by(village.cluster) %>%
#   summarise(watch_p_w = sum(watch_p)
#   watch_


# Export dataset
write_xlsx(data_eccmid_hh, paste0(DirectoryDataOut, "/data_colonisation_eccmid.xlsx"))


###########################################################
# DESCRIPTIVE STATISTICS
###########################################################

# Household size
hist(wash_r0$nmbre_personne_menage)
summary(wash_r0$nmbre_personne_menage)

df_r0 = data_eccmid_hh
# Per household, how many positive
d = df_r0 %>% group_by(village_name, intervention_text, menage_id, esbl_pos) %>%
  summarise(n = n())

samples_per_hh = df_r0 %>% group_by(menage_id) %>%
  summarise(n_samples = n())
  
hh_size = as.data.frame(cbind(df_r0$menage_id,df_r0$nmbre_personne_menage))
names(hh_size) = c("menage_id","hh_size")
hh_size = hh_size[!duplicated(hh_size),]

d = left_join(d, hh_size, by="menage_id") %>%
  left_join(., samples_per_hh) %>% 
  mutate(hh_size = as.numeric(hh_size),
         hh_size_cor = ifelse(hh_size >7, 7, hh_size),
         n = as.numeric(n),
         n_samples = as.numeric(n_samples),
         f_pos = round(n/hh_size,2),
         f_pos_samples_taken = round(n/n_samples,2), # Number of positives over total samples taken
         f_pos_cor = round(n/hh_size_cor,2), # Assuming not all sample individuals are in the R0 database, but max of 7 individuals sampled per household
         f_pos_cor = ifelse(f_pos_cor>1, round(n/hh_size,2), f_pos_cor)
        )
d_pos = d %>% filter(esbl_pos==1)
max(d_pos$f_pos_cor, na.rm=T) # still observations above 1 needs checking


# For the larger households, not all are tested. At one point we stopped at max 5 per household
hist(as.numeric(d_pos$hh_size))

# Plot number of sampled positive ESBLs per household
summary(d_pos$f_pos)
summary(d_pos$f_pos_cor)

median = summary(d_pos$f_pos_cor)[3]
median_i = summary(d_pos$f_pos_cor[d_pos$intervention_text=="intervention"])[3]
median_ni = summary(d_pos$f_pos_cor[d_pos$intervention_text=="contrôle"])[3]
median_t = data.frame(rbind(c("intervention", median_i),c("contrôle", median_ni)))
names(median_t) = c("intervention_text", "Median")
#sorted_clusters <- with(d_pos, reorder(village_name, f_pos_cor, FUN = median))

# Plot boxplot of fraction positive per village per intervention group
bp = ggplot(d_pos, aes(x = village_name, y = f_pos_cor, fill = village_name)) +
  geom_jitter(alpha=0.5) + theme_bw() +
  # geom_hline(yintercept=median, 
  #            color = "red", size=2) +
  geom_boxplot() + theme(legend.position="none",
                         axis.text.x=element_text(size=15, angle=90),
                         axis.title=element_text(size=16,face="bold"),
                         axis.text=element_text(size=14,face="bold"),
                         strip.text=element_text(size=14,face="bold")) +
  ylim(0,1) + 
  facet_wrap(~intervention_text, scales=("free_x")) + 
  geom_abline(data = median_t, aes(intercept = as.numeric(Median), slope = 0),lty="dashed", size=1, col="red")+ 
  labs(title = "Boxplot of ESBL-E positive (%) per village cluster (Baseline)",
       x = "Village",
       y = "% positive")
print(bp)

d_sum = d %>% group_by(village_name,intervention_text) %>%
  summarise(mean = mean(f_pos_cor, na.rm=T),
            median = median(f_pos_cor,na.rm=T),
            q1 = quantile(f_pos_cor,probs=c(0.25), na.rm = T),
            q3 = quantile(f_pos_cor, probs=c(0.75), na.rm = T))
d_sum

# Plot of median fraction positive per village per intervention group
ggplot(d_sum, aes(x = village_name, y = median, col = village_name)) +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=q1,ymax=q3,width=0.5)) +
  facet_wrap(~intervention_text, scales=("free_x")) + 
  labs(title = "Median[IQR] of % positive per village clusters",
       x = "Village",
       y = "% positive (median, IQR)")

# Intervention vs controle groups 
d_pos %>% group_by(intervention_text) %>%
  summarise(mean_cor = mean(f_pos_cor, na.rm=T),
            median_cor = median(f_pos_cor,na.rm=T),
            q1_cor = quantile(f_pos_cor,probs=c(0.25), na.rm = T),
            q3_cor = quantile(f_pos_cor, probs=c(0.75), na.rm = T),
            mean = mean(f_pos, na.rm=T),
            median = median(f_pos,na.rm=T)) # seems rather similar (luckily)

# Plot density intervention vs control
dpi = ggplot(d_pos, aes(x=f_pos_cor, group=intervention_text, fill=intervention_text)) + 
  geom_density(aes(f_pos, ..scaled..)) 
dpi 

ggplot(d_pos, aes(x=f_pos_cor, group=village_name, fill=village_name)) + 
  geom_histogram() + facet_wrap(.~ village_name) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        axis.text=element_text(size=12,face="bold"),
        strip.text=element_text(size=12,face="bold")) +
  xlim(0,1) + 
  labs(x="ESBL-E positive (%) within household", y="number of households")


mean <- d_pos %>% group_by(village_name) %>%
  summarise(mean = mean(f_pos_cor, na.rm=T),)

ggplot(d_pos, aes(x=f_pos_cor, group=village_name, fill = village_name)) + 
  geom_density(aes(f_pos_cor, ..scaled..))+
  labs(x="%positive within hh", y="Density")

dp = ggplot(d_pos, aes(x=f_pos_cor, group=village_name, fill=village_name)) + 
  geom_density(aes(f_pos_cor, ..scaled..))+ 
  facet_wrap(.~village_name)+ geom_vline(data=mean, aes(xintercept=mean))+
  labs(x="%positive within hh", y="Density")
dp

# Number positive per round
prop.table(table(data_eccmid$esbl_pos))
prop.table(table(data_eccmid$r1_esbl_pos, useNA = "always"))
table(data_eccmid$r2_esbl_pos)

# Number of households per round
length(unique(car_r0$household))
length(unique(car_r1$household))
length(unique(car_r2$household))

# Number of household wash
length(unique(wash_r0$menage_id))

# Save plot
pdf(file="./Output/Figures/prevalence_per_village.pdf", width=7, height=4)
print(bp)
dev.off()

pdf(file="./Output/Figures/density_prevalence_per_village.pdf", width=5, height=4)
print(dp)
dev.off()

pdf(file="./Output/Figures/density_prevalence_per_intervention.pdf", width=5, height=4)
print(dpi)
dev.off()

# Export linked data
#write.csv(df_r0,paste0(DirectoryDataOut,"/bf_r0_lab_hh_linked.csv")) 

# missing links
r0_notlink = car_r0[car_r0$found_in_wash==0,] %>% select(menage_id,household,village)
write.csv(df_r0,paste0(DirectoryDataOut,"/need_checking/r0_lab_no_link.csv")) 




###########################################################
# AFTER RECEIVING NEW DATASETS
###########################################################


# Date: 15 April 2024

# Descriptives for which variables to include in analyses
table(HR0_final$esbl_pos)
table(car_r0$germe_c, car_r0$esbl_pos)

# Table 1

wash_r0_table1 = table1(~ factor(intervention_text) + as.numeric(age) +
       + factor(q2_source_princ_saison_seche) 
       + factor(q3_source_princ_saison_pluv)
       + factor(q4_bidon_stock)                    
       + factor(q5a_bidon_ferme_rempli)
       + factor(q5b_bidon_ferme_vide)
       + factor(q5c_bidon_nettoye)
       + factor(q6_traite_eau)
       + factor(q7_type_inst_sanitaire)
       + factor(q9_toilette_partagee)
       + factor(q10_combien_partag)
       + factor(q11_dernier_nettoyage)
       + factor(q12_elimine_selle_enf)
       + factor(q13_vidange_toilette)
       + factor(q14_produit_lavag_main)
       + factor(q15_lave_apr_defec)
       + factor(q16_lave_apr_repas)
       + factor(q17_animaux_menage)
       + factor(q1_diarrhee_prevenu___2)
       + factor(q1_diarrhee_prevenu___3)
       + factor(q1_diarrhee_prevenu___4)          
       + factor(q1_diarrhee_prevenu___5)
       + factor(q1_diarrhee_prevenu___6)
       + factor(q1_diarrhee_prevenu___7)          
       + factor(q1_diarrhee_prevenu___8)
       + factor(q1_diarrhee_prevenu___9)
       + factor(q18_animaux_interieur___1)
       + factor(q18_animaux_interieur___2)
       + factor(q18_animaux_interieur___3)
       + factor(q18_animaux_interieur___4)
       + factor(q18_animaux_interieur___5)
       + factor(q18_animaux_interieur___6)
       + factor(q19_animaux_dehors___1)
       + factor(q19_animaux_dehors___2)
       + factor(q19_animaux_dehors___3)
       + factor(q19_animaux_dehors___4)
       + factor(q19_animaux_dehors___5)
       + factor(q19_animaux_dehors___6)
       + factor(q20_excrement_animaux)
       + factor(q21_animal_malade___1)
       + factor(q21_animal_malade___2)
       + factor(q21_animal_malade___3)
       + factor(q21_animal_malade___4)
       +factor(q21_animal_malade___5)
       + factor(q21_animal_malade___6)
       + factor(eau_assainissement_hygine_complete)| factor(esbl_pos), data=HR0_final)
wash_r0_table1

# By cluster
wash_r0_village = table1(~ factor(intervention_text) + as.numeric(age) +
                          + factor(q2_source_princ_saison_seche) 
                        + factor(q3_source_princ_saison_pluv)
                        + factor(q4_bidon_stock)                    
                        + factor(q5a_bidon_ferme_rempli)
                        + factor(q5b_bidon_ferme_vide)
                        + factor(q5c_bidon_nettoye)
                        + factor(q6_traite_eau)
                        + factor(q7_type_inst_sanitaire)
                        + factor(q9_toilette_partagee)
                        + factor(q10_combien_partag)
                        + factor(q11_dernier_nettoyage)
                        + factor(q12_elimine_selle_enf)
                        + factor(q13_vidange_toilette)
                        + factor(q14_produit_lavag_main)
                        + factor(q15_lave_apr_defec)
                        + factor(q16_lave_apr_repas)
                        + factor(q17_animaux_menage)
                        + factor(q1_diarrhee_prevenu___2)
                        + factor(q1_diarrhee_prevenu___3)
                        + factor(q1_diarrhee_prevenu___4)          
                        + factor(q1_diarrhee_prevenu___5)
                        + factor(q1_diarrhee_prevenu___6)
                        + factor(q1_diarrhee_prevenu___7)          
                        + factor(q1_diarrhee_prevenu___8)
                        + factor(q1_diarrhee_prevenu___9)
                        + factor(q18_animaux_interieur___1)
                        + factor(q18_animaux_interieur___2)
                        + factor(q18_animaux_interieur___3)
                        + factor(q18_animaux_interieur___4)
                        + factor(q18_animaux_interieur___5)
                        + factor(q18_animaux_interieur___6)
                        + factor(q19_animaux_dehors___1)
                        + factor(q19_animaux_dehors___2)
                        + factor(q19_animaux_dehors___3)
                        + factor(q19_animaux_dehors___4)
                        + factor(q19_animaux_dehors___5)
                        + factor(q19_animaux_dehors___6)
                        + factor(q20_excrement_animaux)
                        + factor(q21_animal_malade___1)
                        + factor(q21_animal_malade___2)
                        + factor(q21_animal_malade___3)
                        + factor(q21_animal_malade___4)
                        +factor(q21_animal_malade___5)
                        + factor(q21_animal_malade___6)
                        + factor(eau_assainissement_hygine_complete)| factor(village_name), data=HR0_final)
wash_r0_village

t1flex(wash_r0_table1) %>% 
  save_as_docx(path="./Output/Tables/wash_r0_table1.docx")

# t1flex(wash_r0_village) %>% 
#   save_as_docx(path="./Output/Tables/wash_r0_village.docx",
#                page_size = officer::page_size(width = 21, height = 29.7/2.54, orient = "landscape"))



#####################################################
# ANALYSES FOR ECCMID ABSTRACT
#####################################################

coldata<-readxl::read_xlsx("./Data/BF/clean/data_colonisation_eccmid.xlsx")

# Note that if acquisition occurs in between screening round 0 (baseline) and round 1 (at time of intervention)
# then an acquisition cannot occur between round 1 and round 2 (3 month later) because individual will 
# already be positive at round 1. And if already positive at round 1 individual will not be at risk 
# of a post-intervention acquisition
# Note that field r2_acquisition currently doesn't take into account that individuals are not 
# at risk in round 1 

# Proposed analysi 1 s -  include in analysis only if a negative swab at round one
# i.e r1_esbl_pos =0   and a swab (whether pos or neg) at round 2 (i.e r2_esbl_pos is 1 or 0)
#  outcome is 1 if r2_esbl_pos  
# adjust for clustering at household level (menage_id) and include covariates for age (age), sex (sexe) and 
# intervention (as coded in field intervention_text)
# then as a sensitivity analysis consider only thoe with two initial negative swab

coldata$include_in_analysis1<-coldata$r1_esbl_pos == 0 & (coldata$r2_esbl_pos ==1 |coldata$r2_esbl_pos ==0)

coldata[ , c(4,5,6,7,12,15,16,92)]

coldataforanalysis1<-coldata[coldata$include_in_analysis1==TRUE & !is.na(coldata$include_in_analysis1), c(1,2,4,5,6,7,10,13,16,17,18,19,20,21,22,23,24,93)]
names(coldataforanalysis1)

coldataforanalysis1$intervention<-ifelse(coldataforanalysis1$intervention_text=="intervention",1 , 0)


# first with no random effects or other covariates
logreg0<-glm(r2_esbl_pos ~ intervention , data=coldataforanalysis1,  family = binomial(link = "logit"))
summary(logreg0)
exp(logreg0$coefficients)
exp(confint(logreg0)) # # Effect intervention: OR = 1.14 (CI:0.75 - 1.73)


# then with covariates and random effects 
logreg1<-glmer(r2_esbl_pos ~ age + sexe + intervention + (1|menage_id), data=coldataforanalysis1,  family = binomial(link = "logit"))
summary(logreg1)
exp(fixef(logreg1))
exp(confint(logreg1, method="Wald"))[2:5,]  # Effect intervention: OR = 1.23 (CI:0.74 - 2.01)



# now repeat but with a sensitivity analysis only taking those with first two swabs negative as at risk (as these are those we are most confident of not beng colonised)
coldata$include_in_analysis2<-coldata$esbl_pos==0 & coldata$r1_esbl_pos == 0 & (coldata$r2_esbl_pos ==1 |coldata$r2_esbl_pos ==0)

coldata[ , c(4,5,6,7,12,15,16,93)]


coldataforanalysis2<-coldata[coldata$include_in_analysis2==TRUE & !is.na(coldata$include_in_analysis1), c(1,2,4,5,6,7,10,13,16,17,18,19,20,21,22,23,24,93)]
coldataforanalysis2$intervention<-ifelse(coldataforanalysis2$intervention_text=="intervention",1 , 0)

# first with no random effects or other covariates
logreg0sens1<-glm(r2_esbl_pos ~ intervention , data=coldataforanalysis2,  family = binomial(link = "logit"))
summary(logreg0sens1)
exp(logreg0sens1$coefficients)
exp(confint(logreg0sens1)) # # Effect intervention: OR = 0.90 (CI:0.50 - 1.62)


logreg1sens1<-glmer(r2_esbl_pos ~ age + sexe + intervention + (1|menage_id), data=coldataforanalysis2,  family = binomial(link = "logit"))

# random effects models fails with error "boundary (singular) fit: see help('isSingular')"

# fixed effects with covariates
logreg1sens1a<-glm(r2_esbl_pos ~ age + sexe + intervention , data=coldataforanalysis2,  family = binomial(link = "logit"))
summary(logreg1sens1a)
coef = exp(logreg1sens1a$coefficients)
ci = exp(confint(logreg1sens1a))
cbind(coef,ci) # Effect intervention: OR = 1.02 (CI:0.55 - 1.87)

sum(complete.cases(coldata[, c("esbl_pos", "r1_esbl_pos", "r2_esbl_pos")])) # 892
d = coldata %>% filter(complete.cases(coldata[, c("esbl_pos", "r1_esbl_pos", "r2_esbl_pos")]))
length(unique(d$menage_id))
table(d$esbl_pos)






