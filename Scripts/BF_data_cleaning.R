#####################################################
# DATA CLEANING AND LINKAGE
#####################################################
# This code is cleaning the data and linking the different datasets

# 22 January 2024
# Last update: 

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_manualchange_no_password.xlsx"))
car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_nopassword.xlsx"))


#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-09-20_1301_manual_change_id_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-09-08_1301_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-05-04_manual_change_id_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/bf_esbl_r0.xlsx"), sheet=2)
hh_lab_ids =  read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))

wash_r0 = read_xls(paste0(DirectoryData, "/WP4_WASH_07_09_2023_nopassword.xls"))
wash_r0_stool = wash_r0 %>% filter(!is.na(num_echantillon))

wash_r0_stool_sel = wash_r0_stool %>% select(menage_id,num_echantillon, date_enquete,groupe,
                                             nbre_enf_0_5ans,nbre_menage_conc,
                                             cs_id_individu, dob_age,dob, age, sexe)

villages = read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","intervention_yn")

# Errors in ID (see email Daniel, Daniel is correcting those with Frank)

# Add variables village and household
car_r0$village = substr(car_r0$record_id, start = 1, stop = 2)
#car_r0$household = str_extract(car_r0$record_id_manual_change, "[^-]+")
car_r0$household = str_extract(car_r0$record_id, "[^-]+")
car_r0 = merge(car_r0, villages, by="village")

# Create household variable to link to WP4 survey by adding zero's.
df = data.frame(household = car_r0$household) 
df = df %>%  separate(household, 
                      into = c("text", "num"), 
                      sep = "(?<=[A-Za-z])(?=[0-9])")

df$num_ad = NULL
df$household = car_r0$household

for(i in 1:length(df$num)){
  if(nchar(df$num)[i]==3){
    p = "00000"
    df$num_ad[i] = paste0(p,df$num[i])
  }else if(nchar(df$num)[i]==4){
    p = "0000"
    df$num_ad[i] = paste0(p,df$num[i])
  }else if(nchar(df$num)[i]==5){
    p = "000"
    df$num_ad[i] = paste0(p,df$num[i])
  }
  else if(nchar(df$num)[i]==6){
    p = "00"
    df$num_ad[i] = paste0(p,df$num[i])
  }
}

#This is not yet doing the trick fully as some menage_id have no zero's (see nchar == 8 for some id's)
# Also some nchar == 12 for some in the new df$menage_id which should be 11
df$menage_id = paste0(df$text,df$num_ad)
nchar(df$menage_id)
nchar(wash_r0$menage_id)

car_r0$menage_id = df$menage_id

wash_r0_stool_sel$village = substr(wash_r0_stool_sel$menage_id, start = 1, stop = 2)
wash_r0_stool_sel = merge(wash_r0_stool_sel, villages, by="village")