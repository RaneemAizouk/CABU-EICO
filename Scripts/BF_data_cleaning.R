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

# Lab data
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_manualchange_no_password.xlsx"))
car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_nopassword.xlsx"))
car_r1 = read_xlsx(paste0(DirectoryData, "/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24_nopassword.xlsx"), sheet=1)
car_r2 = read_xlsx(paste0(DirectoryData, "/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24_nopassword.xlsx"), sheet=2)

# Lab ids vs household ids
hh_lab_ids =  read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))

# Household data
wash_r0 = read_xls(paste0(DirectoryData, "/WP4_WASH_07_09_2023_nopassword.xls"))

# Villages (that are the clusters) of CABU-EICO
villages = read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","intervention_yn")


# Check variables between r0 and r1 and r2 - seems we have received the household variables
# of those with a lab sample, not the lab data
which(names(car_r1) %in% names(car_r0))
which(names(car_r1) %in% names(wash_r0))
length(which(names(car_r1) %in% names(wash_r0)))
which(!names(car_r1) %in% names(wash_r0))
names(car_r1[which(names(car_r1) %in% names(car_r0))])
names(car_r1[which(!names(car_r1) %in% names(wash_r0))])
names(car_r1[which(!names(car_r1) %in% names(wash_r0))])

# Add variables village and household
car_r0$village = substr(car_r0$record_id, start = 1, stop = 2)
#car_r0$household = str_extract(car_r0$record_id_manual_change, "[^-]+")
car_r0$household = str_extract(car_r0$record_id, "[^-]+")
car_r0 = merge(car_r0, villages, by="village")

############################################################
# CLEAN LAB DATA
############################################################

# Clean germe; diameters are CLSI guidelines, 
# Brecht discussed with Jan jacobs, there understood ESBL positive is defined as cetriax/cefo <=22 (including I and R, based on CLSI guidelines,
# See in WP6 folder "word file: Interpretation antibiogramme des isolats PORTAGE asymptomatique_ESBL_E. coliKlebsielle.docx), 
# So I understood to use this, instead of esbltest == 1 
# This is then interpreted as, cetriax/cefo resistant following ESBL selective medium)

# Evk 20 September 2023:
# Would be good to do a final check with Yougbare if this is correct or whether
# we can just use all the observations in the car_r0 file (i.e. esbltest == 0 and esbltest == 1). 

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
table(car_r0$esbl_pos==1 & car_r0$testesbl==1) # difference of 7 which need resolving

names(car_r0)

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


#################################################################
# LINK LAB AND HOUSEHOLD DATA
#################################################################

# Create household variable in lab data (car_r0) to link to WASH survey by adding zero's.
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

wash_r0$village = substr(wash_r0$menage_id, start = 1, stop = 2)
wash_r0 = merge(wash_r0, villages, by="village")

# Link lab data with household data
car_r0$found_in_wash[which(car_r0$menage_id %in% wash_r0$menage_id)] = 1
car_r0$found_in_wash[is.na(car_r0$found_in_wash)] = 0
length(unique(car_r0$household[car_r0$found_in_wash==0])) # 28 households can not be found in WASH survey database; could be due to error's in ID or due to zero's that need to be removed, needs checking
unique(car_r0$menage_id[car_r0$found_in_wash==0])
unique(car_r0$household[car_r0$found_in_wash==0])

# Need to adjust the 28 IDs in the household dataset that have 0's lacking
ids_m = as.table(cbind(unique(car_r0$household[car_r0$found_in_wash==0]),
             unique(car_r0$menage_id[car_r0$found_in_wash==0])))
# 11 will be resolved through this method
which(ids_m[,1] %in% unique(wash_r0$menage_id))
ids_m[which(ids_m[,1] %in% unique(wash_r0$menage_id))] # All "CRAS"
unique(df$text[car_r0$found_in_wash==0]) 

# Need to check this (26 January 2024)
p = 1
for(i in ids_m[,1]){
  repl =  wash_r0$menage_id[wash_r0$menage_id==i] 
  wash_r0$menage_id[wash_r0$menage_id==i] = rep(ids_m[p,2], length(repl))
  p = p+1
}

# which left?
car_r0$found_in_wash[which(car_r0$menage_id %in% wash_r0$menage_id)] = 1
car_r0$found_in_wash[is.na(car_r0$found_in_wash)] = 0
length(unique(car_r0$household[car_r0$found_in_wash==0])) # 17 households can not be found in WASH survey database; could be due to error's in ID or due to zero's that need to be removed, needs checking
unique(car_r0$menage_id[car_r0$found_in_wash==0])
unique(car_r0$household[car_r0$found_in_wash==0])

# Need to still correct those, not all household ids are yet corrected to make linkage for all members possible



# Select variables - make database with variables excluding healthcare seeking behaviour survey questions (and related medicine use); as these
# are not 1 observation per household
wash_r0_sel = wash_r0 %>% select(menage_id,village, village_name, intervention_text, intervention_yn, 
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
                                q21_animal_malade___6,eau_assainissement_hygine_comple) %>%
  filter(!is.na(date_enquete))

# Link lab data with household data
# Perform linkage
df_r0 = car_r0%>%
  left_join(wash_r0_sel, by="menage_id", suffix=c("",".y"))%>%
              select(-ends_with(".y"))

table(df_r0$found_in_wash) # 55 households to still link
df_r0$household[which(df_r0$found_in_wash==0)]


###########################################################
# DESCRIPTIVE STATISTICS
###########################################################

# Per household, how many positive
d = df_r0 %>% group_by(village, intervention_text, household, esbl_pos) %>%
  summarise(n = n())

samples_per_hh = df_r0 %>% group_by(household) %>%
  summarise(n_samples = n())
  
hh_size = as.data.frame(cbind(df_r0$household,df_r0$nmbre_personne_menage))
names(hh_size) = c("household","hh_size")
hh_size = hh_size[!duplicated(hh_size),]

d = left_join(d, hh_size, by="household") %>%
  left_join(., samples_per_hh) %>% 
  mutate(hh_size = as.numeric(hh_size),
         hh_size_cor = ifelse(hh_size >5, 5, hh_size),
         n = as.numeric(n),
         n_samples = as.numeric(n_samples),
         f_pos = round(n/hh_size,2),
         f_pos_samples_taken = round(n/n_samples,2), # Number of positives over total samples taken
         f_pos_cor = round(n/hh_size_cor,2) # Assuming not all sample individuals are in the R0 database, but max of 5 individuals sampled per household
        )
d_pos = d %>% filter(esbl_pos==1)

# For the larger households, not all are tested. At one point we stopped at max 5 per household
hist(as.numeric(d_pos$hh_size))

# Plot number of sampled positive ESBLs per household
summary(d_pos$f_pos)
summary(d_pos$f_pos_cor)

sorted_clusters <- with(d_pos, reorder(village, f_pos, FUN = median))

# Plot boxplot of fraction positive per village per intervention group
ggplot(d_pos, aes(x = sorted_clusters, y = f_pos, fill = village)) +
  geom_jitter(alpha=0.5) + 
  geom_boxplot() + 
  facet_wrap(~intervention_text, scales=("free_x")) + 
  labs(title = "Boxplot of % positive per village clusters",
       x = "Village",
       y = "% positive")

d_sum = d %>% group_by(village,intervention_text) %>%
  summarise(mean = mean(f_pos, na.rm=T),
            median = median(f_pos,na.rm=T),
            q1 = quantile(f_pos,probs=c(0.25), na.rm = T),
            q3 = quantile(f_pos, probs=c(0.75), na.rm = T))
d_sum

# Plot of median fraction positive per village per intervention group
ggplot(d_sum, aes(x = village, y = median, col = village)) +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=q1,ymax=q3,width=0.5)) +
  facet_wrap(~intervention_text, scales=("free_x")) + 
  labs(title = "Boxplot of % positive per village clusters",
       x = "Village",
       y = "% positive (median, IQR)")

# Intervention vs controle groups 
d_pos %>% group_by(intervention_text) %>%
  summarise(mean = mean(f_pos, na.rm=T),
            median = median(f_pos,na.rm=T),
            q1 = quantile(f_pos,probs=c(0.25), na.rm = T),
            q3 = quantile(f_pos, probs=c(0.75), na.rm = T)) # seems rather similar (luckily)

# Plot density intervention vs control
ggplot(d_pos, aes(x=f_pos, group=intervention_text, fill=intervention_text)) + 
  geom_density(aes(f_pos, ..scaled..)) 

ggplot(d_pos, aes(x=f_pos, group=village, fill=village)) + 
  geom_histogram() + facet_wrap(.~ village) +
  labs(x="%positive within hh", y="number of households")

mean <- d_pos %>% group_by(village) %>%
  summarise(mean = mean(f_pos, na.rm=T),)

ggplot(d_pos, aes(x=f_pos, group=village, fill = village)) + 
  geom_density(aes(f_pos, ..scaled..))+
  labs(x="%positive within hh", y="Density")

ggplot(d_pos, aes(x=f_pos, group=village, fill=village)) + 
  geom_density(aes(f_pos, ..scaled..))+ 
  facet_wrap(.~village)+ geom_vline(data=mean, aes(xintercept=mean))+
  labs(x="%positive within hh", y="Density")
  # Save plot