#####################################################
# DATA CLEANING AND LINKAGE BURKINA FASO
#####################################################
# This code is cleaning the data and linking the different datasets

# 17 April 2024
# Author: Esther van Kleef

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, msm)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

# Lab data
car_bf = read.csv(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_2024-04-17_1527.csv"), sep=";")

# Lab ids vs household ids
hh_lab_ids =  readxl::read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))
names(hh_lab_ids)
names(car_bf)
names(hh_lab_ids) = c("household", "menage_id", "bras")

# Household data
hh_bf = readxl::read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBWP4_DATA_2024-04-17_1528.xlsx"))

# Villages (that are the clusters) of CABU-EICO
villages = readxl::read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","ajouter")

# Antibiotic use data
abx = read.csv("./Data/BF/clean/watch_acute.csv")


# Add variables village and household
# ALL ROUNDS
####################
car_bf$village = substr(car_bf$record_id, start = 1, stop = 2)
car_bf$household = str_extract(car_bf$record_id, "[^-]+")
car_bf = merge(car_bf, villages, by="village")

############################################################
# CLEAN LAB DATA
############################################################

# Clean germe; diameters are CLSI guidelines, 
# Jan jacobs comm: ESBL positive is defined as cetriax/cefo <=22 (including I and R, based on CLSI guidelines,
# See in WP6 folder "word file: Interpretation antibiogramme des isolats PORTAGE asymptomatique_ESBL_E. coliKlebsielle.docx), 
# So I understood to use this, instead of esbltest == 1 
# This is then interpreted as, cetriax/cefo resistant following ESBL selective medium)

# ALL ROUNDS
####################
unique(car_bf$germe)
car_bf = car_bf %>%
  mutate(germe_c = ifelse(germe %in% c("E.COLI", "E-COLI", "ECOLI", "E.COLI 2", "E.COLI 1", "eE-COLI",
                                       "E-COLI 2","E-CLI","E-CLOI", "E.COLIE","1"),"e.coli", 
                          ifelse(germe %in% c("SALMO", "SALMO SPP", "SALMONELLA SP","SALMONELLA SPP","SALMONELLE SPP","SALMONELLA SSP",
                                              "SALMONELLA", "SELMO"),"salmonella",NA)),
         germe_c = ifelse(morphotyp%in%c(1, NA),germe_c, 
                          paste0(germe_c, "_", morphotyp)),
         esbl_pos = ifelse(diametr_cetriax_or_cefota <= 22, 1, 0),
         date = as.Date(date, format="%d/%m/%Y"),
         date_conserv = as.Date(date, format="%d/%m/%Y")) 



table(car_bf$germe, car_bf$germe_c, useNA= "always") # 7 individuals with no germe indicated, make NA again
car_bf$germe_c[car_bf$germe==""] = NA
table(car_bf$germe, car_bf$germe_c, useNA= "always")

# Number of cases positive
table(car_bf$germe_c, car_bf$esbl_pos)
table(car_bf$esbl_pos, useNA="always") # These are the ESBL positive patients based on cetriax_or_cefota, 2825
table(car_bf$testesbl) # These are the ESBL positive patients based on esbl_pos, 2825
table(car_bf$esbl_pos==1 & car_bf$testesbl==1) # difference; we decided to ignore these differences

# Remove individuals with diametr_cetriax_or_cefota = NA
car_bf = car_bf %>% filter(!is.na(diametr_cetriax_or_cefota)) # 4 removed
table(car_bf$esbl_pos, useNA="always") # These are the ESBL positive patients based on cetriax_or_cefota, 2507

names(car_bf)

#################################################################
# CLEAN HOUSEHOLD DATA
#################################################################
# Good to note is that the household data contains questions from three surveys
# 1) Stool collection survey --> individual-level data (age, sex) of those of whom a stool sample was taken
# 2) WASH survey --> household level data, answered by the household head
# 3) Healthcare utilisation survey --> individual level data
# Here for each household member, the household head is asked about the number of 
# healthcare visits in the last 30 days
# Each of these visits (per provider type) are one line of data
# Therefore there are multiple observations within a household

# with redcap_repeat_instrument, the different survey answers can be recognised
# 1) Stool collection survey --> redcap_repeat_instrument == "formulaire_collecte_de_selles"
# 2) WASH survey --> redcap_repeat_instrument == " " and redcap_event_name == "round_0_arm_1"
# 3) Healthcare utilisation survey --> redcap_repeat_instrument == "visite_structure_sanitaire"
# 4) Healthcare utilisation survey medicines --> redcap_repeat_instrument == "mdicament"


# ALL ROUNDS
################################
hh_bf = hh_bf %>% mutate(
  dob = as.Date(dob, format = "%Y-%m-%d"),
  sexe = factor(sexe, levels=c(1,2), labels=c("Male", "Female")),
  date_enquete = as.Date(date_enquete, format="%Y-%m-%d"),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d"),
  age = tolower(age),
  age = as.numeric(age),
  agegr10 = cut(age, 
                breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf),
                labels = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
                include.lowest = TRUE),
  nmbre_personne_menage = as.numeric(nmbre_personne_menage),
  nbre_enf_0_5ans = as.numeric(nbre_enf_0_5ans),
  nbre_menage_conc = as.numeric(nbre_menage_conc), 
  date_enquete = as.Date(date_enquete, format="%Y-%m-%d"),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d"),
  q1_diarrhee_prevenu___1 = factor(q1_diarrhee_prevenu___1, levels=c(0,1), labels=c("No","Yes")),        
  q1_diarrhee_prevenu___2 = factor(q1_diarrhee_prevenu___2, levels=c(0,1), labels=c("No","Yes")),
  q1_diarrhee_prevenu___3 = factor(q1_diarrhee_prevenu___3, levels=c(0,1), labels=c("No","Yes")),        
  q1_diarrhee_prevenu___4 = factor(q1_diarrhee_prevenu___4, levels=c(0,1), labels=c("No","Yes")),
  q1_diarrhee_prevenu___5 = factor(q1_diarrhee_prevenu___5, levels=c(0,1), labels=c("No","Yes")),         
  q1_diarrhee_prevenu___6 = factor(q1_diarrhee_prevenu___6, levels=c(0,1), labels=c("No","Yes")),
  q1_diarrhee_prevenu___7 = factor(q1_diarrhee_prevenu___7, levels=c(0,1), labels=c("No","Yes")),         
  q1_diarrhee_prevenu___8 = factor(q1_diarrhee_prevenu___8, levels=c(0,1), labels=c("No","Yes")),
  q1_diarrhee_prevenu___9 = factor(q1_diarrhee_prevenu___9, levels=c(0,1), labels=c("No","Yes")), 
  autr_mesur_prev_diarrhe = ifelse(autr_mesur_prev_diarrhe %in% c("certain nouriture","Certain nouriture", "CERTAIN NOURITURE",
                                                                  "Certaine nouriture", "CERTAINE NOURITURE","Certains nouriture",
                                                                  "CERTAINS NOURITURE", "Eviter certain nouriture"), "certain_food",
                                   ifelse(autr_mesur_prev_diarrhe %in% c("LAVAGE DES MAINS AVANT ET APRES LE REPAS"), "hand washing", "")),
  q2_source_princ_saison_seche = factor(q2_source_princ_saison_seche, levels=c(1:9), labels=c("Tap house", 
                                                                                              "Tap concession",
                                                                                              "Tap public/fountain", 
                                                                                              "Borehole",
                                                                                              "improved well (protected)",
                                                                                              "unimproved well (unprotected)",
                                                                                              "rainwater",
                                                                                              "surface water (ponds, dams,rivers,lakes,pits,irrigation canals)",
                                                                                              "bagged water")),
  q3_source_princ_saison_pluv = factor(q3_source_princ_saison_pluv, levels=c(1:10), labels=c("Tap house", 
                                                                                             "Tap concession",
                                                                                             "Tap public/fountain", 
                                                                                             "Borehole",
                                                                                             "improved well (protected)",
                                                                                             "unimproved well(unprotected)",
                                                                                             "rainwater",
                                                                                             "surface water (ponds, dams,rivers,lakes,pits,irrigation canals)",
                                                                                             "bagged water", "bottled water")),  
  q4_bidon_stock = factor(q4_bidon_stock, levels = c(1:3), labels=c("Yes, cans", "Yes, only one large tank", "No")),
  q5a_bidon_ferme_rempli = factor(q5a_bidon_ferme_rempli, levels=c(1:2), labels=c("Yes", "No")),
  q5b_bidon_ferme_vide = factor(q5b_bidon_ferme_vide, levels=c(1:2), labels=c("Yes", "No")),           
  q5c_bidon_nettoye = factor(q5c_bidon_nettoye, levels=c(1,2,3,4,6), labels=c("Yes, with soap", "Yes, with water but no soap","Yes, boiled", "No", "Other")),
  q6_traite_eau = factor(q6_traite_eau, levels = c(1:7), labels= c("No", "Yes,boiling", "Yes,cholinate/add desinfectant", "Yes, filter with cloth",
                                                                   "Yes, filter with filter", "Yes, Solar desinfection (in the sun)","Yes, decant")),               
  q7_type_inst_sanitaire = factor(q7_type_inst_sanitaire, levels = c(3,4,5), labels = c("pit latrine with slab", "pit latrine without slab", "open defecation")),
  q8_autr_lieu_defecation___1 = factor(q8_autr_lieu_defecation___1, levels = c(0,1), labels = c("No","Yes")),      
  q8_autr_lieu_defecation___2 = factor(q8_autr_lieu_defecation___2, levels = c(0,1), labels = c("No","Yes")),
  q8_autr_lieu_defecation___3 = factor(q8_autr_lieu_defecation___3, levels = c(0,1), labels = c("No","Yes")),   
  q8_autr_lieu_defecation___4 = factor(q8_autr_lieu_defecation___4, levels = c(0,1), labels = c("No","Yes")),
  q8_autr_lieu_defecation___5 = factor(q8_autr_lieu_defecation___5, levels = c(0,1), labels = c("No","Yes")),   
  q8_autr_lieu_defecation___6 = factor(q8_autr_lieu_defecation___6, levels = c(0,1), labels = c("No","Yes")),
  q8_autr_lieu_defecation___7 = factor(q8_autr_lieu_defecation___7, levels = c(0,1), labels = c("No","Yes")),  
  q9_toilette_partagee = factor(q9_toilette_partagee, levels=c(1:3), labels=c("Yes, other households (non-public)", "Yes, public", "No")),             
  q10_combien_partag = as.numeric(q10_combien_partag),
  q11_dernier_nettoyage = factor(q11_dernier_nettoyage, levels=c(1:6), labels=c("<24h", ">24h, but <1week", "1-4weeks", ">1month", "Never", "Don't know")), 
  q12_elimine_selle_enf = factor(q12_elimine_selle_enf, levels = c(1:9), labels=c("Child used toilet/latrine", "Thrown/rinsed into toilet/latrine",
                                                                                  "Thrown/rinsed into drainage pit",
                                                                                                   "Trown in garbage", "Buried", "Disposed in open air", 
                                                                                                   "Used as manure", "Other", "NA (no child)")),
  q13_vidange_toilette = factor(q13_vidange_toilette, levels=c(1,2,4,5,7), labels=c("Has not been drained yet", "Don't know", "Removed by service provider and covered in pit",
                                                                                    "Removed by service provider (don't know where)", "Emptied by hh in uncover pit/open ground")),
  q14_produit_lavag_main = factor(q14_produit_lavag_main, levels=c(1:5), labels=c("Yes, soap", "Yes, detergent", "Yes, ash/mud/sand", "No, none available", "No, available but not used")),
  q15_lave_apr_defec = factor(q15_lave_apr_defec, levels=c(1:4), labels=c("Yes, always", "Yes, often", "No or rarerly", "Not sure")),       
  q16_lave_apr_repas = factor(q16_lave_apr_repas, levels=c(1:4), labels=c("Yes, always", "Yes, often", "No or rarerly", "Not sure")),
  q17_animaux_menage = factor(q17_animaux_menage, levels=c(1:4), labels = c("Yes, inside and outside", "Yes, outside next to house","Yes, outside in demarked area", "No")),          
  q18_animaux_interieur___1 = factor(q18_animaux_interieur___1, levels = c(0,1), labels = c("No", "Yes")),
  q18_animaux_interieur___2 = factor(q18_animaux_interieur___2, levels = c(0,1), labels = c("No", "Yes")),       
  q18_animaux_interieur___3 = factor(q18_animaux_interieur___3, levels = c(0,1), labels = c("No", "Yes")),
  q18_animaux_interieur___4 = factor(q18_animaux_interieur___4, levels = c(0,1), labels = c("No", "Yes")),       
  q18_animaux_interieur___5 = factor(q18_animaux_interieur___5, levels = c(0,1), labels = c("No", "Yes")),
  q18_animaux_interieur___6 = factor(q18_animaux_interieur___6, levels = c(0,1), labels = c("No", "Yes")),
  q19_animaux_dehors___1 = factor(q19_animaux_dehors___1, levels = c(0,1), labels = c("No", "Yes")),        
  q19_animaux_dehors___2 = factor(q19_animaux_dehors___2, levels = c(0,1), labels = c("No", "Yes")),
  q19_animaux_dehors___3 = factor(q19_animaux_dehors___3, levels = c(0,1), labels = c("No", "Yes")),        
  q19_animaux_dehors___4 = factor(q19_animaux_dehors___4, levels = c(0,1), labels = c("No", "Yes")),
  q19_animaux_dehors___5 = factor(q19_animaux_dehors___5, levels = c(0,1), labels = c("No", "Yes")),         
  q19_animaux_dehors___6 = factor(q19_animaux_dehors___6, levels = c(0,1), labels = c("No", "Yes")),
  
  q20_excrement_animaux = factor(q20_excrement_animaux, levels = c(1,2,3), labels = c("Yes", "No", "Not possible to determine")),  
  q21_animal_malade___1 = factor(q21_animal_malade___1, levels = c(0,1), labels = c("No", "Yes")),           
  q21_animal_malade___2 = factor(q21_animal_malade___2, levels = c(0,1), labels = c("No", "Yes")),  
  q21_animal_malade___3 = factor(q21_animal_malade___3, levels = c(0,1), labels = c("No", "Yes")),  
  q21_animal_malade___4 = factor(q21_animal_malade___4, levels = c(0,1), labels = c("No", "Yes")),  
  q21_animal_malade___5 = factor(q21_animal_malade___5, levels = c(0,1), labels = c("No", "Yes")),             
  q21_animal_malade___6 = factor(q21_animal_malade___6, levels = c(0,1), labels = c("No", "Yes")),
  eau_assainissement_hygine_complete = factor(eau_assainissement_hygine_complete, levels=c(0,2), labels=c("No", "Yes"))
) %>% select(-dob) %>% rename(
  n.householdmember = "nmbre_personne_menage",
  n.child.0to5 = "nbre_enf_0_5ans",
  n.households.concession = "nbre_menage_conc",
  date.consent = "date_consentement",
  date.stool.collection = "date_recuperation_selle",
  q1.diar.prev.water.pot.covered = "q1_diarrhee_prevenu___1",
  q1.diar.prev.no.finger.in.waterglass = "q1_diarrhee_prevenu___2",
  q1.diar.prev.no.finger.in.waterglass =  "q1_diarrhee_prevenu___2",
  q1.diar.prev.utensil.to.take.water.from.pot = "q1_diarrhee_prevenu___3",    
  q1.diar.prev.cover.food = "q1_diarrhee_prevenu___4",
  q1.diar.prev.boil.water = "q1_diarrhee_prevenu___5",
  q1.diar.prev.filter.water = "q1_diarrhee_prevenu___6",
  q1.diar.prev.other = "q1_diarrhee_prevenu___7",
  q1.diar.prev.cant.be.avoided = "q1_diarrhee_prevenu___8",
  q1.diar.prev.dont.know = "q1_diarrhee_prevenu___9",
  q2.main.water.source.dry= "q2_source_princ_saison_seche",
  q3.main.water.source.rainy= "q3_source_princ_saison_pluv",
  q4.cans.storage.water= "q4_bidon_stock", 
  q5a.cans.storage.water.closed.when.filled= "q5a_bidon_ferme_rempli", 
  q5b.cans.storage.water.closed.wen.empty= "q5b_bidon_ferme_vide",          
  q5c.cans.cleaned.before.reuse= "q5c_bidon_nettoye", 
  q6.treatment.water= "q6_traite_eau", 
  q7.principle.defication= "q7_type_inst_sanitaire", 
  q8.other.defecation.flush.toiled.septic  = "q8_autr_lieu_defecation___1", 
  q8.other.defecation.pit.latrine.ventilation = "q8_autr_lieu_defecation___2",
  q8.other.defecation.pit.latrine.slab = "q8_autr_lieu_defecation___3",   
  q8.other.defecation.pit.latrine.no.slab = "q8_autr_lieu_defecation___4",
  q8.other.defecation.open.defecation = "q8_autr_lieu_defecation___5",   
  q8.other.defecation.other = "q8_autr_lieu_defecation___6",
  q8.other.defecation.none = "q8_autr_lieu_defecation___7",
  q9.shared.toilet= "q9_toilette_partagee",
  q10.n.shared.toilet= "q10_combien_partag",
  q11.toilet.last.cleaned= "q11_dernier_nettoyage",
  q12.disposal.child.stool= "q12_elimine_selle_enf",
  q13.disposal.latrine.pit= "q13_vidange_toilette",
  q14.handwashing.product= "q14_produit_lavag_main",
  q15.handwashing.defecation= "q15_lave_apr_defec",
  q16.handwashing.meals= "q16_lave_apr_repas",
  q17.animals.around.household= "q17_animaux_menage",
  q18.animal.inside.cow = "q18_animaux_interieur___1",
  q18.animal.inside.sheep.goat = "q18_animaux_interieur___2",       
  q18.animal.inside.pig = "q18_animaux_interieur___3",
  q18.animal.inside.donkey.horse = "q18_animaux_interieur___4",       
  q18.animal.inside.chicken.goose.duck = "q18_animaux_interieur___5",
  q18.animal.inside.other = "q18_animaux_interieur___6",
  q19.animal.inside.cow = "q19_animaux_dehors___1",
  q19.animal.inside.sheep.goat = "q19_animaux_dehors___2",       
  q19.animal.inside.pig = "q19_animaux_dehors___3",
  q19.animal.inside.donkey.horse = "q19_animaux_dehors___4",       
  q19.animal.inside.chicken.goose.duck = "q19_animaux_dehors___5",
  q19.animal.inside.other = "q19_animaux_dehors___6",
  q20.animal.excrement.floor= q20_excrement_animaux, 
  q21.when.animal.ill.treatment.with.vet = "q21_animal_malade___1",           
  q21.when.animal.ill.treatment.without.vet= "q21_animal_malade___2",  
  q21.when.animal.ill.bucher = "q21_animal_malade___3",  
  q21.when.animal.ill.eat.meat.at.home= "q21_animal_malade___4",  
  q21.when.animal.ill.burie.dispose = "q21_animal_malade___5",             
  q21.when_animal.ill.autre = "q21_animal_malade___6"
)

table(hh_bf$agegr10)

# Link village and cluster data to household survey
hh_bf$data_row = c(1:nrow(hh_bf)) # This variable we can use for identifying back those individuals that had a sample taken
hh_bf$village = substr(hh_bf$menage_id, start = 1, stop = 2)
hh_bf = merge(hh_bf, villages, by="village")
names(hh_bf)

names(hh_bf) = gsub("_",".",names(hh_bf))
hh_bf = hh_bf %>%
 rename(menage_id = "menage.id",
        village_name = "village.name",
        redcap_event_name = "redcap.event.name")

# Select relevant WASH variables from the household survey
################################################################################################

# variables excluding healthcare seeking behaviour survey questions (and related medicine use); as these
# are not 1 observation per household
hh_bf_wash = hh_bf %>% select(data.row,menage_id,village, village_name, intervention.text,   
                                 redcap_event_name,
                                 date.enquete,
                                 groupe,
                                 n.householdmember, 
                                 n.child.0to5,
                                 n.households.concession,
                                 q1.diar.prev.water.pot.covered,
                                 q1.diar.prev.no.finger.in.waterglass,                    
                                 q1.diar.prev.utensil.to.take.water.from.pot, q1.diar.prev.cover.food,                                  
                                 q1.diar.prev.boil.water, q1.diar.prev.filter.water,                               
                                 q1.diar.prev.other,q1.diar.prev.cant.be.avoided,                            
                                 q1.diar.prev.dont.know, autr.mesur.prev.diarrhe,                                  
                              q2.main.water.source.dry, q3.main.water.source.rainy,                              
                              q4.cans.storage.water, q5a.cans.storage.water.closed.when.filled,                                   
                              q5b.cans.storage.water.closed.wen.empty, q5c.cans.cleaned.before.reuse,                                        
                              q6.treatment.water, q6.autre.traitmen.eau,                                    
                              q7.principle.defication, q7.autr.typ.ins.sanitair,                                 
                              q8.other.defecation.flush.toiled.septic, q8.other.defecation.pit.latrine.ventilation,             
                              q8.other.defecation.pit.latrine.slab, q8.other.defecation.pit.latrine.no.slab,                  
                              q8.other.defecation.open.defecation, q8.other.defecation.other,                                
                              q8.other.defecation.none, q8.autre.preciser,                                        
                              q9.shared.toilet, q10.n.shared.toilet,                                       
                              q11.toilet.last.cleaned, q12.disposal.child.stool,                                    
                              q12.autre.preciser, q13.disposal.latrine.pit,                                     
                              q13.autre.preciser, q14.handwashing.product,                                   
                              q15.handwashing.defecation, q16.handwashing.meals,                                       
                              q17.animals.around.household, q18.animal.inside.cow,                                    
                              q18.animal.inside.sheep.goat, q18.animal.inside.pig,                                    
                              q18.animal.inside.donkey.horse, q18.animal.inside.chicken.goose.duck,                     
                              q18.animal.inside.other,q18.autre.specifie,                                       
                              q19.animal.inside.cow, q19.animal.inside.sheep.goat,                             
                              q19.animal.inside.pig,q19.animal.inside.donkey.horse,                           
                              q19.animal.inside.chicken.goose.duck, q19.animal.inside.other,                                  
                              q19.autre.specifie, q20.animal.excrement.floor,                                    
                              q21.when.animal.ill.treatment.with.vet, q21.when.animal.ill.treatment.without.vet,                
                              q21.when.animal.ill.bucher,q21.when.animal.ill.eat.meat.at.home,                     
                              q21.when.animal.ill.burie.dispose,q21.when.animal.ill.autre,                                
                              eau.assainissement.hygine.complete) %>%
  filter(!is.na(date.enquete)) # Denominator data (i.e. people tested for esbl) for R0


#################################################################
# LINK LAB AND HOUSEHOLD DATA
#################################################################

# Create household variable in lab data (car.bf) to link to WASH survey by adding zero's.

# ALL ROUNDS
###################################################################################

# Link lab data with lab vs hh survey ids (as IDs are in different format)
car_bf = left_join(car_bf, hh_lab_ids, by="household")
# Check if all linked
table(car_bf$bras, useNA= "always")

# Link lab data with hh survey ids
length(unique(car_bf$menage_id)) # 261
length(unique(hh_bf$menage_id)) # 808

# Can all IDs be traced back from lab to wash survey?
car_bf$found_in_hh[which(car_bf$menage_id %in% hh_bf$menage_id)] = 1
car_bf$found_in_hh[is.na(car_bf$found_in_wash)] = 0
length(unique(car_bf$household[car_bf$found_in_wash==0])) # after cleaning, all households can be found in WASH survey database; 
unique(car_bf$menage_id[car_bf$found_in_wash==0])
unique(car_bf$household[car_bf$found_in_wash==0])


# Make dataset with only those household individuals that had a stool sample taken and their individual variables
hh_bf_lab = hh_bf %>% filter(!is.na(cs.id.individu)) %>% # ensure just 1 observation per person of whom individual is esbl positive
  select(data.row, redcap_event_name,cs.id.individu,num.echantillon, menage_id, village, age, agegr10, sexe, date.consent, date.stool.collection)
  
# Merge wash patient characteristics with lab data - NEED TO GET IDs IIN THE SAME FORMAT
which(car_bf$record_id %in% unique(hh_bf_lab$cs.id.individu)) # Have to change the format of both to make sure matching can be done
head(car_bf$record_id ); head(hh_bf_lab$cs.id.individu)
unique(car_bf$record_id)
unique(hh_bf$cs.id.individu)

# SEE HOW TO DO THE LINKAGE
# HH database does no have "-" and puts and "M" before each household member. Also sometimes 01 and sometimes 1 as method of writing
# PROBABLY IF WE TAKE FROM wash_r0 cs_id_individue all "MX" (so M + number after ID), then put that in a seperate column.
# Then we take the household number from the WASH survey (so menage_id), combine these to with a hyphen

# THEN for car_bf we take also menage_Id and digits after the "-", remove the 0's then two can be combined

# Change IDs to the same format

# Create a variable in which we will safe the new formatted and cleaned ID
hh_bf_lab$menage_id_member = NA

# Create a variable that will just store the household member number
hh_bf_lab$member = NA

# Extract the number after M to get household number
hh_bf_lab$member =  gsub(".*M", "", hh_bf_lab$cs.id.individu)
# Remove leading '0's
hh_bf_lab$member = as.character(as.numeric(hh_bf_lab$member))
# Check the one's which are now NA
table(is.na(hh_bf_lab$member))
hh_bf_lab$cs.id.individu[is.na(hh_bf_lab$member)] # 22 individuals have a household member number missing, see if still all lab_ids can be linked when merging with car_bf
#View(wash_r0_lab[is.na(wash_r0_lab$member),])
hh_bf_lab$num.echantillon[is.na(hh_bf_lab$member)] # of 0 we can get them from num_echantillon

# Now create new variable for linking with lab dataset
hh_bf_lab$menage_id_member = paste0(hh_bf_lab$menage_id, "-", hh_bf_lab$member)
hh_bf_lab$menage_id_member
# Make NAs for the one's that still need checking
hh_bf_lab$menage_id_member[is.na(hh_bf_lab$member)] = NA
table(is.na(hh_bf_lab$menage_id_member))


# Check if dubplicates from wash_r0
table(duplicated(hh_bf_lab$menage_id_member)) # 3018 duplicates (as multiple rounds of data)
dups = unique(hh_bf_lab$menage_id_member[duplicated(hh_bf_lab$menage_id_member)]) # 1157
#View(wash_bf_lab[which(wash_bf_lab$menage_id_member%in%dups),])

# Then for remainder keep first record 
hh_bf_lab_de = hh_bf_lab %>% filter(!duplicated(menage_id_member))  # 1237
table(duplicated(hh_bf_lab_de$menage_id_member)) # 

# Now create the same variable for the lab dataset
# Create a variable in which we will safe the new formatted and cleaned ID
car_bf$menage_id_member = NA

# Create a variable that will just store the household member number
car_bf$member = NA

# Extract the number after M to get household number
car_bf$member =  gsub(".*-", "", car_bf$record_id)
# Remove leading '0's
car_bf$member = as.character(as.numeric(car_bf$member))
# Check the one's which are now NA
table(is.na(car_bf$member))
car_bf$record_id[is.na(car_bf$member)] # 43 individuals have a household member number after M
car_bf$member[is.na(car_bf$member)] =  gsub(".*M", "", car_bf$record_id[is.na(car_bf$member)])
table(is.na(car_bf$member))

# Now create new variable for linking with lab dataset
car_bf$menage_id_member = paste0(car_bf$menage_id, "-", car_bf$member)
car_bf$menage_id_member

table(car_bf$germe_c)

# Remove ECC1801 (checked with Franck), which should be ECC01101, and are therefore duplicates so can be removed
which(car_bf$household=="ECC1801")

car_bf = car_bf %>% filter(!household=="ECC1801") 


# Create cohort 
##########################################################################################

# ROUND 0
#####################

car_bf_r0 = car_bf %>% filter(redcap_event_name =="round_0_arm_1")
length(unique(car_bf_r0$menage_id_member))
hh_bf_lab_r0 = hh_bf_lab %>% filter(redcap_event_name =="round_0_arm_1")

# three seperate datasets for e.coli, ecoli_2, ecoli_3 and salmonella so merge goes well
HRe = car_bf %>% filter(germe_c == "e.coli", redcap_event_name =="round_0_arm_1")
HRe2 = car_bf %>% filter(germe_c == "e.coli_2", redcap_event_name =="round_0_arm_1")
HRe3 = car_bf %>% filter(germe_c == "e.coli_3", redcap_event_name =="round_0_arm_1")
HRs = car_bf %>% filter(germe_c == "salmonella", redcap_event_name =="round_0_arm_1")

# check if duplicates
#View(HRe[HRe$menage_id_member %in% c(HRe$menage_id_member[duplicated(HRe$menage_id_member)]), ]) # 3 duplicates
#View(HRe2[HRe2$menage_id_member %in% c(HRe2$menage_id_member[duplicated(HRe2$menage_id_member)]), ]) # 1 dup
#View(HRs[HRs$menage_id_member %in% c(HRs$menage_id_member[duplicated(HRs$menage_id_member)]), ])

# Remove duplicates
HRe = HRe %>% filter(!duplicated(HRe$menage_id_member))
HRe2 = HRe2 %>% filter(!duplicated(HRe2$menage_id_member))

# MERGE individual hh characteristics with car_r0
HRe_l = left_join(HRe,hh_bf_lab, by= c("menage_id_member", "menage_id", "village","redcap_event_name"))
#table(HRe_l$data.row, useNA="always")
#HRe_l$menage_id_member[is.na(HRe_l$data.row)] # Appearantly there are a view IDs with round_0_arm_1 not present in the hh_bf_lab
#hh_bf_lab[hh_bf_lab$menage_id_member%in% HRe_l$menage_id_member[is.na(HRe_l$data.row)],] # Appearantly there are a view IDs with round_0_arm_1 not present in the hh_bf_lab
# but we do want to link with redcap_event_name otherwise we don't get the right date.consent, which we need for model fitting

# For now use an ugly fix 
HRe_l_nl = HRe_l %>% filter(is.na(data.row))
HRe_l_nl = left_join(HRe_l_nl,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""), multiple="last")
HRe_l_nl = HRe_l_nl %>% select(-c(member))
HRe_l_nl = HRe_l_nl %>% mutate(
  date.consent = NA,
  date.stool.collection = NA,
  redcap_event_name="round_0_arm_1") # as the stool collection date is not correct, linked to R1, not R0


HRe_l = HRe_l %>% filter(!is.na(data.row))
HRe_l = rbind(HRe_l, HRe_l_nl)

HRe_w = left_join(HRe_l,hh_bf_wash, by= c("menage_id", "village"), multiple="first") # keeps first matched 
#HRe_w = HRe_w %>% filter(!duplicated(data.row.x)&!(is.na(data.row.x)))
names(HRe_w)

table(HRe$esbl_pos)
table(HRe_w$esbl_pos)

HRe2_l = left_join(HRe2,hh_bf_lab, by= c("menage_id_member", "menage_id", "village","redcap_event_name"))
HRe2_w = left_join(HRe2_l,hh_bf_wash, by= c("menage_id", "village"))

table(HRe2_w$esbl_pos)
table(HRe2$esbl_pos) 

HRs_l = left_join(HRs,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village","redcap_event_name"))
HRs_w = left_join(HRs_l,hh_bf_wash, by= c("menage_id", "village"))
table(HRs_w$esbl_pos)
table(HRs$esbl_pos)

HR0_all = rbind(HRe_w,HRe2_w,HRs_w)
table(HR0_all$esbl_pos)
names(HR0_all)

# remove the duplicate rows
HR0_all = HR0_all %>%
  select(-c(member.y, village_name.y, redcap_event_name.y)) %>%
  rename(member = "member.x",
         village_name = "village_name.x",
         redcap_event_name = "redcap_event_name.x")


# Now add those ones with no E.coli-E or salmonella
testedR0 = hh_bf_lab_de$menage_id_member[hh_bf_lab_de$menage_id_member%in%c(HR0_all$menage_id_member)]
hh_bf_lab_de_R0_no = hh_bf_lab_r0 %>% filter(redcap_event_name == "round_0_arm_1" & !menage_id_member%in%c(testedR0)) 
hh_bf_lab_de_R0_no = left_join(hh_bf_lab_de_R0_no, hh_bf_wash, by= c("menage_id", "village"))
hh_bf_lab_de_R0_no = hh_bf_lab_de_R0_no %>% filter(!duplicated(data.row.x)) 

table(is.na(hh_bf_lab_de_R0_no$menage_id_member)) # # There is one individual that has no ID member, remove
hh_bf_lab_de_R0_no = hh_bf_lab_de_R0_no %>% filter(!is.na(menage_id_member))

# Add column names so datasets of negatives and postives can be rbinded
NOT.names <- names(HR0_all)[!names(HR0_all)%in%names(hh_bf_lab_de_R0_no)]
hh_bf_lab_de_R0_no[,NOT.names] <- NA

# Step 1: Identify the desired order of column names in the reference dataset
desired_order <- names(HR0_all)

# Step 2: Reorder the columns of the dataset based on the desired order
hh_bf_lab_de_R0_no =  hh_bf_lab_de_R0_no[, desired_order]
hh_bf_lab_de_R0_no = hh_bf_lab_de_R0_no[,which(names(hh_bf_lab_de_R0_no) %in% desired_order)] 

HR0_all = rbind(HR0_all,hh_bf_lab_de_R0_no)



# Check if cohort is complete
length(unique(HR0_all$menage_id_member))
table(HR0_all$germe_c, HR0_all$esbl_pos)
table(car_bf_r0$germe_c, car_bf_r0$esbl_pos) # difference of 3 comes from the duplicates


# Keep only relevant variables
names(HR0_all)
# Check which date to take for the analyses
table(HR0_all$date.consent, useNA="always")
table(HR0_all$date.stool.collection, useNA="always")
table(HR0_all$date_conserv, useNA="always")
table(HR0_all$date, useNA="always")
HR0_all%>%filter(is.na(date.consent)) # it seems fairly safe to use the date variable for those which had the stool collection date missing as not with R0 in the hh dataset but now leave as NA
d = HR0_all %>%mutate(
  diff = date - date.consent,
  diff_start = date.consent - min(date.consent, na.rm=T)
)


table((as.numeric(d$diff)))
mean(d$diff,na.rm=T) # could use this to impute the date for those with an unlikely date

table((as.numeric(d$diff_start))) # should not be more than 90-100 days as one collection round took 3 months.

table(HR0_all$date.consent)
#View(HR0_all %>% filter(date.consent>"2023-01-30")) # Just a few have an unlikely date of consent, this may be due to typo. 
# So date.consent has the least number of NAs
max(HR0_all$date, na.rm=T)
max(HR0_all$date.stool.collection, na.rm=T)
max(HR0_all$date.consent, na.rm=T)
#View(HR0_all %>% filter(date.consent == "2023-12-06"))
#View(HR0_all %>% filter(is.na(date.consent)))
HR0_all$date[HR0_all$date.consent>"2023-03-30"]

HR0_all = HR0_all %>% select(-c(redcap_repeat_instrument,redcap_repeat_instance,
                                check_list, id_ecantillon, morphotyp, germe, testesbl,
                                rsultats_antibiogramme_portage_asymptomatique_salm_complete,
                                interpretr_ampici_amoxicil,comment_ampici_amoxic,
                                interpr_amoxi_acid_clavu, comment_amoxi_acid_clavu,
                                interpr_piperaci_tazobac, comment_piperacil_tazobact,
                                interpr_cetriax_or_cefotax, comment_cetriax_cefataxi,
                                interpr_cefepime, comment_cefepime,
                                interpr_meropenem, comment_meropenem,
                                interpr_ertapenem, comment_ertapenem,
                                interpr_gentamycine, comment_gentamycine,
                                interpr_amykacine, comment_amykacine,
                                interpr_pfloxacine, comment_ciprofloxacine,
                                interpr_ciprofloxacine, comment_pfloxacine,
                                interpr_sulfame_trimethop, comment_sulfa_trimethop,
                                groupe, bras, ajouter,
                                found_in_hh, member, num.echantillon,cs.id.individu,
                                autr.mesur.prev.diarrhe,
                                q7.autr.typ.ins.sanitair,
                                q8.autre.preciser,
                                q12.autre.preciser,
                                q18.autre.specifie
                                )) %>%
  rename(data.row.hh.lab = "data.row.x",
         data.row.hh.wash = "data.row.y") %>%
  mutate(esbl_pos = ifelse(is.na(esbl_pos), 0, esbl_pos),
         date = as.Date(date, format = "%d/%m/%Y"),
         date_conserv = as.Date(date_conserv, format ="%d/%m/%Y"),
         intervention.text = factor(intervention.text, levels = c("contrôle", "intervention"), labels=c("control", "intervention")),
         esbl_pos = factor(esbl_pos, levels=c(0,1), labels=c("No", "Yes")),
         r0.esble = ifelse(germe_c %in% c("e.coli","e.coli_2","e.coli_3") & esbl_pos == "Yes", 1, 0),
         r0.salm = ifelse(germe_c %in% c("salmonella"), 1, 0))

#View(HR0_all %>% filter(is.na(record_id)))
HR0_all$redcap_event_name = "round_0_arm_1"
HR0_all = left_join(HR0_all,villages, by="village")
HR0_all = HR0_all %>% select(-c(village_name.x,intervention_text.x, ajouter))%>%
  rename(village_name = "village_name.y",
         intervention_text = "intervention_text.y") %>% 
  mutate(ast_done = 
           ifelse(is.na(record_id), "No", "Yes")) %>%
  select(-c(record_id))
  


# which IDs have a salmonella?
idsalm = HR0_all$menage_id_member[HR0_all$r0.salm == 1]

# which have at least 1 ESBL-E
idesble = HR0_all$menage_id_member[HR0_all$r0.esble == 1]

table(HR0_all$r0.esble, useNA="always")
table(HR0_all$r0.salm)

# Also create a dataset with just one observation per individual
d = HR0_all %>% filter(is.na(germe_c) | germe_c %in% c("e.coli","e.coli_2","e.coli_3","salmonella")) %>% 
  mutate(
    r0.salm = ifelse(menage_id_member%in%idsalm, 1,0),
    r0.esble = ifelse(menage_id_member %in% idesble, 1, 0)
  ) %>% 
  filter(!duplicated(menage_id_member)) 

# Check if correct
table(d$r0.esble) # 679 esbl.e
table(d$esbl_pos) # 679 esbl.e
table(HR0_all$germe_c, HR0_all$esbl_pos) # 675 e.coli with esble and another 97 with second e. coli. There will be overlap between those two

HR0_e = d

HR0_all = HR0_all %>% select(-c(esbl_pos)) # as we created new variable esble
HR0_e = HR0_e %>% select(-c(esbl_pos)) # as we created new variable esble

# Export datasets
write.csv(HR0_all, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_all_r0.csv")) 
write.csv(HR0_e, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r0.csv")) 

rm(d, hh_bf_lab_de_R0_no, HRe, HRe_l, HRe_w, HRe2,HRe2_l, HRe2_w,HRe3, HRs, HRs_l,HRs_w, HRe_l_nl)

# ROUND 1
#############################################################################
car_bf_r1 = car_bf %>% filter(redcap_event_name =="round_1_arm_1")
hh_bf_lab_r1 = hh_bf_lab %>% filter(redcap_event_name =="round_1_arm_1")

length(unique(car_bf_r1$menage_id_member))

# three seperate datasets for e.coli, ecoli_2, ecoli_3 and salmonella so merge goes well
HRe = car_bf %>% filter(germe_c == "e.coli", redcap_event_name =="round_1_arm_1")
HRe2 = car_bf %>% filter(germe_c == "e.coli_2", redcap_event_name =="round_1_arm_1")
HRe3 = car_bf %>% filter(germe_c == "e.coli_3", redcap_event_name =="round_1_arm_1")
HRs = car_bf %>% filter(germe_c == "salmonella", redcap_event_name =="round_1_arm_1")

# check if duplicates
#View(HRe[HRe$menage_id_member %in% c(HRe$menage_id_member[duplicated(HRe$menage_id_member)]), ]) # 4 duplicates
#View(HRe2[HRe2$menage_id_member %in% c(HRe2$menage_id_member[duplicated(HRe2$menage_id_member)]), ]) # no dup
#View(HRs[HRs$menage_id_member %in% c(HRs$menage_id_member[duplicated(HRs$menage_id_member)]), ]) # no dup

# Remove duplicates
HRe = HRe %>% filter(!duplicated(HRe$menage_id_member))

# MERGE indivdual hh characteristics with car_r1
HRe_l = left_join(HRe,hh_bf_lab, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe_w = left_join(HRe_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#HRe_w = HRe_w %>% filter(!duplicated(data.row.x))

table(HRe$esbl_pos)
table(HRe_w$esbl_pos)

HRe2_l = left_join(HRe2,hh_bf_lab, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe2_w = left_join(HRe2_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")

table(HRe2_l$esbl_pos)
table(HRe2_w$esbl_pos) 

HRs_l = left_join(HRs,hh_bf_lab, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"),multiple="first")
HRs_w = left_join(HRs_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")
table(HRs_l$germe_c)
table(HRs_w$germe_c)

HR1_all = rbind(HRe_w,HRe2_w,HRs_w)
table(HR1_all$esbl_pos, useNA="always")

# Now add those ones with no E.coli-E or salmonella
testedR1 = hh_bf_lab_r1$menage_id_member[hh_bf_lab_r1$menage_id_member%in%c(HR1_all$menage_id_member)]

hh_bf_lab_de_R1_no = hh_bf_lab_r1 %>% filter(redcap_event_name == "round_1_arm_1" & !menage_id_member%in%c(testedR1)) 
hh_bf_lab_de_R1_no = left_join(hh_bf_lab_de_R1_no, hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#hh_bf_lab_de_R1_no = hh_bf_lab_de_R1_no %>% filter(!duplicated(data.row.x))


# Add column names so datasets of negatives and postives can be rbinded
NOT.names <- names(HR1_all)[!names(HR1_all)%in%names(hh_bf_lab_de_R1_no)]
hh_bf_lab_de_R1_no[,NOT.names] <- NA

# Step 1: Identify the desired order of column names in the reference dataset
desired_order <- names(HR1_all)

# Step 2: Reorder the columns of the dataset based on the desired order
hh_bf_lab_de_R1_no =  hh_bf_lab_de_R1_no[, desired_order]
hh_bf_lab_de_R1_no = hh_bf_lab_de_R1_no[,which(names(hh_bf_lab_de_R1_no) %in% desired_order)] 

HR1_all = rbind(HR1_all,hh_bf_lab_de_R1_no)

# Check if cohort is complete
length(unique(HR1_all$menage_id_member))
table(HR1_all$germe_c, HR1_all$esbl_pos)
table(car_bf_r1$germe_c, car_bf_r1$esbl_pos) 


# Keep only relevant variables
names(HR1_all)


# Check which date to take for the analyses
table(HR1_all$date.consent, useNA="always")

HR1_all%>%filter(is.na(date.consent)) # it seems fairly safe to use the date variable for those which had the stool collection date missing as not with R0 in the hh dataset but now leave as NA
d = HR1_all %>%mutate(
  diff = date.enquete - date.consent,
  diff_start = date.consent - min(date.consent, na.rm=T)
)

table((as.numeric(d$diff)))
mean(d$diff,na.rm=T) # could use this to impute the date for those with an unlikely date
table((as.numeric(d$diff_start))) # should not be more than 90-100 days as one collection round took 3 months.

table(HR1_all$date.consent)
table(HR1_all$date.enquete)

max(HR1_all$date.consent, na.rm=T)
min(HR1_all$date.consent, na.rm=T) # Something is not right in the hh_bf vs redcap_event_name rounds
#View(HR1_all[HR1_all$date.consent<"2023-03-30",])

HR1_all = HR1_all %>% select(-c(redcap_repeat_instrument,redcap_repeat_instance,
                                check_list, id_ecantillon, morphotyp, germe, testesbl,
                                rsultats_antibiogramme_portage_asymptomatique_salm_complete,
                                interpretr_ampici_amoxicil,comment_ampici_amoxic,
                                interpr_amoxi_acid_clavu, comment_amoxi_acid_clavu,
                                interpr_piperaci_tazobac, comment_piperacil_tazobact,
                                interpr_cetriax_or_cefotax, comment_cetriax_cefataxi,
                                interpr_cefepime, comment_cefepime,
                                interpr_meropenem, comment_meropenem,
                                interpr_ertapenem, comment_ertapenem,
                                interpr_gentamycine, comment_gentamycine,
                                interpr_amykacine, comment_amykacine,
                                interpr_pfloxacine, comment_ciprofloxacine,
                                interpr_ciprofloxacine, comment_pfloxacine,
                                interpr_sulfame_trimethop, comment_sulfa_trimethop,
                                village_name.y,redcap_event_name.x,redcap_event_name.y, groupe, bras, ajouter,
                                found_in_hh, member.x,member.y, num.echantillon,cs.id.individu,
                                autr.mesur.prev.diarrhe,
                                q7.autr.typ.ins.sanitair,
                                q8.autre.preciser,
                                q12.autre.preciser,
                                q18.autre.specifie
)) %>%
  rename(data.row.hh.lab = "data.row.x",
         data.row.hh.wash = "data.row.y",
         village_name = "village_name.x") %>%
  mutate(esbl_pos = ifelse(is.na(esbl_pos), 0, esbl_pos),
         date = as.Date(date, format = "%d/%m/%Y"),
         date_conserv = as.Date(date_conserv, format ="%d/%m/%Y"),
         intervention.text = factor(intervention.text, levels = c("contrôle", "intervention"), labels=c("control", "intervention")),
         esbl_pos = factor(esbl_pos, levels=c(0,1), labels=c("No", "Yes")),
         r1.esble = ifelse(germe_c %in% c("e.coli","e.coli_2","e.coli_3") & esbl_pos == "Yes",1, 0),
         r1.salm = ifelse(germe_c %in% c("salmonella"), 1, 0)
  )

HR1_all$redcap_event_name = "round_1_arm_1"
HR1_all = left_join(HR1_all,villages, by="village")
HR1_all = HR1_all %>% select(-c(village_name.x,intervention_text.x, ajouter))%>%
  rename(village_name = "village_name.y",
         intervention_text = "intervention_text.y")%>% 
  mutate(ast_done = 
           ifelse(is.na(record_id), "No", "Yes"))%>%
  select(-c(record_id))


# which IDs have a salmonella?
idsalm = HR1_all$menage_id_member[HR1_all$r1.salm == 1]

# which have at least 1 ESBL-E
idesble = HR1_all$menage_id_member[HR1_all$r1.esble == 1]

table(HR1_all$r1.esble)
table(HR1_all$r1.salm)

# Also create a dataset with just one observation per individual
d = HR1_all %>% filter(is.na(germe_c) | germe_c %in% c("e.coli","e.coli_2","e.coli_3", "salmonella")) %>% 
  mutate(
    r1.salm = ifelse(menage_id_member%in%idsalm, 1,0),
    r1.esble = ifelse(menage_id_member %in% idesble, 1,0)
  ) %>% 
  filter(!duplicated(menage_id_member)) 

# Check if correct
table(d$r1.esble) # 583 esbl.e
table(d$esbl_pos) # 583 esbl.e
table(HR1_all$germe_c, HR1_all$esbl_pos) # 583 e.coli with esble, overlap with those with a second one

HR1_e = d

HR1_all = HR1_all %>% select(-c(esbl_pos)) # as we created new variable esble
HR1_e = HR1_e %>% select(-c(esbl_pos)) # as we created new variable esble

# Export datasets
write.csv(HR1_all, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_all_r1.csv")) 
write.csv(HR1_e, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r1.csv")) 

rm(d, hh_bf_lab_de_R1_no, HRe, HRe_l, HRe_w, HRe2,HRe2_l, HRe2_w,HRe3, HRs, HRs_l,HRs_w)

# ROUND 2
##############################################################################################
car_bf_r2 = car_bf %>% filter(redcap_event_name =="round_2_arm_1")
hh_bf_lab_r2 = hh_bf_lab %>% filter(redcap_event_name =="round_2_arm_1")

length(unique(car_bf_r2$menage_id_member))

# three seperate datasets for e.coli, ecoli_2, ecoli_3 and salmonella so merge goes well
HRe = car_bf %>% filter(germe_c == "e.coli", redcap_event_name =="round_2_arm_1")
HRe2 = car_bf %>% filter(germe_c == "e.coli_2", redcap_event_name =="round_2_arm_1")
HRe3 = car_bf %>% filter(germe_c == "e.coli_3", redcap_event_name =="round_2_arm_1")
HRs = car_bf %>% filter(germe_c == "salmonella", redcap_event_name =="round_2_arm_1")

# check if duplicates
#View(HRe[HRe$menage_id_member %in% c(HRe$menage_id_member[duplicated(HRe$menage_id_member)]), ]) # 6 duplicates
#View(HRe2[HRe2$menage_id_member %in% c(HRe2$menage_id_member[duplicated(HRe2$menage_id_member)]), ]) # no dup
#View(HRs[HRs$menage_id_member %in% c(HRs$menage_id_member[duplicated(HRs$menage_id_member)]), ]) # no dup

# Remove duplicates
HRe = HRe %>% filter(!duplicated(HRe$menage_id_member))

# MERGE indivdual hh characteristics with car_r1
HRe_l = left_join(HRe,hh_bf_lab, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe_w = left_join(HRe_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#HRe_w = HRe_w %>% filter(!duplicated(data.row.x))

table(HRe$esbl_pos)
table(HRe_w$esbl_pos)

HRe2_l = left_join(HRe2,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe2_w = left_join(HRe2_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")

table(HRe2_w$esbl_pos)
table(HRe2$esbl_pos) 

HRs_l = left_join(HRs,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRs_w = left_join(HRs_l,hh_bf_wash, by= c("menage_id", "village"),multiple="last")
table(HRs_w$germe_c)
table(HRs$germe_c)

HR2_all = rbind(HRe_w,HRe2_w,HRs_w)
table(HR2_all$esbl_pos)

# Now add those ones with no E.coli-E or salmonella
testedR2 = hh_bf_lab_r2$menage_id_member[hh_bf_lab_r2$menage_id_member%in%c(HR2_all$menage_id_member)]

hh_bf_lab_de_R2_no = hh_bf_lab_r2 %>% filter(redcap_event_name == "round_2_arm_1" & !menage_id_member%in%c(testedR2)) 
hh_bf_lab_de_R2_no = left_join(hh_bf_lab_de_R2_no, hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#hh_bf_lab_de_R2_no = hh_bf_lab_de_R2_no %>% filter(!duplicated(data.row.x))


# Add column names so datasets of negatives and postives can be rbinded
NOT.names <- names(HR2_all)[!names(HR2_all)%in%names(hh_bf_lab_de_R2_no)]
hh_bf_lab_de_R2_no[,NOT.names] <- NA

# Step 1: Identify the desired order of column names in the reference dataset
desired_order <- names(HR2_all)

# Step 2: Reorder the columns of the dataset based on the desired order
hh_bf_lab_de_R2_no =  hh_bf_lab_de_R2_no[, desired_order]
hh_bf_lab_de_R2_no = hh_bf_lab_de_R2_no[,which(names(hh_bf_lab_de_R2_no) %in% desired_order)] 

HR2_all = rbind(HR2_all,hh_bf_lab_de_R2_no)

# Check if cohort is complete
length(unique(HR2_all$menage_id_member))
table(HR2_all$germe_c, HR2_all$esbl_pos)
table(car_bf_r2$germe_c, car_bf_r2$esbl_pos) # difference of 6 comes from the duplicates


# Keep only relevant variables
names(HR2_all)

HR2_all = HR2_all %>% select(-c(redcap_repeat_instrument,redcap_repeat_instance,
                                check_list, id_ecantillon, morphotyp, germe, testesbl,
                                rsultats_antibiogramme_portage_asymptomatique_salm_complete,
                                interpretr_ampici_amoxicil,comment_ampici_amoxic,
                                interpr_amoxi_acid_clavu, comment_amoxi_acid_clavu,
                                interpr_piperaci_tazobac, comment_piperacil_tazobact,
                                interpr_cetriax_or_cefotax, comment_cetriax_cefataxi,
                                interpr_cefepime, comment_cefepime,
                                interpr_meropenem, comment_meropenem,
                                interpr_ertapenem, comment_ertapenem,
                                interpr_gentamycine, comment_gentamycine,
                                interpr_amykacine, comment_amykacine,
                                interpr_pfloxacine, comment_ciprofloxacine,
                                interpr_ciprofloxacine, comment_pfloxacine,
                                interpr_sulfame_trimethop, comment_sulfa_trimethop,
                                village_name.y, groupe, bras, ajouter,redcap_event_name.x,redcap_event_name.y,
                                found_in_hh, member.x,member.y, num.echantillon,cs.id.individu,
                                autr.mesur.prev.diarrhe,
                                q7.autr.typ.ins.sanitair,
                                q8.autre.preciser,
                                q12.autre.preciser,
                                q18.autre.specifie
)) %>%
  rename(data.row.hh.lab = "data.row.x",
         data.row.hh.wash = "data.row.y",
         village_name = "village_name.x") %>%
  mutate(esbl_pos = ifelse(is.na(esbl_pos), 0, esbl_pos),
         date = as.Date(date, format = "%d/%m/%Y"),
         date_conserv = as.Date(date_conserv, format ="%d/%m/%Y"),
         intervention.text = factor(intervention.text, levels = c("contrôle", "intervention"), labels=c("control", "intervention")),
         esbl_pos = factor(esbl_pos, levels=c(0,1), labels=c("No", "Yes")),
         r2.esble = ifelse(germe_c %in% c("e.coli","e.coli_2","e.coli_3") & esbl_pos == "Yes",1,0),
         r2.salm = ifelse(germe_c %in% c("salmonella"), 1,0)
  )
HR2_all$redcap_event_name = "round_2_arm_1"
HR2_all = left_join(HR2_all,villages, by="village")
HR2_all = HR2_all %>% select(-c(village_name.x,intervention_text.x, ajouter))%>%
  rename(village_name = "village_name.y",
         intervention_text = "intervention_text.y")%>% 
  mutate(ast_done = 
           ifelse(is.na(record_id), "No", "Yes"))%>%
  select(-c(record_id))

# which IDs have a salmonella?
idsalm = HR2_all$menage_id_member[HR2_all$r2.salm == 1]

# which have at least 1 ESBL-E
idesble = HR2_all$menage_id_member[HR2_all$r2.esble == 1]

table(HR2_all$r2.esble)
table(HR2_all$r2.salm)

# Also create a dataset with just one observation per individual
d = HR2_all %>% filter(is.na(germe_c) | germe_c %in% c("e.coli","e.coli_2","e.coli_3", "salmonella")) %>% 
  mutate(
    r2.salm = ifelse(menage_id_member%in%idsalm, 1,0),
    r2.esble = ifelse(menage_id_member %in% idesble, 1,0)
  ) %>% 
  filter(!duplicated(menage_id_member)) 

# Check if correct
table(d$r2.esble) # 708 esbl.e
table(d$esbl_pos) # 707 esbl.e
table(HR2_all$germe_c, HR2_all$esbl_pos) # 705 e.coli with esble, overlap with those with a second one

HR2_e = d

HR2_all = HR2_all %>% select(-c(esbl_pos)) # as we created new variable esble
HR2_e = HR2_e %>% select(-c(esbl_pos)) # as we created new variable esble

# Export datasets
write.csv(HR2_all, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_all_r2.csv")) 
write.csv(HR2_e, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r2.csv")) 

rm(d, hh_bf_lab_de_R2_no, HRe, HRe_l, HRe_w, HRe2,HRe2_l, HRe2_w,HRe3, HRs, HRs_l,HRs_w)


# ROUND 3
#############################################################################
car_bf_r3 = car_bf %>% filter(redcap_event_name =="round_3_arm_1")
hh_bf_lab_r3 = hh_bf_lab %>% filter(redcap_event_name =="round_3_arm_1")

length(unique(car_bf_r3$menage_id_member))

# three seperate datasets for e.coli, ecoli_2, ecoli_3 and salmonella so merge goes well
HRe = car_bf %>% filter(germe_c == "e.coli", redcap_event_name =="round_3_arm_1")
HRe2 = car_bf %>% filter(germe_c == "e.coli_2", redcap_event_name =="round_3_arm_1")
HRe3 = car_bf %>% filter(germe_c == "e.coli_3", redcap_event_name =="round_3_arm_1")
HRs = car_bf %>% filter(germe_c == "salmonella", redcap_event_name =="round_3_arm_1")

# check if duplicates
#View(HRe[HRe$menage_id_member %in% c(HRe$menage_id_member[duplicated(HRe$menage_id_member)]), ]) # 10 duplicates
#View(HRe2[HRe2$menage_id_member %in% c(HRe2$menage_id_member[duplicated(HRe2$menage_id_member)]), ]) # no dup
#View(HRs[HRs$menage_id_member %in% c(HRs$menage_id_member[duplicated(HRs$menage_id_member)]), ]) # no dup

# Remove duplicates
HRe = HRe %>% filter(!duplicated(HRe$menage_id_member))

# MERGE indivdual hh characteristics with car_r1
HRe_l = left_join(HRe,hh_bf_lab, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe_w = left_join(HRe_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#HRe_w = HRe_w %>% filter(!duplicated(data.row.x))

table(HRe$esbl_pos)
table(HRe_w$esbl_pos)

HRe2_l = left_join(HRe2,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village", "redcap_event_name"), multiple="first")
HRe2_w = left_join(HRe2_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")

table(HRe2_w$esbl_pos)
table(HRe2$esbl_pos) 

HRs_l = left_join(HRs,hh_bf_lab_de, by= c("menage_id_member", "menage_id", "village","redcap_event_name"), multiple="first")
HRs_w = left_join(HRs_l,hh_bf_wash, by= c("menage_id", "village"), multiple="last")
table(HRs_l$germe_c)
table(HRs$germe_c)

HR3_all = rbind(HRe_w,HRe2_w,HRs_w)
table(HR3_all$esbl_pos)

# Now add those ones with no E.coli-E or salmonella
testedR3 = hh_bf_lab_r3$menage_id_member[hh_bf_lab_r3$menage_id_member%in%c(HR3_all$menage_id_member)]

hh_bf_lab_de_R3_no = hh_bf_lab_r3 %>% filter(redcap_event_name == "round_3_arm_1" & !menage_id_member%in%c(testedR3)) 
hh_bf_lab_de_R3_no = left_join(hh_bf_lab_de_R3_no, hh_bf_wash, by= c("menage_id", "village"), multiple="last")
#hh_bf_lab_de_R3_no = hh_bf_lab_de_R3_no %>% filter(!duplicated(data.row.x))


# Add column names so datasets of negatives and postives can be rbinded
NOT.names <- names(HR3_all)[!names(HR3_all)%in%names(hh_bf_lab_de_R3_no)]
hh_bf_lab_de_R3_no[,NOT.names] <- NA

# Step 1: Identify the desired order of column names in the reference dataset
desired_order <- names(HR3_all)

# Step 2: Reorder the columns of the dataset based on the desired order
hh_bf_lab_de_R3_no =  hh_bf_lab_de_R3_no[, desired_order]
hh_bf_lab_de_R3_no = hh_bf_lab_de_R3_no[,which(names(hh_bf_lab_de_R3_no) %in% desired_order)] 

HR3_all = rbind(HR3_all,hh_bf_lab_de_R3_no)

# Check if cohort is complete
length(unique(HR3_all$menage_id_member))
table(HR3_all$germe_c, HR3_all$esbl_pos)
table(car_bf_r3$germe_c, car_bf_r3$esbl_pos) # difference of 5 comes from the duplicates


# Keep only relevant variables
names(HR3_all)

HR3_all = HR3_all %>% select(-c(redcap_repeat_instrument,redcap_repeat_instance,
                                check_list, id_ecantillon, morphotyp, germe, testesbl,
                                rsultats_antibiogramme_portage_asymptomatique_salm_complete,
                                interpretr_ampici_amoxicil,comment_ampici_amoxic,
                                interpr_amoxi_acid_clavu, comment_amoxi_acid_clavu,
                                interpr_piperaci_tazobac, comment_piperacil_tazobact,
                                interpr_cetriax_or_cefotax, comment_cetriax_cefataxi,
                                interpr_cefepime, comment_cefepime,
                                interpr_meropenem, comment_meropenem,
                                interpr_ertapenem, comment_ertapenem,
                                interpr_gentamycine, comment_gentamycine,
                                interpr_amykacine, comment_amykacine,
                                interpr_pfloxacine, comment_ciprofloxacine,
                                interpr_ciprofloxacine, comment_pfloxacine,
                                interpr_sulfame_trimethop, comment_sulfa_trimethop,
                                village_name.y, groupe, bras, ajouter,redcap_event_name.x,redcap_event_name.y,
                                found_in_hh, member.x,member.y, num.echantillon,cs.id.individu,
                                autr.mesur.prev.diarrhe,
                                q7.autr.typ.ins.sanitair,
                                q8.autre.preciser,
                                q12.autre.preciser,
                                q18.autre.specifie
)) %>%
  rename(data.row.hh.lab = "data.row.x",
         data.row.hh.wash = "data.row.y",
         village_name = "village_name.x") %>%
  mutate(esbl_pos = ifelse(is.na(esbl_pos), 0, esbl_pos),
         date = as.Date(date, format = "%d/%m/%Y"),
         date_conserv = as.Date(date_conserv, format ="%d/%m/%Y"),
         intervention.text = factor(intervention.text, levels = c("contrôle", "intervention"), labels=c("control", "intervention")),
         esbl_pos = factor(esbl_pos, levels=c(0,1), labels=c("No", "Yes")),
         r3.esble = ifelse(germe_c %in% c("e.coli","e.coli_2","e.coli_3") & esbl_pos == "Yes", 1,0),
         r3.salm = ifelse(germe_c %in% c("salmonella"), 1, 0)
  ) 
HR3_all$redcap_event_name = "round_3_arm_1"
HR3_all = left_join(HR3_all,villages, by="village")
HR3_all = HR3_all %>% select(-c(village_name.x,intervention_text.x, ajouter))%>%
  rename(village_name = "village_name.y",
         intervention_text = "intervention_text.y")%>% 
  mutate(ast_done = 
           ifelse(is.na(record_id), "No", "Yes"))%>%
  select(-c(record_id))

# which IDs have a salmonella?
idsalm = HR3_all$menage_id_member[HR3_all$r3.salm == 1]

# which have at least 1 ESBL-E
idesble = HR3_all$menage_id_member[HR3_all$r3.esble == 1]

table(HR3_all$r3.esble)
table(HR3_all$r3.salm)

# Also create a dataset with just one observation per individual
d = HR3_all %>% filter(is.na(germe_c) | germe_c %in% c("e.coli","e.coli_2","e.coli_3", "salmonella")) %>% 
  mutate(
    r3.salm = ifelse(menage_id_member%in%idsalm, 1, 0),
    r3.esble = ifelse(menage_id_member %in% idesble, 1,0)
  ) %>% 
  filter(!duplicated(menage_id_member)) 

# Check if correct
table(d$r3.esble) # 554 esbl.e
table(d$esbl_pos) # 554 esbl.e
table(HR3_all$germe_c, HR3_all$esbl_pos) # 554 e.coli with esble, overlap with those with a second one

HR3_e = d

HR3_all = HR3_all %>% select(-c(esbl_pos)) # as we created new variable esble
HR3_e = HR3_e %>% select(-c(esbl_pos)) # as we created new variable esble

# Export datasets
write.csv(HR3_all, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_all_r3.csv")) 
write.csv(HR3_e, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r3.csv")) 

rm(d, hh_bf_lab_de_R3_no, HRe, HRe_l, HRe_w, HRe2,HRe2_l, HRe2_w,HRe3, HRs, HRs_l,HRs_w)

# ALL ROUNDS

# To allow for rbinding, change name of the esbl.e and salmonella columns
d0 = HR0_all
d1 = HR1_all
d2 = HR2_all
d3 = HR3_all

d0 = d0 %>% rename(
  esble = "r0.esble",
  salm = "r0.salm"
)

d1 = d1 %>% rename(
  esble = "r1.esble",
  salm = "r1.salm"
)

d2 = d2 %>% rename(
  esble = "r2.esble",
  salm = "r2.salm"
)

d3 = d3 %>% rename(
  esble = "r3.esble",
  salm = "r3.salm"
)

HR_all_total = rbind(d0,d1,d2,d3)
table(HR_all_total$esble)
table(car_bf$esbl_pos) # Difference likely comes from duplicates that were removed

# Export dataset
write.csv(HR_all_total, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_all_r0123.csv")) 


# For e coli datasets the same
d0 = HR0_e
d1 = HR1_e
d2 = HR2_e
d3 = HR3_e

d0 = d0 %>% rename(
  esble = "r0.esble",
  salm = "r0.salm"
)

d1 = d1 %>% rename(
  esble = "r1.esble",
  salm = "r1.salm"
)

d2 = d2 %>% rename(
  esble = "r2.esble",
  salm = "r2.salm"
)

d3 = d3 %>% rename(
  esble = "r3.esble",
  salm = "r3.salm"
)

HR_e_total = rbind(d0,d1,d2,d3)

# Export dataset
write.csv(HR_e_total, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r0123.csv")) 

# Now create one in wide format
#####################################################################
# For that we will remove the MIC values as these will be different each round
d0 = HR0_e %>% select(-c(redcap_event_name, date, date_conserv, date.consent,date.stool.collection,
                           date_conserv,germe_c, data.row.hh.lab, data.row.hh.wash,date.enquete,
                           diametr_ampici_amoxici, diametr_amoxi_acid_clavu,diametr_piepera_tazobac,                    
                           diametr_cetriax_or_cefota,diametr_cefepime, diametr_meropenem, diametr_ertapenem,                          
                           diametr_gentamycine, diametr_amykacine, diametr_ciprofloxacine, diametr_pfloxacine,                         
                           diametr_sulfamet_trimeth))

d1 = HR1_e %>% select(-c(redcap_event_name, date, date_conserv, date.consent,date.stool.collection,
                           date_conserv,germe_c, data.row.hh.lab, data.row.hh.wash,date.enquete,
                           diametr_ampici_amoxici, diametr_amoxi_acid_clavu,diametr_piepera_tazobac,                    
                           diametr_cetriax_or_cefota,diametr_cefepime, diametr_meropenem, diametr_ertapenem,                          
                           diametr_gentamycine, diametr_amykacine, diametr_ciprofloxacine, diametr_pfloxacine,                         
                           diametr_sulfamet_trimeth))

d2 = HR2_e %>% select(-c(redcap_event_name, date, date_conserv, date.consent,date.stool.collection,
                           date_conserv,germe_c, data.row.hh.lab, data.row.hh.wash,date.enquete,
                           diametr_ampici_amoxici, diametr_amoxi_acid_clavu,diametr_piepera_tazobac,                    
                           diametr_cetriax_or_cefota,diametr_cefepime, diametr_meropenem, diametr_ertapenem,                          
                           diametr_gentamycine, diametr_amykacine, diametr_ciprofloxacine, diametr_pfloxacine,                         
                           diametr_sulfamet_trimeth))

d3 = HR3_e %>% select(-c(redcap_event_name, date, date_conserv, date.consent,date.stool.collection,
                           date_conserv,germe_c, data.row.hh.lab, data.row.hh.wash,date.enquete,
                           diametr_ampici_amoxici, diametr_amoxi_acid_clavu,diametr_piepera_tazobac,                    
                           diametr_cetriax_or_cefota,diametr_cefepime, diametr_meropenem, diametr_ertapenem,                          
                           diametr_gentamycine, diametr_amykacine, diametr_ciprofloxacine, diametr_pfloxacine,                         
                           diametr_sulfamet_trimeth))

HR_e_total_wide = full_join(d0, d1, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""))
# Match again with household survey to get all variables again
HR_e_total_wide = left_join(HR_e_total_wide, HR_e_total, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""),
                            multiple="first")

HR_e_total_wide = full_join(HR_e_total_wide, d2, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""))
# Match again with household survey to get all variables again
HR_e_total_wide = left_join(HR_e_total_wide, HR_e_total, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""),
                            multiple="first")

HR_e_total_wide = full_join(HR_e_total_wide, d3, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""))

# Match again with household survey to get all variables again
HR_e_total_wide = left_join(HR_e_total_wide, HR_e_total, by= c("menage_id_member", "menage_id", "village"), suffix = c("", ""),
                            multiple="first")

HR_e_total_wide = HR_e_total_wide %>% select(-c(redcap_event_name, date, date_conserv, date.consent,date.stool.collection,
                         date_conserv,germe_c, data.row.hh.lab, data.row.hh.wash,date.enquete,
                         diametr_ampici_amoxici, diametr_amoxi_acid_clavu,diametr_piepera_tazobac,                    
                         diametr_cetriax_or_cefota,diametr_cefepime, diametr_meropenem, diametr_ertapenem,                          
                         diametr_gentamycine, diametr_amykacine, diametr_ciprofloxacine, diametr_pfloxacine,                         
                         diametr_sulfamet_trimeth, esble))


table(HR_e_total_wide$r0.esble)
table(HR0_e$r0.esble)

table(HR_e_total_wide$r1.esble)
table(HR1_e$r1.esble)

table(HR_e_total_wide$r2.esble)
table(HR2_e$r2.esble)

table(HR_e_total_wide$r3.esble)
table(HR3_e$r3.esble)

length(unique(HR_e_total_wide$menage_id_member))

# Export dataset
write.csv(HR_e_total_wide, paste0(DirectoryDataOut, "./linked_final/bf_hh_stool_esble_r0123_wide.csv")) 

#################################################################################
# MERGING OF DATASETS FINALSED!
#################################################################################


#################################################################################
# Number of individuals we have complete records of
sum(complete.cases(HR_e_total_wide[, c("r0.esble", "r1.esble", "r2.esble", "r3.esble")])) # 803 individuals


#################################################################################
# CREATE DATASET FOR MSM MODEL
#################################################################################


rm(list=ls())

# Load functions
source("./Scripts/functions/multiplot.R")

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, summarytools, msm,
               survey,tsModel, mgcv)


df = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_esble_r0123_wide.csv")
dfl = read.csv("./Data/BF/clean/linked_final/bf_hh_stool_esble_r0123.csv")

df0 = df %>% filter(!is.na(r0.esble))
dfc = df %>% filter(complete.cases(r1.esble, r2.esble, r3.esble)) # 775 complete cases; however, later
# becomes clear there are individuals with unlikely collection dates. 

# Antibiotic use
abx = read.csv("./Data/BF/clean/watch_acute.csv")
vabx = read.csv("./Data/BF/Raw/village_abx_vs_hh.csv", sep=";") %>% select(village.cluster,village_name)
# Two clusters did not have a pharmacy so no abxuse
abx = left_join(abx,vabx)
abxn = abx %>% filter(site =="Nanoro")

#############################################################
# PREAMBLE
#############################################################

########################
# ABX USE
########################
nsurveys_by_providertype_by_village <- abxn %>%
  group_by(village_name, providertype, round) %>%
  summarise(n_surveys = n(), hcu = mean(hcu))
nsurveys_by_providertype_by_village

pop_by_village = abxn %>% select(village_name,pop_villagecluster) %>%
  filter(!duplicated(village_name))

poststratificationweight <- merge(nsurveys_by_providertype_by_village, pop_by_village, by = "village_name")
poststratificationweight$poststratweight <- (poststratificationweight$pop_villagecluster*poststratificationweight$hcu)/poststratificationweight$n_surveys

poststratificationweight_village = poststratificationweight %>% group_by(village_name,round)%>%
  summarise(weights = sum(poststratweight))

abxn = left_join(abxn, poststratificationweight_village, by=c("village_name", "round"))

# proportion estimates by village and pre/post
surveydesign <- svydesign(id = ~village_name, data = abxn, nest = TRUE)
watchclusterprop <- svyby(~watch, by = ~village_name + round, design = surveydesign, FUN = svymean, na.rm = TRUE)

poststratweightedsurveydesign <- svydesign(id = ~village_name, weights = ~weights, data = abxn, nest = TRUE)


poststratweightedwatchprop <- svyby(~watch, by = ~village_name + round, design = poststratweightedsurveydesign, FUN = svymean, na.rm = TRUE)
poststratweightedabxprop <- svyby(~antibiotic, by = ~village_name + round, design = poststratweightedsurveydesign, FUN = svymean, na.rm = TRUE)

abxn_vw = left_join(poststratweightedwatchprop%>%select(-c(se)),
                    poststratweightedabxprop %>% select(-c(se)))
abxn_vw0 = abxn_vw %>% filter(round=="baseline")

rm(nsurveys_by_providertype_by_village, pop_by_village,poststratweightedsurveydesign,
   poststratificationweight_village, poststratweightedabxprop,poststratweightedwatchprop,
   poststratificationweight, watchclusterprop,surveydesign, abxn,abxn_vw,abx,vabx)

# link with dataframe
dfl = left_join(dfl, abxn_vw0, by="village_name", multiple="first")
df = left_join(df, abxn_vw0, by="village_name", multiple="first")

###############################################################
# DATASET FOR MSM MODEL
###############################################################

# Not all individuals have a stool collection date (i.e. those with no ESBL)
# However, using the consent date for the ones with 
# Check dates
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

ggplot(dfl, aes(x=date.check.collect, group=redcap_event_name)) + geom_histogram() +
  facet_wrap(~ redcap_event_name)

ggplot(dfl, aes(x=date.check.use, group=redcap_event_name)) + geom_histogram() +
  facet_wrap(~ redcap_event_name)

table(dfl$date.check.use, useNA="always") # Date variable can not be used for the time, as has any NAs 

# We can make a minimum date per round so to ensure we only include those with a correct date per round
r0min = min(dfl$date.use[dfl$redcap_event_name=="round_0_arm_1"], na.rm=T) # 
r1min = min(dfl$date.use[dfl$redcap_event_name=="round_1_arm_1"], na.rm=T) # For round 1 not correct as same as round 0
# Should be R0+~90 days
r1min = r0min+90
r2min = min(dfl$date.use[dfl$redcap_event_name=="round_2_arm_1"], na.rm=T)
r3min = min(dfl$date.use[dfl$redcap_event_name=="round_3_arm_1"], na.rm=T)



# START CREATING DATASET BY SELECTING COMPLETE CASES + CASES WITH CORRECT STOOL.COLLECTION/CONSENT.DATE
##############################################################
dfls = dfl %>% filter(time!=0) %>%
  mutate(
    rc = ifelse(time==1 & date.use>=r1min| 
                  time==2 & date.use>=r2min| 
                  time==3 & date.use>=r3min, 1,0)
  )
table(dfls$rc, useNA="always")

dfls = dfls %>% filter(rc!=0) %>% select(-c(rc))

# Which are the individuals with observation in round 1, 2, and 3
complete.cases.id = dfls %>% group_by(menage_id_member) %>%
  summarise(nround = length(menage_id_member)) %>% 
  filter(as.numeric(nround)==3)
complete.cases.id=unique(complete.cases.id$menage_id_member)

data.frame(table(dfls$menage_id_member[dfls$menage_id_member%in%complete.cases.id]))

dfls = dfls %>% filter(menage_id_member %in% complete.cases.id) # 785

d = dfl %>% filter(time==0) %>% group_by(menage_id_member) %>%
  summarise(start0 = min(date.use))

# Create start when round 1 is taken as baseline
d1 = dfls %>% filter(time!=0) %>% group_by(menage_id_member) %>%
  summarise(start1 = min(date.use))

dfls = left_join(dfls,d, by="menage_id_member")
dfls = left_join(dfls,d1, by="menage_id_member")

# Create a variable that calculates the number of days between the collection date and first round
dfls$time.start0 =  dfls$date.use - dfls$start0
dfls$time.start1 =  dfls$date.use - dfls$start1

dfls = dfls[order(dfls$menage_id_member, dfls$time),]


ggplot(dfls, aes(x=time.start1, group=redcap_event_name)) + geom_histogram() +
  facet_wrap(~ redcap_event_name)

# remove those with wrong dates
r2w = dfls[which(dfls$redcap_event_name=="round_2_arm_1" & as.numeric(dfls$time.start1)==0),]
r1w = dfls[which(dfls$redcap_event_name=="round_1_arm_1" & dfls$start1<"2023-01-01"),]

dfls = dfls[-c(which(dfls$redcap_event_name=="round_2_arm_1" & as.numeric(dfls$time.start1)==0)),]
#dfls = dfls[-c(which(dfls$redcap_event_name=="round_0_arm_1" & as.numeric(dfls$time.start)>100)),]

ggplot(dfls, aes(x=time.start1, group=redcap_event_name)) + geom_histogram() +
  facet_wrap(~ redcap_event_name)

dfls = dfls[order(dfls$menage_id_member, dfls$time),]


# start dates per round
start.round = dfls %>% group_by(redcap_event_name) %>%
  summarise(
    start = min(date.use, na.rm=T)
  )
start.round
start.round1 = start.round$start[2] 

# Calculate number of positives per household to add as covariate in the model
hh_pos = dfl %>% filter(menage_id_member%in%complete.cases.id) %>% 
  group_by(time, menage_id) %>%
  summarise(nesbl = sum(esble))
hh_pos = hh_pos[order(hh_pos$menage_id,hh_pos$time),]
hh_pos$nesbl_tminus = lag(hh_pos$nesbl)
hh_pos$nesbl_tminus[hh_pos$time==0] = NA


dfls = left_join(dfls,hh_pos, by=c("menage_id", "time"))
#dfls$nesbl_tminus = lag(dfls$nesbl)

# Add number of tested individual per round per household as denominator for prevalence
d = dfls %>% group_by(menage_id, time) %>%
  summarise(
    n.tested = length(menage_id_member)
  )

d
dfls = left_join(dfls,d, multiple="first")


# Create dataset for analyses
dfls = dfls %>% mutate(
  days = as.numeric(date.use - start1) # use individual first stool collection date, as each village had different start of intervention
) %>% filter(redcap_event_name!="round_0_arm_1" & menage_id_member%in%complete.cases.id)
length(unique(dfls$menage_id_member))
table(dfls$menage_id_member, dfls$time) # ALL IDS HAVE THREE OBSERVATIONS

dfls = dfls[order(dfls$menage_id_member,dfls$time),]
dfls$check.rounds = dfls$time.start1 - lag(dfls$time.start1)
dfls$check.rounds[dfls$redcap_event_name=="round_1_arm_1"] = NA
hist(as.numeric(dfls$check.rounds))

# Remove all those with a negative check.round as then the date.collection in next round is before the previous
table(dfls$check.rounds)
out = which(dfls$check.rounds<1)
dfls = dfls[-c(out),] %>% select(-c(date.check,check.rounds))


# How many unique individuals
length(unique(dfls$menage_id_member)) # 785
dfls = dfls[order(dfls$menage_id_member, dfls$time),]

# Dataset for msm model
dfls = dfls %>% 
  mutate(esble = as.numeric(esble)+1) %>%
  select(-c(date.check.collect, date_conserv, intervention_text,X)) 

reorder_col = c("time", "redcap_event_name","time.start0","time.start1","start0","start1", "days","date.use","month","rainy", "menage_id_member", "menage_id",
                "village_name","intervention.text", "age", "agegr10", "sexe", "n.householdmember","esble","nesbl","nesbl_tminus","n.tested",
                "salm", "ast_done",
                "date","date.stool.collection", "date.consent", "date.check")
nincl = names(dfls)[!names(dfls)%in%reorder_col]
desired_order = c(reorder_col,nincl)
desired_order = desired_order[desired_order%in%names(dfls)]
dfls =  dfls[,desired_order]

# positive per round
dfls %>% group_by(time) %>%
  summarise(
    esblsum = sum(esble-1),
    #acquisitionsum = sum(acquisition,na.rm=T),
    n = length(unique(menage_id_member)),
    prev = esblsum/n
  )

# Also create a dataset with baseline 0 round included
d = dfl %>% filter(menage_id_member%in%complete.cases.id & redcap_event_name=="round_0_arm_1")

d = left_join(d, dfls%>%select(menage_id_member,start0, start1, time.start0,time.start1), multiple="first")
d = d %>% select(-c(X)) %>%
  mutate(days = as.numeric(date.use - start1))

hh_pos = d %>% group_by(time,menage_id) %>%
  summarise(nesbl = sum(esble)
  )
hh_pos = hh_pos[order(hh_pos$menage_id, hh_pos$time),]
hh_pos$nesbl_tminus = NA

d = left_join(d,hh_pos, by=c("menage_id", "time"))
d$esble = d$esble+1

# Add number of tested individual per round per household as denominator for prevalence
ds = d %>%filter(time==0)%>% group_by(menage_id, time) %>%
  summarise(
    n.tested = length(menage_id_member)
  )

d = left_join(d,ds, multiple="first")

table(d$time,d$esble)

names(d)
desired_order

d = d[, desired_order]
dfls0 = rbind(d,dfls)
dfls0 = dfls0[order(dfls0$menage_id_member, dfls0$time),]
table(dfls0$time,dfls0$esble)

dfls0 %>% group_by(time) %>%
  summarise(
    esblsum = sum(esble-1),
    n = length(unique(menage_id_member)),
    prev = esblsum/n
  )

# check the days variable
dfls0 %>% group_by(time) %>%
  ggplot(., aes(x=days)) + geom_histogram() +
  facet_wrap(.~redcap_event_name)

dfls0 %>% group_by(time) %>%
  summarise(median=median(days, na.rm=T),
            q1 = quantile(days,probs=0.25, na.rm=T),
            q3 = quantile(days,probs=0.75,na.rm=T))

# Create new day variable for when fitting MSM with R0, so counting days from R0
dfls0 = dfls0[order(dfls0$menage_id_member,dfls0$time),]
dfls0 = dfls0 %>% group_by(menage_id_member)%>%
  mutate(days2 = lag(days),
         days2 = ifelse(is.na(days2), 0,
                        ifelse(time==1,(days2*-1),
                               ifelse(time%in%c(2,3), days, days2))),
         daysl = lag(days2)) %>%
  group_by(menage_id_member) %>%
  mutate(days2 = ifelse(time==2, days2+daysl,days2),
         daysl = lag(days2))
 
dfls0 = dfls0 %>% group_by(menage_id_member) %>%
   mutate(days2=ifelse(time==3, days2+daysl,days2)
   )    

# Add acquisition and decolonisation variable 
da_lag = dfls0%>% select(time,menage_id_member,esble) %>% group_by(menage_id_member) %>%
  mutate(lagesbl = lag(esble),
         lag2esbl = lag(esble,2))

#da_lag$time = dfls0$time
#da_lag$menage_id_member =dfls0$menage_id_member
#da_lag = da_lag[order(da_lag$menage_id_member,da_lag$time),]
#da_lag$lagesbl[da_lag$time==0] = NA

da_lag$check=da_lag$esble-da_lag$lagesbl
da_lag$check2=da_lag$esble-da_lag$lag2esbl

da_lag$acquisition = ifelse(da_lag$check==1, 1,
                            ifelse(da_lag$lagesbl==2,NA,0))
da_lag$decolonisation = ifelse(da_lag$check==-1, 1,
                               ifelse(da_lag$lagesbl==1, NA,0))
da_lag$acquisition2 = ifelse(da_lag$time==3, da_lag$esble-1, da_lag$acquisition) # For sensitivity analyses, i.e. all individuals at R3 were at risk at R2, assuming so as decolonisation time is 4M on average (paper Bootsma)
da_lag$acquisition3 = ifelse(da_lag$time==2 & da_lag$check==1 & da_lag$check2==1, 0, da_lag$acquisition) # For sensitivity analyses, where all those negative at both R0 and R1 were assumed at risk for R2 (instead of only negative at R1)

da_lag$decolonisation2 = ifelse(da_lag$time==3, da_lag$esble-1, da_lag$decolonisation) # For sensitivity analyses, i.e. all individuals at R3 were at risk at R2, assuming so as decolonisation time is 4M on average (paper Bootsma)
da_lag$decolonisation3 = ifelse(da_lag$time==2 & da_lag$check==-1 & da_lag$check2==-1, 0, da_lag$decolonisation) # For sensitivity analyses, where all those negative at both R0 and R1 were assumed at risk for R2 (instead of only negative at R1)


da_lag = da_lag %>% select(c(time,menage_id_member,acquisition,decolonisation, acquisition2, acquisition3, decolonisation2, decolonisation3))
dfls = left_join(dfls,da_lag)
dfls0 = left_join(dfls0,da_lag)

hh_acq = dfls0 %>% group_by(time,menage_id) %>%
  summarise(nacquisition = sum(acquisition, na.rm=T),
            nacquisition2 = sum(acquisition2, na.rm=T),
            nacquisition3 = sum(acquisition3, na.rm=T),
            ndecolonisation = sum(decolonisation, na.rm=T),
            ndecolonisation2 = sum(decolonisation2, na.rm=T),
            ndecolonisation3 = sum(decolonisation3, na.rm=T)
  )


dfls = left_join(dfls,hh_acq)
dfls0 = left_join(dfls0,hh_acq)

# Acquisitions per round
pa = dfls0 %>% group_by(time) %>%
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

reorder_col = c("time", "redcap_event_name","time.start0","time.start1","start0","start1", "days","days2", "date.use","month", "rainy", "menage_id_member", "menage_id",
                "village_name","intervention.text", "age", "agegr10", "sexe", "n.householdmember","esble","nesbl","nesbl_tminus","n.tested",
                "acquisition","acquisition2", "acquisition3", "nacquisition","nacquisition2","nacquisition3", "decolonisation","decolonisation2", "decolonisation3",
                "ndecolonisation", "ndecolonisation2", "ndecolonisation3",
                "salm", "ast_done",
                "date","date.stool.collection", "date.consent", "date.check")
nincl = names(dfls)[!names(dfls)%in%reorder_col]
desired_order = c(reorder_col,nincl)
dfls0 = dfls0[,desired_order[desired_order!="date.check"]]

desired_order = desired_order[desired_order%in%names(dfls)]
dfls =  dfls[,desired_order]



###############################################################
# DATASET FOR LOGISTIC REGRESSION MODEL
###############################################################

# Combine categories
dfa = df %>%
  mutate(
    q2.main.water.source.dry.c = ifelse(q2.main.water.source.dry %in% c("Tap house", "Tap concession", "Tap public/fountain", "bagged water"),
                                        "tap/bottled",
                                        ifelse(q2.main.water.source.dry %in% c("rainwater", "surface water (ponds, dams,rivers,lakes,pits,irrigation canals)"),
                                               "other",q2.main.water.source.dry)),
    q3.main.water.source.rainy.c = ifelse(q3.main.water.source.rainy %in% c("Tap house", "Tap concession", "Tap public/fountain", "bagged water"),
                                          "tap/bottled",
                                          ifelse(q3.main.water.source.rainy %in% c("rainwater", "surface water (ponds, dams,rivers,lakes,pits,irrigation canals)"),
                                                 "other",q3.main.water.source.rainy)),
    q6.treatment.water.c = ifelse(q6.treatment.water %in%  c("Yes,boiling","Yes,cholinate/add desinfectant",
                                                             "Yes, filter with cloth", "Yes, filter with filter","Yes, Solar desinfection (in the sun)",
                                                             "Yes, decant"), "Yes", "No"),  
    q9.shared.toilet.c = ifelse(q9.shared.toilet %in% c("Yes, public","Yes, other households (non-public)"), "Yes", "No")
  ) %>% filter(menage_id_member %in%complete.cases.id)

# Add number of tested individual per round per household as denominator for prevalence
d = dfa %>% group_by(menage_id) %>%
  summarise(
    r0.ntested = length(menage_id_member[!is.na(r0.esble)]),
    r.ntested = length(menage_id_member),
  )
d

dfa = left_join(dfa, d, multiple="first")
reorder_col = c("redcap_event_name","menage_id_member", "menage_id",
                "village_name","intervention.text", "age", "agegr10", "sexe", "n.householdmember","r0.ntested", "r.ntested",
                "r0.esble","r1.esble","r2.esble","r3.esble","r0.salm","r1.salm","r2.salm","r3.salm",
                "date.use","date","date.stool.collection", "date.consent")
nincl = names(dfa)[!names(dfa)%in%reorder_col]
desired_order = c(reorder_col,nincl)
desired_order = desired_order[desired_order%in%names(dfa)]
dfa =  dfa[,desired_order]

dfls = dfls %>% select(-c(time.start0,time.start1))
dfls0 = dfls0 %>% select(-c(time.start0,time.start1))

table(dfa$r1.esble, useNA="always")
table(dfa$r2.esble, useNA="always")
table(dfa$r3.esble, useNA="always")

# Export dataset
#save(df, file="./Data/BF/clean/use_in_analyses/bf_esbl_wide_all.rda")
save(dfa, file="./Data/BF/clean/use_in_analyses/bf_esbl_wide_completecases.rda")
#load("./Data/BF/clean/use_in_analyses/bf_esbl_wide.rda")

save(dfls0, file="./Data/BF/clean/use_in_analyses/bf_esbl0123_long_completecases.rda")
save(dfls, file="./Data/BF/clean/use_in_analyses/bf_esbl123_long_completecases.rda")
#save(dfl, file="./Data/BF/clean/use_in_analyses/bf_esbl_long_all.rda")

rm(r1w, r2w,hh_pos, start.round,d1, abxn_vw0, dfc)


