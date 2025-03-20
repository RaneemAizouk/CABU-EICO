#################################################################
# CLEAN HOUSEHOLD DATA
#################################################################
# This code is cleaning the household WASH survey data 

# 16 March 2025
# Author: Esther van Kleef

# The CABU-EICO household data contains questions from three surveys
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

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, msm, skimr, wesanderson, patchwork)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

# Villages (that are the clusters) of CABU-EICO
villages = readxl::read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","ajouter")

village_size = readxl::read_xlsx(paste0(DirectoryData, "/Sélection et randomisation grappes CABU_B Nanoro.xlsx"), sheet=1)
names(village_size)
village_size = village_size %>% select(c(`Nom des villages`, ménages, `population active (Mar 2023 update)`))%>%
  rename(
    village_name = `Nom des villages`,
    n.households =  ménages,
    n.population = `population active (Mar 2023 update)`
  ) %>%
  filter(!is.na(village_name))
head(village_size)

villages = left_join(villages, village_size, by="village_name")

# Lab data
#car_bf = read.csv(paste0(DirectoryData, "/householdsurvey/CABUBPortageAsymptom_DATA_2024-04-17_1527.csv"), sep=";")
car_bf = readxl::read_xlsx(paste0(DirectoryData, "/householdsurvey/CABU_WP4_R0-R3.xlsx"))

# Lab ids vs household ids
hh_lab_ids =  readxl::read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))
names(hh_lab_ids)
names(car_bf)
names(hh_lab_ids) = c("household", "menage_id", "bras")

# Add variables village and household to lab data
car_bf$village = substr(car_bf$record_id, start = 1, stop = 2)
car_bf$household = str_extract(car_bf$record_id, "[^-]+")
car_bf = merge(car_bf, villages, by="village")

# SES data
ses = readxl::read_xlsx(paste0(DirectoryData,"/hdss_ses/SES_CABU-EICO hh.xlsx"))
names(ses)
ses = ses %>% select(ménage_id,quintiles) %>%
  rename(
    menage_id = ménage_id,
    ses.quintile = quintiles
  )

ses_pella = readxl::read_xlsx(paste0(DirectoryData,"/hdss_ses/SES_Pella.xlsx")) 
ses_pella = ses_pella %>%
  select(c(menage_id, quintile)) %>%
  rename(
    ses.quintile = quintile
  )

pellahh = unique(ses_pella$menage_id)
# Need to get 
unique(ses_pella$menage_id)

ses_nopella = ses %>%filter(!menage_id%in%pellahh)

d_ses = rbind(ses_nopella, ses_pella)



d_ses = d_ses %>% filter(menage_id%in% unique(ses$menage_id)) %>%
  mutate(
    ses.quintile = factor(ses.quintile, levels=c("lowest", "second", "third", "fourth", "highest", NA)))

table(d_ses$ses.quintile, useNA="always") # 27 still not linked
table(ses$ses.quintile, useNA="always")

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
table(car_bf$esbl_pos, useNA="always") # These are the ESBL positive patients based on cetriax_or_cefota, 2842
table(car_bf$testesbl) # These are the ESBL positive patients based on esbl_pos, 2842
table(car_bf$esbl_pos==1 & car_bf$testesbl==1) # difference; we decided to ignore these differences

# Remove individuals with diametr_cetriax_or_cefota = NA
car_bf = car_bf %>% filter(!is.na(diametr_cetriax_or_cefota)) # 4 removed
table(car_bf$esbl_pos, useNA="always") # These are the ESBL positive patients based on cetriax_or_cefota, 2507; Update: 16 March 2025 this is 2842

names(car_bf)

# Household data
#hh_bf = readxl::read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUBWP4_DATA_2024-04-17_1528.xlsx"))
hh_bf = readxl::read_xlsx(paste0(DirectoryData, "/householdsurvey/CABUWASH_DATA_2024-06-27_R0&R3.xlsx"))


#-------------------------------------------------------------------------------
# HOUSEHOLD SURVEY - WASH INDICATORS
#-------------------------------------------------------------------------------

# BURKINA FASO - ALL ROUNDS
################################
d = hh_bf %>% mutate(
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
  agegr = cut(age,
              breaks = c(0, 5, 18, 50, Inf),
              labels = c("0-4", "5-17", "18-49", "50+"),
              levels = c("0-4", "5-17", "18-49", "50+"),
              include.lowest = TRUE),
  nmbre_personne_menage = as.numeric(nmbre_personne_menage),
  nbre_enf_0_5ans = as.numeric(nbre_enf_0_5ans),
  nbre_menage_conc = as.numeric(nbre_menage_conc), 
  date_enquete = as.Date(date_enquete, format="%Y-%m-%d"),
  date_consentement = as.Date(date_consentement, format="%Y-%m-%d"),
  date_recuperation_selle = as.Date(date_recuperation_selle, format="%Y-%m-%d"),
  # q1_diarrhee_prevenu___1 = factor(q1_diarrhee_prevenu___1, levels=c(0,1), labels=c("No","Yes")),        
  # q1_diarrhee_prevenu___2 = factor(q1_diarrhee_prevenu___2, levels=c(0,1), labels=c("No","Yes")),
  # q1_diarrhee_prevenu___3 = factor(q1_diarrhee_prevenu___3, levels=c(0,1), labels=c("No","Yes")),        
  # q1_diarrhee_prevenu___4 = factor(q1_diarrhee_prevenu___4, levels=c(0,1), labels=c("No","Yes")),
  # q1_diarrhee_prevenu___5 = factor(q1_diarrhee_prevenu___5, levels=c(0,1), labels=c("No","Yes")),         
  # q1_diarrhee_prevenu___6 = factor(q1_diarrhee_prevenu___6, levels=c(0,1), labels=c("No","Yes")),
  # q1_diarrhee_prevenu___7 = factor(q1_diarrhee_prevenu___7, levels=c(0,1), labels=c("No","Yes")),         
  # q1_diarrhee_prevenu___8 = factor(q1_diarrhee_prevenu___8, levels=c(0,1), labels=c("No","Yes")),
  # q1_diarrhee_prevenu___9 = factor(q1_diarrhee_prevenu___9, levels=c(0,1), labels=c("No","Yes")), 
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
                                                                                             "unimproved well (unprotected)",
                                                                                             "rainwater",
                                                                                             "surface water (ponds, dams,rivers,lakes,pits,irrigation canals)",
                                                                                             "bagged water", "bottled water")),  
  q4_bidon_stock = factor(q4_bidon_stock, levels = c(1:3), labels=c("Yes, cans", "Yes, only one large tank", "No")),
  # q5a_bidon_ferme_rempli = factor(q5a_bidon_ferme_rempli, levels=c(1:2), labels=c("Yes", "No")),
  # q5b_bidon_ferme_vide = factor(q5b_bidon_ferme_vide, levels=c(1:2), labels=c("Yes", "No")),           
  q5c_bidon_nettoye = factor(q5c_bidon_nettoye, levels=c(1,2,3,4,6), labels=c("Yes, with soap", "Yes, with water but no soap","Yes, boiled", "No", "Other")),
  q6_traite_eau = factor(q6_traite_eau, levels = c(1:7), labels= c("No", "Yes,boiling", "Yes,cholinate/add desinfectant", "Yes, filter with cloth",
                                                                   "Yes, filter with filter", "Yes, Solar desinfection (in the sun)","Yes, decant")),               
  q7_type_inst_sanitaire = factor(q7_type_inst_sanitaire, levels = c(3,4,5), labels = c("pit latrine with slab", "pit latrine without slab", "open defecation")),
  # q8_autr_lieu_defecation___1 = factor(q8_autr_lieu_defecation___1, levels = c(0,1), labels = c("No","Yes")),      
  # q8_autr_lieu_defecation___2 = factor(q8_autr_lieu_defecation___2, levels = c(0,1), labels = c("No","Yes")),
  # q8_autr_lieu_defecation___3 = factor(q8_autr_lieu_defecation___3, levels = c(0,1), labels = c("No","Yes")),   
  # q8_autr_lieu_defecation___4 = factor(q8_autr_lieu_defecation___4, levels = c(0,1), labels = c("No","Yes")),
  # q8_autr_lieu_defecation___5 = factor(q8_autr_lieu_defecation___5, levels = c(0,1), labels = c("No","Yes")),   
  # q8_autr_lieu_defecation___6 = factor(q8_autr_lieu_defecation___6, levels = c(0,1), labels = c("No","Yes")),
  # q8_autr_lieu_defecation___7 = factor(q8_autr_lieu_defecation___7, levels = c(0,1), labels = c("No","Yes")),  
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
  q20_excrement_animaux = factor(q20_excrement_animaux, levels=c(1:3), labels=c("Yes", "No", NA)),
  # q18_animaux_interieur___1 = factor(q18_animaux_interieur___1, levels = c(0,1), labels = c("No", "Yes")),
  # q18_animaux_interieur___2 = factor(q18_animaux_interieur___2, levels = c(0,1), labels = c("No", "Yes")),       
  # q18_animaux_interieur___3 = factor(q18_animaux_interieur___3, levels = c(0,1), labels = c("No", "Yes")),
  # q18_animaux_interieur___4 = factor(q18_animaux_interieur___4, levels = c(0,1), labels = c("No", "Yes")),       
  # q18_animaux_interieur___5 = factor(q18_animaux_interieur___5, levels = c(0,1), labels = c("No", "Yes")),
  # q18_animaux_interieur___6 = factor(q18_animaux_interieur___6, levels = c(0,1), labels = c("No", "Yes")),
  # q19_animaux_dehors___1 = factor(q19_animaux_dehors___1, levels = c(0,1), labels = c("No", "Yes")),        
  # q19_animaux_dehors___2 = factor(q19_animaux_dehors___2, levels = c(0,1), labels = c("No", "Yes")),
  # q19_animaux_dehors___3 = factor(q19_animaux_dehors___3, levels = c(0,1), labels = c("No", "Yes")),        
  # q19_animaux_dehors___4 = factor(q19_animaux_dehors___4, levels = c(0,1), labels = c("No", "Yes")),
  # q19_animaux_dehors___5 = factor(q19_animaux_dehors___5, levels = c(0,1), labels = c("No", "Yes")),         
  # q19_animaux_dehors___6 = factor(q19_animaux_dehors___6, levels = c(0,1), labels = c("No", "Yes")),
  # 
  # q20_excrement_animaux = factor(q20_excrement_animaux, levels = c(1,2,3), labels = c("Yes", "No", "Not possible to determine")),  
  # q21_animal_malade___1 = factor(q21_animal_malade___1, levels = c(0,1), labels = c("No", "Yes")),           
  # q21_animal_malade___2 = factor(q21_animal_malade___2, levels = c(0,1), labels = c("No", "Yes")),  
  # q21_animal_malade___3 = factor(q21_animal_malade___3, levels = c(0,1), labels = c("No", "Yes")),  
  # q21_animal_malade___4 = factor(q21_animal_malade___4, levels = c(0,1), labels = c("No", "Yes")),  
  # q21_animal_malade___5 = factor(q21_animal_malade___5, levels = c(0,1), labels = c("No", "Yes")),             
  # q21_animal_malade___6 = factor(q21_animal_malade___6, levels = c(0,1), labels = c("No", "Yes")),
  #eau_assainissement_hygine_complete = factor(eau_assainissement_hygine_complete, levels=c(0,2), labels=c("No", "Yes")),
  piam_q1 = as.numeric(piam_q1),
  piam_q2a = as.numeric(piam_q2a),
  piam_q2b = factor(piam_q2b, levels=c(1,2), labels=c("Male", "Female")),
  piam_q6a = as.numeric(piam_q6a),
  piam_q6b = as.numeric(piam_q6b),
  piam_q6c = as.numeric(piam_q6c),
  piam_q6d = as.numeric(piam_q6d),
  piam_q6e = as.numeric(piam_q6e),
  piam_q6f = as.numeric(piam_q6f),
  piam_q6g = as.numeric(piam_q6g),
  piam_q6h = as.numeric(piam_q6h),
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
  q6.traitmen.water.other = "q6_autre_traitmen_eau",
  q7.principle.defication= "q7_type_inst_sanitaire", 
  q7.principle.defication.other = "q7_autr_typ_ins_sanitair",
  q8.other.defecation.flush.toiled.septic  = "q8_autr_lieu_defecation___1", 
  q8.other.defecation.pit.latrine.ventilation = "q8_autr_lieu_defecation___2",
  q8.other.defecation.pit.latrine.slab = "q8_autr_lieu_defecation___3",   
  q8.other.defecation.pit.latrine.no.slab = "q8_autr_lieu_defecation___4",
  q8.other.defecation.open.defecation = "q8_autr_lieu_defecation___5",   
  q8.other.defecation.other = "q8_autr_lieu_defecation___6",
  q8.other.defecation.none = "q8_autr_lieu_defecation___7",
  q8.other.defecation.specified = "q8_autre_preciser",
  q9.shared.toilet= "q9_toilette_partagee",
  q10.n.shared.toilet= "q10_combien_partag",
  q11.toilet.last.cleaned= "q11_dernier_nettoyage",
  q12.disposal.child.stool= "q12_elimine_selle_enf",
  q12.disposal.child.other.specified = "q12_autre_preciser",
  q13.disposal.latrine.pit= "q13_vidange_toilette",
  q13.latrine.pit.other.specified = "q13_autre_preciser",
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
  q18_animal.inside.other.specified = "q18_autre_specifie",
  q19.animal.outside.cow = "q19_animaux_dehors___1",
  q19.animal.outside.sheep.goat = "q19_animaux_dehors___2",       
  q19.animal.outside.pig = "q19_animaux_dehors___3",
  q19.animal.outside.donkey.horse = "q19_animaux_dehors___4",       
  q19.animal.outside.chicken.goose.duck = "q19_animaux_dehors___5",
  q19.animal.outside.other = "q19_animaux_dehors___6",
  q19.animal.outside.other.specified = "q19_autre_specifie",
  q20.animal.excrement.floor= q20_excrement_animaux, 
  q21.when.animal.ill.treatment.with.vet = "q21_animal_malade___1",           
  q21.when.animal.ill.treatment.without.vet= "q21_animal_malade___2",  
  q21.when.animal.ill.sell.bucher = "q21_animal_malade___3",  
  q21.when.animal.ill.slaugter.eat.at.home= "q21_animal_malade___4",  
  q21.when.animal.ill.dies.burie.dispose = "q21_animal_malade___5",             
  q21.when_animal.ill.dies.eat.at.home= "q21_animal_malade___6",
  piam.q1.n.athome.when.survey = "piam_q1",
  piam.q2a.athome.age.random.selected.person = "piam_q2a",
  piam.q2a.athome.sex.random.selected.person = "piam_q2b",
  piam.q3.heard.about.trial.activities = "piam_q3",
  piam.q4.partic.atleastone.trial.activities = "piam_q4",
  piam.q5.perc.key.topic.hwashing.use.soap = "piam_q5___1",
  piam.q5.perc.key.topic.hwashing.when = "piam_q5___2",
  piam.q5.perc.key.topic.abx.what.and.role = "piam_q5___3",
  piam.q5.perc.key.topic.abx.correct.use= "piam_q5___4",
  piam.q5.perc.key.topic.seek.care.no.self.med = "piam_q5___5",
  piam.q5.perc.key.topic.good.sani.latr.wastw= "piam_q5___6",
  piam.q5.perc.key.topic.drink.treat.or.borehole.water= "piam_q5___7",
  piam.q5.perc.key.topic.protect.store.water.clean.container= "piam_q5___8",
  piam.q5.perceived.key.topic.handwashing= "piam_q6a",
  piam.q6.freq.partic.hwashing.public = "piam_q6b",
  piam.q6.freq.partic.abx.use.public = "piam_q6c",
  piam.q6.freq.partic.hwashing.athome = "piam_q6d",
  piam.q6.freq.partic.abx.use.athome= "piam_q6e",
  piam.q6.freq.partic.abx.use.adhoc= "piam_q6f",
  piam.q6.freq.partic.defication.adhoc= "piam_q6g",
  piam.q6.freq.partic.defication.public= "piam_q6h"
)

table(d$agegr10)
table(d$agegr)
names(d)

# Link village and cluster data to household survey
hh_bf$data_row = c(1:nrow(hh_bf))
hh_bf$village = substr(hh_bf$menage_id, start = 1, stop = 2)
d$data_row = c(1:nrow(d)) # This variable we can use for identifying back those individuals that had a sample taken
d$village = substr(d$menage_id, start = 1, stop = 2)
d = merge(d, villages, by="village") 
names(d)

#View(hh_bf[which(!hh_bf$data_row %in% d$data_row),]) # Households that start with SB are not linked. Need to check which village these belong to. 20 March 2025: sent email to Franck.
#unique(villages$village)

names(d) = gsub("_",".",names(d))
d = d %>%
  rename(menage_id = "menage.id",
         village_name = "village.name",
         redcap_event_name = "redcap.event.name")


d = left_join(d, d_ses, by=c("menage_id"))
names(d)
#-------------------------------------------------------------------------------
# DEFINE BINARY WASH INDICATORS
#-------------------------------------------------------------------------------
# 1) Access to/use of safe drinking water
# 2) Cleaning drinking water storage --> if q4 = Yes, cans; Yes, only one large tank, then q5a, q5b and q5c are filled out. Turns out only 3 have no for this question so can just use q5c
# 3) Correct hand washing --> q14 soap available + q15 hand washing after defecation 
# 4) Access/use of improved sanitation facility
# 5) Livestock animal access to the house
# 6) Animal excrement on the floor


# Some background
#---------------------------------------------------------

# Correct handwashing
# Opted for the less strict approach to only define handwashing with soap, detergent, and ash as correct.
# In certain places, including DRCongo, hand washing with powder/liquid is common place and considered correct.
# Handwashing with ash is disputable (upon a quick search, this Cochrane review (Cochrane did not find conclusive evidence),
# but this is promoted under extreme humanitarian conditions (including south sudan), 
# probably less effective than soap. Quickly searching online, this RCT suggests the same. https://bmjopen.bmj.com/content/12/5/e056411
# would suggest we consider combining all forms of handwashing vs not.
# Then in a separate sensitivity analyses we could have handwashing with soap vs all other categories. 
# Why so many handwashing missing?
table(d$q15.handwashing.defecation, useNA="always") # Q15 was only asked if q14 was answered with yes
table(d$q14.handwashing.product, useNA="always") 
table(d$q14.handwashing.product,d$q15.handwashing.defecation, useNA="always")


# Improved sanitation
#1.	Divide improved sanitation in improved vs unimproved. Then for unimproved
# a.	First check Q9 if answer to shared toilet = yes, if yes, then also classify as unimproved
# b.	Pit latrine without slab counts officially as unimproved, not improved (although in practice the difference to exposure risk is likely limited).
# In contrast to the formal classification, one argumentation could be that what matters most is whether there is some form of disposal vs open defecation.
# However, for now, keep to strick definition
table(d$q14.handwashing.product)
table(d$q15.handwashing.defecation)
table(d$q9.shared.toilet)

#dt = rbind(d_wash, d_wash_r2) 

d <- d %>% 
  mutate(main.drinking.water.dry.binary = case_when(
    q2.main.water.source.dry %in%c("Tap house", "Tap concession","Tap public/fountain",
                                   "Borehole", "improved well (protected)","rainwater","bottled water") ~ "Improved",
    q2.main.water.source.dry %in%c("surface water (ponds, dams,rivers,lakes,pits,irrigation canals)", "unimproved well (unprotected)",
                                   "bagged water") ~ "Unimproved",
    TRUE ~ NA),
    main.drinking.water.rainy.binary = case_when(
      q3.main.water.source.rainy %in%c("Tap house", "Tap concession","Tap public/fountain",
                                       "Borehole", "improved well (protected)","rainwater","bottled water") ~ "Improved",
      q3.main.water.source.rainy %in%c("surface water (ponds, dams,rivers,lakes,pits,irrigation canals)", "unimproved well (unprotected)",
                                       "bagged water") ~ "Unimproved",
      TRUE ~ NA),
    cleaning.water.storage.binary = case_when(
      q5c.cans.cleaned.before.reuse %in% c("Yes, with soap","Yes, with water but no soap","Yes, boiled","Other") ~ "Yes",
      q5c.cans.cleaned.before.reuse %in% c("No") ~ "No",
      TRUE ~ NA 
    ),
    correct.handwashing.binary = case_when(
      q15.handwashing.defecation %in% c("Yes, always","Yes, often") & q14.handwashing.product %in% c("Yes, soap", "Yes, detergent","Yes, ash/mud/sand") ~ "Yes",
      q15.handwashing.defecation %in% c("No or rarerly","Not sure") | q15.handwashing.defecation %in% c("Yes, always","Yes, often")
      & q14.handwashing.product %in% c("No, none available", "No, available but not used")| !q14.handwashing.product %in% c("Yes, soap", "Yes, detergent","Yes, ash/mud/sand") ~ "No",
      TRUE ~ NA
    ),
    improved.sanitation.binary = case_when(
      q9.shared.toilet %in% c("Yes, public")|
        q7.principle.defication %in% c("pit latrine without slab","open defecation") ~ "No",
      q7.principle.defication %in% c("pit latrine with slab") &q9.shared.toilet %in% c("No","Yes, other households (non-public)") ~ "Yes",
      TRUE ~ NA
    ),
    livestock.access.house.binary = case_when(
      q17.animals.around.household %in%c("Yes, inside and outside") ~ "Yes",
      q17.animals.around.household %in%c("Yes, outside next to house","Yes, outside in demarked area", "No") ~ "No",
      TRUE ~ NA
    ),
    animal.excrement.floor.binary = factor(q20.animal.excrement.floor, levels=c("No","Yes"))
  )

table(d$main.drinking.water.dry.binary, useNA="always")
table(d$main.drinking.water.rainy.binary, useNA="always")
table(d$cleaning.water.storage.binary, useNA="always")
table(d$correct.handwashing.binary, useNA="always")
d$correct.handwashing.binary = ifelse(is.na(d$q14.handwashing.product), NA, d$correct.handwashing.binary)

table(d$improved.sanitation.binary, useNA="always")
table(d$livestock.access.house.binary, useNA="always")
table(d$animal.excrement.floor.binary, useNA="always")

# Select relevant WASH variables from the household survey
################################################################################################

# variables excluding healthcare seeking behaviour survey questions (and related medicine use); as these
# are not 1 observation per household
d_wash = d %>% select(data.row,menage_id,village, village_name, intervention.text,   
                              redcap_event_name,
                              date.enquete,
                              groupe,
                              n.householdmember, 
                              n.child.0to5,
                              n.households.concession,
                              n.households,
                              n.population,
                              ses.quintile,
                              main.drinking.water.dry.binary,
                              main.drinking.water.rainy.binary, 
                              cleaning.water.storage.binary,
                              correct.handwashing.binary,
                              improved.sanitation.binary,
                              livestock.access.house.binary,
                              animal.excrement.floor.binary,                             
                              q1.diar.prev.water.pot.covered,
                              q1.diar.prev.no.finger.in.waterglass,                    
                              q1.diar.prev.utensil.to.take.water.from.pot, q1.diar.prev.cover.food,                                  
                              q1.diar.prev.boil.water, q1.diar.prev.filter.water,                               
                              q1.diar.prev.other,q1.diar.prev.cant.be.avoided,                            
                              q1.diar.prev.dont.know, autr.mesur.prev.diarrhe,                                  
                              q2.main.water.source.dry, q3.main.water.source.rainy,                              
                              q4.cans.storage.water, q5a.cans.storage.water.closed.when.filled,                                   
                              q5b.cans.storage.water.closed.wen.empty, q5c.cans.cleaned.before.reuse,                                        
                              q6.treatment.water, q6.traitmen.water.other,                                    
                              q7.principle.defication, q7.principle.defication.other,                                 
                              q8.other.defecation.flush.toiled.septic, q8.other.defecation.pit.latrine.ventilation,             
                              q8.other.defecation.pit.latrine.slab, q8.other.defecation.pit.latrine.no.slab,                  
                              q8.other.defecation.open.defecation, q8.other.defecation.other,                                
                              q8.other.defecation.none, q8.other.defecation.specified,                                        
                              q9.shared.toilet, q10.n.shared.toilet,                                       
                              q11.toilet.last.cleaned, q12.disposal.child.stool,                                    
                              q12.disposal.child.other.specified, q13.disposal.latrine.pit,                                     
                              q13.latrine.pit.other.specified, q14.handwashing.product,                                   
                              q15.handwashing.defecation, q16.handwashing.meals,                                       
                              q17.animals.around.household, q18.animal.inside.cow,                                    
                              q18.animal.inside.sheep.goat, q18.animal.inside.pig,                                    
                              q18.animal.inside.donkey.horse, q18.animal.inside.chicken.goose.duck,                     
                              q18.animal.inside.other,q18.animal.inside.other.specified,                                       
                              q19.animal.outside.cow, q19.animal.outside.sheep.goat,                             
                              q19.animal.outside.pig,q19.animal.outside.donkey.horse,                           
                              q19.animal.outside.chicken.goose.duck, q19.animal.outside.other,                                  
                              q19.animal.outside.other.specified, q20.animal.excrement.floor,                                    
                              q21.when.animal.ill.treatment.with.vet, q21.when.animal.ill.treatment.without.vet,                
                              q21.when.animal.ill.sell.bucher,q21.when.animal.ill.slaugter.eat.at.home,                     
                              q21.when.animal.ill.dies.burie.dispose,q21.when.animal.ill.dies.eat.at.home,                                
                              #eau.assainissement.hygine.complete,
                              piam.q1.n.athome.when.survey,
                              piam.q2a.athome.age.random.selected.person,
                              piam.q2a.athome.sex.random.selected.person,
                              piam.q3.heard.about.trial.activities,
                              piam.q4.partic.atleastone.trial.activities,
                              piam.q5.perc.key.topic.hwashing.use.soap,
                              piam.q5.perc.key.topic.hwashing.when,
                              piam.q5.perc.key.topic.abx.what.and.role,
                              piam.q5.perc.key.topic.abx.correct.use,
                              piam.q5.perc.key.topic.seek.care.no.self.med,
                              piam.q5.perc.key.topic.good.sani.latr.wastw,
                              piam.q5.perc.key.topic.drink.treat.or.borehole.water,
                              piam.q5.perc.key.topic.protect.store.water.clean.container,
                              piam.q5.perceived.key.topic.handwashing,
                              piam.q6.freq.partic.hwashing.public,
                              piam.q6.freq.partic.abx.use.public,
                              piam.q6.freq.partic.hwashing.athome,
                              piam.q6.freq.partic.abx.use.athome,
                              piam.q6.freq.partic.abx.use.adhoc,
                              piam.q6.freq.partic.defication.adhoc,
                              piam.q6.freq.partic.defication.public) %>%
  filter(!is.na(date.enquete)) # Denominator data (i.e. people tested for esbl) for R0

d_wash_r2 = d %>% filter(is.na(redcap.repeat.instrument) & redcap_event_name=="round_3_arm_1") %>%
  select(data.row,menage_id,village, village_name, intervention.text,   
         redcap_event_name,
         date.enquete,
         groupe,
         n.householdmember, 
         n.child.0to5,
         n.households.concession,
         n.households,
         n.population,
         ses.quintile,
         main.drinking.water.dry.binary,
         main.drinking.water.rainy.binary, 
         cleaning.water.storage.binary,
         correct.handwashing.binary,
         improved.sanitation.binary,
         livestock.access.house.binary,
         animal.excrement.floor.binary,
         q1.diar.prev.water.pot.covered,
         q1.diar.prev.no.finger.in.waterglass,                    
         q1.diar.prev.utensil.to.take.water.from.pot, q1.diar.prev.cover.food,                                  
         q1.diar.prev.boil.water, q1.diar.prev.filter.water,                               
         q1.diar.prev.other,q1.diar.prev.cant.be.avoided,                            
         q1.diar.prev.dont.know, autr.mesur.prev.diarrhe,                                  
         q2.main.water.source.dry, q3.main.water.source.rainy,                              
         q4.cans.storage.water, q5a.cans.storage.water.closed.when.filled,                                   
         q5b.cans.storage.water.closed.wen.empty, q5c.cans.cleaned.before.reuse,                                        
         q6.treatment.water, q6.traitmen.water.other,                                    
         q7.principle.defication, q7.principle.defication.other,                                 
         q8.other.defecation.flush.toiled.septic, q8.other.defecation.pit.latrine.ventilation,             
         q8.other.defecation.pit.latrine.slab, q8.other.defecation.pit.latrine.no.slab,                  
         q8.other.defecation.open.defecation, q8.other.defecation.other,                                
         q8.other.defecation.none, q8.other.defecation.specified,                                        
         q9.shared.toilet, q10.n.shared.toilet,                                       
         q11.toilet.last.cleaned, q12.disposal.child.stool,                                    
         q12.disposal.child.other.specified, q13.disposal.latrine.pit,                                     
         q13.latrine.pit.other.specified, q14.handwashing.product,                                   
         q15.handwashing.defecation, q16.handwashing.meals,                                       
         q17.animals.around.household, q18.animal.inside.cow,                                    
         q18.animal.inside.sheep.goat, q18.animal.inside.pig,                                    
         q18.animal.inside.donkey.horse, q18.animal.inside.chicken.goose.duck,                     
         q18.animal.inside.other,q18.animal.inside.other.specified,                                       
         q19.animal.outside.cow, q19.animal.outside.sheep.goat,                             
         q19.animal.outside.pig,q19.animal.outside.donkey.horse,                           
         q19.animal.outside.chicken.goose.duck, q19.animal.outside.other,                                  
         q19.animal.outside.other.specified, q20.animal.excrement.floor,                                    
         q21.when.animal.ill.treatment.with.vet, q21.when.animal.ill.treatment.without.vet,                
         q21.when.animal.ill.sell.bucher,q21.when.animal.ill.slaugter.eat.at.home,                     
         q21.when.animal.ill.dies.burie.dispose,q21.when.animal.ill.dies.eat.at.home,                                
         #eau.assainissement.hygine.complete,
         piam.q1.n.athome.when.survey,
         piam.q2a.athome.age.random.selected.person,
         piam.q2a.athome.sex.random.selected.person,
         piam.q3.heard.about.trial.activities,
         piam.q4.partic.atleastone.trial.activities,
         piam.q5.perc.key.topic.hwashing.use.soap,
         piam.q5.perc.key.topic.hwashing.when,
         piam.q5.perc.key.topic.abx.what.and.role,
         piam.q5.perc.key.topic.abx.correct.use,
         piam.q5.perc.key.topic.seek.care.no.self.med,
         piam.q5.perc.key.topic.good.sani.latr.wastw,
         piam.q5.perc.key.topic.drink.treat.or.borehole.water,
         piam.q5.perc.key.topic.protect.store.water.clean.container,
         piam.q5.perceived.key.topic.handwashing,
         piam.q6.freq.partic.hwashing.public,
         piam.q6.freq.partic.abx.use.public,
         piam.q6.freq.partic.hwashing.athome,
         piam.q6.freq.partic.abx.use.athome,
         piam.q6.freq.partic.abx.use.adhoc,
         piam.q6.freq.partic.defication.adhoc,
         piam.q6.freq.partic.defication.public) 
# Denominator data (i.e. people tested for esbl) for r2
names(d)
d_wash_r2 = merge(d_wash_r2, d_wash %>% select(menage_id, n.householdmember,n.child.0to5, n.households.concession), by="menage_id",all.x=T) 
d_wash_r2 = d_wash_r2 %>%
  mutate(n.householdmember.x = n.householdmember.y,
         n.child.0to5.x = n.child.0to5.y,
         n.households.concession.x = n.households.concession.y) %>%
  select(-c(n.householdmember.y,n.child.0to5.y, n.households.concession.y)) %>%
  rename(n.householdmember = "n.householdmember.x",
         n.child.0to5 = "n.child.0to5.x",
         n.households.concession ="n.households.concession.x")

table(d_wash$redcap_event_name)

# Check if any of the 'other' categories needs merging
#---------------------------------------------------------------------
describe(d_wash) # All the 'other' categories have missing observations, so these could be deleted
describe(d_wash_r2) # All the 'other' categories have missing observations, so these could be deleted

d_wash = d_wash %>% select(-c(autr.mesur.prev.diarrhe, q6.traitmen.water.other, q7.principle.defication.other,     
                                       q8.other.defecation.specified, q13.latrine.pit.other.specified, q18.animal.inside.other.specified,
                                      q19.animal.outside.other.specified))

d_wash_r2 = d_wash_r2 %>% select(-c(autr.mesur.prev.diarrhe, q6.traitmen.water.other, q7.principle.defication.other,     
                                      q8.other.defecation.specified, q13.latrine.pit.other.specified, q18.animal.inside.other.specified,
                                      q19.animal.outside.other.specified))

skim(d_wash)
skim(d_wash_r2)

# CHECK FOR DUPLICATES
# R1
length(unique(d_wash$menage_id))
sort(table(d_wash$menage_id, useNA="always")) # 1 duplicate, i.e. EKA00002004
#View(d_wash[d_wash$menage_id=="EKA00002004",]) # Answers not the same, but remove 1
d_wash = d_wash %>% filter(!duplicated(menage_id)) # 808 unique households
 
names(hh_bf)

# R2
length(unique(d_wash_r2$menage_id))
sort(table(d_wash_r2$menage_id, useNA="always")) # 1 duplicate, i.e. EKA00002004
#View(d_wash_r2[d_wash_r2$menage_id=="EKA00002004",]) # All answers the same, except for hh size; remove 1

d_wash_r2 = d_wash_r2 %>% filter(!duplicated(menage_id)) # 776 unique households


#-------------------------------------------------------------------------------
# HOUSEHOLD SURVEY - INDIVIDUAL CHARACTERISTICS
#-------------------------------------------------------------------------------

# ALL ROUNDS
###################################################################################

# Link lab data with lab vs hh survey ids (as IDs are in different format)
car_bf = left_join(car_bf, hh_lab_ids, by="household")
# Check if all linked
table(car_bf$bras, useNA= "always")

# Link lab data with hh survey ids
length(unique(car_bf$menage_id)) # 261
length(unique(d$menage_id)) # 808

# Can all IDs be traced back from lab to wash survey?
car_bf$found_in_hh[which(car_bf$menage_id %in% hh_bf$menage_id)] = 1
car_bf$found_in_hh[is.na(car_bf$found_in_wash)] = 0
length(unique(car_bf$household[car_bf$found_in_wash==0])) # after cleaning, all households can be found in WASH survey database; 
unique(car_bf$menage_id[car_bf$found_in_wash==0])
unique(car_bf$household[car_bf$found_in_wash==0])


# Make dataset with only those household individuals that had a stool sample taken and their individual variables
d_lab = d %>% filter(!is.na(cs.id.individu)) %>% # ensure just 1 observation per person of whom individual is esbl positive
  select(data.row, redcap_event_name,cs.id.individu,num.echantillon, menage_id, village, age, agegr10, agegr, sexe, date.consent, date.stool.collection)

# Merge wash patient characteristics with lab data - NEED TO GET IDs IIN THE SAME FORMAT
which(car_bf$record_id %in% unique(d_lab$cs.id.individu)) # Have to change the format of both to make sure matching can be done
head(car_bf$record_id ); head(d_lab$cs.id.individu)
unique(car_bf$record_id)
unique(d$cs.id.individu)

#--------------------------------------------------------------------------------------------
# ENABLE LINKAGE BETWEEN LAB DATA AND HH DATA BY ADAPTING THE IDs TO THE SAME FORMAT
#--------------------------------------------------------------------------------------------
# HH database does no have "-" and puts and "M" before each household member. Also sometimes 01 and sometimes 1 as method of writing
# PROBABLY IF WE TAKE FROM wash_r0 cs_id_individue all "MX" (so M + number after ID), then put that in a seperate column.
# Then we take the household number from the WASH survey (so menage_id), combine these to with a hyphen

# THEN for car_bf we take also menage_Id and digits after the "-", remove the 0's then two can be combined

# Change IDs to the same format

# Create a variable in which we will safe the new formatted and cleaned ID
d_lab$menage_id_member = NA

# Create a variable that will just store the household member number
d_lab$member = NA

# Extract the number after M to get household number
d_lab$member =  gsub(".*M", "", d_lab$cs.id.individu)
# Remove leading '0's
d_lab$member = as.character(as.numeric(d_lab$member))
# Check the one's which are now NA
table(is.na(d_lab$member))
d_lab$cs.id.individu[is.na(d_lab$member)] # 22 individuals have a household member number missing, see if still all lab_ids can be linked when merging with car_bf
#View(wash_r0_lab[is.na(wash_r0_lab$member),])
d_lab$num.echantillon[is.na(d_lab$member)] # of 0 we can get them from num_echantillon; however it seems also here, member number is missing

# Now create new variable for linking with lab dataset
d_lab$menage_id_member = paste0(d_lab$menage_id, "-", d_lab$member)
d_lab$menage_id_member
# Make NAs for the one's that still need checking
d_lab$menage_id_member[is.na(d_lab$member)] = NA
table(is.na(d_lab$menage_id_member))


# Check if duplicates from wash_r0
table(duplicated(d_lab$menage_id_member)) # 3018 duplicates (as multiple rounds of data)
dups = unique(d_lab$menage_id_member[duplicated(d_lab$menage_id_member)]) # 1157
#View(wash_bf_lab[which(wash_bf_lab$menage_id_member%in%dups),])

# Reorder variables
new_order <- c("data.row", "redcap_event_name", "village", "cs.id.individu", 
               "num.echantillon", "menage_id_member", "member", "menage_id",
               "date.consent", "date.stool.collection",
               "age", "agegr10","agegr", "sexe")

# Check if all columns exist
setdiff(new_order, names(d_lab))  # Shows any missing columns

d_lab <- d_lab[, new_order]  

# Sort out NAs
d_lab$cs.id.individu[is.na(d_lab$member)]
#"SCB0220203"  "SCB0220201"  "SCB02202002" "SCB0220204"  "SCB0220205"  "SCB0220206"  "SCB0240201"  "SCB0240202"  "SCB0240203" 
# "SCB0240204"  "SCB0240205"  "SCB0220206"  "SCB0240207"  "SCB0240208"  "SCB0240209"  "SCB0240210"  "SCB0240211"  "SCB0240212" 
# "SCB0310101"  "SCB0310102"  "SCB0310103"  "SCB0310104"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220203"] = "SCB00002202-3"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220201"] = "SCB00002202-1"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB02202002"] = "SCB00002202-2"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220204"] = "SCB00002202-4"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220205"] = "SCB00002202-5"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220206"] = "SCB00002202-6"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240201"] = "SCB00002402-1"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240202"] = "SCB00002402-2"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240203"] = "SCB00002402-3"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240204"] = "SCB00002402-4"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240205"] = "SCB00002402-5"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0220206"&d_lab$num.echantillon=="SCB0240206"] = "SCB00002402-6"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240207"] = "SCB00002402-7"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240208"] = "SCB00002402-8"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240209"] = "SCB00002402-9"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240210"] = "SCB00002402-10"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240211"] = "SCB00002402-11"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0240212"] = "SCB00002402-12"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0310101"] = "SCB00003101-1"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0310102"] = "SCB00003101-2"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0310103"] = "SCB00003101-3"
d_lab$menage_id_member[d_lab$cs.id.individu=="SCB0310104"] = "SCB00003101-4"

d_lab$member[d_lab$data.row=="4112"] = 2
d_lab$menage_id_member[d_lab$data.row=="4112"] = "ELB00002903-2"
d_lab$member[d_lab$data.row=="3660"] = 2
d_lab$menage_id_member[d_lab$data.row=="3660"] = "EKE00047020-2"
d_lab$member[d_lab$data.row=="2261"] = 1
d_lab$menage_id_member[d_lab$data.row=="2261"] = "EEI00000107-1"
d_lab$member[d_lab$data.row=="1903"] = 3
d_lab$menage_id_member[d_lab$data.row=="1903"] = "EDN00002003-3"


table(d_lab$redcap_event_name) # These are the denominator number of individuals tested per round, i.e. 1221, 1014, 1017, 1003
length(unique(d_lab$menage_id_member[d_lab$redcap_event_name=="round_0_arm_1"])) # 1214
length(unique(d_lab$menage_id_member[d_lab$redcap_event_name=="round_1_arm_1"])) # 1010
length(unique(d_lab$menage_id_member[d_lab$redcap_event_name=="round_2_arm_1"])) # 1011
length(unique(d_lab$menage_id_member[d_lab$redcap_event_name=="round_3_arm_1"])) # 1003

# Lab data only round 1
d_lab_r1 = d_lab %>% filter(redcap_event_name=="round_0_arm_1")
d_lab_r2 = d_lab %>% filter(redcap_event_name=="round_1_arm_1")
d_lab_r3 = d_lab %>% filter(redcap_event_name=="round_2_arm_1")
d_lab_r4 = d_lab %>% filter(redcap_event_name=="round_3_arm_1")

# Remove duplicates - round 1
dup = d_lab_r1$menage_id_member[duplicated(d_lab_r1$menage_id_member)] 
#View(d_lab_r1%>%filter(menage_id_member%in%dup)) # Genuine duplicates
#table(duplicated(d_lab_r1$menage_id_member))# 7 duplicates, these are all genuine duplicates and can be removed 
rm.data.row = d_lab_r1$data.row[duplicated(d_lab_r1$menage_id_member)]
d_lab_r1 = d_lab_r1 %>%filter(!data.row %in%rm.data.row)

# Remove duplicates - round 2
dup2 = d_lab_r2$menage_id_member[duplicated(d_lab_r2$menage_id_member)] 
length(unique(dup2)) # 4 duplicates
#View(d_lab_r2%>%filter(menage_id_member%in%dup2)) # CRAS0101-1 and SJB00003702-2 not duplicates 
#View(d_lab_r1 %>%filter(menage_id_member%in%dup2)) # EAD00000903-4; SJB00003702-2 seem to be different persons in round 0 than in round 1?
#View(d_lab_r3 %>%filter(menage_id_member%in%dup2)) # Checked,CRAS0101-1 = Female 10-19; EAD00000903-4 = Male 0-9 in 3 out of 4 rounds; SJB00003702-2 = Female 20-29 in 3/4 rounds
#View(d_lab_r4 %>%filter(menage_id_member%in%dup2))

d_lab_r2 = d_lab_r2 %>% filter(!data.row %in%c(345,1012,4185,6776))
d_lab_r2$age[d_lab_r2$data.row==6775] = 24
d_lab_r2$agegr10[d_lab_r2$data.row==6775] = "20-29"
d_lab_r2$agegr[d_lab_r2$data.row==6775] = "18-49"
d_lab_r2$sexe[d_lab_r2$data.row==6775] = "Female"

# Remove duplicates - round 3
dup3 = d_lab_r3$menage_id_member[duplicated(d_lab_r3$menage_id_member)] 
length(unique(dup3)) # 2
#View(d_lab_r3 %>%filter(menage_id_member%in%dup3)) # Duplicates
rm.data.row.r3 =  d_lab_r3$data.row[duplicated(d_lab_r3$menage_id_member)]
d_lab_r3 = d_lab_r3 %>%filter(!data.row%in%rm.data.row.r3)

# Do the same for d_lab
d_lab = d_lab %>% filter(!data.row%in%c(rm.data.row,rm.data.row.r3))
d_lab = d_lab %>% filter(!data.row %in%c(345,1012,4185,6776))
d_lab$age[d_lab$data.row==6775] = 24
d_lab$agegr10[d_lab$data.row==6775] = "20-29"
d_lab$agegr[d_lab$data.row==6775] = "18-49"
d_lab$sexe[d_lab$data.row==6775] = "Female"

# Then for remainder keep first record 
d_lab_de = d_lab %>% filter(!duplicated(menage_id_member))  # 1237
table(duplicated(d_lab_de$menage_id_member)) # 
table(d_lab_de$redcap_event_name) # 1191 individuals in round 0, new ones: 21 (round 1); 18 (round 2); 7 (round 3) 
                                  # Update: after converting about IDs this is
                                  # 1205 individuals in round 0, new ones: 8 (round 1); 17 (round 2); 7 (round 3) 


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

# Reorder columns
# Define the variables to move to the start
start_vars <- c("record_id","household",  "redcap_event_name",                                          
                "redcap_repeat_instrument","redcap_repeat_instance","village", "village_name", "intervention_text")

# Define the variables to move to the end
middle_vars <- c("id_ecantillon","menage_id_member","member", "menage_id", "germe","morphotyp", "germe_c", "testesbl", "esbl_pos")

# Get all remaining variables (excluding those in start_vars and end_vars)
end_vars <- setdiff(names(car_bf), c(start_vars, middle_vars))

# Define the new column order
new_order <- c(start_vars, middle_vars, end_vars)
new_order

# Reorder the dataset
car_bf <- car_bf[, new_order]

table(car_bf$esbl_pos)
table(car_bf$testesbl)

# Remove unnecessary variables
names(car_bf)
car_bf <- car_bf %>% select(-c(ajouter, bras, found_in_hh, redcap_repeat_instrument, redcap_repeat_instance, check_list))


#-----------------------------------------------------------
# LINK INDIVIDUAL WITH STOOL DATA
#-----------------------------------------------------------
describe(car_bf)
car_bf = car_bf %>% select(-c(
  "interpretr_ampici_amoxicil", "comment_ampici_amoxic", "interpr_amoxi_acid_clavu", "comment_amoxi_acid_clavu",
  "interpr_piperaci_tazobac", "comment_piperacil_tazobact", "interpr_cetriax_or_cefotax", "interpr_cefepime",
  "comment_cefepime", "comment_meropenem", "comment_ertapenem", "interpr_gentamycine",
  "comment_gentamycine", "interpr_amykacine", "comment_amykacine", "interpr_ciprofloxacine",
  "comment_ciprofloxacine", "interpr_pfloxacine", "comment_pfloxacine", "interpr_sulfame_trimethop",
  "comment_sulfa_trimethop","comment_cetriax_cefataxi"))

d_lab_car_bf_r1 = left_join(d_lab_r1 %>% select(-c(village, redcap_event_name, data.row)), car_bf %>% filter(redcap_event_name=="round_0_arm_1"))
d_lab_car_bf_r2 = left_join(d_lab_r2 %>% select(-c(village, redcap_event_name, data.row)), car_bf %>% filter(redcap_event_name=="round_1_arm_1")) 
d_lab_car_bf_r3 = left_join(d_lab_r3 %>% select(-c(village, redcap_event_name, data.row)), car_bf %>% filter(redcap_event_name=="round_2_arm_1"))
d_lab_car_bf_r4 = left_join(d_lab_r4 %>% select(-c(village, redcap_event_name, data.row)), car_bf %>% filter(redcap_event_name=="round_3_arm_1"))

table(d_lab$redcap_event_name)

length(unique(d_lab_car_bf_r1$menage_id_member)) # 1214
length(unique(d_lab_car_bf_r2$menage_id_member)) # 1010
length(unique(d_lab_car_bf_r3$menage_id_member)) # 1015
length(unique(d_lab_car_bf_r4$menage_id_member)) # 1003

# Create a dataframe with households linked to villages so these can be augmented for those that were not in the lab database
hh_villages = d_wash %>%
  select(village, village_name, intervention.text, menage_id) %>%
  distinct()

hh_villages = hh_villages %>%
  rename(
    intervention_text = "intervention.text"
  )

#anti_join(d_lab_car_bf_r1, hh_villages, by = "menage_id") %>%
#  count(menage_id)

# ROUND 1
d_lab_car_bf_r1 = left_join(d_lab_car_bf_r1 %>%select(-c(village, village_name, intervention_text)), hh_villages, by = "menage_id")
d_lab_car_bf_r1$redcap_event_name = "round_0_arm_1"
d_lab_car_bf_r1$esbl_pos = ifelse(is.na(d_lab_car_bf_r1$esbl_pos), 0,d_lab_car_bf_r1$esbl_pos) 
  
table(d_lab_car_bf_r1$esbl_pos, useNA="always")

# ROUND 2
d_lab_car_bf_r2 = left_join(d_lab_car_bf_r2 %>%select(-c(village, village_name,intervention_text)), hh_villages, by = "menage_id")
d_lab_car_bf_r2$redcap_event_name = "round_1_arm_1"
d_lab_car_bf_r2$esbl_pos = ifelse(is.na(d_lab_car_bf_r2$esbl_pos), 0,d_lab_car_bf_r2$esbl_pos) 

table(d_lab_car_bf_r2$esbl_pos, useNA="always")

# ROUND 3
d_lab_car_bf_r3 = left_join(d_lab_car_bf_r3 %>%select(-c(village, village_name, intervention_text)), hh_villages, by = "menage_id")
d_lab_car_bf_r3$redcap_event_name = "round_2_arm_1"
d_lab_car_bf_r3$esbl_pos = ifelse(is.na(d_lab_car_bf_r3$esbl_pos), 0,d_lab_car_bf_r3$esbl_pos) 

table(d_lab_car_bf_r3$esbl_pos, useNA="always")

# ROUND 4
d_lab_car_bf_r4 = left_join(d_lab_car_bf_r4 %>%select(-c(village, village_name, intervention_text)), hh_villages, by = "menage_id")
d_lab_car_bf_r4$redcap_event_name = "round_3_arm_1"
d_lab_car_bf_r4$esbl_pos = ifelse(is.na(d_lab_car_bf_r4$esbl_pos), 0,d_lab_car_bf_r4$esbl_pos) 

table(d_lab_car_bf_r4$esbl_pos, useNA="always")

d_lab_car_bf_all = rbind(d_lab_car_bf_r1,d_lab_car_bf_r2, d_lab_car_bf_r3, d_lab_car_bf_r4)
missing.age = d_lab_car_bf_all$menage_id_member[is.na(d_lab_car_bf_all$age)]
length(missing.age)

# Identify missing ages and replace
dt <- d_lab_car_bf_all %>%
  group_by(menage_id_member) %>%
  mutate(
    age = ifelse(is.na(age), first(na.omit(age)), age),
    agegr10 = ifelse(is.na(agegr10), first(na.omit(as.character(agegr10))), as.character(agegr10)),
    agegr = ifelse(is.na(agegr), first(na.omit(as.character(agegr))), as.character(agegr))
  ) %>%
  ungroup()

table(is.na(dt$age)) # Replaced 111-38 = 73
table(is.na(d_lab_car_bf_all$age))

# Check
summary(dt$age, rm.na=T)
summary(d_lab_car_bf_all$age, rm.na=T) # Same averages

# Augment collection dates when missing
dt = dt %>% mutate(
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
  rainy = ifelse(month%in%c("06","07","08","09"),"yes", "no")) %>% 
  select(-c(date.check, date.check.collect, date.check.use))

#ggplot(dt, aes(x=date.check.collect, group=redcap_event_name)) + geom_histogram() +
#  facet_wrap(~ redcap_event_name)

# CLEAN DATES FURTHER
table(dt$date.use, useNA="always")
dt$date.use[is.na(dt$date.use)] = dt$date[is.na(dt$date.use)] 

dt <- dt %>%
  group_by(menage_id_member) %>%
  mutate(
    first_round_date = min(date.use[redcap_event_name == "round_0_arm_1"], na.rm = TRUE)
  ) %>%
  ungroup()

# Identify individuals with missing first-round dates
dt <- dt %>%
  mutate(
    first_round_missing = is.na(first_round_date) | first_round_date == Inf
  )

table(dt$first_round_missing)

# Summary
cat("Number of individuals who received a first-round date from a household member:", sum(dt$first_round_missing==T), "\n")

# Fill missing first_round_date using household members (menage_id)
#view(table(dt$menage_id, dt$first_round_date))

dt <- dt %>%
  group_by(menage_id) %>%
  mutate(
    first_round_date = ifelse(first_round_missing,
                              max(first_round_date[first_round_date!=Inf], na.rm=T),  # Borrow from household members
                              first_round_date)
  ) %>%
  ungroup()
table(dt$first_round_date, useNA="always")

dt <- dt %>%
  mutate(
    first_round_missing = first_round_date == -Inf
  )
table(dt$first_round_missing)

dt <- dt %>%
  group_by(village_name) %>%
  mutate(
    first_round_date = ifelse(first_round_missing,
                              max(first_round_date[first_round_date!=-Inf], na.rm=T),  # Borrow from household members
                              first_round_date)
  ) %>%
  ungroup()
table(dt$first_round_date)

# Ensure no Inf values remain
#dt$first_round_date[dt$first_round_date == Inf] <- NA

# Convert to Date format (if numeric)
dt <- dt %>%
  mutate(
    first_round_date = as.Date(first_round_date, origin = "1970-01-01"),
    date.use = as.Date(date.use)  # Ensure date.use is also in Date format
  )

# Calculate the number of days since the first round
dt <- dt %>%
  mutate(
    days_since_first_round = as.numeric(date.use - first_round_date)  # Time difference in days
  )

# Remove helper column
dt <- dt %>% select(-first_round_missing)


# CHECK
ggplot(dt, aes(x = days_since_first_round)) +
  geom_histogram(binwidth = 10, fill = "black", alpha = 0.8) +
  facet_wrap(~redcap_event_name) +
  theme_minimal() +
  labs(title = "Distribution of Days Since First Round by Event",
       x = "Days Since First Round",
       y = "Count")

dt %>%
  group_by(redcap_event_name, time) %>%
  summarise(
    mean_days_since_first = mean(days_since_first_round, na.rm = TRUE),
    median_days_since_first = median(days_since_first_round, na.rm = TRUE),
    min_days = min(days_since_first_round, na.rm = TRUE),
    max_days = max(days_since_first_round, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(time) %>%
  print()

# All looks good for round 0
r0min = min(dt$first_round_date, na.rm = TRUE)  # Start of Round 0
r1min = r0min + 90  # Expected minimum for Round 1 (~3 months after R0)
r2min = r0min + 180 # Expected minimum for Round 2 (~6 months after R0)
r3min = r0min + 365 # Expected minimum for Round 3 (~12 months after R0)

# Ensure correct sequencing by filtering cases that meet the minimum required date
dt = dt %>%
  mutate(
    rc = case_when(
      time == 0 ~ 1,  # Round 0 should always be included
      time == 1 & date.use >= r1min ~ 1,  # Round 1 should be after r1min
      time == 2 & date.use >= r2min ~ 1,  # Round 2 should be after r2min
      time == 3 & date.use >= r3min ~ 1,  # Round 3 should be after r3min
      TRUE ~ 0  # Mark as incorrect if date.use is before the expected minimum
    )
  )# %>%
  #filter(rc == 1)  # Keep only correctly ordered cases

table(dt$rc) # 34 not correct
#View(dt%>% filter(rc==0))


dt = dt %>%
  mutate(date.use = ifelse(rc==0, date, date.use))

table(dt$date.use, useNA = "always") # 13 left with NA 

# Fill missing date.use values using a household member from the same round
dt <- dt %>%
  group_by(menage_id, redcap_event_name) %>%  # Household level within the same round
  mutate(
    date.use = ifelse(is.na(date.use), max(date.use, na.rm = TRUE), date.use)  # Borrow from household members
  ) %>%
  ungroup()

# Fill with village date
dt <- dt %>%
  group_by(village, redcap_event_name) %>%  # Household level within the same round
  mutate(
    date.use = ifelse(date.use==-Inf, max(date.use, na.rm = TRUE), date.use)  # Borrow from household members
  ) %>%
  ungroup()

# Convert to Date format (if numeric)
dt <- dt %>%
  mutate(
    date.use = as.Date(date.use)  # Ensure date.use is also in Date format
  )

table(dt$date.use, useNA = "always") # 2 left with Inf

# Calculate again the number of days since the first round
dt <- dt %>%
  mutate(
    days_since_first_round = as.numeric(date.use - first_round_date)  # Time difference in days
  )

# CHECK
ggplot(dt, aes(x = days_since_first_round)) +
  geom_histogram(binwidth = 10, fill = "black", alpha = 0.8) +
  facet_wrap(~redcap_event_name) +
  theme_minimal() +
  labs(title = "Distribution of Days Since First Round by Event",
       x = "Days Since First Round",
       y = "Count")

dt %>%
  group_by(redcap_event_name, time) %>%
  summarise(
    mean_days_since_first = mean(days_since_first_round, na.rm = TRUE),
    median_days_since_first = median(days_since_first_round, na.rm = TRUE),
    min_days = min(days_since_first_round, na.rm = TRUE),
    max_days = max(days_since_first_round, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(time) %>%
  print()

# ONE STILL NOT CORRECT, THAT IS OF EHL00002504-1, HERE THE FIRST ROUND DATE IS INCORRECT
dt$date.use[dt$menage_id_member=="EHL00002504-1"&dt$redcap_event_name=="round_0_arm_1"] = as.Date("2022-11-06")  
dt$first_round_date[dt$menage_id_member=="EHL00002504-1"&dt$redcap_event_name=="round_0_arm_1"] = as.Date("2022-11-06")  
dt$first_round_date[dt$menage_id_member=="EHL00002504-1"] = as.Date("2022-11-06")  

# Calculate again the number of days since the first round
dt <- dt %>%
  mutate(
    days_since_first_round = as.numeric(date.use - first_round_date)  # Time difference in days
  )

#View(dt%>%filter(days_since_first_round==0&redcap_event_name=="round_1_arm_1"))

# TWO STILL NOT CORRECT, THAT IS OF EDL00000401-1 and EDL00000401-2, HERE THE SECOND ROUND DATE IS INCORRECT
# EDL00000401-1; EDL00000401-2
dt$date.use[dt$menage_id_member=="EDL00000401-1"&dt$redcap_event_name=="round_1_arm_1"] = as.Date("2023-02-24") 
dt$date.use[dt$menage_id_member=="EDL00000401-2"&dt$redcap_event_name=="round_1_arm_1"] = as.Date("2023-02-24") 

# Calculate again the number of days since the first round
dt <- dt %>%
  mutate(
    days_since_first_round = as.numeric(date.use - first_round_date)  # Time difference in days
  )

dt %>%
  group_by(redcap_event_name, time) %>%
  summarise(
    mean_days_since_first = mean(days_since_first_round, na.rm = TRUE),
    median_days_since_first = median(days_since_first_round, na.rm = TRUE),
    min_days = min(days_since_first_round, na.rm = TRUE),
    max_days = max(days_since_first_round, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(time) %>%
  print()

# YES!

# NOW CREATE DATASET FOR MSM MODEL
dt_filtered <- dt %>%
  dplyr::filter(!germe_c %in% c("salmonella")) %>%  # Keep only e.coli records; for those without ESBL-E, assume they all have e.coli
  group_by(menage_id_member, redcap_event_name) %>%
  summarise(
    time = first(time),
    menage_id = first(menage_id),
    village = first(village),
    village_name = first(village_name),
    intervention_text = first(intervention_text),
    age = first(age[!is.na(age)]),
    agegr10 = first(agegr10[!is.na(agegr10)]),
    agegr = first(agegr[!is.na(agegr)]),
    sexe = first(sexe[!is.na(sexe)]),
    month = first(month),
    rainy = first(rainy),
    date.consent = first(date.consent[!is.na(date.consent)]),
    date.stool.collection = first(date.stool.collection[!is.na(date.stool.collection)]),
    date = first(date[!is.na(date)]),
    date_conserv = first(date_conserv[!is.na(date_conserv)]),
    date.use = first(date.use[!is.na(date.use)]),
    days = first(days_since_first_round),
    germe_c = "e.coli",
    esble = ifelse(any(esbl_pos == 1, na.rm = TRUE), 1, 0),  # Set esble = 1 if any observation is positive
    .groups = "drop"
  )

dt %>% 
  filter(germe_c%in%c(NA,"e.coli", "e.coli_2"))%>%
  group_by(redcap_event_name)%>%
  summarise(
  n = length(unique(menage_id_member))
  )

table(dt_filtered$redcap_event_name)

dt %>% 
  filter(germe_c%in%c(NA,"e.coli", "e.coli_2"))%>%
  group_by(redcap_event_name, esbl_pos)%>%
  summarise(
    n = length(unique(menage_id_member))
  )

dt_filtered %>% 
  group_by(redcap_event_name,esble)%>%
  summarise(
    n = length(unique(menage_id_member))
  )


# Filtered to those with complete cases
numround = dt_filtered %>% 
  group_by(menage_id_member)%>%
  summarise(
    n = n())
table(numround$n) # 747 with complete data

com = numround %>% filter(n>3)

# Add number of tested individual per round per household as denominator for prevalence
ntested = dt_filtered %>% group_by(menage_id, time) %>%
  summarise(
    n.tested = length(menage_id_member)
  )

dt_filtered = left_join(dt_filtered, ntested)

dt_filtered_complete = dt_filtered%>%filter(menage_id_member%in%com$menage_id_member)

# NOW LINK DATA WITH WASH HOUSEHOLD DATA
#-------------------------------------------------------------------------------------
d_wash_b = rbind(d_wash, d_wash_r2)

# All stool samples
dt_wash = left_join(dt, d_wash_b%>%select(-c(village,village_name, intervention.text, groupe)))

# MSM data stool samples
dt_filtered_wash = left_join(dt_filtered, d_wash_b%>%select(-c(village,village_name, intervention.text, groupe)))
dt_filtered_complete_wash = left_join(dt_filtered_complete, d_wash_b%>%select(-c(village,village_name, intervention.text, groupe)))

#-----------------------------------------------------------
# LINK INDIVIDUAL WITH HOUSEHOLD DATA
#-----------------------------------------------------------

# Of note: the lab data contains only the ones positive for ESBL and/or with AST conducted.
# Therefore, we need to link individuals in the household survey with individual characteristics 
# With the lab survey to get the denominator of those tested.
# So these datasets contain information on for each individual part of the cohort somewhere (across all rounds, 1237 unique
# individuals in total), what their WASH exposure was in the baseline round and in the post-intervention round

d_lab_wash = left_join(d_lab_de%>%
                         select(-c(village, redcap_event_name, data.row)), d_wash) # 1237

d_lab_wash_r2 = left_join(d_lab_de%>%
                            select(-c(village, redcap_event_name, data.row)), d_wash_r2)

d_lab_wash_b = rbind(d_lab_wash, d_lab_wash_r2)  

# LINKAGE AND EXPORT
#---------------------------------------------------------
# Export WASH dataset for all households
write.csv(d_wash_b, paste0(DirectoryDataOut, "./FINAL_FOR_SHARING/Household_WASH_BF.csv")) 

# Export WASH dataset for those with stool sample
write.csv(d_lab_wash_b, paste0(DirectoryDataOut, "./FINAL_FOR_SHARING/Household_individual_WASH_BF.csv")) 

# Export cleaned and linked stool collection
write.csv(dt_wash, paste0(DirectoryDataOut, "./FINAL_FOR_SHARING/Household_stool_WASH_BF.csv")) 

# Export cleaned and linked data for MSM model
dfls0 = dt_filtered_wash
dfls0 = dfls0 %>%
  rename(
    intervention.text = intervention_text
  ) %>%
  mutate(
    intervention.text = ifelse(intervention.text=="contrôle", "control", intervention.text)
  )

df = d%>%select(c(menage_id, n.householdmember, n.child.0to5, n.households.concession))%>%filter(!duplicated(menage_id))

dfls0 = left_join(dfls0 %>%select(-c(n.householdmember, n.households.concession)),df) 

save(dfls0, file=paste0(DirectoryDataOut, "./use_in_analyses/bf_esbl0123_long_all.rda")) 

dfls0complete = dt_filtered_complete_wash
dfls0complete = dfls0complete %>%
  rename(
    intervention.text = intervention_text
    
  ) %>%
  mutate(
    intervention.text = ifelse(intervention.text=="contrôle", "control", intervention.text)
  )


save(dfls0complete, file=paste0(DirectoryDataOut, "./use_in_analyses/bf_esbl0123_long_completecases.rda")) 


# Spit out typos in menage_id pella SES
typosSES = data.frame(menage_id = ses_pella$menage_id[which(!ses_pella$menage_id%in%c(d$menage_id))], check="typo_ses")
missingSES = data.frame(menage_id = unique(d$menage_id[which(!d$menage_id%in%c(ses_pella$menage_id) & d$village_name=="PELLA")]),
                        check="notinses")

checkses = rbind(typosSES, missingSES)


write.csv(checkses, paste0(DirectoryDataOut, "./need_checking/ses_pella_ids_to_check.csv")) 

#---------------------------------------------------------
# DESCRIPTIVES
#----------------------------------------------------------
# R0 characteristics 
table_wash = table1(~ n.householdmember +
                      n.child.0to5 +
                      n.households.concession +
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
                          factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text), data=d_wash)
table_wash

# R2 characteristics 
table_wash_r2 = table1(~ n.householdmember +
                      n.child.0to5 +
                      n.households.concession +
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
                      factor(q21.when.animal.ill.dies.eat.at.home)| factor(intervention.text), data=d_wash_r2)
table_wash_r2

# PLOT INTERVENTION TIMES AND SAMPLING TIMES
# Colours
palette <- wes_palette("Darjeeling1", n = 5)
palette2 <- wes_palette("BottleRocket2", n = 1)
palette3 <- wes_palette("GrandBudapest1", n = 2)[2]
palette4 <- wes_palette("BottleRocket2", n = 2)[2]
palette5 = c(palette3, palette[2],palette2,palette[5], palette[4],palette4)


# INTERVENTION DATES
# Round 1: 13 Feb - 26 April 2023
# Round 2: 8 May - 13 June 2023
# Round 3: 21 Aug - 3 Nov 2023

village_order <- dfls0complete %>%
  filter(redcap_event_name == "round_0_arm_1") %>%
  group_by(village_name) %>%
  dplyr::summarize(earliest_date = min(date.use, na.rm = TRUE), .groups = "drop") %>%
  arrange(earliest_date)
village_order
print(n=22,village_order)
unique(dfls0$village_name)

# Reorder village_name in the original dataframe
dfls0complete <- dfls0complete %>%
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
p = ggplot(dfls0complete %>% filter(intervention.text=="intervention"), aes(x = date.use, fill = redcap_event_name)) +
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

p2 = ggplot(dfls0complete %>% filter(intervention.text == "control"), aes(x = date.use, fill = redcap_event_name)) +
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

ggsave(file = "./Output/Figures/BF/samplingtimes_villages_UPDATED.pdf", 
       plot = combined_plot, width = 32, height = 15)  



