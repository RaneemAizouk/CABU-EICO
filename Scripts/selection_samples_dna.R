#####################################################
# DATA EXPLORATION
#####################################################
# This code is selecting households per cluster to follow up over time for sequencing.
# The code is selecting from households  with at least 1 positive ESBL only
# There are only few households with no ESBLs (max 3 per cluster)


# 1 September 2023
# Last update: 2 November 2023
# Update on 5 February: realised that "ajoute" parameter in the bf_villages_cabu.xlsx was used for intervention_yn;
# incorrectly. Does not matter for sample selection as for all clusters, 2 or 3 households are selected for sequencing

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr)

# SET DIRECTORY
DirectoryData <- "./Data/BF/Raw"
DirectoryDataOut <- "./Data/BF/clean"

#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_manualchange_no_password.xlsx"))
car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-10-17_nopassword.xlsx"))
car_r1 = read_xlsx(paste0(DirectoryData, "/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24_nopassword.xlsx"), sheet=1)
car_r2 = read_xlsx(paste0(DirectoryData, "/CABUBPortageAsymptom_DATA_R1_R2_2024-01-24_nopassword.xlsx"), sheet=2)


#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-09-20_1301_manual_change_id_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-09-08_1301_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/CABUBPortageAsymptom_DATA_2023-05-04_manual_change_id_nopassword.xlsx"))
#car_r0 = read_xlsx(paste0(DirectoryData,"/bf_esbl_r0.xlsx"), sheet=2)
hh_lab_ids =  read_xlsx(paste0(DirectoryData,"/Correspondande-Code_Lab-ID_Menage.xlsx"))
names(hh_lab_ids)
names(car_r0)
names(hh_lab_ids) = c("household", "menage_id", "bras")

wash_r0 = read_xls(paste0(DirectoryData, "/WP4_WASH_07_09_2023_nopassword.xls"))
wash_r0_stool = wash_r0 %>% filter(!is.na(num_echantillon))

wash_r0_stool_sel = wash_r0_stool %>% select(menage_id,num_echantillon, date_enquete,groupe,
                                             nbre_enf_0_5ans,nbre_menage_conc,
                                             cs_id_individu, dob_age,dob, age, sexe)

villages = read_xlsx(paste0(DirectoryData, "/bf_villages_cabu.xlsx"))
names(villages) = c("village", "village_name","intervention_text","ajoute")

# Villages where no provider was present; thus no provider intervention
no_prov = c("zimidin", "GOULOURE", "MOGODIN")
villages$provider_present = ifelse(villages$village_name%in%no_prov, "no", "yes")

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

# which households with ESBL
car_r0$found_in_wash[which(car_r0$menage_id %in% wash_r0_stool_sel$menage_id)] = 1
car_r0$found_in_wash[is.na(car_r0$found_in_wash)] = 0
length(unique(car_r0$household[car_r0$found_in_wash==0])) # 28 households can not be found in WASH survey database; could be due to error's in ID or due to zero's that need to be removed, needs checking
unique(car_r0$menage_id[car_r0$found_in_wash==0])
unique(car_r0$household[car_r0$found_in_wash==0])


# Number of positives per household
hist(table(car_r0$household))
hhs = data.frame(table(car_r0$household))
names(hhs) = c("household", "hhsize") 

# Household sizes
hist(table(wash_r0_stool_sel$menage_id))

# Number of individuals per household 
hhs = wash_r0_stool_sel %>% 
  group_by(menage_id) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
summary(hhs$n)

hh_at_least_one_esbl = (car_r0 %>% group_by(village) %>%
  summarise(length(unique(household)))) # EA had 23 households in older dataset. I thought we would collect data from 12 households per cluster? 
                                        # Update: CABUBPortageAsymptom_DATA_2023-09-08_1301_manual_change_id_nopassword.xlsx contains updated IDs
                                        # reduced to <=13 per cluster    
                                        # Further updated on 20 09 2023: in CABUBPortageAsymptom_DATA_2023-09-20_1301_manual_change_id_nopassword.xlsx with correction of ids in lines 121-148
hh_per_cluster =  (wash_r0_stool_sel %>% group_by(village) %>%
summarise(length(unique(menage_id))))


hh_total = merge(hh_at_least_one_esbl, hh_per_cluster)

names(hh_total) = c("village", "hh_num_sampled", "hh_num_cluster")

# Number of households that did not have an ESBL
hh_total$hh_num_diff = hh_total$hh_num_cluster - hh_total$hh_num_sampled
sum(hh_total$hh_num_diff[hh_total$hh_num_diff>0]) # 24 households no ESBL 
sum(hh_total$hh_num_diff[hh_total$hh_num_diff>0])/sum(hh_total$hh_num_cluster) # = 5% no ESBL
print(hh_total) 

# The minus mean that for some of the clusters there are still some errors in the sample IDs of the 
# 20 sEPTEMBER 2023: 
# BASED ON BELOW CORRECTIONS NEEDED HAVE UPDATED CABUBPortageAsymptom_DATA_2023-09-08_1301_manual_change_id_nopassword.xlsx WITH CABUBPortageAsymptom_DATA_2023-09-20_1301_manual_change_id_nopassword.xlsx),
# carriage database after manually checking (("record_id_manual_change") as more household IDs 
# than there are households included. 
# WAITING FOR CONFIRMATION Daniel/Frank (20 September 2023)

# Check which ones still not correct (20 September 2023, updated these in )
# EA
# unique(car_r0$household[car_r0$village=="EA"]) # "EAA4801" should be "EAA04801" 
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="EA"])
# 
# # EF
# unique(car_r0$household[car_r0$village=="EF"]) # "EFJ2703" should be "EFJ02703" 
#  # "EFF05402" should be "EFF05102"
#  # "EFF0106" does not exist in wash_r0, correct?
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="EF"])
# 
# # EK
# unique(car_r0$household[car_r0$village=="EK"]) # "EKE01703" should be "EKA01703"
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="EK"])
# 
# # ET
# unique(car_r0$household[car_r0$village=="ET"]) # "ETE00201" should be "ETD01201"
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="ET"])
# 
# # SC
# unique(car_r0$household[car_r0$village=="SC"]) # "SCA0101" should be "SCA101"
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="SC"])
# 
# # SJ
# unique(car_r0$household[car_r0$village=="SJ"]) # #"SJA02002" should be #"SJG02002"  
# # "SJG00601" do not exist in wash_r0
# # "SJG02501" should be "SJG02701"
# unique(wash_r0_stool_sel$menage_id[wash_r0_stool_sel$village=="SJ"])


# Clean germe; below diameters in cancelled out code are EUCAST guidelines (this needs to be CLSI guidelines, 
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
        # esbl_pos = ifelse(testesbl, 1, 0))

        
#sapply(s, function(x) table(x))  

# Number of cases positive
table(car_r0$esbl_pos) # These are the ESBL positive patients based on cetriax_or_cefota, 780
table(car_r0$germe_c, car_r0$esbl_pos)

table(car_r0$esbl_pos==1 & car_r0$testesbl==1)

# Select just esbl positives
d = car_r0 %>% filter(esbl_pos == 1)

# Number of individuals with two morphotypes
d %>% 
  group_by(germe_c) %>%
  summarise(n = n()) #%>%
  #mutate(freq = n / sum(n))


# Number of individuals per village positive (however no village populations)
vill = d %>%
  group_by(village) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

hist(vill$n)
summary(vill$n)
vill$n/(12*5) # rough estimate of fraction ESBL positive
hist(vill$n/(12*5))

d %>%
  group_by(intervention_text) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

d %>%
  group_by(germe_c, intervention_text) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# Number of individuals per household positive (however no household size)
hh = d %>% 
  group_by(household) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) 
hist(hh$n)
hist((hh$n)/5) # rough estimate of fraction positive (we don't have denominator data here)

summary(hh$n) # median of 2 indivuals positive per household

table(d$morphotyp, useNA='always')/length(d$morphotyp)
table(car_r0$morphotyp,car_r0$esbl_pos, useNA='always')

#####################################################################
# Selection of households
#####################################################################

# So 671 unique individuals with an E. coli ESBL over 2010 sampled indvidiuals = 56%
# 1209 is based on earlier presentation of Brecht, needs to be verified (C:\Users\evankleef\OneDrive - ITG\Documenten\GitHub\ITG\Other_projects\AMR\JPI-AMR\Data\BF\CABU_carriageresults.pptx)

# For BF 750 can be sequenced in total
# 1.	First round ESBL E. coli positive = 55.5% (n=671/1210)
# 2.	First round ESBL E. coli with 2 morphotypes = 7.9% (n= 96/1210)
# 3.	First round ESBL E. coli (1 + 2 morphotypes) = 63.4% (n=(96+671)/1210)
# 4.	Max number of complete households that can be sampled per round if 1 isolate per individual 
# -	750/4 = 187 per round
# -	187/(5 members per hh*0.565) = 66 
# -	66, i.e. 33 per intervention arm
# -	With 11 clusters per arm, that is 3 households per cluster

# If we include 2 isolates per individual
# -	187/(5 members per hh*0.646) = 58 
# -	58, i.e. 29 per intervention arm
# -	With 11 clusters per arm, that is 2 to 3 households per cluster
# - 14 clusters with 3 households, 8 clusters with 2 households

# Random selection of households per cluster
# Data to use
# Spit out data of only positive individuals
dat = car_r0 %>% filter(esbl_pos == 1 & germe_c %in% c("e.coli", "e.coli_2"))

dat %>% group_by(village) %>%
  summarise(length(unique(household))) # Number of households per cluster should be the same is in wash_r0 

seed = 120
set.seed(seed)

# First randomly select clusters with 3 vs 2 villages to include 
#dat_three_hh_i = sample(unique(dat$village[dat$intervention_yn=="Oui"]), 7)
dat_three_hh_i = sample(unique(dat$village[dat$ajoute=="Oui"]), 7) # Used the wrong variable here
# But does not matter as for all villages, 2 or 3 households are included

set.seed(seed)
#dat_three_hh_ni = sample(unique(dat$village[dat$intervention_yn=="Non"]), 7)
dat_three_hh_ni = sample(unique(dat$village[dat$ajoute=="Non"]), 7)

dat_three_hh = c(dat_three_hh_i,dat_three_hh_ni)
dat$three_hh = ifelse(dat$village %in%dat_three_hh, 1, 0)
table(dat$village,dat$three_hh)

final_hh_include = NULL
i = 1

for(c in unique(dat$village)){
  d = dat %>% filter(village == c)
  households = unique(d$household)
  print(paste0("cluster = ", i))
  print(paste0("number of household = ", length(households)))
  print(households)
  set.seed(seed)
    if(c %in% dat_three_hh){
      hh_include = sample(households,3)  
    }else{
      hh_include = sample(households,2)  
    }
  print(paste0("households included =", length(hh_include)))
  print(paste0("household ids included ="))
  print(hh_include)
  final_hh_include = c(final_hh_include, hh_include)
    i = i+1
  }
   
final_hh_include     

dat$dna_extraction = ifelse(dat$household %in%final_hh_include, 1, 0)
table(dat$dna_extraction)
table(dat$germe_c,dat$dna_extraction)

# see how many with intervention how many without
table(dat$intervention_text, dat$dna_extraction) # This seems pretty balanced?

ch_dna = dat %>% group_by(intervention_text, dna_extraction)%>%
  summarise(n=n())
ch_dna

dat %>% filter(dna_extraction == 1) %>%
  group_by(intervention_text)%>%
  summarise(n=length(unique(household))) # So 30 households in controle groups; 28 in intervention groups

check = dat %>% filter(dna_extraction == 1) %>%
  group_by(intervention_text, village)%>%
  summarise(n=length(unique(household)))

# Add variable to original data
car_r0$dna_extraction = ifelse(car_r0$household %in%final_hh_include & car_r0$esbl_pos==1 & car_r0$germe_c %in% c("e.coli", "e.coli_2"), 1, 0)
table(car_r0$dna_extraction)
table(car_r0$germe_c,car_r0$dna_extraction)


# Samples to use for extraction
dat_extract = dat %>%filter(household %in% final_hh_include)
#dat_extract_ids = dat_extract %>% select(household,record_id,record_id_manual_change, id_ecantillon)
dat_extract_ids = dat_extract %>% select(household,record_id, id_ecantillon)

# checks
table(dat_extract$morphotyp) # 22 individuals with 2 extracts
length(unique(dat_extract$village)) # 22 villages
length(unique(dat_extract$household)) # 58 households


# Number of individuals per household
sort(table(dat_extract$household))
par(mfrow=c(1,2))
hist(table(dat_extract$household), breaks=unique(table(dat_extract$household)), 
     main="Selected households", xlab= "number of household members")
hist(table(dat$household), breaks=unique(table(dat$household)), main="All households",
     xlab= "number of household members")

# Check if any selected with no WASH link
check = unique(car_r0$menage_id[car_r0$found_in_wash==0])

######################################################################
# UPDATE # Update 5 February 2024: 
# Can use Correspondande-Code_Lab-ID_Menage.xlsx to link wash survey now
# Does not matter for DNA extractions as was merely used to check
# if IDs could be identified back 
# Will add household 2 id here, so not to confuse
######################################################################

# Link ids of dna extracts with hh survey ids
dat_extract_ids = left_join(dat_extract_ids, hh_lab_ids, by="household")
# Check if all linked
table(dat_extract_ids$bras, useNA= "always")
table(dat_extract_ids$menage_id)

dat_extract = left_join(dat_extract, hh_lab_ids, by="household")
# Check if all linked
table(dat_extract$bras, useNA= "always")
dat_extract = dat_extract %>%
  rename(menage_id_selfmade = menage_id.x,
         menage_id = menage_id.y)

# check if ids can be identified in R1 and R2
# Unique ids R1 
r0_ids = unique(dat_extract$menage_id)
r1_ids = unique(car_r1$menage_id)
r2_ids = unique(car_r2$menage_id)

table(r1_ids%in% r0_ids) # 2 households can not be found back
table(r2_ids%in% r0_ids) # 2 households can not be found back
# Have to wait if household members can also be found back

# Save data
file <- paste0(DirectoryDataOut, "/bf_ecoli_esbl.xlsx")
write_xlsx(car_r0, file)

# Save inclusion list
file <- paste0(DirectoryDataOut, "/bf_sample_ids_dna.xlsx")
write_xlsx(dat_extract_ids, file)
