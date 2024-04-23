#####################################################
# DESCRIPTIVES AND RANDOM-EFFECTS MODEL
#####################################################
# This code is cleaning the data and linking the different datasets

# 18 April 2024
# Author: Esther van Kleef

rm(list=ls())

# load package
pacman::p_load(readxl, writexl, lubridate, zoo, ggplot2, tidyverse, Hmisc, stringr,lme4,reshape2, 
               openxlsx, table1, flextable, magrittr, officer, summarytools, msm,
               survey)


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
  date.check.use = date - date.use,
  time = ifelse(redcap_event_name =="round_0_arm_1", 0,
                ifelse(redcap_event_name == "round_1_arm_1", 1,
                       ifelse(redcap_event_name == "round_2_arm_1", 2,3))))

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

dfls = dfls %>% filter(menage_id_member %in% complete.cases.id)

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
hh_pos = dfls %>% group_by(time, menage_id) %>%
  summarise(nesbl = sum(esble)
  )
hh_pos = hh_pos[order(hh_pos$menage_id, hh_pos$time),]

dfls = left_join(dfls,hh_pos, by=c("menage_id", "time"))
dfls$nesbl_tminus = lag(dfls$nesbl)



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
length(unique(dfls$menage_id_member)) # 754
dfls = dfls[order(dfls$menage_id_member, dfls$time),]

# Dataset for msm model
dfls = dfls %>% 
  mutate(esble = as.numeric(esble)+1) %>%
  select(-c(date.check.collect, date_conserv, intervention_text,X)) 

reorder_col = c("time", "redcap_event_name","time.start0","time.start1","start0","start1", "days","date.use","menage_id_member", "menage_id",
                "village_name","intervention.text", "age", "agegr10", "sexe", "n.householdmember","esble","nesbl","nesbl_tminus",
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
    n = length(unique(menage_id_member)),
    prev = esblsum/n
  )


# Also create a dataset with baseline 0 round included
d = dfl %>% filter(menage_id_member%in%complete.cases.id & redcap_event_name=="round_0_arm_1")

d = left_join(d, dfls%>%select(menage_id_member,start0, start1, time.start0,time.start1), multiple="first")
d = d %>% select(-c(X)) %>%
  mutate(days = as.numeric(date.use - start1))

hh_pos = d %>% group_by(time, menage_id) %>%
  summarise(nesbl = sum(esble)
  )
hh_pos = hh_pos[order(hh_pos$menage_id, hh_pos$time),]

d = left_join(d,hh_pos, by=c("menage_id", "time"))
d$nesbl_tminus = lag(d$nesbl)
names(d)
desired_order

d = d[, desired_order]
dfls0 = rbind(d,dfls)
dfls0 = dfls0[order(dfls0$menage_id_member, dfls0$time),]

# check the days variable
dfls0 %>% group_by(time) %>%
  ggplot(., aes(x=days)) + geom_histogram() +
  facet_wrap(.~redcap_event_name)

dfls0 %>% group_by(time) %>%
  summarise(median=median(days, na.rm=T),
            q1 = quantile(days,probs=0.25, na.rm=T),
            q3 = quantile(days,probs=0.75,na.rm=T))

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

reorder_col = c("redcap_event_name","menage_id_member", "menage_id",
                "village_name","intervention.text", "age", "agegr10", "sexe", "n.householdmember",
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
save(dfa, file="./Data/BF/clean/use_in_analyses/bf_esbl_wide.rda")
#load("./Data/BF/clean/use_in_analyses/bf_esbl_wide.rda")

save(dfls0, file="./Data/BF/clean/use_in_analyses/bf_esbl0123_long.rda")
save(dfls, file="./Data/BF/clean/use_in_analyses/bf_esbl123_long.rda")

#############################################################
# DESCRIPTIVES
#############################################################
table(df$r0.esble)
dfa0 = dfa %>% filter(!is.na(r0.esble))


# Table 1

# trial
wash_r0_table1 = table1(~ age + sexe + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          antibiotic + 
                          watch + 
                          q1.diar.prev.water.pot.covered +
                          q1.diar.prev.no.finger.in.waterglass +  
                          q1.diar.prev.utensil.to.take.water.from.pot + 
                          q1.diar.prev.cover.food +
                          q1.diar.prev.boil.water + 
                          q1.diar.prev.filter.water + 
                          q1.diar.prev.other+q1.diar.prev.cant.be.avoided + 
                          q1.diar.prev.dont.know + 
                          q2.main.water.source.dry.c + 
                          q3.main.water.source.rainy.c + 
                          q4.cans.storage.water + 
                          q5a.cans.storage.water.closed.when.filled +
                          q5b.cans.storage.water.closed.wen.empty + 
                          q5c.cans.cleaned.before.reuse +
                          q6.treatment.water.c + 
                          q7.principle.defication + 
                          q8.other.defecation.flush.toiled.septic + 
                          q8.other.defecation.pit.latrine.ventilation + 
                          q8.other.defecation.pit.latrine.slab + 
                          q8.other.defecation.pit.latrine.no.slab +
                          q8.other.defecation.open.defecation + 
                          q8.other.defecation.other +                                
                          q8.other.defecation.none + 
                          q9.shared.toilet.c + 
                          q10.n.shared.toilet +
                          q11.toilet.last.cleaned + 
                          q12.disposal.child.stool + 
                          q13.disposal.latrine.pit +
                          q14.handwashing.product+ 
                          q15.handwashing.defecation + 
                          q16.handwashing.meals +  
                          q17.animals.around.household + 
                          q18.animal.inside.cow +  
                          q18.animal.inside.sheep.goat + 
                          q18.animal.inside.pig + 
                          q18.animal.inside.donkey.horse + 
                          q18.animal.inside.chicken.goose.duck +  
                          q18.animal.inside.other +  
                          q19.animal.inside.cow + q19.animal.inside.sheep.goat +
                          q19.animal.inside.pig +
                          q19.animal.inside.donkey.horse +
                          q19.animal.inside.chicken.goose.duck + 
                          q19.animal.inside.other +   
                          q20.animal.excrement.floor +    
                          q21.when.animal.ill.treatment.with.vet + 
                          q21.when.animal.ill.treatment.without.vet +
                          q21.when.animal.ill.bucher +
                          q21.when.animal.ill.eat.meat.at.home +     
                          q21.when.animal.ill.burie.dispose +
                          q21.when.animal.ill.autre +     
                          eau.assainissement.hygine.complete| factor(intervention.text), data=dfa0)
wash_r0_table1

t1flex(wash_r0_table1) %>% 
  save_as_docx(path="./Output/Tables/wash_r0_table1.docx")


# Effect of the intervention M3?
wash_r2_table1 = table1(~ factor(r2.esble) + age + sexe + 
                          antibiotic + 
                          watch  + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          q1.diar.prev.water.pot.covered +
                          village_name | factor(intervention.text), data=dfa)
wash_r2_table1


# Effect of the intervention M9?
wash_r3_table1 = table1(~ factor(r3.esble) + age + sexe + 
                          n.householdmember +
                          n.child.0to5 +
                          n.households.concession +
                          q1.diar.prev.water.pot.covered +
                          village_name | factor(intervention.text), data=dfa)
wash_r3_table1



###############################################################
# FIGURES
###############################################################

# Antibiotic use per village


# Per household, how many positive
d0 = df0 %>% group_by(village_name, intervention_text, menage_id, r0.esble) %>%
  summarise(n = n())

samples_per_hh = df0 %>% group_by(menage_id) %>%
  summarise(n_samples = n())

hh_size = as.data.frame(cbind(df0$menage_id,df0$n.householdmember))
names(hh_size) = c("menage_id","hh_size")
hh_size = hh_size[!duplicated(hh_size),]

d = left_join(d0, hh_size, by="menage_id") %>%
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
d_pos = d %>% filter(r0.esble==1)
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

# Intervention vs controle groups 
d_pos %>% group_by(intervention_text) %>%
  summarise(mean_cor = mean(f_pos_cor, na.rm=T),
            median_cor = median(f_pos_cor,na.rm=T),
            q1_cor = quantile(f_pos_cor,probs=c(0.25), na.rm = T),
            q3_cor = quantile(f_pos_cor, probs=c(0.75), na.rm = T),
            mean = mean(f_pos, na.rm=T),
            median = median(f_pos,na.rm=T)) # seems rather similar (luckily)

# Distribution %ESBL_pos per cluster
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

# Distribution %ESBL_pos per cluster
dp = ggplot(d_pos, aes(x=f_pos_cor, group=village_name, fill=village_name)) + 
  geom_density(aes(f_pos_cor, ..scaled..))+ 
  facet_wrap(.~village_name)+ geom_vline(data=mean, aes(xintercept=mean))+
  labs(x="%positive within hh", y="Density")
dp


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
# then as a sensitivity analysis consider only thoe with two initial negative swab
dfa_early = dfa %>% filter(r1.esble !=0)
dfa_early_sens = dfa %>% filter(r1.esble !=0 & r0.esble!=0)

dfa_late = dfa %>% filter(r2.esble !=0)

# first with no random effects or other covariates
m0 <- glm(r2.esble ~ intervention.text, data=dfa_early,  family = binomial(link = "logit"))
summary(m0)
exp(m0$coefficients)
exp(confint(m0)) # # Effect intervention: OR = 

# then with covariates and random effects 
m1.early <-glmer(r2.esble ~ intervention.text + (1|menage_id), data=dfa_early,  family = binomial(link = "logit"))
summary(m1.early)
exp(fixef(m1.early))
exp(confint(m1.early, method="Wald"))  # Effect intervention: OR = 

# then with covariates and random effects 
m1.late<-glmer(r3.esble ~ intervention.text + (1|menage_id), data=dfa_late,  family = binomial(link = "logit"))
summary(m1.late)
exp(fixef(m1.late))
exp(confint(m1.late, method="Wald"))  # Effect intervention: OR = 


# Sensitivity analyses where individuals are both negative at R0 and R1
# first with no random effects or other covariates
m0.s <- glm(r2.esble ~ intervention.text, data=dfa_early_sens,  family = binomial(link = "logit"))
summary(m0.s)
exp(m0.s$coefficients)
exp(confint(m0.s)) # # Effect intervention: OR =
# 
# # then with covariates and random effects 
m1.s<-glmer(r2.esble ~ intervention.text + (1|menage_id), data=dfa_early_sens,  family = binomial(link = "logit"))
summary(m1.s)
exp(fixef(m1.s))
exp(confint(m1.s, method="Wald"))  # Effect intervention: OR =
# 


# Now fit both together
m2 <- glmer(esble-1 ~ intervention.text + intervention.text*factor(time) + (1|menage_id_member), data=dfls,  family = binomial(link = "logit"))
summary(m2)
exp(fixef(m2))
exp(confint(m2, method="Wald"))  # Effect intervention: OR = 1.13 (0.83 - 1.53); then each time step, the intervention effect is stronger

# Change in log odds = coef for intervention + coef intervention:time
change_lo_m3 = fixef(m2)[2] + fixef(m2)[3] # m3
prob_in_m3 =1/(1+exp(-change_lo_m3))
prob_in_m3 # Change in log odds compared to T = 1 (intervention start)

change_lo_m9 = fixef(m2)[2] + fixef(m2)[4] # m9
prob_in_m9 =1/(1+exp(-change_lo_m9))
prob_in_m9


# Define a function to calculate the change in log odds
calculate_log_odds_change <- function(coef_intervention, coef_interaction, time_factor) {
  # Calculate the change in log odds for each level of time
  change_in_log_odds <- coef_intervention +
    coef_interaction
  
  return(change_in_log_odds)
}

# Resample the dataset and fit the model
bootstrap_results_r2 <- replicate(100, {
  # Sample with replacement
  bootstrap_sample <- sample_n(dfls, size = nrow(dfls), replace = TRUE)
  
  # Fit the model on bootstrap sample
  bootstrap_model <- glmer(esble-1 ~ intervention.text * factor(time) + (1 | menage_id_member), 
                           data = bootstrap_sample,
                           family = binomial)
  
  # Extract coefficients from the model
  coef_intervention <- fixef(bootstrap_model)["intervention.textintervention"]
  coef_interaction <- fixef(bootstrap_model)["intervention.textintervention:factor(time)2"]
  
  # Get the levels of time factor
  time_levels <- levels(factor(bootstrap_sample$time))[2]
  
  # Calculate change in log odds for each level of time
  log_odds_change <- calculate_log_odds_change(coef_intervention, coef_interaction)

  
  return(log_odds_change)
})

bootstrap_results_r2
# Calculate confidence intervals
boot_ci_r2 <- quantile(bootstrap_results_r2, c(0.025, 0.975))


# Resample the dataset and fit the model
bootstrap_results_r3 <- replicate(100, {
  # Sample with replacement
  bootstrap_sample <- sample_n(dfls, size = nrow(dfls), replace = TRUE)
  
  # Fit the model on bootstrap sample
  bootstrap_model <- glmer(esble-1 ~ intervention.text * factor(time) + (1 | menage_id_member), 
                           data = bootstrap_sample,
                           family = binomial)
  
  # Extract coefficients from the model
  coef_intervention <- fixef(bootstrap_model)["intervention.textintervention"]
  coef_interaction <- fixef(bootstrap_model)["intervention.textintervention:factor(time)3"]
  
  # Get the levels of time factor
  time_levels <- levels(factor(bootstrap_sample$time))[3]
  
  # Calculate change in log odds for each level of time
  log_odds_change <- calculate_log_odds_change(coef_intervention, coef_interaction)
  
  
  return(log_odds_change)
})

bootstrap_results_r3

# Calculate confidence intervals
boot_ci_r3 <- quantile(bootstrap_results_r3, c(0.025, 0.975))
boot_ci_r3

##############################################
# MSM model
##############################################

# Define the transition intensity parameters
lambda <- 0.02  # Transition intensity from state 1 to state 2 (acquisition)
rho <- 1/(5*30)  # Transition intensity from state 2 to state 1 (decolonisation)

# Create the Q matrix
Q <- matrix(0, nrow = 2, ncol = 2)
Q[1, 1] <- 0
Q[2, 2] <- 0
Q[1, 2] <- lambda
Q[2, 1] <- rho
diag(Q) <- -rowSums(Q)

rownames(Q) <- colnames(Q) <- c("no_esbl", "esbl")
crudeinits.msm(esble ~ days, menage_id_member, data=dfls, qmatrix=Q)


# Fit the multi-state model
msm_fit <- msm(
  formula = esble ~ days,          # Model formula
  covariates = list("1-2" = ~factor(intervention.text)+nesbl_tminus, "2-1" = ~factor(intervention.text)),
  subject = menage_id_member,    # ID variable
  data = dfls,                     # Dataset
  qmatrix = Q,                     # Transition intensity matrix
  gen.inits = F
 # method = "BFGS"                  # Optimization method
)

# Summarize the fitted model
summary(msm_fit)

qmatrix.msm(msm_fit)
pmatrix.msm(msm_fit, t=max(dfls$days))
sojourn.msm(msm_fit)
