# Optimising community antibiotic use and environmental infection control with behavioural interventions in rural Burkina Faso and DR Congo (CABU-EICO)
This repository contains the publicly available data, and code underlying the following study:

Aizouk R. et al., *Targerting community-level drivers of antimicrobial resistance in sub-Saharan Africa: The effect of a community-based intervention bundle on household transmission of ESBL-E in rural Burkina Faso - a cluster-randomised trial*. MedRxiv. doi: https://doi.org/10.64898/2025.12.15.25342269 

# Study registration
https://clinicaltrials.gov/ct2/show/NCT05378880

# Project summary
We conducted a cluster-randomised controlled trial in 22 village clusters in Nanoro district, Burkina Faso. We enrolled 12 randomly selected households per cluster to assess intervention impact on ESBL-E household-transmission. The intervention comprised three rounds at three-month intervals and combined WHO AWaRe–based educational feedback for formal and informal medicine providers with a community-wide WASH and antibiotic-use behaviour change campaign. Consenting household members provided stool samples before, during, and after intervention rollout, alongside a pre–post household WASH survey. We estimated intervention effects on ESBL-E acquisition using Bayesian Markov models. Cox frailty models assessed associations between WASH exposures and acquisition. Finally, using quasi-poisson regression models, we estimated changes in pre-specified WASH indicators pre-/post intervention. 

These are the secondary outcomes from the CABU-EICO trial, primary outcomes, the pre- to post-intervention change in prevalence of dispensing Watch antibiotics, and patient management at medicine stores or health centers are presented in the joint submission to this work. 

# Data
  1) bf_esbl0123_long_all.rda                 = Cleaned feacal sample dataset linked with households survey containing the observations of all individiuals
  2) bf_esbl0123_long_completecases.rda       = Cleaned feacal sample dataset linked with households survey containing the observations of all individiuals with four observations (complete follow up)
  3) Household_WASH_BF.csv                    = Household survey WASH observations measured pre- and post-intervention for all households, including those where no stool samples were taken
  4) Household_stool_WASH_BF.csv              = Household survey WASH observations measured pre- and post-intervention for households where stool samples were taken
  
# Data cleaning scripts  
These R scripts can not show the cleaning process to generate the anonymized dataframes available on this repository.
 1)  0_clean_and_describe_hh_survey.R: Cleans, and links the CABU-EICO household survey stool collection data, and Water access, Sanitation and Hygiene (WASH) survey data.
 2)  0_data_prep_markov_model.R: Anonymises and prepares the cleaned CABU-EICO longitudinal household survey data for the continuous-time multi-state (Markov) modelling framework in stan as well as the Cox modelling analyses.
 3)  0_Data_prep_coxme.R: Prepares the cleaned CABU-EICO longitudinal household survey data for Cox proportional hazards mixed-effects modelling.
    
# Stan model codes and scenario analyses
The following Stan-based model scripts implement the Bayesian continuous-time multi-state (Markov) models used to evaluate the intervention and seasonal effects on ESBL-E transmission dynamics:

1) S2_Two_Step_Sine_NonAdd_NonCol.R: Model scenario 1A: Two-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only. This is the base case scenario
2) S1_One_Step_Sine_NonAdd_NonCol.R: Model scenario 2A: One-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only.
3) S2_Two_Step_Spline_NonAdd_NonCol.R: Model scenario 1B: Two-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only.
4) S1_One_Step_Spline_NonAdd_NonCol.R: Model scenario 2B: One-step intervention effect, flexible cubic B-spline seasonal effect, intervention assumed to affect acquisition only.
5) S2_Two_Step_Sine.R: Model scenario 3A: Two-step intervention, sinusoidal seasonality, intervention affects acquisition and decolonisation.
6) S1_One_Step_Sine.R: Model scenario 4A: One-step intervention, sinusoidal seasonality, intervention affects acquisition and decolonisation.
7) S2_Two_Step_Spline.R: Model scenario 3B: Two-step intervention, flexible cubic B-spline seasonality, intervention affects acquisition and decolonisation.
8) S1_One_Step_Spline.R: Model scenario 4B: One-step intervention, flexible cubic B-spline seasonality, intervention affects acquisition and decolonisation.

# Data analyses scripts
 1)  1_scenario_fit_comparison.R: Computes and compares model diagnostics (including LOO, R-hat, and divergence statistics) for all model scenarios.
 2)  2_analyses_change_WASH.R: Fits the quasi-Poisson regression model to estimate pre-/post intervention changes in household WASH indicators.
 3)  3_CoxModel_WASH.R: Fits Cox proportional hazards mixed-effects models to estimate associations between household WASH indicators and ESBL-E acquisition. 
 9)  4_figures_and_tables.R: Generates the final manuscript figures and tables by combining outputs from descriptive analyses, Markov models, and Cox regression analyses. 

# Related publications
1) [CABU-EICO study protocol – Trials (2024)](https://pubmed.ncbi.nlm.nih.gov/38281023/)
2) Ingelbeen B., Valia D. et al., *Effect of a community-based intervention bundle to improve antibiotic use and patient management in Burkina Faso and DR Congo: a cluster randomized controlled trial* (manuscript — to be added).


