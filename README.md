# Optimising Community Antibiotic Use and Infection Control With Behavioural Interventions in Burkina Faso and DR Congo (CABU-EICO project)
This repository contains the study protocol, publicly available data, and code underlying the following study:

Aizouk R. et al., *Targeting community-level drivers of antimicrobial resistance in sub-Saharan Africa: The effect of a community-based intervention bundle on household transmission of ESBL-E in rural Burkina Faso - a cluster-randomised trial*. MedRxiv. doi: https://doi.org/10.64898/2025.12.15.25342269 

# Study registration
https://clinicaltrials.gov/ct2/show/NCT05378880

# Project summary
We developed and evaluated a behavioural intervention bundle, targeting any community-level healthcare or medicine providers and communities, to optimise antibiotic use and improve hygiene, and hence reduce AMR transmission. After a 6-month local co-development phase, the intervention was implemented over 12 months in 22 clusters (villages or neighbourhoods) within health demographic surveillance sites in Nanoro, Burkina Faso and Kimpese, DR Congo. In a cluster RCT, we compared the 22 intervention with 22 control clusters. 

The *primary outcomes* are the change in Watch antibiotic provision from informal and formal medicine providers and change in patient management, assessed via patient exit interviews and simulated client visits. These are presented in a joint submission to this manuscript.

The *secondary outcomes*  are the changes in AMR acquisition dynamics and hygiene practices, which were assessed using a pre-/post-intervention household survey collecting data on Water access, Sanitation and Hygiene (WASH) exposures, and repeated microbiological (stool) sampling among household members (collected 3-months before, at intervention start, 3 months post, and 9 months post-intervention start. These data were collected and evaluated in the Burkina Faso site using a continious-time markov modelling framework. These are presented in this repository.

# Simulated data
 1)  "Public_data/Simulated/Simulated_data_noseasonality_stan_data.rds : Simulated dataset without seasonality in stan format
 2)  "Public_data/Simulated/Simulated_data_seasonality_stan_data.rds : Simulated dataset with seasonality in stan format
    
# Observed data
  1) "/Public_data/Observed/bf_esbl0123_long_all.rds"                 : Cleaned and anonymised feacal sample dataset linked with the CABU-EICO households survey containing the observations of all individiuals.
  2) "/Public_data/Observed/bf_esbl0123_long_completecases.rds"       : Cleaned and anonymised feacal sample dataset linked with the CABU-EICO households survey containing the observations of all individiuals with four observations (complete follow up).
  3) "/Public_data/Observed/Household_WASH_BF.csv"                    : Cleaned and anonymised CABU-EICO household survey WASH observations measured pre- and post-intervention for all households, including those where no stool samples were taken.
  4) "/Public_data/Observed_data/Household_stool_WASH_BF.csv"              : Cleaned and anonymised CABU-EICO Household survey WASH observations measured pre- and post-intervention for households where stool samples were taken.
  
# Scripts to generate simulated data
1) "Scripts/0_create_simulated_data_no_seasonality.R": Generates a simulated dataset without a seasonal pattern in ESBL-E acquisition
2)  "Scripts/0_create_simulated_data_seasonality.R": Generates a simulated dataset a seasonal pattern in ESBL-E acquisition
   
# Scripts to clean CABU-EICO household survey and stool collection data  
These R scripts clean and generate the publicly available anonymized dataframes available on this repository.

 1)  "/Scripts/1_clean_and_describe_hh_survey.R": Cleans, and links the CABU-EICO household survey stool collection data, and WASH survey data.
 2)  "/Scripts/1_data_prep_markov_model.R": Anonymises and prepares the cleaned CABU-EICO longitudinal household survey data for the continuous-time multi-state (Markov) modelling framework in stan as well as the Cox modelling analyses.
 3)  "/Scripts/1_Data_prep_coxme.R": Prepares the cleaned CABU-EICO longitudinal household survey data for Cox proportional hazards mixed-effects modelling.
    
# Scripts with Stan model codes for all intervention scenarios
The following Stan-based model scripts implement the Bayesian continuous-time multi-state (Markov) models used to evaluate the intervention and seasonal effects on ESBL-E transmission dynamics:

| Script | Scenario | Intervention structure | Seasonal effect | Intervention affects |
|------|---------|------------------------|----------------|----------------------|
| `/Scripts/Final_Stan_Codes/S2_Two_Step_Sine_NonAdd_NonCol.R` | 1A (Base case) | Two-step | Sinusoidal | Acquisition only |
| `/Scripts/Final_Stan_Codes/S1_One_Step_Sine_NonAdd_NonCol.R` | 2A | One-step | Sinusoidal | Acquisition only |
| `/Scripts/Final_Stan_Codes/S2_Two_Step_Spline_NonAdd_NonCol.R` | 1B | Two-step | Sinusoidal | Acquisition only |
| `/Scripts/Final_Stan_Codes/S1_One_Step_Spline_NonAdd_NonCol.R` | 2B | One-step | Cubic B-spline | Acquisition only |
| `/Scripts/Final_Stan_Codes/S2_Two_Step_Sine.R` | 3A | Two-step | Sinusoidal | Acquisition and decolonisation |
| `/Scripts/Final_Stan_Codes/S1_One_Step_Sine.R` | 4A | One-step | Sinusoidal | Acquisition and decolonisation |
| `/Scripts/Final_Stan_Codes/S2_Two_Step_Spline.R` | 3B | Two-step | Cubic B-spline | Acquisition and decolonisation |
| `/Scripts/Final_Stan_Codes/S1_One_Step_Spline.R` | 4B | One-step | Cubic B-spline | Acquisition and decolonisation |



1) "/Scripts/Final_Stan_Codes/S2_Two_Step_Sine_NonAdd_NonCol.R": Model scenario 1A: Two-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only. This is the base case scenario
2) "/Scripts/Final_Stan_Codes/S1_One_Step_Sine_NonAdd_NonCol.R": Model scenario 2A: One-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only.
3) "/Scripts/Final_Stan_Codes/S2_Two_Step_Spline_NonAdd_NonCol.R": Model scenario 1B: Two-step intervention effect, sinusoidal seasonal effect, intervention assumed to affect acquisition only.
4) "/Scripts/Final_Stan_Codes/S1_One_Step_Spline_NonAdd_NonCol.R": Model scenario 2B: One-step intervention effect, flexible cubic B-spline seasonal effect, intervention assumed to affect acquisition only.
5) "/Scripts/Final_Stan_Codes/S2_Two_Step_Sine.R": Model scenario 3A: Two-step intervention, sinusoidal seasonality, intervention affects acquisition and decolonisation.
6) "/Scripts/Final_Stan_Codes/S1_One_Step_Sine.R": Model scenario 4A: One-step intervention, sinusoidal seasonality, intervention affects acquisition and decolonisation.
7) "/Scripts/Final_Stan_Codes/S2_Two_Step_Spline.R": Model scenario 3B: Two-step intervention, flexible cubic B-spline seasonality, intervention affects acquisition and decolonisation.
8) "/Scripts/Final_Stan_Codes/S1_One_Step_Spline.R": Model scenario 4B: One-step intervention, flexible cubic B-spline seasonality, intervention affects acquisition and decolonisation.

# Scripts used to analyse stan model output and further data analyses
These R scripts were used to analyse the markov model output, as well as the change in WASH and association between WASH and ESBL-E acquisition.

 1)  "/Scripts/2_scenario_fit_comparison.R": Computes and compares model diagnostics (including LOO, R-hat, and divergence statistics) for all model scenarios.
 2)  "/Scripts/3_analyses_change_WASH.R": Fits the quasi-Poisson regression model to estimate pre-/post intervention changes in household WASH indicators.
 3)  "/Scripts/4_CoxModel_WASH.R": Fits Cox proportional hazards mixed-effects models to estimate associations between household WASH indicators and ESBL-E acquisition. 
 9)  "/Scripts/5_figures_and_tables.R": Generates the final manuscript figures and tables by combining outputs from descriptive analyses, Markov models, and Cox regression analyses. 

# Related publications
1) [CABU-EICO study protocol – Trials (2024)](https://pubmed.ncbi.nlm.nih.gov/38281023/)
2) Ingelbeen B., Valia D. et al., *Effect of a community-based intervention bundle to improve antibiotic use and patient management in Burkina Faso and DR Congo: a cluster randomized controlled trial* (manuscript — to be added).


