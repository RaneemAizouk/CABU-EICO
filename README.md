**Different databases:** 

Burkina Faso
  1) bf_esbl0123_long_all.rda                 = Cleaned feacal sample dataset linked with households survey containing the observations of all individiuals
  2) bf_esbl0123_long_completecases.rda       = Cleaned feacal sample dataset linked with households survey containing the observations of all individiuals with four observations (complete follow up)
  3) Household_WASH_BF.csv                    = Household survey WASH observations measured pre- and post-intervention for all households, including those where no stool samples were taken
  4) Household_stool_WASH_BF.csv              = Household survey WASH observations measured pre- and post-intervention for households where stool samples were taken
  
     
**Different scripts (script folder):**
1) data_set_preparation_spline.R            = Script that converts the stool sampling data in a stan format for the markov model
2) scenario_fit_comparison.R                = Script that compares fit across the different intervention scenarios. Can be used for observed and simulated data
3) descriptive_WASH.R                       = Script that cleans and describes the households survey results, both at individual and household level
4) model_output_summary.R                   = Script that summarises the main model output. Can be used for all intervention scenarios
5) CoxModel_WASH.R                          = Script that fits a Cox proportional hazard model to the complete stool sampling data (i.e. of those with four observations)
5) analyses_change_WASH.R                   = Script that estimates the change in WASH indicators
6) figures_and_tables.R                     = Script that develops the manuscript figures and tables

**Other (Documentation folder):**
1) Interpretation antibiogramme des isolats PORTAGE asymptomatique_ESBL_E. coliKlebsielle.docx
**Data Analysis scripts:**
1. "/0_clean_and_describe_hh_survey.R" is code to 
2. "/0_data_prep_markov_model.R"
3. "1_descriptive_WASH.R"
4. "2_scenario_fit_comparison.R"
5. "3_model_output_summary.R"
6. "4_analyses_change_WASH.R"
7. "5_Data_prep_coxme.R" is code to prepare longitudinal stool sampling and household survey data for
    Cox mixed-effects modeling (`coxme`)
8. "5_CoxModel_WASH.R" is code to fit Cox proportional hazards mixed-effects models (`coxme`) to estimate associations between household WASH     indicators and ESBL-E acquisition using longitudinal stool sampling data.
9. "6_figures_and_tables.R"
