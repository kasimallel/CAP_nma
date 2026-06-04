These files contain the codes (i.e., functions) used in R for the bayesian NMA.

Functions are called as below:

#----------------------------------------#
#IV. Outcomes for adult population, non-severe:
#----------------------------------------#
#BAYESIAN:
#ITT populations:
res_mort_itt_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "Mortality", sample_size_var = "Overall participants ITT")
res_clin_itt_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")
res_ae_itt_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")
res_sae_itt_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")
res_arae_itt_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
#CEvaluable populations
res_clin_cep_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "Treatment failure PP", sample_size_var = "Overall clinically eligible PP")
#MEvaluable populations
res_clin_mep_a_nonsev <- run_full_nma(dat_arm_a_nonsev, outcome = "Failure to Eradication microbiologically evaluable", sample_size_var = "Microbiologically evaluable PP")



#----------------------------------------#
#V. Outcomes for adult population, severe:
#----------------------------------------#
#BAYESIAN
#ITT populations:
res_mort_itt_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "Mortality", sample_size_var = "Overall participants ITT")
res_clin_itt_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")
res_ae_itt_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")
res_sae_itt_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")
res_arae_itt_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
#CEvaluable populations
res_clin_cep_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "Treatment failure PP", sample_size_var = "Overall clinically eligible PP")
#MEvaluable populations
res_clin_mep_a_sev <- run_full_nma(dat_arm_a_sev, outcome = "Failure to Eradication microbiologically evaluable", sample_size_var = "Microbiologically evaluable PP")



#----------------------------------------#
#VI. Outcomes for children population, non-severe:
#----------------------------------------#
#BAYESIAN
#ITT populations:
res_mort_itt_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "Mortality", sample_size_var = "Overall participants ITT")
res_clin_itt_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")
res_ae_itt_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")
res_sae_itt_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")
res_arae_itt_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
#CEvaluable populations
res_clin_cep_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "Treatment failure PP", sample_size_var = "Overall clinically eligible PP")
#MEvaluable populations
res_clin_mep_c_nonsev <- run_full_nma(dat_arm_c_nonsev, outcome = "Failure to Eradication microbiologically evaluable", sample_size_var = "Microbiologically evaluable PP")



#----------------------------------------#
#VII. Outcomes for children population, severe:
#----------------------------------------#
#BAYESIAN
#ITT populations:
res_mort_itt_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "Mortality", sample_size_var = "Overall participants ITT")
res_clin_itt_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")
res_ae_itt_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")
res_sae_itt_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")
res_arae_itt_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
#CEvaluable populations
res_clin_cep_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "Treatment failure PP", sample_size_var = "Overall clinically eligible PP")
#MEvaluable populations
res_clin_mep_c_sev <- run_full_nma(dat_arm_c_sev, outcome = "Failure to Eradication microbiologically evaluable", sample_size_var = "Microbiologically evaluable PP")


