# Bayesian Network Meta-Analysis (NMA) Functions

This repository contains the R functions used to conduct the Bayesian Network Meta-Analyses (NMA) for the Community-Acquired Pneumonia (CAP) project.

The main function is:

```r
run_full_nma(
  data,
  outcome,
  sample_size_var
)
```

where:

* `data`: arm-level dataset.
* `outcome`: outcome variable to be analysed.
* `sample_size_var`: population denominator corresponding to the outcome population.

---

# Adult Population – Non-Severe CAP

## Intention-to-Treat (ITT) Population

```r
res_mort_itt_a_nonsev  <- run_full_nma(dat_arm_a_nonsev, outcome = "Mortality", sample_size_var = "Overall participants ITT")

res_clin_itt_a_nonsev  <- run_full_nma(dat_arm_a_nonsev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")

res_ae_itt_a_nonsev    <- run_full_nma(dat_arm_a_nonsev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")

res_sae_itt_a_nonsev   <- run_full_nma(dat_arm_a_nonsev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")

res_arae_itt_a_nonsev  <- run_full_nma(dat_arm_a_nonsev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
```

## Clinically Evaluable (CE) Population

```r
res_clin_cep_a_nonsev <- run_full_nma(
  dat_arm_a_nonsev,
  outcome = "Treatment failure PP",
  sample_size_var = "Overall clinically eligible PP"
)
```

## Microbiologically Evaluable (ME) Population

```r
res_clin_mep_a_nonsev <- run_full_nma(
  dat_arm_a_nonsev,
  outcome = "Failure to Eradication microbiologically evaluable",
  sample_size_var = "Microbiologically evaluable PP"
)
```

---

# Adult Population – Severe CAP

## Intention-to-Treat (ITT) Population

```r
res_mort_itt_a_sev  <- run_full_nma(dat_arm_a_sev, outcome = "Mortality", sample_size_var = "Overall participants ITT")

res_clin_itt_a_sev  <- run_full_nma(dat_arm_a_sev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")

res_ae_itt_a_sev    <- run_full_nma(dat_arm_a_sev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")

res_sae_itt_a_sev   <- run_full_nma(dat_arm_a_sev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")

res_arae_itt_a_sev  <- run_full_nma(dat_arm_a_sev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
```

## Clinically Evaluable (CE) Population

```r
res_clin_cep_a_sev <- run_full_nma(
  dat_arm_a_sev,
  outcome = "Treatment failure PP",
  sample_size_var = "Overall clinically eligible PP"
)
```

## Microbiologically Evaluable (ME) Population

```r
res_clin_mep_a_sev <- run_full_nma(
  dat_arm_a_sev,
  outcome = "Failure to Eradication microbiologically evaluable",
  sample_size_var = "Microbiologically evaluable PP"
)
```

---

# Children – Non-Severe CAP

## Intention-to-Treat (ITT) Population

```r
res_mort_itt_c_nonsev  <- run_full_nma(dat_arm_c_nonsev, outcome = "Mortality", sample_size_var = "Overall participants ITT")

res_clin_itt_c_nonsev  <- run_full_nma(dat_arm_c_nonsev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")

res_ae_itt_c_nonsev    <- run_full_nma(dat_arm_c_nonsev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")

res_sae_itt_c_nonsev   <- run_full_nma(dat_arm_c_nonsev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")

res_arae_itt_c_nonsev  <- run_full_nma(dat_arm_c_nonsev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
```

## Clinically Evaluable (CE) Population

```r
res_clin_cep_c_nonsev <- run_full_nma(
  dat_arm_c_nonsev,
  outcome = "Treatment failure PP",
  sample_size_var = "Overall clinically eligible PP"
)
```

## Microbiologically Evaluable (ME) Population

```r
res_clin_mep_c_nonsev <- run_full_nma(
  dat_arm_c_nonsev,
  outcome = "Failure to Eradication microbiologically evaluable",
  sample_size_var = "Microbiologically evaluable PP"
)
```

---

# Children – Severe CAP

## Intention-to-Treat (ITT) Population

```r
res_mort_itt_c_sev  <- run_full_nma(dat_arm_c_sev, outcome = "Mortality", sample_size_var = "Overall participants ITT")

res_clin_itt_c_sev  <- run_full_nma(dat_arm_c_sev, outcome = "Treatment failure ITT", sample_size_var = "Overall participants ITT")

res_ae_itt_c_sev    <- run_full_nma(dat_arm_c_sev, outcome = "Adverse events", sample_size_var = "Overall participants ITT")

res_sae_itt_c_sev   <- run_full_nma(dat_arm_c_sev, outcome = "serious adverse events", sample_size_var = "Overall participants ITT")

res_arae_itt_c_sev  <- run_full_nma(dat_arm_c_sev, outcome = "drug-related adverse events", sample_size_var = "Overall participants ITT")
```

## Clinically Evaluable (CE) Population

```r
res_clin_cep_c_sev <- run_full_nma(
  dat_arm_c_sev,
  outcome = "Treatment failure PP",
  sample_size_var = "Overall clinically eligible PP"
)
```

## Microbiologically Evaluable (ME) Population

```r
res_clin_mep_c_sev <- run_full_nma(
  dat_arm_c_sev,
  outcome = "Failure to Eradication microbiologically evaluable",
  sample_size_var = "Microbiologically evaluable PP"
)
```

---

# Outputs

Each analysis automatically generates:

* Bayesian NMA relative treatment effects
* Forest plots
* SUCRA rankings
* Posterior probability of being best
* Closeness-to-best probabilities
* Bayesian non-inferiority analyses
* Network plots
* Heterogeneity diagnostics
* Consistency diagnostics (node-splitting and inconsistency models)
* MCMC convergence diagnostics
* Direct-comparison homogeneity assessments
* Transitivity assessments
* Excel and image outputs for publication-ready reporting

```
```
