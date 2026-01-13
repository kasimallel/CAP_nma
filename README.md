# Impacts of antibiotic choice on community-acquired pneumonia: A network meta-analysis of primary studies identified through an umbrella review

## Overview
This repository contains the data processing scripts and analytical code used to conduct a network meta-analysis evaluating the comparative effectiveness of antibiotic regimens for the treatment of community-acquired pneumonia (CAP). The analysis is based on primary randomised controlled trials identified through an umbrella review of systematic reviews.

The study synthesises evidence across adult and paediatric populations, with additional stratification by infection severity where data permit.

## Study design
The analysis followed a two-stage evidence synthesis approach:
1. An umbrella review of systematic reviews including randomised controlled trials (RCTs).
2. A network meta-analysis of eligible primary RCTs extracted from the included reviews.

## Outcomes
The primary outcome was clinical cure or treatment failure, evaluated across:
- Intention-to-treat (ITT) populations
- Clinically evaluable populations
- Microbiologically evaluable populations

Secondary outcomes included:
- Adverse events (overall, severe, and drug-related)
- Mortality

## Repository structure
## Methods
- Antibiotic regimens were compared by agent, dose, and duration.
- Analyses were conducted separately for adults and children, and stratified by severe and non-severe infections where reported.
- Risk of bias was assessed using the Cochrane Risk of Bias 2 (RoB-2) tool.
- Network meta-analysis was performed using frequentist or Bayesian methods (as specified in the analysis scripts).

## Reproducibility
All analyses were conducted using reproducible workflows. Scripts are organised to allow stepwise reproduction of data preparation, model fitting, and result generation.

## License
This project is licensed under the MIT License. See the LICENSE file for details.


