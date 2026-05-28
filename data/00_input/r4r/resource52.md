---
id: 52
title: Baysian Phase II Single Arm Clinical Trials
description: >
    This application supports the conduct of simulations for a single-arm phase II clinical trial of a new treatment (T) with a binary outcome variable.

website: https://brpnci.shinyapps.io/Phase2sim_v1/
toolTypes:
  - toolType: analysis_tools/statistical_software
  - toolType: community_research_tools/sample_size_calculator
researchAreas:
  - researchArea: cancer_treatment
researchTypes:
  - researchType: clinical_trials
resourceAccess:
  type: open
docs:
  - doc: dctd
pocs: []
---
The null hypothesis to be tested is that pt = pc against the alternative pt > pc.  If criteria for stopping at interim analysis are not met, additional accrual and outcome are simulated until a total of Nt patients are accrued to T. At the end of the trial, the posterior distribution of pt is updated again and the posterior probability pp_final that pt > pc is recomputed. If pp_final is ≥ 1- ε , then the null hypothesis is rejected, and superiority of T over C is declared. For the simulations, the user needs to specify the observed values for nc, xc, assumed values for pt, Nt and the proportion of total enrollment, after which the interim analysis is to be conducted. δ, ε and the number of trials to simulate are also specified by the user.