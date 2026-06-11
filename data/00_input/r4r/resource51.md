---
id: 51
title: Baysian Phase II Single Arm Clinical Trials
description: >
    This application supports the conduct of simulations for a single arm phase II clinical trial of a new treatment T with time to event outcome variable.

website: https://brpnci.shinyapps.io/Phase2sim_survival_v1/
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
For simulations of single-arm trials for treatment (T),  assume that one interim analysis is conducted after nti events have occurred. By generating nti exponentially distributed event times (with rate parameter log(2)/mt), the parameters for the posterior gamma distribution of λt can be updated as αt = 1+nti and βt = 0.001+Σxti where Σxti is calculated as the sum of the nti simulated event times. For the simulations, users should specify the observed values of nc, mc, assumed value for mt, nt and the proportion of total events after which the interim analysis is to be conducted. h, ε and the number of trials to simulate should also be specified by the user. It is assumed that every subject is followed until event occurrence.