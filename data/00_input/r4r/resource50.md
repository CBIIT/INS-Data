---
id: 50
title: Design, Monitoring, and Analysis of Bayesian Basket Discovery Trials
description: >
    This web-based application implements a new method based on Bayesian principles for the design, monitoring, and analysis of basket trials.

website: https://brpnci.shinyapps.io/BasketTrials/
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
This application currently supports sample size planning and interim analysis for basket trials. The user specifies the number of strata and the prior probability that the drug is active in any particular stratum. Also, since in general there would be uncertainty on whether these strata are completely correlated or independent with regard to the distribution of drug activity, the user also specifies a prior probability corresponding to correlation of drug activity among strata. The more general situation where there are asymmetries among the strata with regard to the prior probability of drug activity will be implemented in the next version of this application.