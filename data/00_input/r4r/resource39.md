---
id: 39
title: Biomarker Stratified Randomized Design (Binary)
description: >
    This program performs stratified design with a prospective analysis plan and binary endpoint.

website: https://brb.nci.nih.gov/brb/samplesize/sdpap_binary.html
toolTypes:
  - toolType: analysis_tools/statistical_software
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
A randomized trial comparing new treatment (T) to control (C) includes both classifier-positive and classifier-negative patients. This program presumes  availability of binary classifier for biomarker predictive of benefit for new treatment. Sample size calculation for three analysis plans are provided. The first analysis plan determines sample size for the overall test comparing T to C for all randomized patients at reduced two-sided level alpha.  The second determines sample size for comparing T to C in the classifier-positive subset at two-sided .05 level.  For the third plan, first test for interaction between size of treatment effect and subset (classifier + or classifier -). If the interaction is non-significant, just compare treatments overall at two-sided significance level .05, otherwise, compare treatments within subsets at two-sided .05 level.