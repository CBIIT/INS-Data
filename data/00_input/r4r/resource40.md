---
id: 40
title: Biomarker Stratified Randomized Design (Time-to-Event)
description: >
    This program performs stratified design with prospective analysis plan and time-to-event endpoint.

website: https://brb.nci.nih.gov/brb/samplesize/sdpap_survival.html
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
A randomized trial comparing new treatment (T) to control (C) includes both classifier-positive and classifier-negative patients. The program presumes availability of time-to-event predictive of benefit for new treatment. Sample size calculation for three analysis plans are provided: A) Determine sample size for overall test comparing T to C for all randomized patients at reduced two-sided level alpha.  B) Determine sample size for comparing T to C in classifier positive subset at two-sided .05 level.  C) First test for interaction between size of treatment effect and subset (classifier + or classifier -). If interaction is non-significant, just compare treatments overall at two-sided significance level .05. Otherwise, compare treatments within subsets at two-sided .05 level.