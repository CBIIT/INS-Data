---
id: 38
title: Sample Size Planning for Developing Classifiers Using High Dimensional Data
description: >
    This progam provides estimates of the sample size required for a training set to develop classifiers using high demensional data.

website: https://brb.nci.nih.gov/brb/samplesize/samplesize4GE.html
toolTypes:
  - toolType: analysis_tools/genomic_analysis
  - toolType: analysis_tools/statistical_software
  - toolType: community_research_tools/sample_size_calculator
researchAreas:
  - researchArea: cancer_omics
  - researchArea: bioinformatics
researchTypes:
  - researchType: basic
  - researchType: translational
resourceAccess:
  type: open
docs:
  - doc: dctd
pocs: []
---
This program provides estimates of the sample size required for a training set in order to ensure the resulting binary classifier has an expected accuracy within a tolerance of the optimal accuracy. Classifier performance should also be assessed. This can be done by cross-validation (resampling) or by applying the classifier to an independent validation set. The sample sizes given by this program do not provide the sample size required for a validation set. The sample size provided here also does not address the precision of a cross-validated estimate of prediction accuracy. The program only considers sample sizes below 300. If more than 300 samples are needed, then an error message is returned indicating this.