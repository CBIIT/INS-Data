---
id: 47
title: Pathological Complete Response (pCR)Trial-Level Surrogate Analysis Software
description: >
    This software can be used to access pCR trial-level surrogate analysis.

website: https://brb.nci.nih.gov/programdownload/pCRsoftware.html
toolTypes:
  - toolType: analysis_tools/r_software
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
This is an R source package that contains the code, data, documentation, and results of the analysis reported in "Assessing pathological complete response as a trial-level surrogate endpoint for early-stage breast cancer" and supplement, by E. L. Korn, M. C. Sachs, and L. M. McShane, Ann Oncol 27: 10-15, 2016. The package is organized as follows:

* DESCRIPTION details the analysis and documents the dependencies;
* The /data subdirectory contains the data files in .csv and .RData formats; the datasets are documented in the /man directory.
* The /inst/doc directory contains the code and results of the analysis.
* The /R directory has the R functions that were used in the analysis.
* The /man directory contains the documentation in Rd format. View the documentation in R by typing ?pcrmeta.