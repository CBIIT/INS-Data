---
id: 154
title: Analytic Technique for Assessment of RNAi by Similarity (ATARiS)
description: >
  ATARiS is a computational method designed to analyze the off-target effects in the data generated from RNAi screens.
  
website: https://www.broadinstitute.org/cancer/ataris
toolTypes:
  - toolType: analysis_tools/genomic_analysis
researchAreas:
  - researchArea: cancer_biology
  - researchArea: cancer_omics
  - researchArea: cancer_treatment
researchTypes:
  - researchType: basic
  - researchType: translational
resourceAccess:
  type: register
docs:
  - doc: ccg
poc:
  - email: ataris@broadinstitute.org
---
RNAi reagents designed to target the same gene often induce different degrees of on-target and off-target gene suppression, resulting in inconsistent phenotypes. To address this, ATARiS tries to identify subsets of its RNAi reagents that produce a significantly similar phenotype across the screened smaples. This approach also computes a consistency score that represents the confidence that its observed phenotypic effects are the result of on-target gene suppression.
