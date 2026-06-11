---
id: 128
title: CERES
description: >
  CERES is a method to infer gene essentiality from genome-wide CRISPR-Cas9 screens in cancer cell lines to correct the copy number effect. This approach decreases the false-positive results while taking into account the anti-proliferative copy-number effect.
  
website: https://depmap.org/ceres/
toolTypes:
  - toolType: lab_tools/assays
  - toolType: analysis_tools/genomic_analysis
  - toolType: lab_tools/cell_lines
researchAreas:
  - researchArea: cancer_biology
  - researchArea: cancer_omics
researchTypes:
  - researchType: basic
  - researchType: translational
resourceAccess:
  type: open
docs:
  - doc: ccg
poc:
  - email: David.Quigley@ucsf.edu
    name:
      firstname: David
      lastname: Quigley
---
Studies have shown that genome-wide CRISPR-Cas9 inactivation of genes that are amplified need different analytical approaches for interpretation of the results. The Cas9 induces double strand breaks which lead to false-positive results. A computational method, CERES was developed for inferring gene essentiality from genome-wide CRISPR-Cas9 screens in cancer cell lines to correct the copy number effect. This approach decreases the false-positive results while taking into account the anti-proliferative copy-number effect.