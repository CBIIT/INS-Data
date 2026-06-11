---
id: 120
title: MethylMix
description: >
  MethylMix is an algorithm to identify hyper- and hypomethylated genes for a disease. 
  
website: https://bioconductor.org/packages/3.1/bioc/html/MethylMix.html
toolTypes:
  - toolType: analysis_tools/genomic_analysis
  - toolType: analysis_tools/r_software
researchAreas:
  - researchArea: cancer_biology
  - researchArea: cancer_omics
  - researchArea: cancer_treatment
researchTypes:
  - researchType: basic
  - researchType: translational
resourceAccess:
  type: open
docs:
  - doc: ccg
poc:
  - email: olivier.gevaert@stanford.edu
    name:
      firstname: Olivier
      lastname: Gevaert
---
MethylMix uses a novel statistic, the differential methylation value (DM-value), to define methylation-driven subgroups. This algorithm could be used to identify differentially and transcriptionally predictive methylated genes within a disease by comparing them with the normal DNA methylation state. Matched gene expression data may also be used to identify functional methylation states by focsuing on methylation changes that effect gene expression.