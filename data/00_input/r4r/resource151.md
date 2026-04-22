---
id: 151
title: Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER)
description: >
  The VIPER algorithm allows computational inference of protein activity on an individual sample from gene expression data. 
  
website: http://califano.c2b2.columbia.edu/viper
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
  type: open
docs:
  - doc: ccg
poc:
  - email: af2202@cumc.columbia.edu
    name:
      firstname: Aristidis
      lastname: Floratos
  - email: ac2248@cumc.columbia.edu
    name:
      firstname: Andrea
      lastname: Califano
---
Methods to measure protein abundance on a proteome-wide scale using arrays or mass spectrometry technologies cover only a fraction of proteins, requiring large amounts of tissue, and do not directly capture protein activity. The VIPER algorithm uses the transcripts most directly affected by protein’s activityto rank relative protein activity on a sample-by-sample  basis by transforming a gene expression matrix into a protein activity matrix.
