---
id: 4
title: Driver-gene Inference by Genetical-Genomics and Information Theory (DIGGIT)
description: >
  DIGGIT is a bioconductor package for identifying genetic variants that lie upstream of master regulators and drive cellular phenotypes.

website: http://www.bioconductor.org/packages/release/bioc/html/diggit.html
toolTypes:
  - toolType: analysis_tools/genomic_analysis
  - toolType: analysis_tools/r_software
researchAreas:
  - researchArea: cancer_biology
  - researchArea: cancer_omics
  - researchArea: cancer_treatment
researchTypes:
  - researchType: translational
docs:
  - doc: ccg
resourceAccess:
  type: open
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
Genomic alterations that contribute to aberrant master regulator (MR) activity must be upstream of the MR, although the specific pathways involved may not be known. Master regulators are transcription factors that control the majority of genes differentially expressed between two molecular phenotypes. The DIGGIT package integrates patient-matched genomic mutation and gene expression data with corresponding gene regulatory networks to identify candidate driver mutations that are upstream of and directly perturb master regulators.