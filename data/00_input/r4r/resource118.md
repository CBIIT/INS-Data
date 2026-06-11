---
id: 118
title: Evaluation of Differential DependencY (EDDY)
description: >
  EDDY is a statistical test for estimating differential dependencies for a set of genes between two conditions. 
  
website: https://biocomputing.tgen.org/software/EDDY/
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
  - email: gspeyer@tgen.org
    name:
      firstname: Gil
      lastname: Speyer
---
EDDY is a statistical test for estimating differential dependencies for a set of genes between two conditions. Dependencies can be represented and assessed graphically for the expression of a gene set within a particular cellular context. EDDY then calculates the divergence between the probability distributions of scored graphs for each condition. Finally, the statistical significance of this divergence is computed.