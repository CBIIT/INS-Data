---
id: 132
title: Gene-wise Prior Bayesian Group Factor Analysis (GBGFA)
description: >
  GBGFA explicitly models gene-centric dependencies when integrating genomic alterations data of the same gene from different platforms (e.g. copy number variation, gene expression and mutation data) to prioritize genes supported by multiple inputs. 
  
website: https://github.com/olganikolova/gbgfa
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
  - email: olga.nikolova@gmail.com
    name:
      firstname: Olga
      lastname: Nikolova
  - email: nikolova@ohsu.edu
---
GBGFA explicitly models gene-centric dependencies when integrating genomic alterations data of the same gene from different platforms (e.g. copy number variation, gene expression and mutation data) to prioritize genes supported by multiple inputs. The multitask approach of this algorithm provides the ability to leverage similarities in the response profiles of drug groups, that are more likely to correspond to true biological effects.
