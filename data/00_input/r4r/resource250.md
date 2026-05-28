---
id: 250
title: Algorithm for the Reconstruction of Accurate Cellular Networks (ARACNe)
description: >
  ARACNe is an algorithm for inferring direct regulatory relationships between transcriptional regulator proteins and target genes.

website: https://califano.c2b2.columbia.edu/aracne
toolTypes:
  - toolType: datasets_databases/genomic_datasets
  - toolType: analysis_tools/genomic_analysis
researchAreas:
  - researchArea: cancer_biology
  - researchArea: cancer_omics
  - researchArea: bioinformatics
researchTypes:
  - researchType: basic
  - researchType: translational
resourceAccess:
  type: open
docs:
  - doc: ocg
  - doc: ccg
poc:
  - email: ac2248@cumc.columbia.edu
    name:
      prefix: Dr.
      firstname: Andrea
      lastname: Califano
    title: Clyde and Helen Wu Professor of Chemical and Systems Biology
---
ARACNe uses microarray expression profiles to reconstruct tissue-specific gene regulatory transcriptional interactions in cellular networks. The method uses mutual information of two random variables to identify genes that are co-expressed, then applies the data processing inequality to filter out interactions that are likely to be indirect. This eliminates the vast majority of false-positive transcriptional interactions typically inferred by pairwise gene-expression correlation analysis. This tool could be used by researchers to determine novel driver genes and drug mechanisms of action.
