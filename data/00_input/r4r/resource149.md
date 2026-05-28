---
id: 149
title: Master Regulator Inference algorithm (MARINa)
description: >
  MARINa is an algorithm that could be used to identfy transcription factors (TFs) that control the transition between two cellular phenotypes.

website: http://califano.c2b2.columbia.edu/marina
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
  - email: af2202@cumc.columbia.edu
    name:
      firstname: Aristidis
      lastname: Floratos
  - email: ac2248@cumc.columbia.edu
    name:
      firstname: Andrea
      lastname: Califano
---
Phenotypic changes eﬀected by pathophysiological events are captured by gene expression proﬁle measurements, determining mRNA abundance on a genome-wide scale in a cellular population. Furthermore, mRNA expression does not constitute a reliable predictor of protein activity, as it fails to capture a variety of post-transcriptional and post-translational events that are involved in its modulation. To negate this problem, MARINa uses the transcriptional targets of each TF as a multiplexed reporter assay to infer the TFs controlling the transition between cellular phenotypes. This task is performed by computing the effect that enrichment of each regulon (i.e., its activated and repressed targets) has on the differentially expressed genes between two phenotypic states.