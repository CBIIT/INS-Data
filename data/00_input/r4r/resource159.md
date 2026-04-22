---
id: 159
title: Screening Bayesian Evaluation and Analysis Method (ScreenBEAM)
description: >
  ScreenBEAM is an algorithm that measures gene-level activity to assess the effect of high-throughput RNA interference (RNAi) or Clustered Regularly Interspaced Short Palindromic Repeats (CRISPR) screens through Bayesian hierarchical modeling. 
  
website: https://github.com/jyyu/ScreenBEAM
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
  - email: yujiyang@gmail.com
    name:
      firstname: Jiyang
      lastname: Yu
---
ScreenBEAM is an algorithm that measures gene-level activity to assess the effect of high-throughput RNAi or CRISPR screens through Bayesian hierarchical modeling. For both RNAi and CRISPR, multiple shRNAs or sgRNAs (respectively) are used to target a single gene. ScreenBEAM analyzes gene-level activity for the whole set of shRNAs or sgRNAs targeting the same gene (multi-probe analysis) instead of analyzing the effect of each individual shRNA or sgRNA on a given gene. This reduces false positive and negative rates of high-throughput RNAi or CRISPR screens. This algorithm can handle both microarray and next generation sequencing data as input.