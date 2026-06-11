---
id: 277
title: Automated Germline Variant Pathogenicity (AutoGVP)
description: >
  Automated Germline Variant Pathogenicity (AutoGVP) integrates ClinVar variant annotation with a modified InterVar classification approach, based on American College of Medical Genetics-Association for Molecular Pathology (ACMG-AMP) guidelines, to generate germline variant classification. Since AutoGVP input only requires a VCF file, it can facilitate large-scale, clinically focused classification of germline sequence variants.
website: https://dceg.cancer.gov/tools/analysis/autogvp
toolTypes:
  - toolType: community_research_tools/risk_assessment
  - toolType: clinical_research_tools/guidelines_protocols
  - toolType: datasets_databases/epidemiologic_data
  - toolType: analysis_tools/genomic_analysis
  - toolType: analysis_tools/r_software
researchAreas:
  - researchArea: screening_detection
  - researchArea: causes_of_cancer
researchTypes:
  - researchType: epidemiologic
resourceAccess:
  type: open
docs:
  - doc: cbiit
poc:
  - email: jung.kim2@nih.gov
    name:
      prefix: Dr.
      firstname: Jung
      lastname: Kim
    title: Staff Scientist
---
Automated Germline Variant Pathogenicity (AutoGVP) integrates ClinVar variant annotation with a modified InterVar classification approach, based on American College of Medical Genetics-Association for Molecular Pathology (ACMG-AMP) guidelines, to generate germline variant classification. Since AutoGVP input only requires a VCF file, it can facilitate large-scale, clinically focused classification of germline sequence variants.   
Resource Full Description: Automated Germline Variant Pathogenicity (AutoGVP) was developed to integrate ClinVar variant annotation with a modified InterVar classification approach based on American College of Medical Genetics-Association for Molecular Pathology (ACMG-AMP) guidelines to output germline variant classification. Since AutoGVP input only requires a VCF file, it can facilitate large-scale, clinically focused classification of germline sequence variants.   
Using automated germline classification reduces hands-on time as well as allows for reproducibility of variant classification.  
AutoGVP is an open source dockerized workflow implemented in R and freely available on GitHub. 