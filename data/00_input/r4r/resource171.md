---
id: 171
title: Standardized Occupation Coding for Computer-assisted Epidemiologic Research (SOCcer)
description: >
  SOCcer is software to code free-text job descriptions into standardized occupation classification codes to assist researchers in incorporating occupational risk factors into their studies.
  
website: https://dceg.cancer.gov/tools/design/soccer-tool
researchAreas:
  - researchArea: causes_of_cancer
researchTypes:
  - researchType: epidemiologic
resourceAccess:
  type: open
docs:
  - doc: dceg
poc:
  - email: NCISOCcerWebAdmin@mail.nih.gov
---
SOCcer is a publicly available application that was developed to assist epidemiological researchers incorporate occupational risk into their studies. The application is not intended to replace expert coders, but rather prioritizes job descriptions that would most benefit from expert coders. Low scoring job descriptions are more likely to require expert review than high scoring job description. The coding is performed using an ensemble classifier, which combines the results of multiple classifiers to produce a single classifier that performs better than any single classifier in the ensemble.