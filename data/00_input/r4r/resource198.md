---
id: 198
title: Interactive Risk Attributable Program (IRAP)
description: >
  IRAP is an interactive computer program that calculates both point estimates and confidence intervals for attributable risk, based on regression models. Modeling consists of three conceptual phases: creating a library of variables, defining the model, and running the model. 
  
website: https://dceg.cancer.gov/tools/risk-assessment/irap
toolTypes:
  - toolType: community_research_tools/risk_assessment
researchAreas:
  - researchArea: causes_of_cancer
researchTypes:
  - researchType: epidemiologic
resourceAccess:
  type: register
  notes: License Agreement
docs:
  - doc: dceg
poc:
  - email: gailm@mail.nih.gov
    name:
      prefix: Dr.
      firstname: Mitchell
      lastname: Gail
---
IRAP keeps a "library" of variables that the user has defined. A variable can take its value directly from the data file or it can be defined in terms of other variables in the library. Typically, the user will create a separate library for each data file that he or she uses. A library can be saved to disk, then loaded whenever the corresponding data file is used.
Once a library of variables has been created (or loaded from removable media), the user may define a model. IRAP supports six data sampling methods:

* Simple Random Sampling
* Stratified Random Sampling
* Frequency Matching
* Individual Matching
* Cohort Analysis
* Cross-Sectional Analysis

After choosing the appropriate sampling method, the user defines a model by specifying the appropriate variable names. A model can contain both discrete and continuous variables. IRAP correctly parameterizes the discrete variables, so the user need not explicitly create "dummy variables."
Once the model has been specified, the user tells IRAP to fit the model. The program asks for the names of the input and output files. It then reads the data and fits the model.
While these three phases are conceptually distinct, they can be interwoven. For example, after fitting a model, the user may decide to combine two levels of exposure into one. It is easy to define a new variable, then rerun the model, substituting the new exposure variable for the old.