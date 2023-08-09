# INS-Data Repository

Welcome to the INS-Data repository for the [Index of NCI Studies (INS)](https://studycatalog.cancer.gov/)! This repository is designed to use a curated list of Key NCI Programs to build a list of all associated extramural grants using the [NIH RePORTER API](https://api.reporter.nih.gov/). These associated grants can then be used to connect research outputs back to Key Programs. 

## Table of Contents

- [Introduction](#introduction)
- [Usage](#usage)
- [Data Structure](#data-structure)

## Introduction

The INS-Data repository workflow follows the general outline below:  
1. **Process Qualtrics CSV**
    - Receive curated CSV of Key Programs from the NCI Office of Data Sharing (ODS). This CSV is an export of survey results from the Qualtrics survey tool. Each Key Program in this export includes Notices of Funding Opportunities (NOFOs) and/or Grant IDs (in long or short form).
    - Save cleaned Key Programs CSV for reference and downstream use
    - The fields expected are defined and can be modified in `config.py`
    - Function(s) defined in `data_preparation.py` module

2. **Get grants data from NIH RePORTER API**
    - For each Key Program, query the NIH RePORTER API to gather a list of all associated extramural grants along with descriptive data for each grant. 
    - The NOFOs (e.g. `RFA-CA-21-038`; `PAR21-346`) and/or grant IDs (e.g. `1 U24 CA274274-01`; `P50CA221745`; `3U24CA055727-26S1`) provided for each Key Program are used as the query. 
    - The following exclusion are also applied within the query:
        - Subprojects are excluded
        - Grants prior to fiscal year 2000 are excluded
        - Grants where NCI is not the administrative agency are excluded
    - Function(s) defined in `nih_reporter_api.py` module

3. **Process grants data**
    - Reformat the data received from the NIH RePORTER API for use within INS. 
        - Remove extraneous fields and rename fields to match the existing INS data model.
        - Flatten nested JSON structures. In particular, the PI, PO, and agency funding fields have this structure. 
        - Format names to standardize capitalization
    - The fields expected are defined and can be modified in `config.py`
    - Function(s) defined in `clean_grants_data.py` module

4. **Save grants data for each Key Program**
    - Save the grants data associated with each Key Program into a TSV. Store those in a versioned folder within the `data/processed/` directory. 
    - Because not all Key Programs have an ID, acronym, etc., each file is saved using a form of the Key Program name where spaces and non-alphanumeric characters are removed. (e.g. `BarrettsEsophagusTranslationalResearchNetworkBETRNet.tsv`)
    - These TSVs are intended to be available for future processing/cleaning modules or directly loaded into the INS database.


## Usage

1. Clone the repo to your local machine.
2. Add or update the Qualtrics CSV received from ODS to the `data/raw/` folder. It should be in the format `qualtrics_output_{version}_{type}.csv` (e.g. `qualtrics_output_2023-07-19_raw.csv`)
3. Update the values for `QUALTRICS_VERSION` and `QUALTRICS_TYPE` in `config.py` to match the Qualtrics CSV as needed.
4. In the command terminal, run `python main.py` from the INS-Data root directory. This will run all steps of the workflow and save output files. 

## Structure

The entire workflow is captured within `main.py` for simplicity of use and reproducibility of data. Functions used within `main.py` are defined in scripts within the `modules/` directory to allow for additional processing steps if needed in the future. 

```
INS-Data
├── data/
│   ├── cleaned/
│   │   ├── key_programs_{version}.csv
│   │   └── ...other versions
│   ├── raw/
│   │   ├── qualtrics_output_{version}_{type}.csv
│   │   └── ...other versions/types
│   └── processed/
│       ├── {version}/
│       │   ├── KeyProgramName01.tsv
│       │   ├── KeyProgramName02.tsv
│       │   └── ...other grant data TSVs for each Key Program
│       └── ...other versions
├── modules/
│   ├── clean_grants_data.py
│   ├── data_preparation.py
│   └── nih_reporter_api.py
├── notebooks/
│   └── Non-production Jupyter notebooks used to test development
├── config.py
└── main.py
```
