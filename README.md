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
    - Add a program.program_id to each grant.
        -Because not all Key Programs have an ID, acronym, etc., IDs are created using a form of the Key Program name where spaces and non-alphanumeric characters are removed. (e.g. `BarrettsEsophagusTranslationalResearchNetworkBETRNet`)
    - Combine grants data from all programs and store as a versioned `project.tsv` within the `data/processed/` directory.

5. **Generate summary statistics**
    - Build reports useful for testing and validation but not intended for ingestion into the site
        - `grantsStatsByProgram.csv` groups grants data by Key Program and aggregates counts of grants, projects, searched values, and earliest fiscal year
        - `sharedProjectsByProgramPair.csv` lists pairs of Key Programs and counts of projects that are associated with both

## Usage

1. Clone the repo to your local machine.
2. Install either [Conda or Miniconda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/download.html#anaconda-or-miniconda). Setup environment and install packages with `conda env create -f environment.yaml` run in terminal from the INS-Data directory.
3. If necessary, update the Qualtrics CSV received from ODS. Rename and place it in the `data/raw/` folder. It should be in the format `qualtrics_output_{version}_{type}.csv` (e.g. `qualtrics_output_2023-07-19_raw.csv`). If the Qualtrics CSV is updated, also update the values for `QUALTRICS_VERSION` and `QUALTRICS_TYPE` in `config.py` to match the Qualtrics CSV as needed.
5. In the command terminal, run `python main.py` from the INS-Data root directory. This will run all steps of the workflow and save output files. 

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
│   │   └── ...other Qualtrics versions/types
│   └── processed/
│       ├── {version}/
│       │   ├── api-gathered-{gathering date}/
│       │   │   ├── project.tsv
│       │   └── ...other data gathered on different dates
│       └── ...other Qualtrics versions
├── modules/
│   ├── clean_grants_data.py
│   ├── data_preparation.py
│   ├── nih_reporter_api.py
│   └── summary_statistics.py
├── notebooks/
│   └── Non-production Jupyter notebooks used during development
├── reports/
│   ├── {version}/
│   │   ├── api-gathered-{gathering date}/
│   │   │   ├── grantStatsByProgram.csv
│   │   │   ├── sharedProjectsByProgramPair.csv
│   │   │   └── ... additional reports as needed
│   │   └── ...other reports for data gathered on different dates
│   └── ...other Qualtrics versions
├── config.py
├── environment.yaml
├── main.py
└── README.md
```
