"""
gather_cedcd_data.py
2025-05-14 ZD

This module defines primary function gather_cedcd_data that will process a CSV
of cohorts from the Center for Epidemiology Descriptive Cohort Database (CEDCD).
The intermediate cedcd_datasets.csv will then be finalized with 
`package_output_data.py` into a TSV ready for loading within INS.

Note: The "gather" filename is used for consistency with other modules that handle
input data even though this module isn't pulling from automated external sources.

"""

import os
import sys

import pandas as pd
import uuid


# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config



def get_newest_cohort_versions(df:pd.DataFrame):
    """Filter cohorts to keep only the newest version of each unique title.

    Each new CEDCD entry receives a new dataset_id that counts upward. 
    In order to update cohort informaton, CEDCD adds a new entry and redirects
    the old id(s) url to the newest one. In the export for INS, all entries 
    are present and need to be filtered. 

    Args:
        df: Dataframe of CEDCD cohort metadata
    
    Returns:
        Deduplicated dataframe with newest version of each cohort
    
    """

    # Sort newest versions first
    df_sorted = df.sort_values(by='dataset_id', ascending=False)
    
    # Drop rows with duplicate titles but keep the first (newest)
    df_filtered = df_sorted.drop_duplicates(subset='dataset_title', keep='first')

    return df_filtered



def gather_cedcd_data():
    """Process a CSV of CEDCD Cohort metadata into an intermediate CSV for INS.

    Input and output filepaths are handled in config.py
    """

    print(f"\n---\nDATASETS - CEDCD:\n"
          f"Processing CEDCD cohort data...\n---\n")

    # Define paths from config
    input_csv = config.CEDCD_INPUT_CSV
    output_csv = config.CEDCD_INTERMED_CSV

    # Load input CSV received from CEDCD team
    df = pd.read_csv(input_csv)
    print(f"Loaded CEDCD cohort data from {input_csv}")

    # Keep only the newest version of each dataset title
    df = get_newest_cohort_versions(df)

    # Remove leading commas found in some PI values
    df['principal_investigators'] = df['principal_investigators'].str.lstrip(', ')

    # Build CEDCD url by combining base with CEDCD ID
    base_url = 'https://cedcd.nci.nih.gov/cohort?id='
    df['dataset_source_url'] = base_url + df['dataset_id'].astype(str)

    # Add dataset uuid
    df['dataset_uuid'] = df.apply(lambda row: uuid.uuid4(), axis=1)

    # Rename provided columns where names need to be reused
    df.rename(columns={'primary_disease': 'related_diseases'}, inplace = True)

    # Add hard-coded values
    df['type'] = 'cedcd_dataset'
    df['dataset_source_repo'] = 'CEDCD'
    df['primary_disease'] = 'Unspecified'

    # Add empty dataset columns to avoid downstream issues
    df['GPA'] = ''
    df['limitations_for_reuse'] = ''
    df['dataset_pmid'] = ''
    df['funding_source'] = ''
    df['assay_method'] = ''
    df['release_date'] = ''
    df['sample_count'] = ''
    df['related_genes'] = ''

    # Drop extra columns to avoid downstream issues
    df.drop(columns=['cohort_acronym'], inplace=True)

    # Make output directories if they don't exist
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    # Save output
    df.to_csv(output_csv, index=False)

    print(f"\n\nSuccess! CEDCD cohorts saved to {output_csv}.\n"
          f"Total CEDCD dataset records:    {len(df)}")



# Run module as a standalone script when called directly
if __name__ == "__main__":

    print(f"Running {os.path.basename(__file__)} as standalone module...")

    # Run main function for data validation file generation
    gather_cedcd_data()