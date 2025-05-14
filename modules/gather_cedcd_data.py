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


# Append the project's root directory to the Python path
# This allows for importing config when running as part of main.py or alone
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config


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

    # Transform step placeholders

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