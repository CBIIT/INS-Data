"""
main.py
2023-07-26 ZD

Main function for the INS-Data project meant to be run as a single command. 

Inputs required:
- Curated Qualtrics CSV of Key Programs
- iCite database download for publication details
- Reviewed and corrected list of invalid NOFOs/Awards (Optional on subsequent runs)
- All other inputs are gathered via API

Ouputs generated:
- data/01_intermediate/
    - CSVs for all gathered data types, including programs, publications, and 
    projects. These will be passed through data packaging before INS ingestion.
    - Additional CSV outputs and checkpoint files gathered during the process. 
    These are useful for troubleshooting and validation.

- data/02_output/
    - Finalized TSVs for all gathered data types. These are ready for ingestion
    into INS

- reports/:
    - Summary statistic report csvs with high-levelgrant data
    - Report of publications removed from publication.csv with reasons
"""

import os
import sys

import pandas as pd

import config
from modules.data_preparation import load_and_clean_programs
from modules.gather_grants_data import gather_grants_data
from modules.summary_statistics import get_summary_statistics
from modules.gather_publication_data import gather_publication_data


def main():
    """Main function for the INS Data Gathering pipeline."""

    # STEP 1: PROGRAMS
    print(f"Loading and processing {config.QUALTRICS_CSV_PATH}...")

    # Load and clean Key Programs
    continue_bool, key_programs_df = load_and_clean_programs(
                                    csv_filepath = config.QUALTRICS_CSV_PATH, 
                                    col_dict = config.QUALTRICS_COLS)

    # If issues are found, user is prompted to stop or continue
    if continue_bool == True:
        # Export cleaned Key Programs file
        key_programs_df.to_csv(config.CLEANED_KEY_PROGRAMS_CSV, index=False)
        print(f"Success! Saved {config.CLEANED_KEY_PROGRAMS_CSV}.")

    else:
        sys.exit("Process ended by user. Cleaned Key Programs csv not saved. "
            "Make manual edits to input CSV and retry.")

    # STEP 2: GRANTS
    # Gather, format, and save grants data for each Key Program
    all_cleaned_grants = gather_grants_data(key_programs_df, print_meta=False)

    # STEP 3: STATS
    # Build and save reports describing the programs and grants data
    get_summary_statistics(all_cleaned_grants)

    # STEP 4: PUBLICATIONS
    # Gather, process, and save publication data
    gather_publication_data(all_cleaned_grants, print_meta=False)


if __name__ == "__main__":
    main()
