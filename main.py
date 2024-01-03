"""
main.py
2023-07-26 ZD

Main function for INS-Data gathering processes intended to be run as a 
single command. 

Inputs required:
- Curated Qualtrics CSV of Key Programs
- iCite database download for publication details
- Reviewed and corrected list of invalid NOFOs/Awards (Optional on subsequent runs)
- All other inputs are gathered via API

Ouputs generated:
- data/01_intermediate/
    - CSVs for all gathered data types, including programs, publications, 
    grants, and projects. These will be passed through data packaging before 
    INS ingestion.
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
from modules.gather_program_data import gather_program_data
from modules.gather_grant_data import gather_grant_data
from modules.summary_statistics import get_summary_statistics
from modules.gather_publication_data import gather_publication_data
from modules.package_output_data import package_output_data


def main():
    """Main function for the INS Data Gathering pipeline."""

    # STEP 1: PROGRAMS
    # Gather program data. Load, clean, and validate from curated file.
    programs_df = gather_program_data(config.QUALTRICS_CSV_PATH)

    # STEP 2: GRANTS
    # Gather, format, and save grants data for each Key Program
    all_cleaned_grants = gather_grant_data(programs_df, print_meta=False)

    # STEP 3: STATS
    # Build and save reports describing the programs and grants data
    get_summary_statistics(all_cleaned_grants)

    # STEP 4: PUBLICATIONS
    # Gather, process, and save publication data
    gather_publication_data(all_cleaned_grants, print_meta=False)

    # STEP 5: PACKAGE
    # Final packaging steps to store output files as TSVs
    package_output_data()


if __name__ == "__main__":
    main()
