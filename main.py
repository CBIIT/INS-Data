# main.py
# 2023-07-26 ZD
#
# Main function for the INS-Data project. This will accept a Qualtrics CSV of 
# Key Programs as input and generate the following outputs:
#   - data/cleaned/: 
#       - Clean `key_programs_{date}.tsv`
#   - data/processed/:
#       - Clean 'project.tsv' containing data for grants associated with key
#         programs
#   - reports/:
#       - Summary statistic report csvs with high-level grant data

import os
import pandas as pd
import config

from modules.data_preparation import load_and_clean_programs
from modules.nih_reporter_api import get_nih_reporter_grants
from modules.clean_grants_data import clean_grants_data
from modules.summary_statistics import get_summary_statistics

def main():


    # STEP 1: Data Prep - Load and clean Key Programs data
    print(f"Loading and processing {config.QUALTRICS_CSV_PATH}...")

    # Load and clean Key Programs
    key_programs_df = load_and_clean_programs(
        csv_filepath = config.QUALTRICS_CSV_PATH, 
        col_dict = config.QUALTRICS_COLS)
        # TODO: add detection of bad NOFO lists (e.g. ;; or PA18-;91)
        # TODO: Add detection for column shifts in data entry 
        # (e.g. Awards listed in NOFO column when no NOFO present)

    # Export cleaned Key Programs file
    key_programs_df.to_csv(config.CLEANED_KEY_PROGRAMS_CSV, index=False)
    print(f"Success! Saved {config.CLEANED_KEY_PROGRAMS_CSV}.")


    # STEP 2: NIH RePORTER API - Get grants info for each Key Program
    print(f"---\nGathering, cleaning, and saving grants data from NIH "
          f"Reporter API for each Key Program...")

    # Create empty df to fill with grants
    all_cleaned_grants = pd.DataFrame()

    # Initialize a counter for total records
    total_records_count = 0

    # Iterate through each Key Program to get name, nofos, and awards
    for index, row in key_programs_df.iterrows():
        program_name = row['program_name']
        nofo_str = row['nofo']
        award_str = row['award']

        print(f"---\n{program_name}")
        
        # Convert semicolon-separated award string into listed values
        if isinstance(award_str, str): 
            award_list = [x.strip() for x in award_str.split(';')]
        else: award_list = []
        # Convert semicolon-separated NOFO string into listed values
        if isinstance(nofo_str, str): 
            nofo_list = [x.strip() for x in nofo_str.split(';')]
        else: nofo_list = []        

        # If neither NOFO nor Awards provided
        if pd.isna(nofo_str) and pd.isna(award_str):
            print(f"No NOFOs or Awards defined. Skipping.")
            continue

        # Gather grants data from NIH RePORTER
        award_grants_data = get_nih_reporter_grants(award_list, 'award', 
                                                    print_meta=True)
        nofo_grants_data = get_nih_reporter_grants(nofo_list, 'nofo', 
                                                   print_meta=True)
        
        # Combine and report
        grants_data = award_grants_data + nofo_grants_data
        total_records_count = total_records_count + len(grants_data)
        print(f"Complete. {len(grants_data)} grants found.")
        if award_grants_data and nofo_grants_data:
            print(f"Both Awards and NOFOs provided. Results may include "
                  f"duplicate grants gathered by both NOFO and Award "
                  f"searches.")
            continue

        # Skip remaining steps if no grants data gathered
        if len(grants_data) == 0:
            print(f"NOTICE: No grants found for {program_name}")
            continue

        # STEP 3: Data Cleaning - Clean and format grants data

        # Run cleaning steps
        cleaned_grants_data = clean_grants_data(grants_data)

        # Add program id column from alphanumeric program name
        program_id = ''.join(filter(str.isalnum, program_name))
        cleaned_grants_data[config.PROGRAM_ID_FIELDNAME] = program_id

    # STEP 4: Combine and save cleaned grants data into single tsv

        all_cleaned_grants = pd.concat([all_cleaned_grants,
                                        cleaned_grants_data])   


    # Define versioned output directory using config.py
    clean_grants_directory = config.PROCESSED_DIR
    if not os.path.exists(clean_grants_directory):
        os.makedirs(clean_grants_directory)
    project_filename = os.path.join(clean_grants_directory,
                                    config.PROJECTS_OUTPUT_FILENAME)

    # Sort by program and project for consistency
    all_cleaned_grants.sort_values(by=[config.PROGRAM_ID_FIELDNAME,
                                       config.API_FIELD_RENAMER.values()[0]
                                       ], inplace=True, ignore_index=True)

    # Export to csv
    all_cleaned_grants.to_csv(project_filename, sep = '\t', index=False)

    
    print(f"---\nSuccess! NIH RePORTER API data gathered, cleaned, and saved.\n"
          f"{total_records_count} grants gathered across all Awards and NOFOs.\n"
          f"{len(all_cleaned_grants)} NCI-funded grants retained for INS. \n"
          f"Results can be found in {clean_grants_directory}.\n---") 

        # TODO: Save all the printed console output to versioned txt 


    # STEP 5: Stats Report - Generate and output summary stats
    
    # Run summary statistic module that processes and exports reports as csv
    get_summary_statistics(all_cleaned_grants)
    print(f"Summary statistic reports for grants data successfully generated. \n"
          f"Results can be found in {config.REPORTS_DIR}.\n---")


if __name__ == "__main__":
    main()
